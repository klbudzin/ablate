#include "pressureTankInlet.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "utilities/mathUtilities.hpp"

ablate::boundarySolver::lodi::PressureTankInlet::PressureTankInlet(std::shared_ptr<eos::EOS> eos, PetscReal tankPressure, PetscReal tankTemperature, PetscReal valveArea,
                                                                   PetscReal Rgas, PetscReal gamma, std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling)
    : LODIBoundary(std::move(eos), std::move(pressureGradientScaling)), tankPressure(tankPressure), tankTemperature(tankTemperature),
      valveArea(valveArea), gammaInlet(gamma ? gamma : 1.4), RInlet(Rgas? Rgas : 287.0), tankDensity(tankPressure/RInlet/tankTemperature) {}

void ablate::boundarySolver::lodi::PressureTankInlet::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    ablate::boundarySolver::lodi::LODIBoundary::Setup(bSolver);
    bSolver.RegisterFunction(PressureTankInletFunction, this, fieldNames, fieldNames, {});
    bSolver.RegisterFunction(UpdateConservedVariables, this, fieldNames, {finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD});
}

PetscErrorCode ablate::boundarySolver::lodi::PressureTankInlet::PressureTankInletFunction(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                            const PetscFVCellGeom *boundaryCell, const PetscInt uOff[], const PetscScalar *boundaryValues,
                                                                                            const PetscScalar *stencilValues[], const PetscInt aOff[], const PetscScalar *auxValues,
                                                                                            const PetscScalar *stencilAuxValues[], PetscInt stencilSize, const PetscInt stencil[],
                                                                                            const PetscScalar stencilWeights[], const PetscInt sOff[], PetscScalar source[], void *ctx) {
    PetscFunctionBeginUser;
    auto Inlet = (PressureTankInlet *)ctx;

    // Compute the transformation matrix
    PetscReal transformationMatrix[3][3];
    utilities::MathUtilities::ComputeTransformationMatrix(dim, fg->normal, transformationMatrix);

    // Compute the pressure/values on the boundary
    PetscReal boundaryDensity;
    PetscReal boundaryTemperature;
    PetscReal boundaryVel[3];
    PetscReal boundaryNormalVelocity = 0.0;
    PetscReal boundarySpeedOfSound;
    PetscReal boundaryPressure;

    // Get the velocity and pressure on the surface
    {
        boundaryDensity = boundaryValues[uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            boundaryVel[d] = boundaryValues[uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / boundaryDensity;
            boundaryNormalVelocity += boundaryVel[d] * fg->normal[d];
        }
        PetscCall(Inlet->computeTemperature.function(boundaryValues, &boundaryTemperature, Inlet->computeTemperature.context.get()));
        PetscCall(Inlet->computeSpeedOfSound.function(boundaryValues, boundaryTemperature, &boundarySpeedOfSound, Inlet->computeSpeedOfSound.context.get()));
        PetscCall(Inlet->computePressureFromTemperature.function(boundaryValues, boundaryTemperature, &boundaryPressure, Inlet->computePressureFromTemperature.context.get()));
    }

    // Map the boundary velocity into the normal coord system
    PetscReal boundaryVelNormCord[3];
    utilities::MathUtilities::Multiply(dim, transformationMatrix, boundaryVel, boundaryVelNormCord);

    // Compute each stencil point
    std::vector<PetscReal> stencilDensity(stencilSize);
    std::vector<std::vector<PetscReal>> stencilVel(stencilSize, std::vector<PetscReal>(dim));
    std::vector<PetscReal> stencilNormalVelocity(stencilSize);
    std::vector<PetscReal> stencilPressure(stencilSize);

    for (PetscInt s = 0; s < stencilSize; s++) {
        stencilDensity[s] = stencilValues[s][uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            stencilVel[s][d] = stencilValues[s][uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / stencilDensity[s];
            stencilNormalVelocity[s] += stencilVel[s][d] * fg->normal[d];
        }
        PetscCall(Inlet->computePressure.function(stencilValues[s], &stencilPressure[s], Inlet->computePressure.context.get()));
    }

    // Interpolate the normal velocity gradient to the surface
    PetscScalar dVeldNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryNormalVelocity, stencilSize, &stencilNormalVelocity[0], stencilWeights, dVeldNorm);
    PetscScalar dPdNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryPressure, stencilSize, &stencilPressure[0], stencilWeights, dPdNorm);

    PetscReal boundaryCp, boundaryCv;
    Inlet->computeSpecificHeatConstantPressure.function(boundaryValues, boundaryTemperature, &boundaryCp, Inlet->computeSpecificHeatConstantPressure.context.get());
    Inlet->computeSpecificHeatConstantVolume.function(boundaryValues, boundaryTemperature, &boundaryCv, Inlet->computeSpecificHeatConstantVolume.context.get());

    // Compute the enthalpy
    PetscReal boundarySensibleEnthalpy;
    Inlet->computeSensibleEnthalpyFunction.function(boundaryValues, boundaryTemperature, &boundarySensibleEnthalpy, Inlet->computeSensibleEnthalpyFunction.context.get());

    // get_vel_and_c_prims(PGS, velwall[0], C, Cp, Cv, velnprm, Cprm);
    PetscReal velNormPrim, speedOfSoundPrim;
    Inlet->GetVelAndCPrims(boundaryNormalVelocity, boundarySpeedOfSound, boundaryCp, boundaryCv, velNormPrim, speedOfSoundPrim);

    // get_eigenvalues
    std::vector<PetscReal> lambda(Inlet->nEqs);
    Inlet->GetEigenValues(boundaryNormalVelocity, boundarySpeedOfSound, velNormPrim, speedOfSoundPrim, &lambda[0]);

    // Compute alpha2
    PetscReal alpha2 = 1.0;
    if (Inlet->pressureGradientScaling) {
        alpha2 = PetscSqr(Inlet->pressureGradientScaling->GetAlpha());
    }

    // Get scriptL
    std::vector<PetscReal> scriptL(Inlet->nEqs);
    scriptL[1 + dim] = lambda[1 + dim] * (dPdNorm - boundaryDensity * dVeldNorm * alpha2 * (velNormPrim - boundaryNormalVelocity - speedOfSoundPrim));  // Outgoing Accoustic Wave
    PetscScalar F;
    F = boundaryCp/boundaryCv*boundaryNormalVelocity*speedOfSoundPrim;
    F = F/(boundarySpeedOfSound*boundarySpeedOfSound/alpha2 + boundaryCp/boundaryCv*boundaryNormalVelocity*(velNormPrim-boundaryNormalVelocity));
    F = (1+F)/(1-F);
    scriptL[0] = F*scriptL[1 + dim]; //incoming acoustic wave assuming constant mass flux
//    scriptL[1] = 0; // Entropy wave, 0 because of isentropic assumption with BC derivation
    //Below is isothermal assumption
    scriptL[1] = 0.5 * (boundaryCp/boundaryCv-1) * (scriptL[1+dim] + scriptL[0]) -
                 0.5 * (boundaryCp/boundaryCv+1) * (scriptL[0] - scriptL[1+dim]) * (velNormPrim-boundaryNormalVelocity)/speedOfSoundPrim ;
    for (int d = 1; d < dim; d++) {
        scriptL[1 + d] = 0.e+0;  // Tangential velocities
    }
    // Species
    for (int ns = 0; ns < Inlet->nSpecEqs; ns++) {
        scriptL[2 + dim + ns] = 0.e+0;
    }
    // Extra variables
    for (int ne = 0; ne < Inlet->nEvEqs; ne++) {
        scriptL[2 + dim + Inlet->nSpecEqs + ne] = 0.e+0;
    }

    // Directly compute the source terms, note that this may be problem in the future with multiple source terms on the same boundary cell
    Inlet->GetmdFdn(sOff,
                             boundaryVelNormCord,
                             boundaryDensity,
                             boundaryTemperature,
                             boundaryCp,
                             boundaryCv,
                             boundarySpeedOfSound,
                             boundarySensibleEnthalpy,
                             velNormPrim,
                             speedOfSoundPrim,
                             boundaryValues,
                             uOff,
                             scriptL.data(),
                             transformationMatrix,
                             source);

    PetscFunctionReturn(0);
}
//Function to ensure energy at inlet coincides with the isentropic values predicted from the current state of the tank
PetscErrorCode ablate::boundarySolver::lodi::PressureTankInlet::UpdateConservedVariables(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg, const PetscFVCellGeom *boundaryCell,
                                                                           const PetscInt *uOff, PetscScalar *boundaryValues, const PetscScalar *stencilValues, const PetscInt *aOff,
                                                                           PetscScalar *auxValues, const PetscScalar *stencilAuxValues, void *ctx) {
    PetscFunctionBeginUser;
    auto Inlet = (PressureTankInlet *)ctx;
    //Assuming Temperature is correct at the boundary, update values using the isentropic relations and choked/unchoked mass flux
    PetscReal boundaryDensity = boundaryValues[uOff[Inlet->eulerId] + RHO];
    PetscReal boundaryPressure = Inlet->tankPressure * pow(boundaryDensity/Inlet->tankDensity,Inlet->gammaInlet);
    PetscReal boundaryTemperature = Inlet->tankTemperature * pow(boundaryDensity/Inlet->tankDensity,Inlet->gammaInlet-1);

    // determine if mass flux is choked or not
    PetscReal massFlux, vel, KE = 0;
    if (boundaryPressure/Inlet->tankPressure <= 0.5238) //Fix later to actual value based on gamma....TODO
    { //Flow is choked
        massFlux = Inlet->tankPressure/sqrt(Inlet->tankTemperature)*sqrt(Inlet->gammaInlet/Inlet->RInlet) *
                   pow((Inlet->gammaInlet+1)/2, -(Inlet->gammaInlet+1)/2/(Inlet->gammaInlet-1));
    } else {
        throw std::runtime_error("The Mass flux has become unchoked In the pressure tank inlet! Consider adding the unchocked part now");
    }

    // set boundary mass flux and KE
    for (PetscInt d = 0; d < dim; d++) {
        //The normal should be pointing outside of the domain
        //So the velocity should be -massflux \hat{n}_i
        vel = -massFlux * fg->normal[d];
        boundaryValues[uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] = vel;
        KE += PetscSqr(vel/boundaryDensity);
    }
    KE=KE/2.;
    // set total energy = h_valve + KE - P_valve/\rho_valve
    PetscReal boundarySensibleEnthalpy;
    Inlet->computeSensibleEnthalpyFunction.function(boundaryValues, boundaryTemperature, &boundarySensibleEnthalpy, Inlet->computeSensibleEnthalpyFunction.context.get());
    boundaryValues[uOff[Inlet->eulerId] + RHOE] = boundaryDensity*(KE + boundarySensibleEnthalpy) - boundaryPressure;
    PetscFunctionReturn(0);
}

////Function to ensure energy at inlet coincides with the isentropic values predicted from the current state of the tank
//PetscErrorCode ablate::boundarySolver::lodi::PressureTankInlet::UpdateConservedVariables(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg, const PetscFVCellGeom *boundaryCell,
//                                                                           const PetscInt *uOff, PetscScalar *boundaryValues, const PetscScalar *stencilValues, const PetscInt *aOff,
//                                                                           PetscScalar *auxValues, const PetscScalar *stencilAuxValues, void *ctx) {
//    PetscFunctionBeginUser;
//    auto Inlet = (PressureTankInlet *)ctx;
//    PetscReal boundaryTemperature = auxValues[aOff[0]];
//    //Assuming Temperature is correct at the boundary, update values using the isentropic relations and choked/unchoked mass flux
//    PetscReal boundaryCp, boundaryCv;
//    Inlet->computeSpecificHeatConstantPressure.function(boundaryValues, boundaryTemperature, &boundaryCp, Inlet->computeSpecificHeatConstantPressure.context.get());
//    Inlet->computeSpecificHeatConstantVolume.function(boundaryValues, boundaryTemperature, &boundaryCv, Inlet->computeSpecificHeatConstantVolume.context.get());
//    PetscReal gamma = boundaryCp/boundaryCv;
//    PetscReal R = boundaryCp*(gamma-1)/gamma;
//    PetscReal tankDensity = Inlet->tankPressure/R/Inlet->tankTemperature;
//    PetscReal boundaryDensity = tankDensity*pow(boundaryTemperature/Inlet->tankTemperature,1/(gamma-1));
//    PetscReal boundaryPressure = boundaryDensity*R*boundaryTemperature;
//
//    // determine if mass flux is choked or not
//    PetscReal massFlux, vel, KE = 0;
//    if (boundaryPressure/Inlet->tankPressure <= 0.5238) //Fix later to actual value based on gamma....TODO
//    { //Flow is choked
//        massFlux = Inlet->tankPressure/sqrt(Inlet->tankTemperature)*sqrt(gamma/R) *
//                   pow((gamma+1)/2, -(gamma+1)/2/(gamma-1));
//    } else {
//        throw std::runtime_error("The Mass flux has become unchoked In the pressure tank inlet! Consider adding the unchocked part now");
//    }
//    boundaryValues[uOff[Inlet->eulerID] + finiteVolume::CompressibleFlowFields::RHO] = boundaryDensity;
//    // set boundary mass flux and KE
//    for (PetscInt d = 0; d < dim; d++) {
//        //The normal should be pointing outside of the domain
//        //So the velocity should be -massflux \hat{n}_i
//        vel = -massFlux * fg->normal[d];
//        boundaryValues[uOff[Inlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] = vel;
//        KE += PetscSqr(vel/boundaryDensity);
//    }
//    // set total energy = h_valve + KE - P_valve/\rho_valve
//    PetscReal boundarySensibleEnthalpy;
//    Inlet->computeSensibleEnthalpyFunction.function(boundaryValues, boundaryTemperature, &boundarySensibleEnthalpy, Inlet->computeSensibleEnthalpyFunction.context.get());
//    boundaryValues[uOff[Inlet->eulerId] + RHOE] = boundaryDensity*(KE + boundarySensibleEnthalpy) - boundaryPressure;
//    PetscFunctionReturn(0);
//}

#include "registrar.hpp"
REGISTER(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::lodi::PressureTankInlet, "Enforces a stagnation enthalpy driven inlet",
         ARG(ablate::eos::EOS, "eos", "The EOS describing the flow field at the wall"),
         ARG(PetscReal, "tankPressure", "The pressure of the pressurized tank"),
         ARG(PetscReal, "tankTemperature", "The temperature of the pressurized tank"),
         ARG(PetscReal, "valveArea", "The stagnation enthalpy that is driving the inlet"),
         OPT(PetscReal, "gamma", "The specific heat ratio expected at the valve"),
         OPT(PetscReal, "Rgas", "The specific heat ratio expected at the valve"),
         OPT(ablate::finiteVolume::processes::PressureGradientScaling, "pgs", "Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows"));
