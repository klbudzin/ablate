#include "stagnationEnthalpy.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "utilities/mathUtilities.hpp"

ablate::boundarySolver::lodi::StagnationEnthalpy::StagnationEnthalpy(std::shared_ptr<eos::EOS> eos, PetscReal stagEnthalpy, std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling)
    : LODIBoundary(std::move(eos), std::move(pressureGradientScaling)), stagnationEnthalpy(stagEnthalpy) {}

void ablate::boundarySolver::lodi::StagnationEnthalpy::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    ablate::boundarySolver::lodi::LODIBoundary::Setup(bSolver);
    bSolver.RegisterFunction(StagnationEnthalpyFunction, this, fieldNames, fieldNames, {});
    bSolver.RegisterFunction(UpdateEnergy, this, {finiteVolume::CompressibleFlowFields::EULER_FIELD}, {finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD});
}

PetscErrorCode ablate::boundarySolver::lodi::StagnationEnthalpy::StagnationEnthalpyFunction(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                            const PetscFVCellGeom *boundaryCell, const PetscInt uOff[], const PetscScalar *boundaryValues,
                                                                                            const PetscScalar *stencilValues[], const PetscInt aOff[], const PetscScalar *auxValues,
                                                                                            const PetscScalar *stencilAuxValues[], PetscInt stencilSize, const PetscInt stencil[],
                                                                                            const PetscScalar stencilWeights[], const PetscInt sOff[], PetscScalar source[], void *ctx) {
    PetscFunctionBeginUser;
    auto stagInlet = (StagnationEnthalpy *)ctx;

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
        boundaryDensity = boundaryValues[uOff[stagInlet->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            boundaryVel[d] = boundaryValues[uOff[stagInlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / boundaryDensity;
            boundaryNormalVelocity += boundaryVel[d] * fg->normal[d];
        }
        PetscCall(stagInlet->computeTemperature.function(boundaryValues, &boundaryTemperature, stagInlet->computeTemperature.context.get()));
        PetscCall(stagInlet->computeSpeedOfSound.function(boundaryValues, boundaryTemperature, &boundarySpeedOfSound, stagInlet->computeSpeedOfSound.context.get()));
        PetscCall(stagInlet->computePressureFromTemperature.function(boundaryValues, boundaryTemperature, &boundaryPressure, stagInlet->computePressureFromTemperature.context.get()));
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
        stencilDensity[s] = stencilValues[s][uOff[stagInlet->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            stencilVel[s][d] = stencilValues[s][uOff[stagInlet->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / stencilDensity[s];
            stencilNormalVelocity[s] += stencilVel[s][d] * fg->normal[d];
        }
        PetscCall(stagInlet->computePressure.function(stencilValues[s], &stencilPressure[s], stagInlet->computePressure.context.get()));
    }

    // Interpolate the normal velocity gradient to the surface
    PetscScalar dVeldNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryNormalVelocity, stencilSize, &stencilNormalVelocity[0], stencilWeights, dVeldNorm);
    PetscScalar dPdNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryPressure, stencilSize, &stencilPressure[0], stencilWeights, dPdNorm);

    PetscReal boundaryCp, boundaryCv;
    stagInlet->computeSpecificHeatConstantPressure.function(boundaryValues, boundaryTemperature, &boundaryCp, stagInlet->computeSpecificHeatConstantPressure.context.get());
    stagInlet->computeSpecificHeatConstantVolume.function(boundaryValues, boundaryTemperature, &boundaryCv, stagInlet->computeSpecificHeatConstantVolume.context.get());

    // Compute the enthalpy
    PetscReal boundarySensibleEnthalpy;
    stagInlet->computeSensibleEnthalpyFunction.function(boundaryValues, boundaryTemperature, &boundarySensibleEnthalpy, stagInlet->computeSensibleEnthalpyFunction.context.get());

    // get_vel_and_c_prims(PGS, velwall[0], C, Cp, Cv, velnprm, Cprm);
    PetscReal velNormPrim, speedOfSoundPrim;
    stagInlet->GetVelAndCPrims(boundaryNormalVelocity, boundarySpeedOfSound, boundaryCp, boundaryCv, velNormPrim, speedOfSoundPrim);

    // get_eigenvalues
    std::vector<PetscReal> lambda(stagInlet->nEqs);
    stagInlet->GetEigenValues(boundaryNormalVelocity, boundarySpeedOfSound, velNormPrim, speedOfSoundPrim, &lambda[0]);

    // Compute alpha2
    PetscReal alpha2 = 1.0;
    if (stagInlet->pressureGradientScaling) {
        alpha2 = PetscSqr(stagInlet->pressureGradientScaling->GetAlpha());
    }

    // Get scriptL
    std::vector<PetscReal> scriptL(stagInlet->nEqs);
    scriptL[1 + dim] = lambda[1 + dim] * (dPdNorm - boundaryDensity * dVeldNorm * alpha2 * (velNormPrim - boundaryNormalVelocity - speedOfSoundPrim));  // Outgoing Accoustic Wave
    PetscScalar term1, term2;
    term1 = (boundaryCp / boundaryCv + 1)/(boundaryCp / boundaryCv - 1) * (velNormPrim -boundaryNormalVelocity) / (speedOfSoundPrim);
    term2 = boundaryNormalVelocity/(speedOfSoundPrim*alpha2);
    scriptL[0] = -scriptL[1 + dim]*(1+term1+term2)/(1-term1-term2);  // Incoming acoustic wave
    scriptL[1] = 0; // Entropy wave, 0 because of isentropic assumption with BC derivation
    for (int d = 1; d < dim; d++) {
        scriptL[1 + d] = 0.e+0;  // Tangential velocities
    }
    // Species
    for (int ns = 0; ns < stagInlet->nSpecEqs; ns++) {
        scriptL[2 + dim + ns] = 0.e+0;
    }
    // Extra variables
    for (int ne = 0; ne < stagInlet->nEvEqs; ne++) {
        scriptL[2 + dim + stagInlet->nSpecEqs + ne] = 0.e+0;
    }

    // Directly compute the source terms, note that this may be problem in the future with multiple source terms on the same boundary cell
    stagInlet->GetmdFdn(sOff,
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
//Function to ensure energy at inlet coincides with the correct stagnation enthalpy
PetscErrorCode ablate::boundarySolver::lodi::StagnationEnthalpy::UpdateEnergy(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg, const PetscFVCellGeom *boundaryCell,
                                                                           const PetscInt *uOff, PetscScalar *boundaryValues, const PetscScalar *stencilValues, const PetscInt *aOff,
                                                                           PetscScalar *auxValues, const PetscScalar *stencilAuxValues, void *ctx) {
    PetscFunctionBeginUser;
    auto stagInlet = (StagnationEnthalpy *)ctx;
    PetscScalar boundaryPressure;
    PetscScalar boundaryDensity = boundaryValues[uOff[0] + RHO];
    PetscCall(stagInlet->computePressureFromTemperature.function(boundaryValues, auxValues[aOff[0]], &boundaryPressure, stagInlet->computePressureFromTemperature.context.get()));
    //Total energy = h_stag - P/\rho |_{boundary}
    boundaryValues[uOff[0] + RHOE] = boundaryDensity*stagInlet->stagnationEnthalpy - boundaryPressure ;
    PetscFunctionReturn(0);
}

#include "registrar.hpp"
REGISTER(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::lodi::StagnationEnthalpy, "Enforces a stagnation enthalpy driven inlet",
         ARG(ablate::eos::EOS, "eos", "The EOS describing the flow field at the wall"),
         ARG(PetscReal, "stagnationEnthalpy", "The stagnation enthalpy that is driving the inlet"),
         OPT(ablate::finiteVolume::processes::PressureGradientScaling, "pgs", "Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows"));
