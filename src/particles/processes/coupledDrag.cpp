#include "coupledDrag.hpp"
#include "particles/burningParticleSolver.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "eos/transport/transportModel.hpp"
#include "utilities/constants.hpp"

ablate::particles::processes::CoupledDrag::CoupledDrag(std::shared_ptr<ablate::eos::transport::TransportModel> transportModel, const std::vector<PetscReal>& gravity)
    : gravity(gravity) {
    if (!transportModel)
        throw std::invalid_argument("A transport model is needed to use the coupled drag model.");
    else
        this->viscosityFunction = transportModel->GetTransportTemperatureFunction(eos::transport::TransportProperty::Viscosity, {});
}


void ablate::particles::processes::CoupledDrag::ComputeRHS(PetscReal time, ablate::particles::accessors::SwarmAccessor& swarmAccessor, ablate::particles::accessors::RhsAccessor& rhsAccessor,
                ablate::particles::accessors::EulerianAccessor& eulerianAccessor) {

    //Actual Source terms (Not an analytical solution)

    //Parcel/particle Properties
    auto partVelocity = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleVelocity];
    auto partDensity = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleDensity];
    auto avgDiameter = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleDiameter];
    const auto Nparticles = swarmAccessor.GetNumberParticles();

    //Grab the gas temperature and density
    const auto dim = eulerianAccessor.GetDimensions();
    auto gasEuler = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];
    auto gasTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
    auto gasVelocity = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::VELOCITY_FIELD];


    //Setup the particle RHS
    auto particleVelocityRHS = rhsAccessor[ablate::particles::BurningParticleSolver::ParticleVelocity];
    auto particleCoordsRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleCoordinates];
    //declare some variables
    PetscReal gasMu,ReP, Cd, tauRelax, gasDensity, velDiffMag =0;
    PetscReal vd[dim];
    for(PetscInt p = 0; p < Nparticles; p++) {
        gasDensity = gasEuler(p,ablate::finiteVolume::CompressibleFlowFields::RHO);
        for( PetscInt d=0; d < dim; d++) {
            vd[d] = gasVelocity(p,d) - partVelocity(p,d);
            velDiffMag += vd[d]*vd[d];
        }
        this->viscosityFunction.function(NULL, gasTemperature(p), &gasMu, this->viscosityFunction.context.get());
        velDiffMag = PetscSqrtReal(velDiffMag);
        ReP = velDiffMag * avgDiameter(p)* gasDensity / gasMu;
        Cd = 24. * (1. + 0.15 * PetscPowReal(ReP, 0.687) ) / ( ReP + ablate::utilities::Constants::tiny) +
                0.42/(1.+42500./PetscPowReal(ReP + ablate::utilities::Constants::tiny, 1.16));
        tauRelax = PetscMin(4./3. * partDensity(p) / gasDensity * avgDiameter(p) / (Cd*velDiffMag), ablate::utilities::Constants::large);
        //Add sources to rhs
        for( PetscInt d=0; d < dim; d++) {
            particleCoordsRHS(p,d) += partVelocity(p,d);
            particleVelocityRHS(p,d) += gravity[d] + vd[d]/tauRelax;
        }
    }
}

void ablate::particles::processes::CoupledDrag::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, ablate::particles::accessors::SwarmAccessor& swarmAccessorPreStep,
                                                                                  ablate::particles::accessors::SwarmAccessor& swarmAccessorPostStep,
                                                                                  ablate::particles::accessors::EulerianSourceAccessor& eulerianSourceAccessor) {
    //We need to communicate the momentum that was transfer to/from the eulerian field to the particle back to the eulerian field rhoEnergy
    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    auto source = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];

    //Evaluate the timestep
    PetscReal dt = endTime-startTime;

    //Grab Next and Old Swarm Temperatures (Assume Cp and mass does not change ever)
    auto parcelMassNew = swarmAccessorPostStep[ablate::particles::BurningParticleSolver::ParticleMass];
    auto parcelMassOld = swarmAccessorPreStep[ablate::particles::BurningParticleSolver::ParticleMass];
    auto particleVelocityOld = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleVelocity];
    auto particleVelocityNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleVelocity];

    // Compute the function
    for (PetscInt p = 0; p < np; ++p) {
        for (PetscInt d = 0; d < ndims; ++d) {
            source(p,ablate::finiteVolume::CompressibleFlowFields::RHOU+d) -= (parcelMassNew(p)*particleVelocityNew(p,d) - parcelMassOld(p)*particleVelocityOld(p,d)) - (parcelMassNew(p)+parcelMassOld(p))*dt/2.*gravity[d];
        }
    }
}

#include "registrar.hpp"
REGISTER(ablate::particles::processes::Process, ablate::particles::processes::CoupledDrag, "Couples Drag to fluid",
         ARG(ablate::eos::transport::TransportModel, "transport", "The transport model used to get the fluid viscosity"),
         ARG(std::vector<PetscReal>, "gravity", "The gravitational vector"));
