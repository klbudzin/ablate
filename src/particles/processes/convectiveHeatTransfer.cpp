#include "convectiveHeatTransfer.hpp"
//#include "particles/particleSolver.hpp"
#include "particles/burningParticleSolver.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"

#include <utility>

ablate::particles::processes::ConvectiveHeatTransfer::ConvectiveHeatTransfer(PetscReal h)
    : convCoefficient(h) {

}


void ablate::particles::processes::ConvectiveHeatTransfer::ComputeRHS(PetscReal time, ablate::particles::accessors::SwarmAccessor& swarmAccessor, ablate::particles::accessors::RhsAccessor& rhsAccessor,
                ablate::particles::accessors::EulerianAccessor& eulerianAccessor) {
    // $\pder{T_Parcel}{t} = 1/(Mass_{parcel}*cp)*h(T_{E}-T_{parcel}) SA_{tot} $
    //SA_tot needs avgDiameter and particlesperparcel
    auto avgDiameter = swarmAccessor[ablate::particles::ParticleSolver::ParticleDiameter];
    auto particlesPerParcel = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleNPP];

    //Parcel Properties
    auto cp = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleCP];
    auto massParcel = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleMass];
    auto Tp = swarmAccessor[ablate::particles::BurningParticleSolver::ParticleTemperature];

    //Grab the eulerian Temperature
    auto eulerianTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];

    //Setup the particle RHS
    const auto Nparticles = swarmAccessor.GetNumberParticles();
    auto particleTemperatureRHS = rhsAccessor[ablate::particles::BurningParticleSolver::ParticleTemperature];
    PetscReal SAtot;
    for(PetscInt p = 0; p < Nparticles; p++) {
        SAtot = particlesPerParcel(p)*PETSC_PI*std::pow(avgDiameter(p),2);
        particleTemperatureRHS(p) = convCoefficient*(eulerianTemperature(p)-Tp(p))*SAtot/(massParcel(p)*cp(p));
    }
}

void ablate::particles::processes::ConvectiveHeatTransfer::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, ablate::particles::accessors::SwarmAccessor& swarmAccessorPreStep,
                                                                                  ablate::particles::accessors::SwarmAccessor& swarmAccessorPostStep,
                                                                                  ablate::particles::accessors::EulerianSourceAccessor& eulerianSourceAccessor) {
    //We need to communicate the energy that was transfer to/from the eulerian field to the particle back to the eulerian field rhoEnergy
    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    //The Energy field is in the eulerian field (rhoE is the second index (1))
    auto source = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];

    //Grab Next and Old Swarm Temperatures (Assume Cp and mass does not change ever)
    auto parcelTempNew = swarmAccessorPostStep[ablate::particles::BurningParticleSolver::ParticleTemperature];
    auto parcelTempOld = swarmAccessorPreStep[ablate::particles::BurningParticleSolver::ParticleTemperature];
    auto parcelMass = swarmAccessorPostStep[ablate::particles::BurningParticleSolver::ParticleMass];
    auto parcelCp = swarmAccessorPostStep[ablate::particles::BurningParticleSolver::ParticleCP];
    PetscReal sourceVal;
    // Compute the function
    for (PetscInt p = 0; p < np; ++p) {
        sourceVal = parcelMass(p)*parcelCp(p)*(parcelTempNew(p) - parcelTempOld(p));
        //The 1 is to access the rhoE source field
        source(p,1) -= sourceVal;
    }
}

#include "registrar.hpp"
REGISTER(ablate::particles::processes::Process, ablate::particles::processes::ConvectiveHeatTransfer, "Convectively cools are heats solid particles (Assumes no change in Mass or Cp)",
         ARG(PetscReal, "convectionCoefficient", "The convection coefficient"));
