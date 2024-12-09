#include "twoZoneBurn.hpp"
#include "utilities/constants.hpp"

ablate::particles::processes::burningModel::TwoZoneBurn::TwoZoneBurn(PetscReal convectionCoeff,
       PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx, PetscReal Lv, PetscReal heatOfCombustion,
       const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
       PetscReal extinguishmentOxygenMassFraction, std::shared_ptr<eos::EOS> eosIn) :
       BurningProcess(std::move(eosIn), {}, massFractionsProducts, ignitionTemperature, nuOx, Lv, heatOfCombustion, extinguishmentOxygenMassFraction),
       burnRate(burnRate), convectionCoeff(convectionCoeff)
       { }


void ablate::particles::processes::burningModel::TwoZoneBurn::ComputeRHS(PetscReal time, accessors::SwarmAccessor &swarmAccessor, accessors::RhsAccessor &rhsAccessor, accessors::EulerianAccessor &eulerianAccessor)
{
    //Grab wanted fields from the eulerian field accessor
    //In this case we only need the Oxygen and temperature fields
    auto farFieldTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
//    auto farFieldSpecies = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD];

    //Grab The fields of the particle that we will be changing
    auto massRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleMass];
    //Burning particles
    auto tempRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleTemperature];

    //Grab the Particle Fields we need
    auto partDensity = swarmAccessor[ablate::particles::ParticleSolver::ParticleDensity];
    auto particlesPerParcel = swarmAccessor[ablate::particles::ParticleSolver::ParticleNPP];
    auto partDiameterAvg = swarmAccessor[ablate::particles::ParticleSolver::ParticleDiameter];
    auto partBurning = swarmAccessor[ablate::particles::ParticleSolver::ParticleBurning];
    auto partTemperature = swarmAccessor[ablate::particles::ParticleSolver::ParticleTemperature];
    auto partCp = swarmAccessor[ablate::particles::ParticleSolver::ParticleCP];
    auto parcelMass = swarmAccessor[ablate::particles::ParticleSolver::ParticleMass];
    auto numParticles = swarmAccessor.GetNumberParticles();
    PetscReal SAtot;
    //Calculate Each Particles RHS's
    for( auto np = 0; np < numParticles; np++) {
        //Ideally this is a farField Value, but for this test case I'm not letting the particles burn until the particle itself reaches
        // the ignition temperature, Particle may be burning now then when it was updated in the decode call, this will be called twice
        if (partBurning(np))
            massRHS(np) -= particlesPerParcel(np)*(partDensity(np)*PETSC_PI/4*partDiameterAvg(np)*burnRate);
        //If it's not burning, we are just heating with our shitty convection model
        else {
            SAtot = particlesPerParcel(np)*PETSC_PI*std::pow(partDiameterAvg(np),2);
            tempRHS(np) += convectionCoeff*(farFieldTemperature(np)-partTemperature(np))*SAtot/(parcelMass(np)*partCp(np));
            }
    }
}


void ablate::particles::processes::burningModel::TwoZoneBurn::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor &swarmAccessorPreStep, accessors::SwarmAccessor &swarmAccessorPostStep, accessors::EulerianSourceAccessor &eulerianSourceAccessor)
{
    //We need to communicate the energy that was transfer to/from the eulerian field to the particle back to the eulerian field rhoEnergy
    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    //The mass and energy  are in the eulerian field (rhoE is the second index (1), rho is the first index)
    auto sourceEuler = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];
    auto sourceSpecies = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD];

    //Grab Next and Old Swarm Temperatures (Assume Cp and mass does not change ever)
    auto parcelTempNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelTempOld = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelMassNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    auto parcelMassOld = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    //Determine if the particle was burning or not
    auto parcelBurning = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleBurning];

    auto parcelCp = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleCP];
    PetscReal mdot;
    //Product Species will be the same for all particles so calculate them here
    //The below should work since the function should be a constant function and not dependent on space or time

    //Compute the source terms due to burning/ heat transfer
    for (PetscInt p = 0; p < np; ++p) {
        if (parcelBurning(p)) {
            // Start by adding in the mass source terms
            mdot = parcelMassOld(p)-parcelMassNew(p); //mass going into farfield
            sourceEuler(p,ablate::finiteVolume::CompressibleFlowFields::RHO) += mdot;
            // now onto the species mass
            //Now to take care of species terms
            for (PetscInt ns = 0; ns < numberSpecies; ns ++)
                sourceSpecies(p,ns) += mdot*(1+properties.nuOx)*properties.massFractionsProducts[ns]; //add the products
            sourceSpecies(p,oxygenOffset) -= properties.nuOx*mdot; //remove the oxygen
            //add in the burning energy
            sourceEuler(p,1) -= mdot*properties.HC;
        }
        else {
            //Since this is assumed no mass loss unless its burning, add in the particle heating/cooling coupling
            sourceEuler(p,1) += parcelMassNew(p)*parcelCp(p)*(parcelTempNew(p)-parcelTempOld(p));
        }
    }
}


#include "registrar.hpp"
REGISTER(ablate::particles::processes::Process, ablate::particles::processes::burningModel::TwoZoneBurn, "Example of an no evaporation/Burning Model",
         ARG(PetscReal, "convectionCoefficient", "The convection Coefficient for the simple heating mode"),
         OPT(PetscReal, "IgnitionTemperature", "The temperature of the far field where it is said to *ignite*"),
         OPT(PetscReal, "burnRate", "The constant burning rate K for the d^2 law (defaults to 0)"),
         OPT(PetscReal, "nuOx", "stoichiometric oxidizer mass coefficient (defaults to 0)"),
         OPT(PetscReal, "Lv", "LatentHeat of Vaporization (defaults to 1e5)"),
         OPT(PetscReal, "Hc", "heat of combustion per unit mass of fuel (defaults to 1e6"),
         ARG(ablate::mathFunctions::FieldFunction, "productsMassFractions", "The product species functions"),
         OPT(PetscReal, "extinguishmentOxygenMassFraction", "The minimum oxygen mass fraction needed for burning (defaults to 0.1)"),
         ARG(ablate::eos::EOS, "eos", "the eos used to compute various gas properties") );
