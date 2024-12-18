#ifndef ABLATELIBRARY_PARTICLE_BURNINGPROCESS_HPP
#define ABLATELIBRARY_PARTICLE_BURNINGPROCESS_HPP

#include "particles/accessors/eulerianSourceAccessor.hpp"
#include "particles/processes/coupledProcess.hpp"
#include "eos/zerork.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "particles/particleSolver.hpp"

namespace ablate::particles::processes{

struct DropletMaterialProperties {
    //Holds droplet material properties such as, it's vapor/ burning product species, oxidizer mass coefficient, laten heat, and heat of combustion
    PetscReal* massFractionsProducts;
    PetscReal* massFractionsVapor;
    const PetscReal nuOx; //Oxidizer stoichiometric mass coefficient
    const PetscReal Lv; //Latent Heat of Vaporization
    const PetscReal HC; //Heat of Combustion per unit mass of fuel, i.e. particle

    /** The Constructor
     *
     * @param massFractionsVapor
     * @param massFractionsProducts
     * @param nuOx
     * @param LatentHeat
     */
    DropletMaterialProperties(PetscInt Nspecies, const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsVaporIn,
                      const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProductsIn,
                      PetscReal nuOx, PetscReal LatentHeat, PetscReal heatOfCombustion): nuOx(nuOx), Lv(LatentHeat), HC(heatOfCombustion)
        {
            massFractionsVapor = new PetscReal[Nspecies];
            massFractionsProducts = new PetscReal[Nspecies];
            //update the new vapor/products arrays to their set values
            if (massFractionsVaporIn) {
                mathFunctions::PetscFunction func = massFractionsVaporIn->GetFieldFunction()->GetPetscFunction();
                func(0,0.0,nullptr,Nspecies, massFractionsVapor, massFractionsVaporIn->GetFieldFunction()->GetContext());
            } else {
                // If not defined, it better not be used! Set all species to 0 except final numerical dilute species
                for(int ns = 0; ns < (Nspecies); ns++) massFractionsVapor[ns] = 0;
                massFractionsVapor[Nspecies-1] = 1;
            }
            if (massFractionsProductsIn) {
                mathFunctions::PetscFunction func = massFractionsProductsIn->GetFieldFunction()->GetPetscFunction();
                func(0,0.0,nullptr,Nspecies, massFractionsProducts, massFractionsProductsIn->GetFieldFunction()->GetContext());
            } else {
                // If not defined, it better not be used! Set all species to 0 except final numerical dilute species
                for(int ns = 0; ns < (Nspecies); ns++) massFractionsVapor[ns] = 0;
                massFractionsVapor[Nspecies-1] = 1;
            }
        }
    /**
     * Empty Default constructor
     */
//    DropletMaterialProperties() = default;
    ~DropletMaterialProperties() = default;
};

/**
 * interface used for coupled processes to compute the eulerian source terms
 */
class BurningProcess : public CoupledProcess {
    protected:
        //eos needed to keep track of species and do any particular needed burning calculations (cp etc.)
        std::shared_ptr<eos::zerorkEOS> eos;
        PetscInt numberSpecies = 0;
        PetscInt oxygenOffset = -1;// offset for the oxygen species

        //Keep track of droplet material properties
        DropletMaterialProperties properties;
        //Different thresholds for ignition/burning
        PetscReal ignitionTemperature; //Temperature threshold for when the particle ignites in air
        PetscReal extinguishmentOxygenMassFraction; //The mass fraction of oxygen threshold for when the particle extinguishes

        virtual void IsParticleBurning(PetscReal YO2inf, PetscReal TCheck, double *burning){
            if (burning && (YO2inf < extinguishmentOxygenMassFraction) ) burning[0] = PETSC_FALSE;
            else if ((TCheck >= ignitionTemperature) && (YO2inf >= extinguishmentOxygenMassFraction)) burning[0] = PETSC_TRUE;
        }

    public:

        BurningProcess(std::shared_ptr<eos::zerorkEOS> eosIn,
                     const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsVapor,
                     const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
                     PetscReal ignitionTemperature, PetscReal nuOx, PetscReal LatentHeat,
                     PetscReal heatOfCombustion, PetscReal extinguishmentOxygenMassFraction)
                     : eos(std::move(eosIn)),
                       properties(eos->GetSpeciesVariables().size(), massFractionsVapor, massFractionsProducts, nuOx, LatentHeat, heatOfCombustion),
                       ignitionTemperature(ignitionTemperature ? ignitionTemperature : 0),
                       extinguishmentOxygenMassFraction(extinguishmentOxygenMassFraction ? extinguishmentOxygenMassFraction : 0.1)
       {
            if (!std::dynamic_pointer_cast<eos::zerorkEOS>(eosIn)) {
                throw std::invalid_argument("The twozone model only accepts eos::zerorkEOS as the input sorry in advance, Kolos");
            }

            // determine the offset for the O2 species in the eos, and if there is none error out
            auto speciesList = eos->GetSpeciesVariables();
            numberSpecies = speciesList.size();
            for( PetscInt idx = 0; idx < numberSpecies; idx++ )
                if (speciesList.at(idx) == "O2") {
                    oxygenOffset = idx;
                    continue;
                }

            if (oxygenOffset == -1)
                    throw std::runtime_error("The particle burning model requires an O2 species to work.");
            if (this->extinguishmentOxygenMassFraction < 0 || this->extinguishmentOxygenMassFraction > 1)
                throw std::invalid_argument(" The extinguishment oxygen mass fraction should be between 0 and 1!");
       }

//       BurningProcess() = default;
       ~BurningProcess() = default;

//         /**
//         * Check if the Particle is burning based off the farfield temperature an Oxygen Massfraction
//         */
        virtual void UpdateParticleBurning(accessors::EulerianAccessor& eA, accessors::MutableSwarmAccessor& mA) {
            auto numParticles = mA.GetNumberParticles();
            auto burning = mA[ablate::particles::ParticleSolver::ParticleBurning];
            auto farFieldTemperature = eA[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
            auto farFieldMassFractions = eA[ablate::finiteVolume::CompressibleFlowFields::YI_FIELD];
            for (PetscInt p = 0; p < numParticles; p++){
                IsParticleBurning(farFieldMassFractions(p,oxygenOffset),farFieldTemperature(p), &burning(p));
            }
        }

    /**
     * default to do nothing
     */
    virtual void ComputeRHS(PetscReal time, accessors::SwarmAccessor& swarmAccessor, accessors::RhsAccessor& rhsAccessor, accessors::EulerianAccessor& eulerianAccessor) override {};

    /**
     * default to do nothing (might change to simit eulerian source as default calculator)
     */
    virtual void ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor& swarmAccessorPreStep, accessors::SwarmAccessor& swarmAccessorPostStep,
                               accessors::EulerianSourceAccessor& eulerianSourceAccessor) override {};

    /**
     *  The point of this method is to allow someone to do an analytical solve of the burning process, i.e. this should allow
     *  us the ability to add an analytical solve boolean to the burning particle solver that will bypass any RHS updates in the burning part of the solver
     *  and allow the burning field to analytically update the new mass and temperature of the particle. This currently is not used though
     */
     [[maybe_unused]] virtual void UpdateAuxFields() {};

};

}  // namespace ablate::particles::processes

#endif  // ABLATELIBRARY_PARTICLE_BURNINGPROCESS_HPP