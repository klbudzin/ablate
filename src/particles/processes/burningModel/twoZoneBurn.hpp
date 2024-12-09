#ifndef ABLATELIBRARY_TWOZONEBURN_HPP
#define ABLATELIBRARY_TWOZONEBURN_HPP
#include "particles/processes/burningProcess.hpp"
#include "particles/particleSolver.hpp"

namespace ablate::particles::processes::burningModel {


class TwoZoneBurn : public ablate::particles::processes::BurningProcess
{
    public:
    const PetscReal burnRate; //D^2 Law constant K
    const PetscReal convectionCoeff; //simple convection heating mode
    const PetscReal YiFuel[];

    TwoZoneBurn(PetscReal convectionCoeff, PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx,
                              PetscReal Lv, PetscReal Hc, const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
                              PetscReal extinguishmentOxygenMassFraction,std::shared_ptr<eos::EOS> eos);


    /**
     * Overload default calls to be do nothing calls
     */
    void ComputeRHS(PetscReal time, accessors::SwarmAccessor& swarmAccessor, accessors::RhsAccessor& rhsAccessor, accessors::EulerianAccessor& eulerianAccessor) override;

    /**
     * overload default calls to be do nothing calls (might change this to be default normal compute eulerian source
     */
    void ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor& swarmAccessorPreStep, accessors::SwarmAccessor& swarmAccessorPostStep,
                               accessors::EulerianSourceAccessor& eulerianSourceAccessor) override;

    /** overide the update particle burning method because I want it based on particle temperature for this process
     * *
     */
     void UpdateParticleBurning(accessors::EulerianAccessor& eA, accessors::MutableSwarmAccessor& mA) override {
        auto numParticles = mA.GetNumberParticles();
        auto burning = mA[ablate::particles::ParticleSolver::ParticleBurning];
        auto farFieldMassFractions = eA[ablate::finiteVolume::CompressibleFlowFields::YI_FIELD];
        auto particleTemperature = mA[ablate::particles::ParticleSolver::ParticleTemperature];
        for (PetscInt p = 0; p < numParticles; p++){
            IsParticleBurning(farFieldMassFractions(p,oxygenOffset), particleTemperature(p), &burning(p));
        }
     }
};


}
#endif //ABLATELIBRARY_TWOZONEBURN_HPP
