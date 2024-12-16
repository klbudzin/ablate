#ifndef ABLATELIBRARY_TWOZONEBURN_HPP
#define ABLATELIBRARY_TWOZONEBURN_HPP
#include "particles/processes/burningProcess.hpp"
#include "particles/particleSolver.hpp"
#include "eos/eos.hpp"
#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
//#include "particles/processes/burningModel/liquidFuels/waxFuel.hpp"

namespace ablate::particles::processes::burningModel {


class TwoZoneBurn : public ablate::particles::processes::BurningProcess
{
    public:
    const PetscReal burnRate; //D^2 Law constant K
    const PetscReal convectionCoeff; //simple convection heating mode
    const PetscReal YiFuel[];



//    std::shared_ptr<eos> eos; //TODO convert EOS into zerorkEOS
    private:
    std::shared_ptr<ablate::particles::processes::burningModel::LiquidFuel> liquidFuel;


    struct fuel{

        double rhoSolid = 678;
        double Tsurf = 364;
        double latentheat=36000;
        double heatOfCombustion=5000;
        double MW = 450;
        double k =12;

    };

    struct farField{

        double Temperature = 350;
        double Pressure;
        std::vector<double> Y(int nSpc);
        double Yox;

    };

    struct innerZone{
        double k;
        double Cp;
        double Diff;
        double rho;
        double T;
        double gamma;
        double Le;
        double Lambda;
        double MdotF_D_Mdot1;

    };

    struct outerZone{
        double k;
        double Cp;
        double Diff;
        double rho;
        double T;
        double gamma;
        double Le;
        double Lambda;
        double MdotOX_D_Mdot2;
    };


    TwoZoneBurn(PetscReal convectionCoeff, PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx,
                              PetscReal Lv, PetscReal Hc, const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
                              PetscReal extinguishmentOxygenMassFraction,std::shared_ptr<eos::EOS> eos,
                              const std::string& fuelType);

    void CalcBurnRate();

    void SolveTwoZone(double YFs,ablate::particles::processes::burningModel::TwoZoneBurn::farField farfield,
                                                                               ablate::particles::processes::burningModel::TwoZoneBurn::fuel fuel,
                                                                               ablate::particles::processes::burningModel::TwoZoneBurn::innerZone innerzone,
                                                                               ablate::particles::processes::burningModel::TwoZoneBurn::outerZone outerzone);

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

    private:


     double YFs;
     double MWs;
     double Ts;
     double ql;
     double qlimfac=0.1;
     int nSpc=10;
     const double Pivalue=3.14159;


};



}
#endif //ABLATELIBRARY_TWOZONEBURN_HPP
