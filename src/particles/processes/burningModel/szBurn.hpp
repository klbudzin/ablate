#ifndef ABLATELIBRARY_SZBURN_HPP
#define ABLATELIBRARY_SZBURN_HPP
#include "particles/processes/burningProcess.hpp"
#include "particles/particleSolver.hpp"
#include "eos/eos.hpp"
#include "eos/zerork.hpp"
#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
#include "particles/processes/burningModel/dropletFlame/cxhyozFlame.hpp"


namespace ablate::particles::processes::burningModel {


class SZBurn : public ablate::particles::processes::BurningProcess
{
    public:
    const PetscReal burnRate; //D^2 Law constant K
    const PetscReal convectionCoeff; //simple convection heating mode
//    const PetscReal YiFuel[];

    // the eos used to species the species and compute properties
//    std::shared_ptr<eos::zerorkEOS> eos;
    SZBurn(PetscReal convectionCoeff, PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx,
           PetscReal Lv, PetscReal Hc, const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
           PetscReal extinguishmentOxygenMassFraction,std::shared_ptr<eos::EOS> eos,
           std::string  fuelin);

    private:
    std::shared_ptr<ablate::particles::processes::burningModel::LiquidFuel> liquidFuel;
    std::shared_ptr<ablate::particles::processes::burningModel::CxHyOzFlame> flame;




    struct farFieldProp{
        double rox;
        double Temperature;
        double Pressure = 101325;
        //TODO could set up later to include a vector for chemical equilibrium...
//        std::vector<double> Y(int nSpc);
        double Yox = 1;

    };

    farFieldProp farField;
    const std::string fuelType;



    void CalcBurnRate();

    void SolveSZBurn(double* YFsguess,std::vector<double>* res,ablate::particles::processes::burningModel::SZBurn::farFieldProp farfield);

    void UpdateFarfield(farFieldProp* FarField,double T, double P, double YiO2);

    void TwoZoneTransport(double* T,double* density,double* Yi[],double* Cp,double* k,double* D);


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

     double nSpec;
     double YFs;
     double MWair = 28.96;
     double Dp;
     double Ts;
     double ql;
     double mDot;
     double K;
     double qlimfac=0.1;
     int nSpc=10;
     const double Pivalue=3.14159;

     struct results{
         double F;
         double Yfs_new;
         double Ts;
         double mdot1;
         double rstar;
         double Tf;
     };

     struct TransportConstsstruct{
         inline const static PetscReal pr = 0.707;
         inline const static PetscReal muo = 1.716e-5;
         inline const static PetscReal to = 273.e+0;
         inline const static PetscReal so = 111.e+0;
         inline const static PetscReal sc = 0.707;
     };

     TransportConstsstruct TransportConst;


};



}
#endif //ABLATELIBRARY_SZBURN_HPP
