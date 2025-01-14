#ifndef ABLATELIBRARY_SZBURN_HPP
#define ABLATELIBRARY_SZBURN_HPP
#include "particles/processes/burningProcess.hpp"
#include "particles/particleSolver.hpp"
#include "eos/eos.hpp"
#include "eos/zerork.hpp"
#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
#include "particles/processes/burningModel/dropletFlame/cxhyozFlame.hpp"
#include <vector>
#include <functional>


namespace ablate::particles::processes::burningModel {



class SZBurn : public ablate::particles::processes::BurningProcess
{
    public:
    PetscReal K; //D^2 Law constant K
    const PetscReal convectionCoeff; //simple convection heating mode
    PetscReal Td;
    double Dp;
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


    double MWair = 28.96;
    double mDot;

    struct Constants{
        const double qlimfac=0.1;
    };
    Constants dropletConstants;

    struct resultsStruct{
        double F;
        double Yfs_new;
        double Ts;
        double mdot;
        double rstar; //need this for later...
        double Tf;
        double K;
    };
    resultsStruct result;

    struct TransportConstStruct{
        const double pr = 0.707;
        const double muo = 1.716e-5;
        const double to = 273.e+0;
        const double so = 111.e+0;
        const double sc = 0.707;
    };
    TransportConstStruct TransportConst;

    struct farFieldProp{
        double rox;
        double Temperature;
        double Pressure = 101325;
        double Yox = 1;
        double kg;
        double Cpg;

    };
    farFieldProp farField;


    std::string fuelType;


    double BisectionMethodSolve(std::function<double(double)> func, double a, double b, double tol, int maxIter);

    // store the name of species used in the ode solver
    const std::vector<std::string> speciesNames;
    std::vector<int> speciesOffSet;

    void CalcBurnRate();

    double SolveSZBurn(double YFsguess);

    void CalcEvapRate();

    double SolveSZEvap(double YFsguess);

    void UpdateFarfield(double T, double P, double YiO2);

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
#endif //ABLATELIBRARY_SZBURN_HPP
