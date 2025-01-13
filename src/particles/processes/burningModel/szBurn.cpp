#include "szBurn.hpp"
#include "utilities/constants.hpp"
#include "particles/processes/burningModel/liquidFuels/waxFuel.hpp"
#include "particles/processes/burningModel/dropletFlame/waxFlame.hpp"
//#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
#include <math.h>

#include <utility>

ablate::particles::processes::burningModel::SZBurn::SZBurn(PetscReal convectionCoeff,
       PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx, PetscReal Lv, PetscReal heatOfCombustion,
       const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
       PetscReal extinguishmentOxygenMassFraction, std::shared_ptr<eos::EOS> eosIn,std::string  fuelin) :
       BurningProcess(std::move(eosIn), {}, massFractionsProducts, ignitionTemperature, nuOx, Lv, heatOfCombustion, extinguishmentOxygenMassFraction),
       burnRate(burnRate), convectionCoeff(convectionCoeff),speciesNames{fuelin, "O2", "N2", "CO2", "H2O"},speciesOffSet(5, -1)

       {

    //one can easily add more fuels here
    //TODO CHANGE THIS to C32H66 after done with debugging!!!!!!!!!
    if (fuelin == "CH4") {
        fuelType="wax";
        liquidFuel = std::make_shared<ablate::particles::processes::burningModel::waxFuel>();
        flame = std::make_shared<ablate::particles::processes::burningModel::waxFlame>();
    } else {
        throw std::invalid_argument("No liquid fuel, flamelet of the fuel with formula: " + fuelin);
    }

    //initialize farfield //TODO make these as optional inputs
    farField.Temperature = 350.0;
    farField.Pressure = 101325.0;
    farField.Yox = 1.0;


    auto speciesList = eos->GetSpeciesVariables();
    numberSpecies = speciesList.size();
    for (int idxsp = 0; idxsp < 5;idxsp++) {
        for (PetscInt idx = 0; idx < numberSpecies; idx++)
            if (speciesList.at(idx) == speciesNames[idxsp]) {
                speciesOffSet[idxsp] = idx;
                continue;
            }

    }

    //Check if all species are found in the mechanism;
    int cnt = count(speciesOffSet.begin(), speciesOffSet.end(), -1);
    if (cnt > 0)
        throw std::runtime_error("The Shvab-Zeldovich model requires fuelFormula,O2,N2,CO2,H2O to work.");

    //Initialize surface temperature to 300
    Ts=300;

}


void ablate::particles::processes::burningModel::SZBurn::ComputeRHS(PetscReal time, accessors::SwarmAccessor &swarmAccessor, accessors::RhsAccessor &rhsAccessor, accessors::EulerianAccessor &eulerianAccessor)
{
    auto partDiam2 = swarmAccessor[ablate::particles::ParticleSolver::ParticleDiameter];
    auto partVel2 = swarmAccessor[ablate::particles::ParticleSolver::ParticleVelocity];
    auto partDens2 = swarmAccessor[ablate::particles::ParticleSolver::ParticleDensity];
    auto velocityRhs2 = rhsAccessor[ablate::particles::ParticleSolver::ParticleVelocity];

    //Grab wanted fields from the eulerian field accessor
    //In this case we only need the Oxygen temperature and pressure fields
    auto farFieldTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
    auto farFieldPressure = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::PRESSURE_FIELD];
    auto farFieldSpecies = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::YI_FIELD];
    //TODO get the actual O2 mass fractions;
    double YiO2=1;
    UpdateFarfield(&farField, farFieldTemperature.values[0],farFieldPressure.values[0],YiO2);


    //TODO figure this thing out...
    //Burning particles
    auto tempRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleTemperature];
    //Grab The fields of the particle that we will be changing
    auto massRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleMass];
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
        CalcBurnRate();
//        double mdot

//        if (partBurning(np))
//            massRHS(np) -= particlesPerParcel(np)*(partDensity(np)*PETSC_PI/4*partDiameterAvg(np)*burnRate);
//        //If it's not burning, we are just heating with our shitty convection model
//        else {
//            SAtot = particlesPerParcel(np)*PETSC_PI*std::pow(partDiameterAvg(np),2);
//            tempRHS(np) += convectionCoeff*(farFieldTemperature(np)-partTemperature(np))*SAtot/(parcelMass(np)*partCp(np));
//            }

    }
}

void ablate::particles::processes::burningModel::SZBurn::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor &swarmAccessorPreStep, accessors::SwarmAccessor &swarmAccessorPostStep, accessors::EulerianSourceAccessor &eulerianSourceAccessor)
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

    PetscReal mdot=mDot;
    //Product Species will be the same for all particles so calculate them here
    //The below should work since the function should be a constant function and not dependent on space or time

    //Compute the source terms due to burning/ heat transfer
    for (PetscInt p = 0; p < np; ++p) {


        if (parcelBurning(p)) {
            //Calculate the droplet burn rate
            CalcBurnRate();

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

    void ablate::particles::processes::burningModel::SZBurn::CalcBurnRate() {

        double Yguess=1;

        SolveSZBurn(&Yguess,&result,&farField);

    //Calculate burnrate, Energy source,
    burnRate = result.K;
}

void ablate::particles::processes::burningModel::SZBurn::UpdateFarfield(ablate::particles::processes::burningModel::SZBurn::farFieldProp* farfield,
                                                                        double Tfar, double Pfar, double YiO2far) {
    //Update properties
    farfield->Temperature=Tfar;
    farfield->Pressure=Pfar;
    farfield->Yox=YiO2far;

    double nN2dnO2=((1-YiO2far)*28)/(YiO2far*32);

    //rox = (x+y/4)*(MWO2+3.76*MWN2)/MWF
    farfield->rox = (flame->x+flame->y/4)*(flame->MWO2+nN2dnO2*flame->MWN2)/flame->MWFuel;

    //Calculate gas transport (conductivity)

    double mug = TransportConst.muo * PetscSqrtReal(Tfar / TransportConst.to) * (Tfar / TransportConst.to) * (TransportConst.to + TransportConst.so) / (Tfar + TransportConst.so);

//    *Cp = eos->mech->getMassCpFromTY();
    farfield->kg = mug * liquidFuel->fuelProperties.Cp / TransportConst.pr;


};


void ablate::particles::processes::burningModel::SZBurn::SolveSZBurn(double* YFsguess,resultsStruct* results,ablate::particles::processes::burningModel::SZBurn::farFieldProp* farfield){

    double YFs_old = *YFsguess;
    double MWs = 1/(YFs_old/liquidFuel->fuelProperties.MW + (1-YFs_old)/MWair);
    double Pvap = YFs_old * farfield->Pressure * MWs / liquidFuel->fuelProperties.MW;

    //Bound the vapor pressure
    Pvap = std::min(Pvap,farfield->Pressure);
    Pvap = std::max(Pvap, 1E-20);

    //calculate the surface temperature
    double Tsnew = 290;
    liquidFuel->Tvap(&Tsnew,&Pvap);

    double ql=0; //TODO implement a model for heat transfer into surface
//    mDot =0; // TODO Need this to calculate heatflux

    double BoxT = (farfield->Yox*flame->Hc/farfield->rox+liquidFuel->fuelProperties.Cp*(farfield->Temperature-Tsnew))/liquidFuel->fuelProperties.Hvap;

    double YFs_new = BoxT/(1+BoxT);

    YFs_new = std::min(YFs_new, 1.);
    YFs_new = std::max(YFs_new, 0.);
    results->Yfs_new=YFs_new;

    results->F=YFs_new-YFs_old;
    results->K=8*farfield->kg/(liquidFuel->fuelProperties.Cp*liquidFuel->fuelProperties.rhol)* log(1+BoxT);
//    results->mdot=
}




#include "registrar.hpp"

REGISTER(ablate::particles::processes::Process, ablate::particles::processes::burningModel::SZBurn, "Example of an no evaporation/Burning Model",
         OPT(PetscReal, "convectionCoefficient", "The convection Coefficient for the simple heating mode"),
         OPT(PetscReal, "IgnitionTemperature", "The temperature of the far field where it is said to *ignite*"),
         OPT(PetscReal, "burnRate", "The constant burning rate K for the d^2 law (defaults to 0)"),
         OPT(PetscReal, "nuOx", "stoichiometric oxidizer mass coefficient (defaults to 0)"),
         OPT(PetscReal, "Lv", "LatentHeat of Vaporization (defaults to 1e5)"),
         OPT(PetscReal, "Hc", "heat of combustion per unit mass of fuel (defaults to 1e6)"),
         OPT(ablate::mathFunctions::FieldFunction, "productsMassFractions", "The product species functions"),
         OPT(PetscReal, "extinguishmentOxygenMassFraction", "The minimum oxygen mass fraction needed for burning (defaults to 0.1)"),
         ARG(ablate::eos::EOS, "eos", "the eos used to compute various gas properties"),
         ARG(std::string, "fuelFormula", "Fomula of the fuel just as it is in the mechanism, curently implmented: 'C32H66'(wax)") );
