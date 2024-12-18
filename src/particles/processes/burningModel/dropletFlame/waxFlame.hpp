#ifndef ABLATELIBRARY_WAXFLAME_H
#define ABLATELIBRARY_WAXFLAME_H
#include "cxhyozFlame.hpp"
#include "math.h"
namespace ablate::particles::processes::burningModel {


class waxFlame : public ablate::particles::processes::burningModel::CxHyOzFlame{

   public:
    // Constructor with initializer list
    waxFlame()     // Initialize Antoine constants (A, B, C)
    {
        x=31;
        y=64;
        z=0;


        double YO2_inf=0.23;

        nuO2 = (x+y/4.-z/2.)*MWO2/MWFuel;
        nuCO2 = x*MWCO2/MWFuel;
        nuH2O = 0.5*y*MWH2O/MWFuel;


        //TODO calculate these properties with Canetra for wax
        //Adiabatic Flame temperature
        Tad=2700;
        //Heat of combustion
        Hc=42E6;


    }

    public:
     const std::string fuelname = "wax";

};
}
#endif  // ABLATELIBRARY_WAXFLAME_H
