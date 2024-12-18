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
        nuC32H66=1;
        nuCO2=1;
        nuH2O=1;
        nuO2=1;

        Tad=2700;

        Hc=42E6;
    }

    public:
     const std::string fuelname = "wax";

};
}
#endif  // ABLATELIBRARY_WAXFLAME_H
