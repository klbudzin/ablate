#ifndef ABLATELIBRARY_WAXFUEL_H
#define ABLATELIBRARY_WAXFUEL_H
#include "liquidFuel.hpp"
#include "math.h"
namespace ablate::particles::processes::burningModel {


class waxFuel : public ablate::particles::processes::burningModel::LiquidFuel{

   public:
    // Constructor with initializer list
    waxFuel() : fuelProperties{ 0.33E6, 41.1E6,450, 800, 0.12,3330},  // Fuel properties (Hvap, Hc, MW, rhol, kl, Cp)
          AntoineConstants{7.1356, 2276.9, 75.9}       // Initialize Antoine constants (A, B, C)
    {}
    // Fuel properties are consistent with simit
    // Cp evaluated at 1000K from simit


    public:

     const std::string fuelname = "wax";

    const ablate::particles::processes::burningModel::LiquidFuel::fuelPropeties fuelProperties;

    const ablate::particles::processes::burningModel::LiquidFuel::AntoineConstants AntoineConstants;

    void Tvap(double* Tvap, double* Pvap) {

        double P_hgmm=*Pvap/133.322;			// converting Pin (Pa) into Hgmm
        double Tresult =(AntoineConstants.B/(AntoineConstants.A-log10(P_hgmm)))-AntoineConstants.C;    	//input pressures in Hgmm and returns Temperature in C
        *Tvap=Tresult+273.15;		//converting the result (Celsius) into Kelvin //const.Tref is 298.15 which is room temperature

    }
    /******************************************************************************/
   void Pvap(double* Tvap,double* Pvap) {

        double T_celsius=*Tvap-273.15;		//converting Tin (K) into C
        double Presult = (pow(10,AntoineConstants.A - AntoineConstants.B/(T_celsius+AntoineConstants.C)));					//input Temperature in C and returns Pressure in Hgmm
        *Pvap=Presult*133.322;		//converting the result(Hgmm) into Pa

    }
};
}
#endif  // ABLATELIBRARY_WAXFUEL_H
