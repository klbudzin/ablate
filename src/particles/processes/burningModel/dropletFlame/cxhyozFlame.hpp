#ifndef ABLATELIBRARY_CXHYOZFLAME_H
#define ABLATELIBRARY_CXHYOZFLAME_H
#include <string>


/**
* for a simple hydrocarbon reaction with air of the form:
*
* CxHyOz +  (x+y/4-z/2) (02 + 3.76 N2)  -->  y/2 H2O + x CO2 + (x+y/4-z/2) 3.76 N2
*/
namespace ablate::particles::processes::burningModel {

class CxHyOzFlame {
   public:
    CxHyOzFlame(){};

    double x;
    double y;
    double z;

    const double MWO2=32;
    const double MWN2=28;
    const double MWH2O=18;
    const double MWCO2=44;
    const double MWFuel=450; //TODO Decide if its C31 or C32 and correct it...

    double nuO2;
    double nuC32H66;
    double nuCO2;
    double nuH2O;
    double N2_D_O2;

    //Adiabatic Flame temperature
    double Tad;

    //Heat of combustion
    double Hc;

    //
    double Mdoti_D_Mdot2_O2;
    double Mdoti_D_Mdot2_CO2;
    double Mdoti_D_Mdot2_H2O;

    //mass fractions in the fuel stream (before ignition)
    double YiZ1_PMO2;
    double YiZ1_PMN2;
    double YiZ1_PMFuel;

    // "fuel" stream (Z=1.) after reactions in limit as time --> infinity.
    double YiZ1_infFuel;
    double YiZ1_infCO2;
    double YiZ1_infH2O;
    double YiZ1_infN2;


    // products as time --> infinity (Zst)
    double YiZst_infCO2;
    double YiZst_infH2O;
    double YiZst_infN2;


    const std::string fuelname;                       // Name of the fuel

    virtual ~CxHyOzFlame() = default;


    virtual void FlameUpdate(double* YiO2, double* Hc, double* Tad);


};

}
#endif  // ABLATELIBRARY_CXHYOZFLAME_H
