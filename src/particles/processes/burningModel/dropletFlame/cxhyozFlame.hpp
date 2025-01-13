#ifndef ABLATELIBRARY_CXHYOZFLAME_H
#define ABLATELIBRARY_CXHYOZFLAME_H
#include <string>

/**
* for a simple hydrocarbon reaction with air of the form:
*
* CxHyOz +  (x+y/4-z/2) (O2 + 3.76 N2)  -->  y/2 H2O + x CO2 + (x+y/4-z/2) 3.76 N2
 */
namespace ablate::particles::processes::burningModel {

class CxHyOzFlame {
   public:
    CxHyOzFlame(){};

    double x;
    double y;
    double z;

    const double MWO2 = 32;
    const double MWN2 = 28;
    const double MWH2O = 18;
    const double MWCO2 = 44;
    const double MWFuel = 450; // TODO Decide if it's C31 or C32 and correct it...

    double nuO2;
    double nuC32H66;
    double nuCO2;
    double nuH2O;
    double N2_D_O2;

    // Adiabatic Flame temperature
    double Tad;

    // Heat of combustion
    double Hc;

    //
    double Mdoti_D_Mdot2_O2;
    double Mdoti_D_Mdot2_CO2;
    double Mdoti_D_Mdot2_H2O;

    // Mass fractions in the fuel stream (before ignition)
    double YiZ1_PMO2;
    double YiZ1_PMN2;
    double YiZ1_PMFuel;

    // "Fuel" stream (Z=1) after reactions in limit as time --> infinity.
    double YiZ1_infFuel;
    double YiZ1_infCO2;
    double YiZ1_infH2O;
    double YiZ1_infN2;

    // Products as time --> infinity (Zst)
    double YiZst_infCO2;
    double YiZst_infH2O;
    double YiZst_infN2;

    const std::string fuelname; // Name of the fuel

    virtual ~CxHyOzFlame() = default;

//    virtual void FlameUpdate(double* YiO2, double* Hcin, double* Tadin) {
//        // Implementation of the FlameUpdate method
//        // For demonstration purposes, let's assume we perform some calculations here
//        // You may need to adjust this implementation based on your specific requirements
//
//        // Example calculation (this is just a placeholder):
//        *YiO2 = nuO2;  // Update YiO2 with the nuO2 value
////        *Hc = this->Hc;  // Update Hc with the heat of combustion
////        *Tad = this->Tad;  // Update Tad with the adiabatic flame temperature
//
//        // Add any additional calculations or logic as needed
//    }
};

}  // namespace ablate::particles::processes::burningModel

#endif  // ABLATELIBRARY_CXHYOZFLAME_H