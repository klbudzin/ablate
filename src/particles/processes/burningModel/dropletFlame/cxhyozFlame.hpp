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

    virtual ~CxHyOzFlame() = default;

    double nuO2;
    double nuC32H66;
    double nuCO2;
    double nuH2O;
    double N2_D_O2;

    //Adiabatic Flame temperature
    double Tad;

    //Heat of combustion
    double Hc;

    const std::string fuelname;                       // Name of the fuel




    virtual void FlameUpdate(double* YiO2, double* Hc, double* Tad);


};

}
#endif  // ABLATELIBRARY_CXHYOZFLAME_H
