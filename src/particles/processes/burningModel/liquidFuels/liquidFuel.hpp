#ifndef ABLATELIBRARY_LIQUIDFUEL_H
#define ABLATELIBRARY_LIQUIDFUEL_H
#include <string>


/**
 * interface used for liquid fuels required for the twozone model
 */
namespace ablate::particles::processes::burningModel {

class LiquidFuel {
   public:

    virtual ~LiquidFuel() = default;

    struct fuelPropeties {
    double Hvap;
    double Hc;
    double MW;
    double rhol;
    double kl;
    };

    struct AntoineConstants{
        double A;
        double B;
        double C;
    };

    const std::string fuelname;                       // Name of the fuel
//    const fuelPropeties fuelProperties;              // Fuel properties
//    const AntoineConstants AntoineConstants;         // Antoine constants



    virtual void Tvap(double* Tvap, double* Pvap);

    virtual void Pvap(double* Tvap,double* Pvap);
};

}
#endif  // ABLATELIBRARY_LIQUIDFUEL_H
