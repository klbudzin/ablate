#ifndef ABLATELIBRARY_LIQUIDFUEL_H
#define ABLATELIBRARY_LIQUIDFUEL_H

#include <string>

/**
 * Interface used for liquid fuels required for the two-zone model
 */
namespace ablate::particles::processes::burningModel {

class LiquidFuel {
   public:
    virtual ~LiquidFuel() = default;

    struct FuelProperties {
        double Hvap;
        double Hc;
        double MW;
        double rhol;
        double kl;
        double Cp;
    };

    struct AntoineConstants {
        double A;
        double B;
        double C;
    };

    const std::string fuelname;  // Name of the fuel
    const FuelProperties fuelProperties;  // Fuel properties
    const AntoineConstants antoineConstants;  // Antoine constants

    // Constructor for LiquidFuel
    LiquidFuel(const std::string& name, const FuelProperties& properties, const AntoineConstants& constants)
        : fuelname(name), fuelProperties(properties), antoineConstants(constants) {}

    virtual void Tvap(double* Tvap, double* Pvap) = 0;
    virtual void Pvap(double* Tvap, double* Pvap) = 0;
};

}  // namespace ablate::particles::processes::burningModel

#endif  // ABLATELIBRARY_LIQUIDFUEL_H