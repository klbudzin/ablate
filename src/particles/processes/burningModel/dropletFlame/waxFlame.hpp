#ifndef ABLATELIBRARY_WAXFLAME_H
#define ABLATELIBRARY_WAXFLAME_H

#include "cxhyozFlame.hpp"
#include <cmath>

namespace ablate::particles::processes::burningModel {

class waxFlame : public ablate::particles::processes::burningModel::CxHyOzFlame {
   public:
    waxFlame() {
        x = 32;
        y = 66;
        z = 0;

        double YO2_inf = 0.23;

        nuO2 = (x + y / 4.0 - z / 2.0) * MWO2 / MWFuel;
        nuCO2 = x * MWCO2 / MWFuel;
        nuH2O = 0.5 * y * MWH2O / MWFuel;

        Tad = 2700;
        Hc = 41.1E6;
    }

    const std::string fuelname = "wax";
};

} // namespace ablate::particles::processes::burningModel

#endif // ABLATELIBRARY_WAXFLAME_H