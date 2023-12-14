#ifndef ABLATELIBRARY_SURFACERADIATION_HPP
#define ABLATELIBRARY_SURFACERADIATION_HPP

#include "domain/range.hpp"
#include "domain/reverseRange.hpp"
#include "radiation.hpp"
#include "utilities/constants.hpp"

namespace ablate::radiation {

class SurfaceRadiation : public ablate::radiation::Radiation {
   protected:
    //! used to look up from the face id to range index
    ablate::domain::ReverseRange indexLookup;

   public:
    SurfaceRadiation(const std::string& solverId, const std::shared_ptr<domain::Region>& region, const PetscInt raynumber, std::shared_ptr<eos::radiationProperties::RadiationModel> radiationModelIn,
                     std::shared_ptr<ablate::monitors::logs::Log> = {});
    ~SurfaceRadiation();

    void Initialize(const ablate::domain::Range& cellRange, ablate::domain::SubDomain& subDomain) override;
    /**
     * Computes the normal component for this ray
     * @param normal
     * @param iCell
     * @param nphi
     * @param ntheta
     * @return
     */
    PetscReal SurfaceComponent(const PetscReal normal[], PetscInt iCell, PetscInt nphi, PetscInt ntheta) override;

    /**
     * Represents the name of the class for logging and other utilities
     * @return
     */
    static inline std::string GetClassType() { return "SurfaceRadiation"; }

    /**
     * Compute total intensity (pre computed gains + current loss) with
     * @param faceId the current face id
     * @param temperature the temperature of the face
     * @param emissivity the emissivity of the surface
     * @return
     */
    virtual inline void GetSurfaceIntensity(PetscReal* intensity, PetscInt faceId, PetscReal temperature, PetscReal emissivity = 1.0, PetscReal absorptivity = 1.0) {
        // add in precomputed gains
        for (int i = 0; i < (int)absorptivityFunction.propertySize; ++i) {  // Compute the losses
            // Emitted From Surface = \epsilon \sigma T^4
            PetscReal netIntensity = -emissivity*ablate::utilities::Constants::sbc * temperature * temperature * temperature * temperature;

            //Absorbed at the surface = q''_Rad * \alpha
            netIntensity += absorptivity*evaluatedGains[absorptivityFunction.propertySize * indexLookup.GetAbsoluteIndex(faceId) + i];

            // This line just ensures that the netIntenisty isn't going to extreme values, If this clipping is really needed, it probably should be fixed somewhere else -klb
            intensity[i] = abs(netIntensity) > ablate::utilities::Constants::large ? ablate::utilities::Constants::large * PetscSignReal(netIntensity) : netIntensity;
        }
    }
};
}  // namespace ablate::radiation
#endif  // ABLATELIBRARY_SURFACERADIATION_HPP
