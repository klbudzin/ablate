#ifndef ABLATELIBRARY_CONVECTIVEHEATTRANSFER_HPP
#define ABLATELIBRARY_CONVECTIVEHEATTRANSFER_HPP

#include <string>
#include "coupledProcess.hpp"

namespace ablate::particles::processes {

class ConvectiveHeatTransfer : public CoupledProcess {
   private:
    // convective coefficient in h (T_{eularian}-T_{particle}) *SA_{tot}
    const PetscReal convCoefficient;

   public:
    /**
     * Adds an arbitrary source function for each particle to the eulerian field
     * @param coupledFieldName the name of the eulerian coupled field
     * @param sourceFunction the function to compute the source
     */
    ConvectiveHeatTransfer(PetscReal h);

    /**
     * RHS function for the particle
     */
    void ComputeRHS(PetscReal time, accessors::SwarmAccessor& swarmAccessor, accessors::RhsAccessor& rhsAccessor, accessors::EulerianAccessor& eulerianAccessor) override;

    /**
     * Add the arbitrary source to the eulerianSourceAccessor
     * @param startTime
     * @param endTime
     * @param swarmAccessorPreStep
     * @param swarmAccessorPostStep
     * @param eulerianSourceAccessor
     */
    void ComputeEulerianSource(PetscReal startTime, PetscReal endTime, accessors::SwarmAccessor& swarmAccessorPreStep, accessors::SwarmAccessor& swarmAccessorPostStep,
                               accessors::EulerianSourceAccessor& eulerianSourceAccessor) override;
};

}  // namespace ablate::particles::processes
#endif  // ABLATELIBRARY_CONVECTIVEHEATTRANSFER_HPP
