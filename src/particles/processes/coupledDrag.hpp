#ifndef ABLATELIBRARY_COUPLEDDRAG_HPP
#define ABLATELIBRARY_COUPLEDDRAG_HPP

#include <string>
#include "coupledProcess.hpp"
#include "eos/eos.hpp"
#include "eos/transport/transportModel.hpp"
namespace ablate::particles::processes {

class CoupledDrag : public CoupledProcess {
   private:
    //Need the function to calculate viscosity from local temperature.
    eos::ThermodynamicTemperatureFunction viscosityFunction;
    const std::vector<PetscReal> gravity;

   public:
    /**
     * Adds an arbitrary source function for each particle to the eulerian field
     * @param coupledFieldName the name of the eulerian coupled field
     * @param sourceFunction the function to compute the source
     */
    CoupledDrag(std::shared_ptr<ablate::eos::transport::TransportModel> transportModel, const std::vector<PetscReal>& gravity);

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
    void ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor& swarmAccessorPreStep, accessors::SwarmAccessor& swarmAccessorPostStep,
                               accessors::EulerianSourceAccessor& eulerianSourceAccessor) override;
};

}  // namespace ablate::particles::processes
#endif  // ABLATELIBRARY_CONVECTIVEHEATTRANSFER_HPP
