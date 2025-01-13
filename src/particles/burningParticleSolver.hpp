#ifndef ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP
#define ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP

#include "coupledParticleSolver.hpp"
#include "processes/coupledProcess.hpp"
#include "solver/rhsFunction.hpp"
#include "particles/processes/burningProcess.hpp"
namespace ablate::particles {

class BurningParticleSolver : public CoupledParticleSolver{
    //minimum diameter a particle can be before it is extinguished
    const PetscReal minimumDiameter;
    std::shared_ptr<processes::BurningProcess> burningModel;

   public:
    /**
     * default constructor
     * @param solverId
     * @param options
     * @param fields
     * @param processes
     * @param initializer
     * @param fieldInitialization
     * @param exactSolutions
     * @param coupledFields the fields to couple to the flow solver.  If not specified all solution fields will be coupled
     */
    BurningParticleSolver(std::string solverId, std::shared_ptr<domain::Region>, std::shared_ptr<parameters::Parameters> options, std::vector<FieldDescription> fields,
                          std::vector<std::shared_ptr<processes::Process>> processes, std::shared_ptr<initializers::Initializer> initializer,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, PetscReal minimumDiameterIn,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions = {}, const std::vector<std::string>& coupledFields = {});

    /**
     * shared pointer version of the constructor
     * @param solverId
     * @param options
     * @param fields
     * @param processes
     * @param initializer
     * @param fieldInitialization
     * @param exactSolutions
     * @param coupledFields the fields to couple to the flow solver.  If not specified all solution fields will be coupled
     */
    BurningParticleSolver(std::string solverId, std::shared_ptr<domain::Region>, std::shared_ptr<parameters::Parameters> options, const std::vector<std::shared_ptr<FieldDescription>>& fields,
                          std::vector<std::shared_ptr<processes::Process>> processes, std::shared_ptr<initializers::Initializer> initializer,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, PetscReal minimumDiameterIn,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions = {}, const std::vector<std::string>& = {});

    /**
    * Check if the Particle Diameter is small enough not to just assume extinguihsed
    */
    inline PetscBool CheckMinimumDiameterLimit(PetscReal particleDiameter) const { return (particleDiameter < minimumDiameter) ? PETSC_TRUE : PETSC_FALSE;}

    void DecodeSolverAuxVariables(double dt) override;

    //! cleanup any petsc objects
    ~BurningParticleSolver() override = default;

};

}  // namespace ablate::particles
#endif  // ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP
