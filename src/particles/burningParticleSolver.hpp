#ifndef ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP
#define ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP

#include "coupledParticleSolver.hpp"
#include "processes/coupledProcess.hpp"
#include "solver/rhsFunction.hpp"

namespace ablate::particles {

/**
 * This is an extension of the particle solver that allows for fully coupled simulations.
 * The class implements the RHSFunction to allow inserting source terms back to main TS/flowfield
 */
class BurningParticleSolver : public CoupledParticleSolver{

   public:
    inline static const char ParticleCP[] = "ParticleSpecificHeat";
    inline static const char ParticleTemperature[] = "ParticleTemperature";
    inline static const char ParticleNPP[] = "ParticlesPerParcel";
    inline static const char ParticleMass[] = "ParcelMass";
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
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions = {},
                          const std::vector<std::string>& coupledFields = {});

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
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions = {},
                          const std::vector<std::string>& = {});

    void DecodeSolverAuxVariables() override;

    //! cleanup any petsc objects
    ~BurningParticleSolver() override = default;

};

}  // namespace ablate::particles
#endif  // ABLATELIBRARY_BURNINGPARTICLESOLVER_HPP
