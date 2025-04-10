#ifndef ABLATELIBRARY_BURNINGPARTICELARBITRARYEJECTOR_HPP
#define ABLATELIBRARY_BURNINGPARTICELARBITRARYEJECTOR_HPP

#include "boundarySolver/boundaryProcess.hpp"
#include "particles/particleSolver.hpp"

//Simple example of a particle Ejector where only one partical per boundary face can ever be ejected per time step
//This class has no use other than showing a tutorial on how a particle ejector can work -klb 2/27/24

namespace ablate::boundarySolver::particles {

class BurningParticleArbitraryEjector : public BoundaryProcess  {
   private:
    inline const static std::string particleEjectorID = "BurningParticle_Arbitrary_Ejector";
     /**
     * Keep the shared pointer to the provided particle Solver
     */
    std::shared_ptr<ablate::particles::ParticleSolver> particleModel = nullptr;

     //! Store the boundary type
    const BoundarySolver::BoundarySourceType boundarySourceType = BoundarySolver::BoundarySourceType::Flux;

    //The math function to be used to compute whether a particle ejected or not
    const std::shared_ptr<ablate::mathFunctions::MathFunction> mathFunction;

    //A limiting value to decide whether to eject a droplet or not (i.e. mathFunction(x) > Limiting Math Value -> Eject droplet)
    const PetscReal limitingMathValue;
    //Offset factor for placing particle in domain
    const PetscReal offset;
    //Stores the current time of the solver
    PetscReal currentTime;

    //Values for the particles coming off (Constant In this case)
    PetscReal* updateValues;
    std::vector<ablate::particles::Field> updateFields;

    public:
        //Constructor
        explicit BurningParticleArbitraryEjector(std::shared_ptr<ablate::particles::ParticleSolver> particleSolver,
                                          std::shared_ptr<ablate::mathFunctions::MathFunction> mathFunc,
                                          PetscReal LimitingValue, PetscReal NormalOffset, PetscReal ParcelMass,
                                          PetscReal PPP, PetscReal Cp, PetscReal T, PetscReal rho);

        //Inherited Methods for a typical Boundary process
        void Setup(ablate::boundarySolver::BoundarySolver &bSolver) override;
       // void Initialize(ablate::boundarySolver::BoundarySolver &bSolver) override;

        // RHS source function that will be hooked into to compute droplet size and eject droplets
        static PetscErrorCode computeAndEjectDroplets(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                             const PetscFVCellGeom *boundaryCell, const PetscInt *uOff, const PetscScalar *boundaryValues,
                                             const PetscScalar **stencilValues, const PetscInt *aOff, const PetscScalar *auxValues,
                                             const PetscScalar **stencilAuxValues, PetscInt stencilSize, const PetscInt *stencil,
                                             const PetscScalar *stencilWeights, const PetscInt *sOff, PetscScalar *source, void *ctx);

};
}
#endif //ABLATELIBRARY_BURNINGPARTICLEARBITRARYEJECTOR_HPP
