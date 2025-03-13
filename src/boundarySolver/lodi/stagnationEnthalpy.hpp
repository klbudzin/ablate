#ifndef ABLATELIBRARY_STAGNATIONENTHALPY_HPP
#define ABLATELIBRARY_STAGNATIONENTHALPY_HPP

#include "lodiBoundary.hpp"
namespace ablate::boundarySolver::lodi {

class StagnationEnthalpy : public LODIBoundary {
   public:
    explicit StagnationEnthalpy(std::shared_ptr<eos::EOS> eos, PetscReal stagEnthalpy, std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling = {});

    void Setup(ablate::boundarySolver::BoundarySolver& bSolver) override;

    static PetscErrorCode StagnationEnthalpyFunction(PetscInt dim, const boundarySolver::BoundarySolver::BoundaryFVFaceGeom* fg, const PetscFVCellGeom* boundaryCell, const PetscInt uOff[],
                                                 const PetscScalar* boundaryValues, const PetscScalar* stencilValues[], const PetscInt aOff[], const PetscScalar* auxValues,
                                                 const PetscScalar* stencilAuxValues[], PetscInt stencilSize, const PetscInt stencil[], const PetscScalar stencilWeights[], const PetscInt sOff[],
                                                 PetscScalar source[], void* ctx);

   private:
    static PetscErrorCode UpdateEnergy(PetscInt dim, const BoundarySolver::BoundaryFVFaceGeom* fg, const PetscFVCellGeom* boundaryCell, const PetscInt uOff[], PetscScalar* boundaryValues,
                                        const PetscScalar* stencilValues, const PetscInt aOff[], PetscScalar* auxValues, const PetscScalar* stencilAuxValues, void* ctx);
    PetscReal stagnationEnthalpy;
};

}  // namespace ablate::boundarySolver::lodi
#endif  // ABLATELIBRARY_STAGNATIONENTHALPY_HPP
