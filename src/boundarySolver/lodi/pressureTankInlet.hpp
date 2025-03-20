#ifndef ABLATELIBRARY_PRESSURETANKINLET_HPP
#define ABLATELIBRARY_PRESSURETANKINLET_HPP

#include "lodiBoundary.hpp"
namespace ablate::boundarySolver::lodi {

class PressureTankInlet : public LODIBoundary {
   public:
    explicit PressureTankInlet(std::shared_ptr<eos::EOS> eos, PetscReal tankPressure, PetscReal tankTemperature, PetscReal valveArea,
                               PetscReal Rgas, PetscReal gamma, std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling = {});

    void Setup(ablate::boundarySolver::BoundarySolver& bSolver) override;

    static PetscErrorCode PressureTankInletFunction(PetscInt dim, const boundarySolver::BoundarySolver::BoundaryFVFaceGeom* fg, const PetscFVCellGeom* boundaryCell, const PetscInt uOff[],
                                                 const PetscScalar* boundaryValues, const PetscScalar* stencilValues[], const PetscInt aOff[], const PetscScalar* auxValues,
                                                 const PetscScalar* stencilAuxValues[], PetscInt stencilSize, const PetscInt stencil[], const PetscScalar stencilWeights[], const PetscInt sOff[],
                                                 PetscScalar source[], void* ctx);

   private:
    static PetscErrorCode UpdateConservedVariables(PetscInt dim, const BoundarySolver::BoundaryFVFaceGeom* fg, const PetscFVCellGeom* boundaryCell, const PetscInt uOff[], PetscScalar* boundaryValues,
                                        const PetscScalar* stencilValues, const PetscInt aOff[], PetscScalar* auxValues, const PetscScalar* stencilAuxValues, void* ctx);
    PetscReal tankPressure;
    PetscReal tankTemperature;
    PetscReal valveArea;
    PetscReal gammaInlet;
    PetscReal RInlet;
    PetscReal tankDensity;
};

}  // namespace ablate::boundarySolver::lodi
#endif  // ABLATELIBRARY_PRESSURETANKINLET_HPP
