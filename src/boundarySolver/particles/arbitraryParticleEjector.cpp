#include "arbitraryParticleEjector.hpp"

ablate::boundarySolver::particles::ArbitraryParticleEjector::ArbitraryParticleEjector(std::shared_ptr<ablate::particles::ParticleSolver> particleSolver, std::shared_ptr<ablate::mathFunctions::MathFunction> mathFunc, PetscReal LimitingValue, PetscReal NormalOffset)
    : particleModel(std::move(particleSolver)),
      mathFunction(std::move(mathFunc)),
      limitingMathValue(LimitingValue ? LimitingValue : 0.5),
      offset(NormalOffset ? NormalOffset : 1e-5){
          this->currentTime = 0;
          this->numParticleFields = particleModel->getFields().size();
    }

void ablate::boundarySolver::particles::ArbitraryParticleEjector::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    //Register the compute and eject droplets function to be called from the boundarysolver
    bSolver.RegisterFunction(computeAndEjectDroplets, this, std::vector<std::string>{}, std::vector<std::string>{}, std::vector<std::string>{}, boundarySourceType);
    //Could also register a preRHS function below if we needed to update properties here before computing (in this case we update the current time)
        bSolver.RegisterPreStep([this](auto ts, auto &solver) {
        PetscFunctionBeginUser;
        PetscCall(TSGetTime(ts, &(this->currentTime)));
        PetscFunctionReturn(0);
    });
}




PetscErrorCode ablate::boundarySolver::particles::ArbitraryParticleEjector::computeAndEjectDroplets(PetscInt _dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                         const PetscFVCellGeom *boundaryCell, const PetscInt *uOff, const PetscScalar *boundaryValues,
                                                                                         const PetscScalar **stencilValues, const PetscInt *aOff, const PetscScalar *auxValues,
                                                                                         const PetscScalar **stencilAuxValues, PetscInt stencilSize, const PetscInt *stencil,
                                                                                         const PetscScalar *stencilWeights, const PetscInt *sOff, PetscScalar *source, void *ctx) {
    PetscFunctionBeginUser;
    // Grab the ejector and particleSolver objects from the context
    auto particleEjector = (ablate::boundarySolver::particles::ArbitraryParticleEjector *)ctx;

    //Grab the Face Coords, for below
    auto faceCoords = fg->centroid;
    auto faceNorm =fg->normal;
    if (particleEjector->mathFunction->Eval(faceCoords,_dim, particleEjector->currentTime) > particleEjector->limitingMathValue) {
        //We add a particle in to the ejector by grabbing the new particl vector
        auto &newParticlesVec = particleEjector->particleModel->getNewParticleVector();
        //Create New coords pointer for this face
        PetscReal* particleCoords = new PetscReal[_dim];
        //Offset the particle coordinates so they spawn inside the domain
        for(auto i = 0; i < _dim; i ++)
            particleCoords[i] = -faceNorm[i]*particleEjector->offset + faceCoords[i];
        //Create a pointer for the fields data
        PetscInt numParticles = particleEjector->numParticleFields;
        PetscReal* randomFieldData = new PetscReal[numParticles];
        for(auto i = 0; i < particleEjector->numParticleFields; i++)
            randomFieldData[i] = particleEjector->limitingMathValue;
        newParticlesVec.push_back(particleCoords);
    }
    PetscFunctionReturn(0);
}

#include "registrar.hpp"
REGISTER(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::particles::ArbitraryParticleEjector, "Ejects a particle from the boundary face every set time",
         ARG(ablate::particles::ParticleSolver, "particleSolver", "The particle solver holding the particle swarm dm"),
         ARG(ablate::mathFunctions::MathFunction, "function", "Function to compute every time step to see if we should eject a particle"),
         OPT(PetscReal, "limitingValue", "Value used to determine whether a particle is ejection from the math function"),
         OPT(PetscReal, "offset", "Offset from the boundary the particles are placed")
         );