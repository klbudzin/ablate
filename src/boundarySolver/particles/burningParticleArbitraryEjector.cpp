#include "burningParticleArbitraryEjector.hpp"

ablate::boundarySolver::particles::BurningParticleArbitraryEjector::BurningParticleArbitraryEjector(
        std::shared_ptr<ablate::particles::ParticleSolver> particleSolver, std::shared_ptr<ablate::mathFunctions::MathFunction> mathFunc,
        PetscReal LimitingValue, PetscReal NormalOffset, PetscReal ParcelMass, PetscReal PPP, PetscReal Cp,
        PetscReal T, PetscReal rho)
    : particleModel(std::move(particleSolver)), mathFunction(std::move(mathFunc)), limitingMathValue(LimitingValue ? LimitingValue : 0.5),
      offset(NormalOffset ? NormalOffset : 1e-5)
      {
          this->currentTime = 0;
          //Set the field values
          this->updateValues = new PetscReal[5]{ParcelMass,PPP,Cp,T,rho};
    }

void ablate::boundarySolver::particles::BurningParticleArbitraryEjector::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    //Create the updateFields Vector
    updateFields.push_back(particleModel->GetField(ablate::particles::ParticleSolver::ParticleMass));
    updateFields.push_back(particleModel->GetField(ablate::particles::ParticleSolver::ParticleNPP));
    updateFields.push_back(particleModel->GetField(ablate::particles::ParticleSolver::ParticleCP));
    updateFields.push_back(particleModel->GetField(ablate::particles::ParticleSolver::ParticleTemperature));
    updateFields.push_back(particleModel->GetField(ablate::particles::ParticleSolver::ParticleDensity));

    //Register the compute and eject droplets function to be called from the boundarysolver
    bSolver.RegisterFunction(computeAndEjectDroplets, this, std::vector<std::string>{}, std::vector<std::string>{}, std::vector<std::string>{}, boundarySourceType);
    //Could also register a preRHS function below if we needed to update properties here before computing (in this case we update the current time)
    bSolver.RegisterPreStep([this](auto ts, auto &solver) {
        PetscFunctionBeginUser;
        PetscCall(TSGetTime(ts, &(this->currentTime)));
        PetscFunctionReturn(0);
    });
}




PetscErrorCode ablate::boundarySolver::particles::BurningParticleArbitraryEjector::computeAndEjectDroplets(PetscInt _dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                         const PetscFVCellGeom *boundaryCell, const PetscInt *uOff, const PetscScalar *boundaryValues,
                                                                                         const PetscScalar **stencilValues, const PetscInt *aOff, const PetscScalar *auxValues,
                                                                                         const PetscScalar **stencilAuxValues, PetscInt stencilSize, const PetscInt *stencil,
                                                                                         const PetscScalar *stencilWeights, const PetscInt *sOff, PetscScalar *source, void *ctx) {
    PetscFunctionBeginUser;
    // Grab the ejector and particleSolver objects from the context
    auto particleEjector = (ablate::boundarySolver::particles::BurningParticleArbitraryEjector *)ctx;
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
        //All the field data should be set to at least initialize a burning particles
        newParticlesVec.push_back(ablate::particles::Particle(particleCoords,5,particleEjector->updateValues,&particleEjector->updateFields));
    }
    PetscFunctionReturn(0);
}

#include "registrar.hpp"
REGISTER(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::particles::BurningParticleArbitraryEjector, "Ejects a particle from the boundary face every set time",
         ARG(ablate::particles::ParticleSolver, "particleSolver", "The particle solver holding the particle swarm dm"),
         ARG(ablate::mathFunctions::MathFunction, "function", "Function to compute every time step to see if we should eject a particle"),
         OPT(PetscReal, "limitingValue", "Value used to determine whether a particle is ejection from the math function"),
         OPT(PetscReal, "offset", "Offset from the boundary the particles are placed"),
         OPT(PetscReal, "parcelMass", "Mass of the parcel"),
         OPT(PetscReal, "particlesPerParcel", "Number of Particles per Parcel"),
         OPT(PetscReal, "Cp", "Droplet specific heat"),
         OPT(PetscReal, "temperature", "Particle Temperature"),
         OPT(PetscReal, "density", "density of the particle")
         );