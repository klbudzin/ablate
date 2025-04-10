#include "burningParticleSolver.hpp"
#include "utilities/vectorUtilities.hpp"
#include "accessors/mutableSwarmAccessor.hpp"

    //The whole point of defining this specific class serperate from the coupled particle solver is to ensure that certain fields wanted for all burning particles are already defined
    //Thus the constructor will just call the coupled particle solver initializer with preset Fields already defined
    ablate::particles::BurningParticleSolver::BurningParticleSolver(std::string solverId, std::shared_ptr<domain::Region> region, std::shared_ptr<parameters::Parameters> options, std::vector<FieldDescription> additionalFields,
                          std::vector<std::shared_ptr<processes::Process>> processes, std::shared_ptr<initializers::Initializer> initializer,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, PetscReal minimumDiameterIn, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions,
                          const std::vector<std::string>& coupledFields) : CoupledParticleSolver(std::move(solverId), std::move(region), std::move(options),
        //Create the Default Fields for a Burning Particle, Coordinates are already defulated in for the particle Solver
        std::vector<FieldDescription>{   ablate::particles::FieldDescription( ParticleMass, domain::FieldLocation::SOL ),
         ablate::particles::FieldDescription( ParticleTemperature, domain::FieldLocation::SOL ),
         ablate::particles::FieldDescription( ParticleVelocity, domain::FieldLocation::SOL, std::vector<std::string>{"vel" + domain::FieldDescription::DIMENSION} ),
         ablate::particles::FieldDescription( ParticleDiameter, domain::FieldLocation::AUX ),
         ablate::particles::FieldDescription( ParticleDensity, domain::FieldLocation::AUX ),
         ablate::particles::FieldDescription( ParticleNPP, domain::FieldLocation::AUX),
         ablate::particles::FieldDescription( ParticleCP, domain::FieldLocation::AUX),
         ablate::particles::FieldDescription( ParticleBurning, domain::FieldLocation::AUX)
        }, std::move(processes), std::move(initializer), std::move(fieldInitialization), std::move(exactSolutions),
        std::move(coupledFields)), minimumDiameter(minimumDiameterIn ? minimumDiameterIn : 1e-5)
    {
        // Make sure there is a burning model and filter it out
        {
            std::vector<std::shared_ptr<processes::BurningProcess>> temporary;
            temporary = ablate::utilities::VectorUtilities::Filter<processes::BurningProcess>(coupledProcesses);
            if (!temporary.empty())  burningModel = temporary.at(0);
            if (temporary.size() > 1)
                throw std::invalid_argument("The Burning Particle solver expects only 1 burning process in the processes");
        }
    }

    ablate::particles::BurningParticleSolver::BurningParticleSolver(std::string solverId, std::shared_ptr<domain::Region> region, std::shared_ptr<parameters::Parameters> options, const std::vector<std::shared_ptr<FieldDescription>>& additionalFields,
                          std::vector<std::shared_ptr<processes::Process>> processes, std::shared_ptr<initializers::Initializer> initializer,
                          std::vector<std::shared_ptr<mathFunctions::FieldFunction>> fieldInitialization, PetscReal minimumDiameterIn, std::vector<std::shared_ptr<mathFunctions::FieldFunction>> exactSolutions,
                          const std::vector<std::string>& coupledFields)
                          : BurningParticleSolver(std::move(solverId), std::move(region), std::move(options), ablate::utilities::VectorUtilities::Copy(additionalFields),
                                                  std::move(processes), std::move(initializer), std::move(fieldInitialization), minimumDiameterIn, std::move(exactSolutions), coupledFields)
                          {}

    void ablate::particles::BurningParticleSolver::CheckForRemovedParticles() {
        /** Get the swarm cell DM, cell Id and coordinate fields */
//        DMSwarmCellDM cellDm;
//        const char *cellid;
//        PetscInt* swarm_index;
//        DMSwarmGetCellDMActive(swarmDm, &cellDm) >> utilities::PetscUtilities::checkError;
//        DMSwarmCellDMGetCellID(cellDm, &cellid) >> utilities::PetscUtilities::checkError;
//        DMSwarmGetField(swarmDm, cellid, nullptr, nullptr, (void**)&swarm_index) >> utilities::PetscUtilities::checkError;

        // determine if we should cachePointData
        auto cachePointData = processes.size() != 1;
        Vec packedSolutionVec;
        DMSwarmCreateGlobalVectorFromField(swarmDm, PackedSolution, &packedSolutionVec) >> utilities::PetscUtilities::checkError;
        {
             accessors::SwarmAccessor swarmAccessor(cachePointData, swarmDm, fieldsMap, packedSolutionVec);
             auto diameter = swarmAccessor[ParticleDiameter];
             const auto np = swarmAccessor.GetNumberParticles();
              for(PetscInt p = 0; p < np; p++) {
//                  if(!subDomain->InRegion(swarm_index[np]) || diameter(p) <= 0){
                    if(diameter(p) <= 0){
                      this->RemovingParticleIndices.push_back(p);
                  }
              }
        }
        DMSwarmDestroyGlobalVectorFromField(swarmDm, PackedSolution, &packedSolutionVec) >> utilities::PetscUtilities::checkError;
//        DMSwarmRestoreField(swarmDm, cellid, nullptr, nullptr, (void**)&swarm_index) >> utilities::PetscUtilities::checkError;
    }

    void ablate::particles::BurningParticleSolver::DecodeSolverAuxVariables(double dt) {
        // determine if we should cachePointData
        auto cachePointData = processes.size() != 1;
        Vec packedSolutionVec;

        // extract the vectors again
        DMSwarmCreateLocalVectorFromField(swarmDm, PackedSolution, &packedSolutionVec) >> utilities::PetscUtilities::checkError;
        {
            // use a mutable swarm accessor to access and change the right fields
            accessors::MutableSwarmAccessor swarmAccessor(cachePointData, swarmDm, fieldsMap, packedSolutionVec);
            // We need the eulerian field to update the burning state of the particle
            accessors::EulerianAccessor eulerianAccessor(cachePointData, subDomain, swarmAccessor, timeFinal);

            //Allow the Burning Model to handel updating whether a particle is burning or not
            if(burningModel) {
                // can use this updateAuxFields call to do analytical time stepped solutions, should be able to have the macrostepper
                // pass in startTime end Time here as well, can even have an analytical toggle
                this->burningModel->UpdateAuxFields();
                this->burningModel->UpdateParticleBurning(eulerianAccessor, swarmAccessor);
            }

            //Now check the minimum diameter and mark particles to be removed that have burned fully
            const auto np = swarmAccessor.GetNumberParticles();
            auto parcelMass = swarmAccessor[ParticleMass];
            auto Npp = swarmAccessor[ParticleNPP];
            auto density = swarmAccessor[ParticleDensity];
            auto diameter = swarmAccessor[ParticleDiameter];

            for(PetscInt p = 0; p < np; p++) {
                //If the Mass went negative, force it and the diameter to 0
                parcelMass(p) = PetscMax(0,parcelMass(p));
                diameter(p) = std::pow(6/PETSC_PI*parcelMass(p)/Npp(p)/density(p),1./3.);
                //Check if the diameter has gone lower than the minimum diameter, If it has Send the particle to 0 Mass
                if ( CheckMinimumDiameterLimit(diameter(p)) ) {
                    parcelMass(p) = 0;
                    diameter(p) = 0;
                }
            }
        }
        DMSwarmDestroyLocalVectorFromField(swarmDm, PackedSolution, &packedSolutionVec) >> utilities::PetscUtilities::checkError;
    }








#include "registrar.hpp"
REGISTER(ablate::solver::Solver, ablate::particles::BurningParticleSolver, "Burning lagrangian particle solver", ARG(std::string, "id", "the name of the particle solver"),
         OPT(ablate::domain::Region, "region", "the region to apply this solver.  Default is entire domain"),
         OPT(ablate::parameters::Parameters, "options", "the options passed to PETSC for the flow"),
         OPT(std::vector<ablate::particles::FieldDescription>, "fields", "any additional fields beside Mass, Velocity, Temperature, Diameter, Density, SpecificHeat, and Num particles per parcel"),
         ARG(std::vector<ablate::particles::processes::Process>, "processes", "the processes used to describe the particle source terms"),
         ARG(ablate::particles::initializers::Initializer, "initializer", "the initial particle setup methods"),
         OPT(std::vector<ablate::mathFunctions::FieldFunction>, "fieldInitialization", "the initial particle fields values"),
         OPT(PetscReal, "minimumDiameter", "minimum Diameter a particle can be before it is marked for removal"),
         OPT(std::vector<ablate::mathFunctions::FieldFunction>, "exactSolutions", "particle fields (SOL) exact solutions"));