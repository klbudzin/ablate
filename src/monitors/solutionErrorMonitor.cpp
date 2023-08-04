#include "solutionErrorMonitor.hpp"
#include <monitors/logs/stdOut.hpp>
#include "mathFunctions/mathFunction.hpp"
#include "utilities/petscUtilities.hpp"

ablate::monitors::SolutionErrorMonitor::SolutionErrorMonitor(ablate::monitors::SolutionErrorMonitor::Scope errorScope, ablate::utilities::MathUtilities::Norm normType,
                                                             const std::shared_ptr<logs::Log>& logIn)
    : errorScope(errorScope), normType(normType), log(logIn ? logIn : std::make_shared<logs::StdOut>()) {}

PetscErrorCode ablate::monitors::SolutionErrorMonitor::MonitorError(TS ts, PetscInt step, PetscReal crtime, Vec u, void* ctx) {
    PetscFunctionBeginUser;

    DM dm;
    PetscDS ds;
    PetscCall(TSGetDM(ts, &dm));
    PetscCall(DMGetDS(dm, &ds));

    // Check for the number of DS, this should be relaxed
    PetscInt numberDS;
    PetscCall(DMGetNumDS(dm, &numberDS));
    if (numberDS > 1) {
        SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "This monitor only supports a single DS in a DM");
    }

    // Get the number of fields
    PetscInt numberOfFields;
    PetscCall(PetscDSGetNumFields(ds, &numberOfFields));
    PetscInt* numberComponentsPerField;
    PetscCall(PetscDSGetComponents(ds, &numberComponentsPerField));

    auto* errorMonitor = (SolutionErrorMonitor*)ctx;

    // if this is the first time step init the log
    if (!errorMonitor->log->Initialized()) {
        errorMonitor->log->Initialize(PetscObjectComm((PetscObject)dm));
    }

    std::vector<PetscReal> ferrors;
    try {
        ferrors = errorMonitor->ComputeError(ts, crtime, u);
    } catch (std::exception& exception) {
        SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_LIB, "%s", exception.what());
    }

    // get the error type
    std::stringstream errorTypeStream;
    errorTypeStream << errorMonitor->normType;
    auto errorTypeName = errorTypeStream.str();
    // Change the output depending upon type
    switch (errorMonitor->errorScope) {
        case Scope::VECTOR:
            errorMonitor->log->Printf("Timestep: %04d time = %-8.4g \t %s", (int)step, (double)crtime, errorTypeName.c_str());
            errorMonitor->log->Print("error", ferrors, "%2.3g");
            errorMonitor->log->Print("\n");
            break;
        case Scope::COMPONENT: {
            errorMonitor->log->Printf("Timestep: %04d time = %-8.4g \t %s error:\n", (int)step, (double)crtime, errorTypeName.c_str());
            PetscInt fieldOffset = 0;
            for (PetscInt f = 0; f < numberOfFields; f++) {
                PetscObject field;
                PetscCall(DMGetField(dm, f, NULL, &field));
                const char* name;
                PetscObjectGetName((PetscObject)field, &name);

                errorMonitor->log->Print("\t ");
                errorMonitor->log->Print(name, numberComponentsPerField[f], &ferrors[fieldOffset], "%2.3g");
                errorMonitor->log->Print("\n");
                fieldOffset += numberComponentsPerField[f];
            }
        } break;
        default: {
            SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_LIB, "Unknown error scope");
        }
    }

    PetscFunctionReturn(0);
}

std::vector<PetscReal> ablate::monitors::SolutionErrorMonitor::ComputeError(TS ts, PetscReal time, Vec u) {
    DM dm;
    PetscDS ds;
    TSGetDM(ts, &dm) >> utilities::PetscUtilities::checkError;
    DMGetDS(dm, &ds) >> utilities::PetscUtilities::checkError;

    // Get the number of fields
    PetscInt numberOfFields;
    PetscDSGetNumFields(ds, &numberOfFields) >> utilities::PetscUtilities::checkError;
    PetscInt* numberComponentsPerField;
    PetscDSGetComponents(ds, &numberComponentsPerField) >> utilities::PetscUtilities::checkError;

    // compute the total number of components
    PetscInt totalComponents = 0;

    // Get the exact funcs and context
    std::vector<ablate::mathFunctions::PetscFunction> exactFuncs(numberOfFields);
    std::vector<void*> exactCtxs(numberOfFields);
    for (auto f = 0; f < numberOfFields; ++f) {
        PetscDSGetExactSolution(ds, f, &exactFuncs[f], &exactCtxs[f]) >> utilities::PetscUtilities::checkError;
        if (!exactFuncs[f]) {
            throw std::invalid_argument("The exact solution has not set");
        }
        totalComponents += numberComponentsPerField[f];
    }

    // Create an vector to hold the exact solution
    Vec exactVec;
    VecDuplicate(u, &exactVec) >> utilities::PetscUtilities::checkError;
    DMProjectFunction(dm, time, &exactFuncs[0], &exactCtxs[0], INSERT_ALL_VALUES, exactVec) >> utilities::PetscUtilities::checkError;

    // If we treat this as a single vector or multiple components change how this is done
    totalComponents = errorScope == Scope::VECTOR ? 1 : totalComponents;

    // Update the block size
    VecSetBlockSize(exactVec, totalComponents) >> utilities::PetscUtilities::checkError;

    // Compute the l2 errors
    std::vector<PetscReal> ferrors(totalComponents);

    // compute the norm
    ablate::utilities::MathUtilities::ComputeNorm(normType, u, exactVec, ferrors.data()) >> utilities::PetscUtilities::checkError;

    VecDestroy(&exactVec) >> utilities::PetscUtilities::checkError;
    return ferrors;
}

std::ostream& ablate::monitors::operator<<(std::ostream& os, const ablate::monitors::SolutionErrorMonitor::Scope& v) {
    switch (v) {
        case SolutionErrorMonitor::Scope::VECTOR:
            return os << "vector";
        case SolutionErrorMonitor::Scope::COMPONENT:
            return os << "component";
        default:
            return os;
    }
}

std::istream& ablate::monitors::operator>>(std::istream& is, ablate::monitors::SolutionErrorMonitor::Scope& v) {
    std::string enumString;
    is >> enumString;

    if (enumString == "vector") {
        v = SolutionErrorMonitor::Scope::VECTOR;
    } else if (enumString == "component") {
        v = SolutionErrorMonitor::Scope::COMPONENT;
    } else {
        throw std::invalid_argument("Unknown Scope type " + enumString);
    }
    return is;
}

#include "registrar.hpp"
REGISTER(ablate::monitors::Monitor, ablate::monitors::SolutionErrorMonitor, "Computes and reports the error every time step",
         ENUM(ablate::monitors::SolutionErrorMonitor::Scope, "scope", "how the error should be calculated ('vector', 'component')"),
         ENUM(ablate::utilities::MathUtilities::Norm, "type", "norm type ('l1','l1_norm','l2', 'linf', 'l2_norm')"),
         OPT(ablate::monitors::logs::Log, "log", "where to record log (default is stdout)"));
