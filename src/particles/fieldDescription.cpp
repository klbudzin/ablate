#include "fieldDescription.hpp"

#include <utility>
#include "utilities/petscUtilities.hpp"
#include "domain/fieldDescription.hpp"
#include <regex>

ablate::particles::FieldDescription::FieldDescription(std::string name, ablate::domain::FieldLocation type, const std::vector<std::string>& componentsIn, PetscDataType dataTypeIn)
    : name(std::move(name)),
      components(componentsIn.empty() ? std::vector<std::string>{"_"} : componentsIn),
      location(type),
      dataType(dataTypeIn == PETSC_DATATYPE_UNKNOWN ? PETSC_REAL : dataTypeIn) {}

void ablate::particles::FieldDescription::DecompressComponents(PetscInt ndims) {
    //If there is a dimension component, convert it the field description to contain that many dimensions -klb
    for (std::size_t c = 0; c < components.size(); c++) {
        if (components[c].find(ablate::domain::FieldDescription::DIMENSION) != std::string::npos) {
            auto baseName = components[c];

            // Delete this component
            components.erase(components.begin() + c);

            for (PetscInt d = ndims - 1; d >= 0; d--) {
                auto newName = std::regex_replace(baseName, std::regex(ablate::domain::FieldDescription::DIMENSION), std::to_string(d));
                components.insert(components.begin() + c, newName);
            }
        }
    }
}
using namespace ablate::utilities;
#include "registrar.hpp"
REGISTER_DEFAULT(ablate::particles::FieldDescription, ablate::particles::FieldDescription, "Describes a single field for particles", ARG(std::string, "name", "the name of the field"),
                 ARG(EnumWrapper<ablate::domain::FieldLocation>, "location", "if it is a solution (SOL) or auxiliary (aux) field"),
                 OPT(std::vector<std::string>, "components", "the components in this field. (default is a single component)"),
                 OPT(EnumWrapper<PetscDataType>, "dataType", "Possible PETSc data type (default is PETSC_REAL) "));
