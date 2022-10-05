target_sources(ablateLibrary
        PRIVATE
        perfectGas.cpp
        stiffenedGas.cpp
        tChem.cpp
        tChemSoot.cpp
        chemTab.cpp
        chemistryModel.cpp
        twoPhase.cpp

        PUBLIC
        eos.hpp
        perfectGas.hpp
        tChem.hpp
        tChemSoot.hpp
        stiffenedGas.hpp
        chemTab.hpp
        chemistryModel.hpp
        twoPhase.hpp
        )

add_subdirectory(transport)
add_subdirectory(tChem)
add_subdirectory(tChemSoot)
add_subdirectory(radiationProperties)
