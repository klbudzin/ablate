target_sources(ablateLibrary
        PRIVATE
        perfectGas.cpp
        stiffenedGas.cpp
        tChem.cpp
        tChemSoot.cpp
        PUBLIC
        eos.hpp
        perfectGas.hpp
        tChem.hpp
        tChemSoot.hpp
        stiffenedGas.hpp
        )

add_subdirectory(transport)
add_subdirectory(tChem)
add_subdirectory(radiationProperties)