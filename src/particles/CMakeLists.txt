target_sources(ablateLibrary
        PRIVATE
        particleSolver.cpp
        fieldDescription.cpp
        coupledParticleSolver.cpp
        burningParticleSolver.cpp

        PUBLIC
        field.hpp
        fieldDescription.hpp
        particleSolver.hpp
        coupledParticleSolver.hpp
        burningParticleSolver.hpp
        )

add_subdirectory(initializers)
add_subdirectory(processes)
add_subdirectory(accessors)

