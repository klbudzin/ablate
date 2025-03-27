target_sources(ablateLibrary
        PRIVATE
        lodiBoundary.cpp
        isothermalWall.cpp
        openBoundary.cpp
        inlet.cpp
        stagnationEnthalpy.cpp
        pressureTankInlet.cpp
        massFluxInlet.cpp

        PUBLIC
        lodiBoundary.hpp
        isothermalWall.hpp
        openBoundary.hpp
        inlet.hpp
        stagnationEnthalpy.hpp
        pressureTankInlet.hpp
        massFluxInlet.hpp
        )

