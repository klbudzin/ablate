target_sources(ablateLibrary
        PRIVATE
        tracer.cpp
        inertial.cpp
        arbitraryEulerianSource.cpp
        convectiveHeatTransfer.cpp

        PUBLIC
        process.hpp
        tracer.hpp
        inertial.hpp
        coupledProcess.hpp
        arbitraryEulerianSource.hpp
        convectiveHeatTransfer.hpp
        )

add_subdirectory(drag)