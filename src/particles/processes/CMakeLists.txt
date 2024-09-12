target_sources(ablateLibrary
        PRIVATE
        tracer.cpp
        inertial.cpp
        arbitraryEulerianSource.cpp
        convectiveHeatTransfer.cpp
        coupledDrag.cpp

        PUBLIC
        process.hpp
        tracer.hpp
        inertial.hpp
        coupledProcess.hpp
        arbitraryEulerianSource.hpp
        convectiveHeatTransfer.hpp
        coupledDrag.hpp
        )

add_subdirectory(drag)