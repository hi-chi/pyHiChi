cmake_minimum_required(VERSION 3.1)

project(pyHiChi)

include_directories(include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ../Modules/TightFocusing/include
    ${FFT_INCLUDES})

if (USE_OMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(pyHiChi_source
    include/pyGrid.h
    src/pyHiChi.cpp)

pybind11_add_module(pyHiChi ${pyHiChi_source})

if (UNIX)
    target_link_libraries(pyHiChi PUBLIC -Wl,--start-group ${FFT_LIBS} -Wl,--end-group)
else()
    target_link_libraries(pyHiChi PUBLIC ${FFT_LIBS})
endif()
