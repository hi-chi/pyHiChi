cmake_minimum_required(VERSION 3.1)

project(simulation)

include_directories(include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
)

add_library(simulation_interface INTERFACE)
set(SIMULATION_INCLUDE_DIR include)
set(SIMULATION_HEADER_DIR ${SIMULATION_INCLUDE_DIR})

set(simulation_headers
    ${SIMULATION_HEADER_DIR}/Simulation.h)
add_custom_target(simulation SOURCES ${simulation_headers})

if (USE_MPI)
    find_package(MPI)
    if (MPI_FOUND)
        message(STATUS "MPI FOUND")
        add_definitions(-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK)
        include_directories(simulation_mpi ${MPI_INCLUDE_PATH})
    else()
        message(STATUS "MPI NOT FOUND")
    endif()
endif()

if (USE_OMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()
add_executable(simulation_mpi src/Simulation.cpp)