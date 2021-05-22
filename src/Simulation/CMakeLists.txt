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