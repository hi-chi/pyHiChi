cmake_minimum_required(VERSION 3.1)

add_subdirectory(Core)
add_subdirectory(Simulation)
add_subdirectory(Modules/Particle)
add_subdirectory(Modules/Field)
set(ADIOS2_DIR ${CMAKE_SOURCE_DIR}/../ADIOS2/build)
find_package(ADIOS2 REQUIRED)
add_subdirectory(pyHiChi)

if (USE_TESTS)
   add_subdirectory(UnitTests)
endif()