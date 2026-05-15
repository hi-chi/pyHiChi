cmake_minimum_required(VERSION 3.1)

add_subdirectory(Core)
add_subdirectory(Modules/Particle)
add_subdirectory(Modules/Field)
add_subdirectory(pyHiChi)

if (USE_TESTS)
   add_subdirectory(UnitTests)
   add_subdirectory(MPI_Tests)
endif()

if (USE_PTESTS)
   add_subdirectory(PerformanceTests)
endif()

add_subdirectory(MPI)