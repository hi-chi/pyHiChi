cmake_minimum_required(VERSION 3.1)

add_subdirectory(Core)
add_subdirectory(Modules/Particle)
add_subdirectory(Modules/Field)

if (USE_PYTHON_INTERFACE)
   add_subdirectory(pyHiChi)
endif()

if (USE_TESTS)
   add_subdirectory(UnitTests)
endif()

if (USE_PTESTS)
   add_subdirectory(PerformanceTests)
endif()