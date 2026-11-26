cmake_minimum_required(VERSION 4.2)

project(particleModules)

set(PARTICLEMODULES_INCLUDE_DIRS include ../../Core/include)
set(PARTICLEMODULES_HEADER_DIR include)

file(GLOB_RECURSE particleModules_headers CONFIGURE_DEPENDS "${PARTICLEMODULES_HEADER_DIR}/*.h*")

add_library(particleModules INTERFACE ${particleModules_headers})

target_include_directories(particleModules INTERFACE ${CORE_INCLUDE_DIRS})
if (USE_MPI)
	target_link_libraries(particleModules INTERFACE MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(particleModules INTERFACE OpenMP::OpenMP_CXX)
endif()
