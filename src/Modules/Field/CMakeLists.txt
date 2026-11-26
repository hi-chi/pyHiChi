cmake_minimum_required(VERSION 4.2)

project(fieldModules)

set(FIELDMODULES_INCLUDE_DIRS include ../../Core/include ${FFT_INCLUDES})
set(FIELDMODULES_HEADER_DIR include)

file(GLOB_RECURSE fieldModules_headers CONFIGURE_DEPENDS "${FIELDMODULES_HEADER_DIR}/*.h*")

add_library(fieldModules INTERFACE ${fieldModules_headers})

target_include_directories(fieldModules INTERFACE ${CORE_INCLUDE_DIRS})
if (USE_MPI)
	target_link_libraries(fieldModules INTERFACE MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(fieldModules INTERFACE OpenMP::OpenMP_CXX)
endif()
