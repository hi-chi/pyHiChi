cmake_minimum_required(VERSION 4.2)

project(core)

set(CORE_INCLUDE_DIRS include ${FFT_INCLUDES})
set(CORE_HEADER_DIR include)

file(GLOB_RECURSE core_headers CONFIGURE_DEPENDS "${CORE_HEADER_DIR}/*.h*")

add_library(core INTERFACE ${core_headers})

target_include_directories(core INTERFACE ${CORE_INCLUDE_DIRS})
if (USE_MPI)
	target_link_libraries(core INTERFACE MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(core INTERFACE OpenMP::OpenMP_CXX)
endif()
