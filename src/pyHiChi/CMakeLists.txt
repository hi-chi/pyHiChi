cmake_minimum_required(VERSION 4.2)

project(pyHiChi)

set(PYHICHI_INCLUDE_DIRS
	include
	../Core/include
	../Modules/Particle/include
	../Modules/Field/include
)
set(PYHICHI_HEADER_DIR include)
set(PYHICHI_SOURCE_DIR src)

file(GLOB_RECURSE pyHiChi_includes CONFIGURE_DEPENDS "${PYHICHI_HEADER_DIR}/*.h*")
file(GLOB_RECURSE pyHiChi_sources CONFIGURE_DEPENDS "${PYHICHI_SOURCE_DIR}/*.c*")

set(pyHiChi_source ${pyHiChi_sources} ${pyHiChi_includes})

pybind11_add_module(pyHiChi ${pyHiChi_source})

target_include_directories(pyHiChi PUBLIC ${PYHICHI_INCLUDE_DIRS})
target_link_libraries(pyHiChi PUBLIC ${FFT_TARGET})
if (USE_MPI)
    target_link_libraries(pyHiChi PUBLIC MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(pyHiChi PUBLIC OpenMP::OpenMP_CXX)
endif()
