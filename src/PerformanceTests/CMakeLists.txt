cmake_minimum_required(VERSION 4.2)

project(ptests)

set(CMAKE_CXX_STANDARD 17)

set(PTESTS_INCLUDE_DIRS
	include
	../../3rdparty/google/googletest/googletest/include
	../Core/include
	../Modules/Particle/include
	../Modules/Field/include
)
set(PTESTS_HEADER_DIR include)
set(PTESTS_SOURCE_DIR src)

file(GLOB_RECURSE ptests_includes CONFIGURE_DEPENDS "${PTESTS_HEADER_DIR}/*.h*")
file(GLOB_RECURSE ptests_sources CONFIGURE_DEPENDS "${PTESTS_SOURCE_DIR}/*.c*")

add_executable(ptests ${ptests_sources} ${ptests_includes})

target_include_directories(ptests PUBLIC ${PTESTS_INCLUDE_DIRS})
target_link_libraries(ptests PUBLIC benchmark ${FFT_TARGET})
if (USE_MPI)
	target_link_libraries(ptests PUBLIC MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(ptests PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(benchmark PROPERTIES FOLDER UnitTests)
set_target_properties(ptests PROPERTIES FOLDER UnitTests)
