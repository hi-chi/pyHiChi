cmake_minimum_required(VERSION 4.2)

project(tests)

set(CMAKE_CXX_STANDARD 17)

set(TESTS_INCLUDE_DIRS
	include
	../../3rdparty/google/googletest/googletest/include
	../Core/include
	../Modules/Particle/include
	../Modules/Field/include
)
set(TESTS_HEADER_DIR include)
set(TESTS_SOURCE_DIR src)

file(GLOB_RECURSE tests_includes CONFIGURE_DEPENDS "${TESTS_HEADER_DIR}/*.h*")
file(GLOB_RECURSE tests_sources CONFIGURE_DEPENDS "${TESTS_SOURCE_DIR}/*.c*")

add_executable(tests ${tests_sources} ${tests_includes})

target_include_directories(tests PUBLIC ${TESTS_INCLUDE_DIRS})
target_link_libraries(tests PUBLIC gtest ${FFT_TARGET})
if (USE_MPI)
	target_link_libraries(tests PUBLIC MPI::MPI_CXX)
endif()
if (USE_OMP)
    target_link_libraries(tests PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(gtest PROPERTIES FOLDER UnitTests)
set_target_properties(tests PROPERTIES FOLDER UnitTests)
