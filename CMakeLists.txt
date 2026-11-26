cmake_minimum_required(VERSION 4.2)

set(CMAKE_CXX_STANDARD 17)

option(USE_MKL OFF)
option(USE_FFTW OFF)
option(USE_OMP ON)

project(hiChi)

include(cmake/patch_libstdcxx_noexcept.cmake)

add_compile_options(
	"$<$<AND:$<COMPILE_LANG_AND_ID:C,IntelLLVM>,$<PLATFORM_ID:Windows>>:/fp:precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>,$<PLATFORM_ID:Windows>>:/fp:precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:C,IntelLLVM>,$<NOT:$<PLATFORM_ID:Windows>>>:-fp-model=precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>,$<NOT:$<PLATFORM_ID:Windows>>>:-fp-model=precise>"
)

if (USE_MPI)
	add_definitions(-D__USE_MPI__)
	find_package(MPI REQUIRED)
endif()

if (USE_OMP)
	add_definitions(-D__USE_OMP__)
	find_package(OpenMP REQUIRED)
endif()

if (USE_MKL OR USE_FFTW)
	add_definitions(-D__USE_FFT__)
	include(cmake/find_fft.cmake)
	if(USE_MKL)
		find_mkl_fft()
	else()
		find_fftw()
	endif()
endif()


set(PYBIND11_FINDPYTHON NEW)
add_subdirectory(3rdparty/pybind11)

if (USE_TESTS)
	set(GTEST_HAS_PTHREAD, 0)
	add_subdirectory(3rdparty/google/googletest)
endif()

if (USE_PTESTS)
	add_subdirectory(3rdparty/google/benchmark)
endif()

add_subdirectory(src)
