cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

option(USE_MKL OFF)
option(USE_FFTW OFF)
option(USE_OMP ON)

project(hiChi)

if (USE_OMP)
    add_definitions(-D__USE_OMP__)
endif()

include(cmake/functions.cmake)
link_fft_libs()

if (USE_MKL OR USE_FFTW)
    add_definitions(-D__USE_FFT__)
endif()


add_subdirectory(3rdparty/pybind11)

add_subdirectory(src)