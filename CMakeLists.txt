cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

option(USE_MKL OFF)
option(USE_FFTW OFF)
option(USE_OMP ON)

project(hiChi)

if (USE_OMP)
    add_definitions(-D__USE_OMP__)

    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

include(cmake/functions.cmake)
link_fft_libs()

if (USE_MKL OR USE_FFTW)
    add_definitions(-D__USE_FFT__)
endif()


add_subdirectory(3rdparty/pybind11)

if (USE_PTESTS)
   add_subdirectory(3rdparty/google/benchmark)
endif()

add_subdirectory(src)