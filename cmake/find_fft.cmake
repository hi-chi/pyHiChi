####################### FFTW ############################
function(find_fftw)

    file(TO_CMAKE_PATH ${FFTW_DIR} FFTW_DIR)

    if (NOT EXISTS ${FFTW_DIR})
        message(FATAL_ERROR "Cannot find fftw install directory '${FFTW_DIR}'")
    endif()

    if (EXISTS "${FFTW_DIR}/include")
        set(FFTW_INCLUDE_DIR "${FFTW_DIR}/include")
    else()
        message(FATAL_ERROR "Cannot find 'include' directory in '${FFTW_DIR}' directory")
    endif()

    find_library(
        FFTW3_LIB
        NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}
        PATHS "${FFTW_DIR}/lib" "${FFTW_DIR}/lib64"
    )
    if (NOT FFTW3_LIB)
        message(FATAL_ERROR "Cannot find fftw3 library in '${FFTW_DIR}/lib' or '${FFTW_DIR}/lib64' directory")
    endif()

    set(FFTW_LIBS)
    if (USE_OMP)
        find_library(
            FFTW3_OMP_LIB
            NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX}
            PATHS "${FFTW_DIR}/lib" "${FFTW_DIR}/lib64"
        )
        if (NOT FFTW3_OMP_LIB)
            message(FATAL_ERROR "Cannot find fftw3_omp library in '${FFTW_DIR}/lib' or '${FFTW_DIR}/lib64' directory")
        endif()
        list(APPEND FFTW_LIBS ${FFTW3_OMP_LIB})
    endif()
    list(APPEND FFTW_LIBS ${FFTW3_LIB})

    if (NOT TARGET fftw_target)
        add_library(fftw_target INTERFACE)
        target_include_directories(fftw_target INTERFACE "${FFTW_INCLUDE_DIR}")
        target_link_libraries(fftw_target INTERFACE ${FFTW_LIBS})
    endif()

    set(FFT_TARGET fftw_target PARENT_SCOPE)
    message(STATUS "Using FFTW: ${FFTW_LIBS}")

endfunction()
####################### END FFTW ########################


####################### MKL #############################
function(find_mkl_fft)

    find_package(MKL CONFIG REQUIRED)

    if (NOT TARGET mkl_target)
        add_library(mkl_target INTERFACE)
        target_link_libraries(mkl_target INTERFACE MKL::MKL)

        # add FFTW headers
        set(MKL_FFTW_INCLUDES "")
        if (DEFINED ENV{MKLROOT} AND EXISTS "$ENV{MKLROOT}/include/fftw/fftw3.h")
            set(MKL_FFTW_INCLUDES "$ENV{MKLROOT}/include/fftw")
        endif()

        if (MKL_FFTW_INCLUDES)
            target_include_directories(mkl_target INTERFACE "${MKL_FFTW_INCLUDES}")
        else()
            message(FATAL_ERROR "MKL found, but fftw3.h is missing. MKL's FFTW interface will not be usable")
        endif()
    endif()

    set(FFT_TARGET mkl_target PARENT_SCOPE)
    message(STATUS "Using MKL")

endfunction()
####################### END MKL #########################
