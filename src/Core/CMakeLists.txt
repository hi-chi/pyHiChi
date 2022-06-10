cmake_minimum_required(VERSION 3.1)

project(core)

include_directories(include
    ${FFT_INCLUDES})

add_library(core_interface INTERFACE)
set(CORE_INCLUDE_DIR include)	
set(CORE_HEADER_DIR ${CORE_INCLUDE_DIR})

set(core_headers
    ${CORE_HEADER_DIR}/Allocators.h
    ${CORE_HEADER_DIR}/AnalyticalField.h
    ${CORE_HEADER_DIR}/Constants.h
    ${CORE_HEADER_DIR}/Enums.h
    ${CORE_HEADER_DIR}/Dimension.h
    ${CORE_HEADER_DIR}/Ensemble.h
    ${CORE_HEADER_DIR}/FieldValue.h
    ${CORE_HEADER_DIR}/FormFactor.h
    ${CORE_HEADER_DIR}/FourierTransform.h
    ${CORE_HEADER_DIR}/FP.h
    ${CORE_HEADER_DIR}/Grid.h
    ${CORE_HEADER_DIR}/GridTypes.h
    ${CORE_HEADER_DIR}/Particle.h
    ${CORE_HEADER_DIR}/ParticleArray.h
    ${CORE_HEADER_DIR}/ParticleTraits.h
    ${CORE_HEADER_DIR}/ParticleTypes.h
    ${CORE_HEADER_DIR}/ScalarField.h
    ${CORE_HEADER_DIR}/Vectors.h
    ${CORE_HEADER_DIR}/VectorsProxy.h
    ${CORE_HEADER_DIR}/macros.h)

if (FFT_INCLUDES)
    set(core_headers ${core_headers} 
        ${FFT_INCLUDES}/fftw3.h)
endif()
add_custom_target(core SOURCES ${core_headers})
