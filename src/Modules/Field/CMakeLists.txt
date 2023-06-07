cmake_minimum_required(VERSION 3.1)

project(fieldModules)


include_directories(include
    ../../Core/include)

add_library(fieldModules_interface INTERFACE)
set(FIELDMODULES_INCLUDE_DIR include)	
set(FIELDMODULES_HEADER_DIR ${FIELDMODULES_INCLUDE_DIR})
set(fieldModules_headers
    ${FIELDMODULES_HEADER_DIR}/AnalyticalFieldSolver.h
    ${FIELDMODULES_HEADER_DIR}/Mapping.h
    ${FIELDMODULES_HEADER_DIR}/Field.h
    ${FIELDMODULES_HEADER_DIR}/FieldConfiguration.h
    ${FIELDMODULES_HEADER_DIR}/Fdtd.h
    ${FIELDMODULES_HEADER_DIR}/FieldBoundaryCondition.h
    ${FIELDMODULES_HEADER_DIR}/FieldBoundaryConditionFdtd.h
    ${FIELDMODULES_HEADER_DIR}/FieldBoundaryConditionSpectral.h
    ${FIELDMODULES_HEADER_DIR}/FieldGenerator.h
    ${FIELDMODULES_HEADER_DIR}/FieldGeneratorFdtd.h
    ${FIELDMODULES_HEADER_DIR}/FieldGeneratorSpectral.h
    ${FIELDMODULES_HEADER_DIR}/FieldSolver.h
    ${FIELDMODULES_HEADER_DIR}/Pml.h
    ${FIELDMODULES_HEADER_DIR}/PmlFdtd.h
    ${FIELDMODULES_HEADER_DIR}/PmlPstd.h
    ${FIELDMODULES_HEADER_DIR}/PmlPsatd.h
    ${FIELDMODULES_HEADER_DIR}/PmlSpectralTimeStraggered.h
    ${FIELDMODULES_HEADER_DIR}/Pstd.h
    ${FIELDMODULES_HEADER_DIR}/PsatdTimeStaggered.h
    ${FIELDMODULES_HEADER_DIR}/Psatd.h)
add_custom_target(fieldModules SOURCES ${fieldModules_headers})