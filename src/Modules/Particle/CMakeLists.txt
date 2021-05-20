cmake_minimum_required(VERSION 3.1)

project(particleModules)

include_directories(include
    ../../Core/include)

add_library(particleModules_interface INTERFACE)
set(PARTICLEMODULES_INCLUDE_DIR include)	
set(PARTICLEMODULES_HEADER_DIR ${PARTICLEMODULES_INCLUDE_DIR})
set(particleModules_headers
    ${PARTICLEMODULES_HEADER_DIR}/Pusher.h
    ${PARTICLEMODULES_HEADER_DIR}/QED_AEG.h
    ${PARTICLEMODULES_HEADER_DIR}/Species.h
    ${PARTICLEMODULES_HEADER_DIR}/synchrotron.h

    ${PARTICLEMODULES_HEADER_DIR}/Merging.h
    ${PARTICLEMODULES_HEADER_DIR}/Thinning.h
)

set(particleModulesSYCL_headers
    ${PARTICLEMODULES_HEADER_DIR}/sycl/PusherSYCL.h
)

add_custom_target(particleModules SOURCES ${particleModules_headers} ${particleModulesSYCL_headers})

source_group("headers"   FILES ${particleModules_headers})
source_group("sycl_headers"     FILES ${particleModulesSYCL_headers})