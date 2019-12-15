cmake_minimum_required(VERSION 3.1)

project(tightFocusing)

include_directories(include
    ../../TightFocusing/include)

add_library(tightFocusingModules_interface INTERFACE)
set(TIGHTFOCUSINGMODULES_INCLUDE_DIR include)	
set(TIGHTFOCUSINGMODULES_HEADER_DIR ${TIGHTFOCUSINGMODULES_INCLUDE_DIR})
file(GLOB tightFocusingModules_headers include/*.h)

add_custom_target(tightFocusing SOURCES ${tightFocusingModules_headers})