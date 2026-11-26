cmake_minimum_required(VERSION 4.2)

set(CMAKE_CXX_STANDARD 17)

option(USE_OMP ON)

project(fftw)

add_compile_options(
	"$<$<AND:$<COMPILE_LANG_AND_ID:C,IntelLLVM>,$<PLATFORM_ID:Windows>>:/fp:precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>,$<PLATFORM_ID:Windows>>:/fp:precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:C,IntelLLVM>,$<NOT:$<PLATFORM_ID:Windows>>>:-fp-model=precise>"
	"$<$<AND:$<COMPILE_LANG_AND_ID:CXX,IntelLLVM>,$<NOT:$<PLATFORM_ID:Windows>>>:-fp-model=precise>"
)

if (USE_OMP)
	add_definitions(-D__USE_OMP__)
	find_package(OpenMP REQUIRED)
endif()

# Rewrites the cmake version of FFTW
# This is correct and does not generate build problems
set(FFTW_CMAKE_PATH_FILE "${CMAKE_BINARY_DIR}/patch_fftw.cmake")
file(WRITE ${FFTW_CMAKE_PATH_FILE}
"   file(READ \"\${FILE}\" TEXT)\n\
    string(REGEX REPLACE \"cmake_minimum_required[^\\n]*\" \"cmake_minimum_required(VERSION 4.2)\" TEXT \"\${TEXT}\")\n\
    file(WRITE \"\${FILE}\" \"\${TEXT}\")
"
)

set(FFTW_VERSION 3.3.10)
set(FFTW_DIR "${CMAKE_SOURCE_DIR}")
set(FFTW_INSTALL_DIR "${FFTW_DIR}/fftw-install")
set(FFTW_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

include(ExternalProject)
ExternalProject_Add(fftw
		URL "http://fftw.org/fftw-${FFTW_VERSION}.tar.gz"
		PREFIX ${FFTW_BUILD_DIR}
		DOWNLOAD_DIR ${FFTW_BUILD_DIR}
		SOURCE_DIR ${FFTW_DIR}/fftw
		INSTALL_DIR ${FFTW_INSTALL_DIR}
		TMP_DIR ${FFTW_BUILD_DIR}/tmp
		STAMP_DIR ${FFTW_BUILD_DIR}/stamp
		BINARY_DIR ${FFTW_BUILD_DIR}/build
        BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release -j
		PATCH_COMMAND ${CMAKE_COMMAND} -DFILE=<SOURCE_DIR>/CMakeLists.txt -P ${FFTW_CMAKE_PATH_FILE}
		CMAKE_ARGS
			"-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
			"-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
			"-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
			"-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
			"-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
			"-DBUILD_SHARED_LIBS=OFF"
			"-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
			"-DENABLE_AVX2=ON"
			"-DCMAKE_INSTALL_PREFIX=${FFTW_INSTALL_DIR}" 
			"-DENABLE_OPENMP=${USE_OMP}"
			"-DBUILD_TESTS=OFF"
)
install(DIRECTORY "${FFTW_DIR}" DESTINATION .)

message(STATUS "FFTW install directory is '${FFTW_INSTALL_DIR}'")
