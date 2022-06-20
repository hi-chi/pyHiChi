####################### FFTW ############################
function(link_fftw)
	
	file(TO_CMAKE_PATH ${FFTW_DIR} FFTW_DIR)
		
	if (NOT EXISTS ${FFTW_DIR})		
		message(FATAL_ERROR "Cannot find fftw install directory '${FFTW_DIR}'")
	else()
	
		if (EXISTS "${FFTW_DIR}/include")
			set(FFT_INCLUDES ${FFTW_DIR}/include)
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
		
		if(USE_OMP)
			find_library(
				FFTW3_OMP_LIB
				NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX}
				PATHS "${FFTW_DIR}/lib" "${FFTW_DIR}/lib64"
			)
			if (NOT FFTW3_OMP_LIB)
				message(FATAL_ERROR "Cannot find fftw3_omp library in '${FFTW_DIR}/lib' or '${FFTW_DIR}/lib64' directory")
			endif()
		endif()
		
		set(FFT_LIBS ${FFTW3_OMP_LIB} ${FFTW3_LIB})
	
		set(FFT_INCLUDES ${FFT_INCLUDES} PARENT_SCOPE)
		set(FFT_LIBS ${FFT_LIBS} PARENT_SCOPE)
	
		message(STATUS "Using FFTW: ${FFT_LIBS}")
	
	endif()
	
endfunction()
####################### END FFTW ########################

####################### MKL #############################
function(link_mkl_fft)
	
	list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
	find_package(MKL)
	if (MKL_FOUND)
		set(FFT_INCLUDES ${MKL_INCLUDE_DIRS} PARENT_SCOPE)
		set(FFT_LIBS ${MKL_LIBRARIES} PARENT_SCOPE)
		message(STATUS "Using MKL")
	else()
		message(FATAL_ERROR "Cannot find MKL")
	endif()
	
endfunction()
####################### END MKL #########################
