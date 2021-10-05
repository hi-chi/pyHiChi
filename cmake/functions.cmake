function(link_fft_libs)
	
	####################### MKL #############################
	
	if(USE_MKL)
		list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
		find_package(MKL)
		if (MKL_FOUND)
			set(FFT_INCLUDES ${MKL_INCLUDE_DIRS} PARENT_SCOPE)
			set(FFT_LIBS ${MKL_LIBRARIES} PARENT_SCOPE)
			message(STATUS "Using MKL")
		else()
			message(WARNING "Cannot find MKL")
		endif()
	endif()
	
	####################### END MKL #########################
	
	####################### FFTW ############################
	
	if (USE_FFTW)
		
		if (NOT FFTW_DIR)
			
			set(FFTW_VERSION 3.3.8)
			set(INSTALL_DIR "${CMAKE_BINARY_DIR}/3rdparty")
			include(ExternalProject)
			ExternalProject_Add(project_fftw
					URL "http://fftw.org/fftw-${FFTW_VERSION}.tar.gz"
					PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
					CMAKE_ARGS
						"-DCMAKE_CFLAGS=${CMAKE_C_FLAGS}"
						"-DBUILD_SHARED_LIBS=OFF"
						"-DCMAKE_POSITION_INDEPENDENT_CODE=ON"
						"-DENABLE_AVX2=ON"
						"-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/fftw" 
						"-DENABLE_OPENMP=${USE_OMP}"
						"-DBUILD_TESTS=OFF"
			)
			install(DIRECTORY "${INSTALL_DIR}" DESTINATION .)
			set(FFTW_DIR ${INSTALL_DIR}/fftw)
			
			set(FFT_INCLUDES ${FFTW_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
			
			set(FFTW3_LIB ${FFTW_DIR}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
			if(USE_OMP)
				set(FFTW3_OMP_LIB ${FFTW_DIR}/${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX})
			endif()
		
		else()
		
			if (EXISTS "${FFTW_DIR}/include")
				set(FFT_INCLUDES ${FFTW_DIR}/include)
			else()
				message(WARNING "Cannot find 'include' directory in '${FFTW_DIR}' directory")
			endif()
			
			find_library(
				FFTW3_LIB
				NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}
				PATHS "${FFTW_DIR}/lib" "${FFTW_DIR}/lib64"
			)
			if (NOT FFTW3_LIB)
				message(WARNING "Cannot find fftw3 library in '${FFTW_DIR}/lib' or '${FFTW_DIR}/lib64' directory")
			endif()
			
			if(USE_OMP)
				find_library(
					FFTW3_OMP_LIB
					NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX}
					PATHS "${FFTW_DIR}/lib" "${FFTW_DIR}/lib64"
				)
				if (NOT FFTW3_OMP_LIB)
					message(WARNING "Cannot find fftw3_omp library in '${FFTW_DIR}/lib' or '${FFTW_DIR}/lib64' directory")
				endif()
			endif()
		
		endif()
		
		set(FFT_LIBS ${FFTW3_OMP_LIB} ${FFTW3_LIB})
		
		set(FFT_INCLUDES ${FFT_INCLUDES} PARENT_SCOPE)
		set(FFT_LIBS ${FFT_LIBS} PARENT_SCOPE)
		
		message(STATUS "Using FFTW: ${FFT_LIBS}")
		
	endif()
	
	####################### END FFTW ########################
	
endfunction()
