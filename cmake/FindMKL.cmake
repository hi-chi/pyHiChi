#  Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - MKL libraries
#
#  The environment variable MKLROOT is used
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()

# If already in cache, be silent
if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES)
    set (MKL_FIND_QUIETLY TRUE)
endif()

set(MKLROOT $ENV{MKLROOT})
    
find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS ${MKLROOT}/include)
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR} ${MKL_INCLUDE_DIR}/fftw)    
set(MKL_LIBRARIES_PATHS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)

if (USE_OMP)

    function(find_libomp5_win32)
        find_library(IOMP5_LIB libiomp5md${CMAKE_STATIC_LIBRARY_SUFFIX} $ENV{INTEL}/compiler/lib/intel64)
        if (IOMP5_LIB)
            set(IOMP5 ${IOMP5_LIB} PARENT_SCOPE)
            message(STATUS "Found libiomp5")
        else()
            message(WARNING "Cannot find libiomp5")
        endif()
    endfunction()

    # static linking is preferred
    find_library( MKL_CDFT_CORE_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_cdft_core${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_cdft_core PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_INTEL_ILP64_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_intel_ilp64${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_intel_ilp64 PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_INTEL_THREAD_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_intel_thread${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_intel_thread PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_CORE_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_core${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_core PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_BLACS_INTELMPI_ILP64_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_blacs_intelmpi_ilp64${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_blacs_intelmpi_ilp64 PATHS ${MKL_LIBRARIES_PATHS} )
    
    if (MKL_INCLUDE_DIRS AND
         MKL_CDFT_CORE_LIB AND
         MKL_INTEL_ILP64_LIB AND
         MKL_INTEL_THREAD_LIB AND 
         MKL_CORE_LIB AND
         MKL_BLACS_INTELMPI_ILP64_LIB)
         
        set(MKL_LIBRARIES ${MKL_CDFT_CORE_LIB}
                ${MKL_INTEL_ILP64_LIB} 
                ${MKL_INTEL_THREAD_LIB} 
                ${MKL_CORE_LIB} 
                ${MKL_BLACS_INTELMPI_ILP64_LIB})
        
        # ICC
        if (${CMAKE_CXX_COMPILER} MATCHES "icc.*$" OR ${CMAKE_CXX_COMPILER} MATCHES "icl.*$")
            if (WIN32)
                find_libomp5_win32()
                set(MKL_LIBRARIES ${MKL_LIBRARIES} ${IOMP5})
            elseif(UNIX)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -liomp5 -lpthread -lm -ldl")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
            endif()
        endif()
        
        # GCC
        if (${CMAKE_CXX_COMPILER} MATCHES "gcc.*$")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp -lpthread -lm -ldl")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
        endif()
        
        # MSVC
        if (MSVC AND NOT ${CMAKE_CXX_COMPILER} MATCHES "icl.*$") 
            find_libomp5_win32()
            set(MKL_LIBRARIES ${MKL_LIBRARIES} ${IOMP5})
        endif()
           
        # Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
        # all listed variables are TRUE.
        INCLUDE(FindPackageHandleStandardArgs)
        FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS)
        
        MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES)
        
    else()
        message(WARNING "Cannot find some necessary mkl components")
        message(STATUS "MKL_INCLUDE_DIRS=${MKL_INCLUDE_DIRS}")
        message(STATUS "MKL_CDFT_CORE_LIB=${MKL_CDFT_CORE_LIB}")
        message(STATUS "MKL_INTEL_ILP64_LIB=${MKL_INTEL_ILP64_LIB}")
        message(STATUS "MKL_INTEL_THREAD_LIB=${MKL_INTEL_SEQUENTIAL_LIB}")
        message(STATUS "MKL_CORE_LIB=${MKL_CORE_LIB}")
        message(STATUS "MKL_BLACS_INTELMPI_ILP64_LIB=${MKL_BLACS_INTELMPI_ILP64_LIB}")
    endif()
    
    
else()  # USE_OMP

    # static linking is preferred
    find_library( MKL_CDFT_CORE_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_cdft_core${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_cdft_core PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_INTEL_ILP64_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_intel_ilp64${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_intel_ilp64 PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_INTEL_SEQUENTIAL_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_sequential${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_sequential PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_CORE_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_core${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_core PATHS ${MKL_LIBRARIES_PATHS} )
    find_library( MKL_BLACS_INTELMPI_ILP64_LIB NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}mkl_blacs_intelmpi_ilp64${CMAKE_STATIC_LIBRARY_SUFFIX} mkl_blacs_intelmpi_ilp64 PATHS ${MKL_LIBRARIES_PATHS} )
        
    if (MKL_INCLUDE_DIRS AND
         MKL_CDFT_CORE_LIB AND
         MKL_INTEL_ILP64_LIB AND
         MKL_INTEL_SEQUENTIAL_LIB AND 
         MKL_CORE_LIB AND
         MKL_BLACS_INTELMPI_ILP64_LIB)
         
        set(MKL_LIBRARIES ${MKL_CDFT_CORE_LIB}
                ${MKL_INTEL_ILP64_LIB} 
                ${MKL_INTEL_SEQUENTIAL_LIB} 
                ${MKL_CORE_LIB} 
                ${MKL_BLACS_INTELMPI_ILP64_LIB})
                
        if (UNIX)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -lm -ldl")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
        endif()
           
        # Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
        # all listed variables are TRUE.
        INCLUDE(FindPackageHandleStandardArgs)
        FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS)
        
        MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES)
        
    else()
        message(WARNING "Cannot find some necessary mkl components")
        message(STATUS "MKL_INCLUDE_DIRS=${MKL_INCLUDE_DIRS}")
        message(STATUS "MKL_CDFT_CORE_LIB=${MKL_CDFT_CORE_LIB}")
        message(STATUS "MKL_INTEL_ILP64_LIB=${MKL_INTEL_ILP64_LIB}")
        message(STATUS "MKL_INTEL_SEQUENTIAL_LIB=${MKL_INTEL_SEQUENTIAL_LIB}")
        message(STATUS "MKL_CORE_LIB=${MKL_CORE_LIB}")
        message(STATUS "MKL_BLACS_INTELMPI_ILP64_LIB=${MKL_BLACS_INTELMPI_ILP64_LIB}")
    endif()
    

endif()  # USE_OMP