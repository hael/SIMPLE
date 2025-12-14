# Find and configure all dependencies

set(SIMPLE_LIBRARIES "")

# ============================================================================
# OpenMP
# ============================================================================
if(USE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran C CXX)
    if(OpenMP_Fortran_FOUND)
        add_definitions(-DOPENMP)
        list(APPEND SIMPLE_LIBRARIES OpenMP::OpenMP_Fortran)
        message(STATUS "OpenMP enabled")
    else()
        message(FATAL_ERROR "OpenMP requested but not found")
    endif()
endif()

# ============================================================================
# OpenACC
# ============================================================================
if(USE_OPENACC)
    # OpenACC support for GCC (requires -fopenacc flag)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fopenacc>)
    add_definitions(-DOPENACC)
    message(STATUS "OpenACC enabled (requires GCC 6.0+)")
endif()

# ============================================================================
# Coarrays
# ============================================================================
if(USE_COARRAYS)
    # Coarray support for GFortran
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=single>)
    message(STATUS "Coarray support enabled (single-image mode)")
    message(STATUS "For multi-image coarrays, use: -fcoarray=lib and link with OpenCoarrays")
endif()

# ============================================================================
# MPI
# ============================================================================
if(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran C)
    if(MPI_Fortran_FOUND)
        add_definitions(-DUSING_MPI)
        list(APPEND SIMPLE_LIBRARIES MPI::MPI_Fortran)
        message(STATUS "MPI enabled")
    else()
        message(FATAL_ERROR "MPI requested but not found")
    endif()
endif()

# ============================================================================
# FFTW3
# ============================================================================
find_package(FFTW REQUIRED)
if(FFTW_FOUND)

    # Need single-precision FFTW
    find_library(FFTW_FLOAT_LIBRARY
        NAMES fftw3f libfftw3f libfftw3f-3
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    if(NOT FFTW_FLOAT_LIBRARY)
        message(FATAL_ERROR "FFTW single-precision library (fftw3f) not found")
    endif()
    list(APPEND SIMPLE_LIBRARIES ${FFTW_FLOAT_LIBRARY})

    # Need double-precision FFTW
    find_library( FFTW_DOUBLE_PRECISION_LIBRARIES
        NAMES fftw3 libfftw3 libfftw3-3
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    if(NOT FFTW_DOUBLE_PRECISION_LIBRARIES)
        message(FATAL_ERROR "FFTW double-precision library (fftw3) not found")
    endif()
    list(APPEND SIMPLE_LIBRARIES ${FFTW_DOUBLE_PRECISION_LIBRARIES})
    
    # Also need FFTW with threading support
    find_library(FFTW_THREADS_LIBRARY 
        NAMES fftw3f_threads libfftw3f_threads
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    
    find_library(FFTW_OMP_LIBRARY 
        NAMES fftw3f_omp libfftw3f_omp
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    
    # if(TARGET FFTW::FFTW)
    #     list(APPEND SIMPLE_LIBRARIES FFTW::FFTW)
    # else()
    #     include_directories(${FFTW_INCLUDE_DIRS})
    #     list(APPEND SIMPLE_LIBRARIES ${FFTW_LIBRARIES})
    # endif()
    
    # Add threading library if found
    if(FFTW_THREADS_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${FFTW_THREADS_LIBRARY})
        message(STATUS "FFTW3 threads library found: ${FFTW_THREADS_LIBRARY}")
    elseif(FFTW_OMP_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${FFTW_OMP_LIBRARY})
        message(STATUS "FFTW3 OMP library found: ${FFTW_OMP_LIBRARY}")
    endif()
    
    message(STATUS "FFTW3 found: ${FFTW_LIBRARIES}")
else()
    message(FATAL_ERROR "FFTW3 not found. Install with: brew install fftw (macOS) or apt-get install libfftw3-dev (Linux)")
endif()

# ============================================================================
# TIFF and dependencies
# ============================================================================
if(USE_LIBTIFF)
    find_package(TIFF REQUIRED)
    
    # Find TIFF dependencies
    find_library(JPEG_LIBRARY NAMES jpeg REQUIRED)
    find_library(ZLIB_LIBRARY NAMES z REQUIRED)
    find_library(LZMA_LIBRARY NAMES lzma)
    find_library(ZSTD_LIBRARY NAMES zstd)
    find_library(JBIG_LIBRARY NAMES jbig)
    
    add_definitions(-DUSING_TIFF=1)
    include_directories(${TIFF_INCLUDE_DIR})
    
    list(APPEND SIMPLE_LIBRARIES 
        ${TIFF_LIBRARY}
        ${JPEG_LIBRARY}
        ${ZLIB_LIBRARY}
    )
    
    # Optional TIFF dependencies
    if(LZMA_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${LZMA_LIBRARY})
    endif()
    if(ZSTD_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${ZSTD_LIBRARY})
    endif()
    if(JBIG_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${JBIG_LIBRARY})
    endif()
    
    # Platform-specific TIFF dependencies
    if(NOT APPLE)
        find_library(WEBP_LIBRARY NAMES webp)
        find_library(DEFLATE_LIBRARY NAMES deflate)
        if(WEBP_LIBRARY)
            list(APPEND SIMPLE_LIBRARIES ${WEBP_LIBRARY})
        endif()
        if(DEFLATE_LIBRARY)
            list(APPEND SIMPLE_LIBRARIES ${DEFLATE_LIBRARY})
        endif()
    endif()
    
    message(STATUS "TIFF support enabled")
endif()

# ============================================================================
# Threading
# ============================================================================
find_package(Threads REQUIRED)
list(APPEND SIMPLE_LIBRARIES Threads::Threads)

# Math library (usually needed on Linux)
if(UNIX AND NOT APPLE)
    list(APPEND SIMPLE_LIBRARIES m)
endif()

# Export the library list for use in other CMakeLists
set(SIMPLE_LIBRARIES ${SIMPLE_LIBRARIES} PARENT_SCOPE)

message(STATUS "All dependencies found")