# cmake/Dependencies.cmake
#
# Finds and configures external dependencies.
# Populates SIMPLE_LIBRARIES for linking by the SIMPLE library and executables.

set(SIMPLE_LIBRARIES "")

# ------------------------------------------------------------------------------
# OpenMP (required by default)
# ------------------------------------------------------------------------------

if(USE_OPENMP)
    # CMake 3.25 has good OpenMP Fortran support with GCC 14+
    find_package(OpenMP COMPONENTS Fortran C CXX)
    if(OpenMP_Fortran_FOUND AND OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
        message(STATUS "OpenMP: using CMake OpenMP targets")
        add_compile_definitions(OPENMP)
        list(APPEND SIMPLE_LIBRARIES
            OpenMP::OpenMP_Fortran
            OpenMP::OpenMP_C
            OpenMP::OpenMP_CXX
        )
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        # Fallback: manual flags for GNU if FindOpenMP fails for some reason
        message(WARNING
            "OpenMP: CMake's FindOpenMP failed; falling back to -fopenmp for GNU.\n"
            "Check your toolchain / CMake version if this persists."
        )
        string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
        string(APPEND CMAKE_C_FLAGS       " -fopenmp")
        string(APPEND CMAKE_CXX_FLAGS     " -fopenmp")
        add_compile_definitions(OPENMP)
        # libgomp is usually pulled automatically by -fopenmp
    else()
        message(FATAL_ERROR
            "OpenMP requested but could not be configured.\n"
            "Compiler is not GNU and FindOpenMP failed."
        )
    endif()
endif()

# ------------------------------------------------------------------------------
# OpenACC (optional, Fortran only, GNU)
# ------------------------------------------------------------------------------
if(USE_OPENACC)
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "OpenACC currently only wired for GFortran.")
    endif()
    add_compile_definitions(OPENACC)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fopenacc>)
    message(STATUS "OpenACC enabled (-fopenacc)")
endif()

# ------------------------------------------------------------------------------
# Coarrays (optional, single-image coarrays via GFortran)
# ------------------------------------------------------------------------------
if(USE_COARRAYS)
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Coarrays currently only wired for GFortran.")
    endif()
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=single>)
    message(STATUS "Coarray support enabled (single-image).")
    message(STATUS "For multi-image coarrays, use -fcoarray=lib and link OpenCoarrays manually.")
endif()

# ------------------------------------------------------------------------------
# MPI (optional)
# ------------------------------------------------------------------------------
if(USE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran C)
    if(MPI_Fortran_FOUND)
        add_compile_definitions(USING_MPI)
        list(APPEND SIMPLE_LIBRARIES MPI::MPI_Fortran)
        message(STATUS "MPI enabled")
    else()
        message(FATAL_ERROR "MPI requested but not found")
    endif()
endif()

# ------------------------------------------------------------------------------
# FFTW3 (required)
#   Uses custom FindFFTW.cmake (FFTW::FFTW imported target)
# ------------------------------------------------------------------------------
find_package(FFTW REQUIRED)
if(FFTW_FOUND)
    # Single-precision FFTW
    find_library(FFTW_FLOAT_LIBRARY
        NAMES fftw3f libfftw3f libfftw3f-3
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    if(NOT FFTW_FLOAT_LIBRARY)
        message(FATAL_ERROR "FFTW single-precision library (fftw3f) not found")
    endif()
    list(APPEND SIMPLE_LIBRARIES ${FFTW_FLOAT_LIBRARY})
    # Double-precision FFTW
    find_library(FFTW_DOUBLE_PRECISION_LIBRARIES
        NAMES fftw3 libfftw3 libfftw3-3
        HINTS ${FFTW_LIBRARY_DIRS}
        PATH_SUFFIXES lib lib64
    )
    if(NOT FFTW_DOUBLE_PRECISION_LIBRARIES)
        message(FATAL_ERROR "FFTW double-precision library (fftw3) not found")
    endif()
    list(APPEND SIMPLE_LIBRARIES ${FFTW_DOUBLE_PRECISION_LIBRARIES})
    # Threaded FFTW (optional)
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
    if(FFTW_THREADS_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${FFTW_THREADS_LIBRARY})
        message(STATUS "FFTW3 threads library found: ${FFTW_THREADS_LIBRARY}")
    elseif(FFTW_OMP_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${FFTW_OMP_LIBRARY})
        message(STATUS "FFTW3 OMP library found: ${FFTW_OMP_LIBRARY}")
    endif()
    message(STATUS "FFTW3 base libraries: ${FFTW_LIBRARIES}")
else()
    message(FATAL_ERROR "FFTW3 not found. Install e.g. libfftw3-dev / fftw.")
endif()

# ------------------------------------------------------------------------------
# TIFF and dependencies (optional)
# ------------------------------------------------------------------------------
if(USE_LIBTIFF)
    find_package(TIFF REQUIRED)
    find_library(JPEG_LIBRARY NAMES jpeg REQUIRED)
    find_library(ZLIB_LIBRARY NAMES z REQUIRED)
    find_library(LZMA_LIBRARY NAMES lzma)
    find_library(ZSTD_LIBRARY NAMES zstd)
    find_library(JBIG_LIBRARY NAMES jbig)
    add_compile_definitions(USING_TIFF=1)
    include_directories(${TIFF_INCLUDE_DIR})
    list(APPEND SIMPLE_LIBRARIES
        ${TIFF_LIBRARY}
        ${JPEG_LIBRARY}
        ${ZLIB_LIBRARY}
    )
    if(LZMA_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${LZMA_LIBRARY})
    endif()
    if(ZSTD_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${ZSTD_LIBRARY})
    endif()
    if(JBIG_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${JBIG_LIBRARY})
    endif()
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

# ------------------------------------------------------------------------------
# Threads + libm
# ------------------------------------------------------------------------------
find_package(Threads REQUIRED)
list(APPEND SIMPLE_LIBRARIES Threads::Threads)
if(UNIX AND NOT APPLE)
    list(APPEND SIMPLE_LIBRARIES m)
endif()

# ------------------------------------------------------------------------------
# Export to parent
# ------------------------------------------------------------------------------
set(SIMPLE_LIBRARIES "${SIMPLE_LIBRARIES}" PARENT_SCOPE)
message(STATUS "All dependencies configured")
