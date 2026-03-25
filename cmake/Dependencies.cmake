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
# Coarrays (optional, multi-image coarrays via GFortran)
# ------------------------------------------------------------------------------
if(USE_COARRAYS)
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Coarrays currently only wired for GFortran.")
    endif()
    find_package(OpenCoarrays REQUIRED)
    find_library(COARRAYS_LIBRARY NAMES caf_mpi libcaf_mpi 
        PATHS /mnt/nasapps/development/modules_elmlund/opencoarrays/2.10.3/lib64
    REQUIRED)
    list(APPEND SIMPLE_LIBRARIES ${COARRAYS_LIBRARY})
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fcoarray=lib>)
    message(STATUS "Coarray support enabled (multi-image).")
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
find_package(FFTW) # (required)
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
    message(FATAL_ERROR
        "================================================================================\n"
        "  FFTW3 REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE requires FFTW3 development headers to compile.\n"
        "  Install fftw3-devel via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libfftw3-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install fftw3-devel\n"
        "    or: sudo dnf install fftw3-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S fftw\n"
        "  macOS (Homebrew):\n"
        "    brew install fftw\n"
        "  macOS (MacPorts):\n"
        "    sudo port install fftw\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
endif()

# ------------------------------------------------------------------------------
# TIFF and dependencies (required) 
# ------------------------------------------------------------------------------

find_package(TIFF) # (required)
if(TIFF_FOUND)
    find_package(JPEG) # (required by TIFF)
    if(NOT JPEG_FOUND)
        message(FATAL_ERROR
        "================================================================================\n"
        "  JPEG LIBRARY REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE requires JPEG development headers to compile.\n"
        "  Install libjpeg-devel via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libjpeg-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install libjpeg-devel\n"
        "    or: sudo dnf install libjpeg-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S libjpeg\n"
        "  macOS (Homebrew):\n"
        "    brew install libjpeg\n"
        "  macOS (MacPorts):\n"
        "    sudo port install jpeg\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
    endif()
    find_package(ZLIB) # (required by TIFF)
    if(NOT ZLIB_FOUND)
        message(FATAL_ERROR
            "================================================================================\n"
            "  ZLIB LIBRARY REQUIRED but NOT FOUND\n"
            "================================================================================\n"
            "  SIMPLE requires ZLIB development headers to compile.\n"
            "  Install zlib-devel via your system package manager:\n"
            "  Linux (Debian/Ubuntu):\n"
            "    sudo apt-get update\n"
            "    sudo apt-get install zlib1g-dev\n"
            "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
            "    sudo yum install zlib-devel\n"
            "    or: sudo dnf install zlib-devel\n"
            "  Linux (Arch):\n"
            "    sudo pacman -S zlib\n"
            "  macOS (Homebrew):\n"
            "    brew install zlib\n"
            "  macOS (MacPorts):\n"
            "    sudo port install zlib\n"
            "  After installation, reconfigure CMake:\n"
            "    rm -rf build/\n"
            "    cmake -B build\n"
            "================================================================================\n"
        )
    endif()
    find_library(LZMA_LIBRARY NAMES lzma)  # (optional, used by TIFF if available)
    find_library(ZSTD_LIBRARY NAMES zstd)  # (optional, used by TIFF if available)
    find_library(JBIG_LIBRARY NAMES jbig)  # (optional, used by TIFF if available)
    if(NOT APPLE)
        find_library(WEBP_LIBRARY NAMES webp) # (optional, used by TIFF if available)
        find_library(DEFLATE_LIBRARY NAMES deflate) # (optional, used by TIFF if available)
    endif()
    add_compile_definitions(USING_TIFF=1)
    include_directories(${TIFF_INCLUDE_DIR})
    list(APPEND SIMPLE_LIBRARIES
        ${TIFF_LIBRARY}
        ${JPEG_LIBRARY}
        ${ZLIB_LIBRARY}
    )
    if(LZMA_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${LZMA_LIBRARY})
        message(STATUS "  TIFF optional dep: lzma found: ${LZMA_LIBRARY}")
    else()
        message(STATUS "  TIFF optional dep: lzma NOT found")
    endif()
    if(ZSTD_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${ZSTD_LIBRARY})
        message(STATUS "  TIFF optional dep: zstd found: ${ZSTD_LIBRARY}")
    else()
        message(STATUS "  TIFF optional dep: zstd NOT found")
    endif()
    if(JBIG_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${JBIG_LIBRARY})
        message(STATUS "  TIFF optional dep: jbig found: ${JBIG_LIBRARY}")
    else()
        message(STATUS "  TIFF optional dep: jbig NOT found")
    endif()
    if(WEBP_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${WEBP_LIBRARY})
        message(STATUS "  TIFF optional dep: webp found: ${WEBP_LIBRARY}")
    else()
        message(STATUS "  TIFF optional dep: webp NOT found (or Apple)")
    endif()
    if(DEFLATE_LIBRARY)
        list(APPEND SIMPLE_LIBRARIES ${DEFLATE_LIBRARY})
        message(STATUS "  TIFF optional dep: deflate found: ${DEFLATE_LIBRARY}")
    else()
        message(STATUS "  TIFF optional dep: deflate NOT found (or Apple)")
    endif()
else()
    message(FATAL_ERROR
        "================================================================================\n"
        "  LIBTIFF REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE requires LIBTIFF development headers to compile.\n"
        "  Install libtiff-devel via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libtiff-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install libtiff-devel\n"
        "    or: sudo dnf install libtiff-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S libtiff\n"
        "  macOS (Homebrew):\n"
        "    brew install libtiff\n"
        "  macOS (MacPorts):\n"
        "    sudo port install libtiff\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
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
# librt - only available on unix
# ------------------------------------------------------------------------------
if(UNIX AND NOT APPLE)
    find_package(LibRt REQUIRED)
    list(APPEND SIMPLE_LIBRARIES LIBRT::LIBRT)
endif()

# ------------------------------------------------------------------------------
# libCurl (REQUIRED)
# ------------------------------------------------------------------------------
find_package(CURL)
if(NOT CURL_FOUND)
    message(FATAL_ERROR
        "================================================================================\n"
        "  libcurl REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE requires libcurl development headers to compile.\n"
        "  Install libcurl-devel via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libcurl4-openssl-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install libcurl-devel\n"
        "    or: sudo dnf install libcurl-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S curl\n"
        "  macOS (Homebrew):\n"
        "    brew install curl\n"
        "  macOS (MacPorts):\n"
        "    sudo port install curl\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
endif()
list(APPEND SIMPLE_LIBRARIES CURL::libcurl)
message(STATUS "libcurl: FOUND")

# ------------------------------------------------------------------------------
# Final dependency summary
# ------------------------------------------------------------------------------
message(STATUS "=== SIMPLE Dependency Summary ===")

message(STATUS " FFTW3:               YES")
message(STATUS " JPEG_LIBRARY:        YES")
message(STATUS " LIBCURL:             YES")
message(STATUS " LIBTIFF:             YES")
if(USE_MPI)
    message(STATUS " MPI:                 YES")
else()
    message(STATUS " MPI:                 NO")
endif()
message(STATUS " OpenMP:              YES")
message(STATUS " ZLIB_LIBRARY:        YES")
message(STATUS "==================================")
