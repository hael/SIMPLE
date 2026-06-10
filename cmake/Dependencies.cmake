# cmake/Dependencies.cmake
#
# Finds and configures external dependencies.
# Populates SIMPLE_LIBRARIES for linking by the SIMPLE library and executables.


set(SIMPLE_LIBRARIES "")
set(SIMPLE_DEPENDENCY_TARGETS "")

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
    # OpenMP device offloading & some CUDA support
    if( USE_OPENMP_OFFLOAD )
        if(APPLE)
            set(USE_OPENMP_OFFLOAD OFF)
            message(WARNING
                "OpenMP offload requested but not supported on Apple platforms.\n"
                "Falling back to standard OpenMP."
            )
        else()
            find_package(CUDAToolkit REQUIRED)
            if(NOT CUDAToolkit_Fortran_FOUND)
                list(APPEND SIMPLE_LIBRARIES
                    CUDA::cudart
                    CUDA::cufft
                    CUDA::cublas)
            else()
                message(FATAL_ERROR
                    "The CUDAToolkit library was not found.\n"
                )
            endif()
            string(APPEND CMAKE_Fortran_FLAGS " -foffload=nvptx-none")
            add_compile_definitions(USE_OPENMP_OFFLOAD)
        endif()
    endif()
else()
    set(USE_OPENMP_OFFLOAD OFF)
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
    add_compile_definitions(USE_COARRAYS)
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
# BLAS / LAPACK / ARPACK (required)
# ------------------------------------------------------------------------------
# SIMPLE uses the standard LP64 Fortran ABI for these libraries. If you point
# CMake at an ILP64 vendor build, integer sizes will not match the current code.

set(SIMPLE_OPENBLAS_VERSION "0.3.27" CACHE STRING "OpenBLAS version used for the bundled BLAS/LAPACK build")
set(SIMPLE_ARPACK_NG_VERSION "3.9.1" CACHE STRING "ARPACK-NG version used for the bundled ARPACK build")
set(SIMPLE_NUMERICS_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/_deps/numerics/install"
    CACHE PATH "Install prefix for BLAS/LAPACK/ARPACK built in the CMake build folder"
)

set(SIMPLE_BLAS_STATUS "system")
set(SIMPLE_LAPACK_STATUS "system")
set(SIMPLE_ARPACK_STATUS "system")
set(SIMPLE_NUMERICS_MISSING "")
set(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS OFF)
set(SIMPLE_NUMERICS_BUILD_ARPACK OFF)
set(SIMPLE_NUMERICS_BUNDLED_LINK_DIRS "")

find_package(BLAS)
if(BLAS_FOUND)
    message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
    set(SIMPLE_BLAS_STATUS "system")
else()
    message(STATUS "BLAS libraries: NOT FOUND")
    list(APPEND SIMPLE_NUMERICS_MISSING "BLAS")
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
    message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
    set(SIMPLE_LAPACK_STATUS "system")
else()
    message(STATUS "LAPACK libraries: NOT FOUND")
    list(APPEND SIMPLE_NUMERICS_MISSING "LAPACK")
endif()

find_library(ARPACK_LIBRARY
    NAMES arpack arpack-ng libarpack
    HINTS
        ${ARPACK_ROOT}
        $ENV{ARPACK_ROOT}
        $ENV{ARPACK_DIR}
        $ENV{ARPACKDIR}
    PATH_SUFFIXES lib lib64
    PATHS
        /opt/homebrew/opt/arpack
        /usr/local/opt/arpack
        /usr/local
        /usr
        /opt/local
        /opt/homebrew
        /opt
)
if(ARPACK_LIBRARY)
    message(STATUS "ARPACK library found: ${ARPACK_LIBRARY}")
    set(SIMPLE_ARPACK_STATUS "system")
else()
    message(STATUS "ARPACK library: NOT FOUND")
    list(APPEND SIMPLE_NUMERICS_MISSING "ARPACK")
endif()

if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    set(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS ON)
endif()
if(NOT ARPACK_LIBRARY)
    set(SIMPLE_NUMERICS_BUILD_ARPACK ON)
endif()

if(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS OR SIMPLE_NUMERICS_BUILD_ARPACK)
    include(ExternalProject)
    file(MAKE_DIRECTORY
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib64"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/include"
    )
    list(APPEND SIMPLE_NUMERICS_BUNDLED_LINK_DIRS
        "-L${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib"
        "-L${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib64"
    )
    set(SIMPLE_NUMERICS_BUNDLED_PROJECTS "")
    if(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS)
        list(APPEND SIMPLE_NUMERICS_BUNDLED_PROJECTS
            "OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (BLAS/LAPACK)"
        )
    endif()
    if(SIMPLE_NUMERICS_BUILD_ARPACK)
        list(APPEND SIMPLE_NUMERICS_BUNDLED_PROJECTS
            "ARPACK-NG ${SIMPLE_ARPACK_NG_VERSION}"
        )
    endif()
    list(JOIN SIMPLE_NUMERICS_BUNDLED_PROJECTS ", " SIMPLE_NUMERICS_BUNDLED_TEXT)
    list(JOIN SIMPLE_NUMERICS_MISSING ", " SIMPLE_NUMERICS_MISSING_TEXT)
    set(SIMPLE_NUMERICS_REASON
        "  The following required numerical libraries were NOT FOUND: ${SIMPLE_NUMERICS_MISSING_TEXT}\n"
    )
    message(WARNING
        "================================================================================\n"
        "  BLAS / LAPACK / ARPACK REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "${SIMPLE_NUMERICS_REASON}"
        "  SIMPLE requires BLAS, LAPACK and ARPACK development libraries for numerical\n"
        "  routines, dense solvers, and sparse eigenvalue routines.\n"
        "\n"
        "  Install them via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libopenblas-dev liblapack-dev libarpack2-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install openblas-devel lapack-devel arpack-devel\n"
        "    or: sudo dnf install openblas-devel lapack-devel arpack-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S openblas lapack arpack\n"
        "  macOS (Homebrew):\n"
        "    brew install openblas lapack arpack\n"
        "  macOS (MacPorts):\n"
        "    sudo port install OpenBLAS lapack arpack-ng\n"
        "\n"
        "  If installed in a non-standard prefix, pass it to CMake:\n"
        "    cmake -B build -DCMAKE_PREFIX_PATH=<dependency_prefix>\n"
        "\n"
        "  SIMPLE expects the LP64 ABI, where Fortran integer arguments are 32-bit.\n"
        "  Avoid ILP64 BLAS/LAPACK/ARPACK builds unless SIMPLE is changed consistently.\n"
        "\n"
        "  WARNING: SIMPLE will continue by downloading and compiling bundled libraries\n"
        "  inside the CMake build folder during the build step:\n"
        "    ${SIMPLE_NUMERICS_BUNDLED_TEXT}\n"
        "  Build/download location:\n"
        "    ${CMAKE_BINARY_DIR}/_deps\n"
        "  Bundled install prefix:\n"
        "    ${SIMPLE_NUMERICS_INSTALL_PREFIX}\n"
        "\n"
        "  To use system libraries instead, install the missing libraries and reconfigure:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
endif()

if(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS)
    find_program(SIMPLE_MAKE_PROGRAM NAMES gmake make)
    if(NOT SIMPLE_MAKE_PROGRAM)
        message(FATAL_ERROR
            "Building bundled OpenBLAS requires GNU make or make."
        )
    endif()

    set(SIMPLE_OPENBLAS_LIBRARY
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    set(SIMPLE_OPENBLAS_URL
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${SIMPLE_OPENBLAS_VERSION}/OpenBLAS-${SIMPLE_OPENBLAS_VERSION}.tar.gz"
    )
    set(SIMPLE_OPENBLAS_MAKE_ARGS
        "CC=${CMAKE_C_COMPILER}"
        "FC=${CMAKE_Fortran_COMPILER}"
        "USE_THREAD=1"
        "USE_OPENMP=0"
        "NO_SHARED=1"
        "NO_LAPACK=0"
        "INTERFACE64=0"
    )

    ExternalProject_Add(simple_openblas_ep
        URL "${SIMPLE_OPENBLAS_URL}"
        PREFIX "${CMAKE_BINARY_DIR}/_deps/openblas"
        DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/_deps/downloads"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE TRUE
        BUILD_COMMAND "${SIMPLE_MAKE_PROGRAM}" ${SIMPLE_OPENBLAS_MAKE_ARGS} libs netlib
        INSTALL_COMMAND "${SIMPLE_MAKE_PROGRAM}" ${SIMPLE_OPENBLAS_MAKE_ARGS}
            "PREFIX=${SIMPLE_NUMERICS_INSTALL_PREFIX}" install
        BUILD_BYPRODUCTS "${SIMPLE_OPENBLAS_LIBRARY}"
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
    )
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_openblas_ep)

    set(SIMPLE_BLAS_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION}")
    set(SIMPLE_LAPACK_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION}")
    set(SIMPLE_NUMERICS_PROVIDER_LIBRARIES openblas)
else()
    set(SIMPLE_NUMERICS_PROVIDER_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
endif()

if(SIMPLE_NUMERICS_BUILD_ARPACK)
    set(SIMPLE_ARPACK_LIBRARY
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}arpack${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    set(SIMPLE_ARPACK_NG_URL
        "https://github.com/opencollab/arpack-ng/archive/refs/tags/${SIMPLE_ARPACK_NG_VERSION}.tar.gz"
    )
    if(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS)
        set(SIMPLE_ARPACK_BLAS_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
        set(SIMPLE_ARPACK_LAPACK_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
    else()
        set(SIMPLE_ARPACK_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
        set(SIMPLE_ARPACK_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
    endif()
    string(REPLACE ";" "|" SIMPLE_ARPACK_BLAS_LIBRARIES_ARG "${SIMPLE_ARPACK_BLAS_LIBRARIES}")
    string(REPLACE ";" "|" SIMPLE_ARPACK_LAPACK_LIBRARIES_ARG "${SIMPLE_ARPACK_LAPACK_LIBRARIES}")

    set(SIMPLE_ARPACK_CMAKE_ARGS
        "-DCMAKE_INSTALL_PREFIX:PATH=${SIMPLE_NUMERICS_INSTALL_PREFIX}"
        "-DCMAKE_INSTALL_LIBDIR:PATH=lib"
        "-DCMAKE_PREFIX_PATH:PATH=${SIMPLE_NUMERICS_INSTALL_PREFIX}"
        "-DBUILD_SHARED_LIBS:BOOL=OFF"
        "-DEXAMPLES:BOOL=OFF"
        "-DTESTS:BOOL=OFF"
        "-DMPI:BOOL=OFF"
        "-DICB:BOOL=OFF"
        "-DINTERFACE64:BOOL=OFF"
        "-DBLAS_LIBRARIES:STRING=${SIMPLE_ARPACK_BLAS_LIBRARIES_ARG}"
        "-DLAPACK_LIBRARIES:STRING=${SIMPLE_ARPACK_LAPACK_LIBRARIES_ARG}"
        "-DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}"
        "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
    )
    if(CMAKE_BUILD_TYPE)
        list(APPEND SIMPLE_ARPACK_CMAKE_ARGS
            "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        )
    endif()

    ExternalProject_Add(simple_arpack_ep
        URL "${SIMPLE_ARPACK_NG_URL}"
        PREFIX "${CMAKE_BINARY_DIR}/_deps/arpack-ng"
        DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/_deps/downloads"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        LIST_SEPARATOR "|"
        CMAKE_ARGS ${SIMPLE_ARPACK_CMAKE_ARGS}
        BUILD_BYPRODUCTS "${SIMPLE_ARPACK_LIBRARY}"
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
    )
    if(SIMPLE_NUMERICS_USE_BUNDLED_OPENBLAS)
        add_dependencies(simple_arpack_ep simple_openblas_ep)
    endif()
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_arpack_ep)

    set(SIMPLE_ARPACK_STATUS "bundled ARPACK-NG ${SIMPLE_ARPACK_NG_VERSION}")
    set(SIMPLE_ARPACK_LINK_LIBRARY arpack)
else()
    set(SIMPLE_ARPACK_LINK_LIBRARY ${ARPACK_LIBRARY})
endif()

# Link ARPACK before its numerical providers; this matters for static builds.
if(SIMPLE_NUMERICS_BUNDLED_LINK_DIRS)
    list(APPEND SIMPLE_LIBRARIES ${SIMPLE_NUMERICS_BUNDLED_LINK_DIRS})
endif()
list(APPEND SIMPLE_LIBRARIES
    ${SIMPLE_ARPACK_LINK_LIBRARY}
    ${SIMPLE_NUMERICS_PROVIDER_LIBRARIES}
)

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

message(STATUS " ARPACK:              YES (${SIMPLE_ARPACK_STATUS})")
message(STATUS " BLAS:                YES (${SIMPLE_BLAS_STATUS})")
message(STATUS " FFTW3:               YES")
message(STATUS " JPEG_LIBRARY:        YES")
message(STATUS " LAPACK:              YES (${SIMPLE_LAPACK_STATUS})")
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
