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
            # REQUIRED already hard-fails when the toolkit is absent; the previous inner
            # conditional tested CUDAToolkit_Fortran_FOUND, which FindCUDAToolkit never sets,
            # so its FATAL_ERROR branch fired only by accident of an undefined variable.
            find_package(CUDAToolkit REQUIRED)
            list(APPEND SIMPLE_LIBRARIES
                CUDA::cudart
                CUDA::cufft
                CUDA::cublas)
            string(APPEND CMAKE_Fortran_FLAGS " -foffload=nvptx-none")
            add_compile_definitions(USE_OPENMP_OFFLOAD)
        endif()
    endif()
else()
    set(USE_OPENMP_OFFLOAD OFF)
endif()

# ------------------------------------------------------------------------------
# Coarrays (optional, multi-image coarrays via GFortran)
# ------------------------------------------------------------------------------
if(USE_COARRAYS)
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Coarrays currently only wired for GFortran.")
    endif()
    add_compile_definitions(USE_COARRAYS)
    find_package(Coarray)
    if(NOT Coarray_FOUND)
        message(FATAL_ERROR "USE_COARRAYS=ON requires Coarray support (cafrun, caf_mpi/libcaf_mpi, and MPI Fortran).")
    endif()
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

set(SIMPLE_OPENBLAS_VERSION "0.3.34")
set(SIMPLE_OPENBLAS_SHA256 "cd7e129868320cc2d033afa920e31202dfe0b8066a5b66661900ccc0f197dfed")
set(SIMPLE_ARPACK_NG_VERSION "3.9.1")
set(SIMPLE_ARPACK_NG_SHA256 "f6641deb07fa69165b7815de9008af3ea47eb39b2bb97521fbf74c97aba6e844")
set(SIMPLE_OPENBLAS_NUM_PARALLEL 4)

# The OpenMP build of OpenBLAS threads serial-region products and stays on the
# calling thread inside SIMPLE's OpenMP teams. Toolchains that cannot share the
# same OpenMP runtime use the pthreads flavour and SIMPLE's internal team path.
if(USE_OPENMP AND CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(SIMPLE_OPENBLAS_USE_OPENMP 1)
    set(SIMPLE_OPENBLAS_FLAVOUR "OpenMP")
else()
    set(SIMPLE_OPENBLAS_USE_OPENMP 0)
    set(SIMPLE_OPENBLAS_FLAVOUR "pthreads")
endif()

set(SIMPLE_OPENBLAS_BUILD_FLAGS
    "USE_THREAD=1"
    "USE_OPENMP=${SIMPLE_OPENBLAS_USE_OPENMP}"
    # Four simultaneous calculation calls cover the decisive GitHub team. In
    # OpenBLAS 0.3.34 callers above this cap wait for a buffer set; keeping the
    # cap finite avoids the NUM_PARALLEL x runtime-thread buffer footprint.
    "NUM_PARALLEL=${SIMPLE_OPENBLAS_NUM_PARALLEL}"
    "NO_SHARED=1"
    "NO_LAPACK=0"
    "NO_LAPACKE=1"
    "INTERFACE64=0"
)
string(JOIN " " SIMPLE_OPENBLAS_RECIPE
    "OpenBLAS=${SIMPLE_OPENBLAS_VERSION}"
    "system=${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
    "C=${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}-${CMAKE_C_COMPILER}"
    "Fortran=${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER}"
    ${SIMPLE_OPENBLAS_BUILD_FLAGS}
    "ARPACK-NG=${SIMPLE_ARPACK_NG_VERSION}"
    "ARPACK_BUILD_TYPE=Release"
    "ARPACK_MPI=OFF"
    "ARPACK_INTERFACE64=OFF"
)
string(MAKE_C_IDENTIFIER
    "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}-${SIMPLE_OPENBLAS_FLAVOUR}"
    SIMPLE_NUMERICS_CACHE_KEY
)
set(SIMPLE_NUMERICS_CACHE_ROOT "${CMAKE_SOURCE_DIR}/.deps/numerics" CACHE PATH
    "Persistent cache for SIMPLE's pinned OpenBLAS and ARPACK-NG builds")
set(SIMPLE_NUMERICS_INSTALL_PREFIX
    "${SIMPLE_NUMERICS_CACHE_ROOT}/${SIMPLE_NUMERICS_CACHE_KEY}/openblas-${SIMPLE_OPENBLAS_VERSION}-arpack-${SIMPLE_ARPACK_NG_VERSION}"
)
set(SIMPLE_NUMERICS_RECIPE_STAMP "${SIMPLE_NUMERICS_INSTALL_PREFIX}/simple-numerics-recipe.txt")
set(SIMPLE_OPENBLAS_LIBRARY
    "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
set(SIMPLE_ARPACK_LIBRARY
    "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}arpack${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
set(SIMPLE_OPENBLAS_EP_PREFIX "${CMAKE_BINARY_DIR}/_deps/openblas-bundle")
set(SIMPLE_ARPACK_EP_PREFIX "${CMAKE_BINARY_DIR}/_deps/arpack-ng-bundle")

set(SIMPLE_BLAS_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (${SIMPLE_OPENBLAS_FLAVOUR})")
set(SIMPLE_LAPACK_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (${SIMPLE_OPENBLAS_FLAVOUR})")
set(SIMPLE_ARPACK_STATUS "bundled ARPACK-NG ${SIMPLE_ARPACK_NG_VERSION}")

# SIMPLE's numerical provider is deliberately invariant across platforms. The
# build never probes or links a system BLAS, LAPACK, or ARPACK: ARPACK-NG is
# built against this same OpenBLAS so it cannot pull another BLAS into the
# process through a system dependency.

message(STATUS
    "SIMPLE numerical provider: bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} "
    "and ARPACK-NG ${SIMPLE_ARPACK_NG_VERSION} (persistent cache: ${SIMPLE_NUMERICS_INSTALL_PREFIX})"
)

set(SIMPLE_NUMERICS_CACHE_CURRENT FALSE)
if(EXISTS "${SIMPLE_OPENBLAS_LIBRARY}" AND EXISTS "${SIMPLE_ARPACK_LIBRARY}" AND
   EXISTS "${SIMPLE_NUMERICS_RECIPE_STAMP}")
    file(READ "${SIMPLE_NUMERICS_RECIPE_STAMP}" SIMPLE_NUMERICS_CACHED_RECIPE)
    string(STRIP "${SIMPLE_NUMERICS_CACHED_RECIPE}" SIMPLE_NUMERICS_CACHED_RECIPE)
    if(SIMPLE_NUMERICS_CACHED_RECIPE STREQUAL SIMPLE_OPENBLAS_RECIPE)
        set(SIMPLE_NUMERICS_CACHE_CURRENT TRUE)
    endif()
endif()

if(SIMPLE_NUMERICS_CACHE_CURRENT)
    message(STATUS "Reusing cached SIMPLE numerical libraries; no OpenBLAS or ARPACK-NG build is scheduled")
else()
    if(EXISTS "${SIMPLE_NUMERICS_INSTALL_PREFIX}")
        message(STATUS "Cached SIMPLE numerical recipe changed or is incomplete; rebuilding ${SIMPLE_NUMERICS_INSTALL_PREFIX}")
        file(REMOVE_RECURSE "${SIMPLE_NUMERICS_INSTALL_PREFIX}")
    endif()
    # A recipe change can occur inside an existing build tree. Discard only
    # the generated dependency scratch so stale ExternalProject completion
    # stamps cannot suppress the required rebuild.
    file(REMOVE_RECURSE "${SIMPLE_OPENBLAS_EP_PREFIX}" "${SIMPLE_ARPACK_EP_PREFIX}")
    find_program(SIMPLE_MAKE_PROGRAM NAMES gmake make)
    if(NOT SIMPLE_MAKE_PROGRAM)
        message(FATAL_ERROR
            "Building bundled OpenBLAS requires GNU make or make."
        )
    endif()

    set(SIMPLE_NUMERICS_DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/.deps/downloads")
    set(SIMPLE_OPENBLAS_ARCHIVE "${SIMPLE_NUMERICS_DOWNLOAD_DIR}/OpenBLAS-${SIMPLE_OPENBLAS_VERSION}.tar.gz")
    set(SIMPLE_ARPACK_ARCHIVE "${SIMPLE_NUMERICS_DOWNLOAD_DIR}/arpack-ng-${SIMPLE_ARPACK_NG_VERSION}.tar.gz")
    file(MAKE_DIRECTORY
        "${SIMPLE_NUMERICS_DOWNLOAD_DIR}"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/include"
    )

    function(simple_fetch_pinned_archive label url archive sha256)
        if(NOT EXISTS "${archive}")
            message(STATUS "Fetching ${label}: ${url}")
            file(DOWNLOAD "${url}" "${archive}"
                EXPECTED_HASH "SHA256=${sha256}"
                STATUS _download_status
                SHOW_PROGRESS
                TLS_VERIFY ON
            )
            list(GET _download_status 0 _download_error)
            if(NOT _download_error EQUAL 0)
                file(REMOVE "${archive}")
                message(FATAL_ERROR "Could not download ${label}: ${_download_status}")
            endif()
        else()
            file(SHA256 "${archive}" _archive_hash)
            if(NOT _archive_hash STREQUAL sha256)
                message(FATAL_ERROR
                    "Cached ${label} archive has SHA-256 ${_archive_hash}; expected ${sha256}: ${archive}"
                )
            endif()
        endif()
    endfunction()

    simple_fetch_pinned_archive(
        "OpenBLAS ${SIMPLE_OPENBLAS_VERSION}"
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${SIMPLE_OPENBLAS_VERSION}/OpenBLAS-${SIMPLE_OPENBLAS_VERSION}.tar.gz"
        "${SIMPLE_OPENBLAS_ARCHIVE}"
        "${SIMPLE_OPENBLAS_SHA256}"
    )
    simple_fetch_pinned_archive(
        "ARPACK-NG ${SIMPLE_ARPACK_NG_VERSION}"
        "https://github.com/opencollab/arpack-ng/archive/refs/tags/${SIMPLE_ARPACK_NG_VERSION}.tar.gz"
        "${SIMPLE_ARPACK_ARCHIVE}"
        "${SIMPLE_ARPACK_NG_SHA256}"
    )
    file(WRITE "${SIMPLE_NUMERICS_RECIPE_STAMP}" "${SIMPLE_OPENBLAS_RECIPE}\n")
    include(ExternalProject)

    set(SIMPLE_OPENBLAS_MAKE_ARGS
        "CC=${CMAKE_C_COMPILER}"
        "FC=${CMAKE_Fortran_COMPILER}"
        ${SIMPLE_OPENBLAS_BUILD_FLAGS}
    )

    # Scratch is disposable; the installed library and verified archive persist
    # outside build/ and are reused by clean Debug and Release builds.
    ExternalProject_Add(simple_openblas_ep
        URL "${SIMPLE_OPENBLAS_ARCHIVE}"
        URL_HASH "SHA256=${SIMPLE_OPENBLAS_SHA256}"
        PREFIX "${SIMPLE_OPENBLAS_EP_PREFIX}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE TRUE
        BUILD_COMMAND "${SIMPLE_MAKE_PROGRAM}" ${SIMPLE_OPENBLAS_MAKE_ARGS} libs netlib
        INSTALL_COMMAND "${SIMPLE_MAKE_PROGRAM}" ${SIMPLE_OPENBLAS_MAKE_ARGS}
            "PREFIX=${SIMPLE_NUMERICS_INSTALL_PREFIX}" install
        BUILD_BYPRODUCTS "${SIMPLE_OPENBLAS_LIBRARY}"
        LOG_DOWNLOAD TRUE
        LOG_BUILD TRUE
        LOG_INSTALL TRUE
        LOG_OUTPUT_ON_FAILURE TRUE
    )
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_openblas_ep)

    set(SIMPLE_BLAS_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (${SIMPLE_OPENBLAS_FLAVOUR})")
    set(SIMPLE_LAPACK_STATUS "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (${SIMPLE_OPENBLAS_FLAVOUR})")
    set(SIMPLE_ARPACK_BLAS_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
    set(SIMPLE_ARPACK_LAPACK_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
    if(SIMPLE_OPENBLAS_USE_OPENMP)
        # the static OpenMP OpenBLAS needs the OpenMP runtime wherever it is linked
        if(OpenMP_Fortran_LIBRARIES)
            list(APPEND SIMPLE_ARPACK_BLAS_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
            list(APPEND SIMPLE_ARPACK_LAPACK_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
        else()
            list(APPEND SIMPLE_ARPACK_BLAS_LIBRARIES "-fopenmp")
            list(APPEND SIMPLE_ARPACK_LAPACK_LIBRARIES "-fopenmp")
        endif()
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
        "-DCMAKE_BUILD_TYPE:STRING=Release"
    )

    ExternalProject_Add(simple_arpack_ep
        URL "${SIMPLE_ARPACK_ARCHIVE}"
        URL_HASH "SHA256=${SIMPLE_ARPACK_NG_SHA256}"
        PREFIX "${SIMPLE_ARPACK_EP_PREFIX}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        LIST_SEPARATOR "|"
        CMAKE_ARGS ${SIMPLE_ARPACK_CMAKE_ARGS}
        BUILD_BYPRODUCTS "${SIMPLE_ARPACK_LIBRARY}"
        LOG_DOWNLOAD TRUE
        LOG_CONFIGURE TRUE
        LOG_BUILD TRUE
        LOG_INSTALL TRUE
        LOG_OUTPUT_ON_FAILURE TRUE
    )
    add_dependencies(simple_arpack_ep simple_openblas_ep)
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_arpack_ep)

    # Use the archive by full path so an installed libarpack cannot win the
    # linker's name search and reintroduce a system BLAS transitively.
    set(SIMPLE_ARPACK_LINK_LIBRARY "${SIMPLE_ARPACK_LIBRARY}")
endif()

# Use both archives by full path. This prevents an installed library with the
# same name from winning the linker's search and reintroducing a system BLAS.
set(SIMPLE_NUMERICS_PROVIDER_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
set(SIMPLE_ARPACK_LINK_LIBRARY "${SIMPLE_ARPACK_LIBRARY}")

# The bundled OpenMP flavour stays on the calling thread inside OpenMP teams;
# the pthreads flavour does not carry that contract.
if(SIMPLE_OPENBLAS_USE_OPENMP)
    add_compile_definitions(SIMPLE_BLAS_TEAM_SAFE)
    set(SIMPLE_BLAS_TEAM_SAFE ON)
    set(SIMPLE_BLAS_TEAM_STATUS
        "ON (OpenMP OpenBLAS stays on the calling thread; NUM_PARALLEL=${SIMPLE_OPENBLAS_NUM_PARALLEL})"
    )
else()
    set(SIMPLE_BLAS_TEAM_SAFE OFF)
    set(SIMPLE_BLAS_TEAM_STATUS
        "OFF (pthreads OpenBLAS; SIMPLE uses its internal team products)"
    )
endif()

# Link ARPACK before its numerical providers; this matters for static builds.
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
message(STATUS " BLAS in OMP teams:   ${SIMPLE_BLAS_TEAM_STATUS}")
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
