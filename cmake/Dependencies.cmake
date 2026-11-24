# cmake/Dependencies.cmake
#
# Finds and configures external dependencies.
# Populates SIMPLE_LIBRARIES for linking by the SIMPLE library and executables.


set(SIMPLE_LIBRARIES "")
set(SIMPLE_DEPENDENCY_TARGETS "")

# Download a pinned archive into the persistent cache once, or verify the cached copy
# (OpenBLAS, gtk-fortran).
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
# BLAS / LAPACK (required): SIMPLE's pinned OpenBLAS
# ------------------------------------------------------------------------------
# SIMPLE uses the standard LP64 Fortran ABI for these libraries. If you point
# CMake at an ILP64 vendor build, integer sizes will not match the current code.
#
# The build never probes or links a system BLAS or LAPACK: OpenBLAS provides
# both, is built once per toolchain into a cache outside build/ and is linked
# by full path, so an installed library of the same name cannot win the
# linker's search.

set(SIMPLE_OPENBLAS_VERSION "0.3.34")
set(SIMPLE_OPENBLAS_SHA256 "cd7e129868320cc2d033afa920e31202dfe0b8066a5b66661900ccc0f197dfed")

# NUM_THREADS sizes the library independently of the machine that builds it:
# its buffer allocator has 2*NUM_THREADS slots and every thread inside a BLAS
# call holds one (the OpenMP build also reserves one per runtime thread for its
# serial-region products); past the table OpenBLAS warns, corrupts its table
# and crashes, it does not wait. It also caps the threads of a serial-region
# product. Any node with up to this many cores is covered.
set(SIMPLE_OPENBLAS_NUM_THREADS 128)

# The OpenMP build threads serial-region products on SIMPLE's own OpenMP pool
# and stays on the calling thread inside SIMPLE's teams. A C compiler that
# cannot share the Fortran OpenMP runtime gets the sequential build with the
# locked allocator (USE_LOCKING: safe when called from several threads), which
# also stays on the calling thread; its serial-region products run on one core.
if(USE_OPENMP AND CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(SIMPLE_OPENBLAS_FLAVOUR "OpenMP")
    set(SIMPLE_OPENBLAS_THREAD_FLAGS "USE_THREAD=1" "USE_OPENMP=1")
else()
    set(SIMPLE_OPENBLAS_FLAVOUR "sequential")
    set(SIMPLE_OPENBLAS_THREAD_FLAGS "USE_THREAD=0" "USE_OPENMP=0" "USE_LOCKING=1")
endif()

set(SIMPLE_OPENBLAS_BUILD_FLAGS
    ${SIMPLE_OPENBLAS_THREAD_FLAGS}
    "NUM_THREADS=${SIMPLE_OPENBLAS_NUM_THREADS}"
    # one buffer set for threaded calls: SIMPLE never threads two products at once
    "NUM_PARALLEL=1"
    "NO_SHARED=1"
    "NO_LAPACK=0"
    "NO_LAPACKE=1"
    "INTERFACE64=0"
)
string(JOIN " " SIMPLE_OPENBLAS_RECIPE
    "OpenBLAS=${SIMPLE_OPENBLAS_VERSION}"
    "scratch=build-local-v1"
    "system=${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
    "C=${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}-${CMAKE_C_COMPILER}"
    "Fortran=${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER}"
    ${SIMPLE_OPENBLAS_BUILD_FLAGS}
)
string(MAKE_C_IDENTIFIER
    "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}-${SIMPLE_OPENBLAS_FLAVOUR}"
    SIMPLE_NUMERICS_CACHE_KEY
)
set(SIMPLE_NUMERICS_CACHE_ROOT "${CMAKE_SOURCE_DIR}/.deps/numerics" CACHE PATH
    "Persistent cache for SIMPLE's pinned OpenBLAS build")
set(SIMPLE_NUMERICS_INSTALL_PREFIX
    "${SIMPLE_NUMERICS_CACHE_ROOT}/${SIMPLE_NUMERICS_CACHE_KEY}/openblas-${SIMPLE_OPENBLAS_VERSION}"
)
set(SIMPLE_NUMERICS_RECIPE_STAMP "${SIMPLE_NUMERICS_INSTALL_PREFIX}/simple-numerics-recipe.txt")
set(SIMPLE_OPENBLAS_LIBRARY
    "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
set(SIMPLE_OPENBLAS_EP_PREFIX "${CMAKE_BINARY_DIR}/_deps/openblas-bundle")
set(SIMPLE_BLAS_STATUS
    "bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} (${SIMPLE_OPENBLAS_FLAVOUR}, NUM_THREADS=${SIMPLE_OPENBLAS_NUM_THREADS})"
)

message(STATUS
    "SIMPLE numerical provider: bundled OpenBLAS ${SIMPLE_OPENBLAS_VERSION} "
    "(persistent cache: ${SIMPLE_NUMERICS_INSTALL_PREFIX})"
)

set(SIMPLE_NUMERICS_CACHE_CURRENT FALSE)
if(EXISTS "${SIMPLE_OPENBLAS_LIBRARY}" AND EXISTS "${SIMPLE_NUMERICS_RECIPE_STAMP}")
    file(READ "${SIMPLE_NUMERICS_RECIPE_STAMP}" SIMPLE_NUMERICS_CACHED_RECIPE)
    string(STRIP "${SIMPLE_NUMERICS_CACHED_RECIPE}" SIMPLE_NUMERICS_CACHED_RECIPE)
    if(SIMPLE_NUMERICS_CACHED_RECIPE STREQUAL SIMPLE_OPENBLAS_RECIPE)
        set(SIMPLE_NUMERICS_CACHE_CURRENT TRUE)
    endif()
endif()

if(SIMPLE_NUMERICS_CACHE_CURRENT)
    message(STATUS "Reusing the cached OpenBLAS; no OpenBLAS build is scheduled")
else()
    if(EXISTS "${SIMPLE_NUMERICS_INSTALL_PREFIX}")
        message(STATUS "Cached OpenBLAS recipe changed or is incomplete; rebuilding ${SIMPLE_NUMERICS_INSTALL_PREFIX}")
        file(REMOVE_RECURSE "${SIMPLE_NUMERICS_INSTALL_PREFIX}")
    endif()
    # A recipe change can occur inside an existing build tree. Discard only
    # the generated dependency scratch so stale ExternalProject completion
    # stamps cannot suppress the required rebuild.
    file(REMOVE_RECURSE "${SIMPLE_OPENBLAS_EP_PREFIX}")
    find_program(SIMPLE_MAKE_PROGRAM NAMES gmake make)
    if(NOT SIMPLE_MAKE_PROGRAM)
        message(FATAL_ERROR
            "Building bundled OpenBLAS requires GNU make or make."
        )
    endif()

    set(SIMPLE_NUMERICS_DOWNLOAD_DIR "${CMAKE_SOURCE_DIR}/.deps/downloads")
    set(SIMPLE_OPENBLAS_ARCHIVE "${SIMPLE_NUMERICS_DOWNLOAD_DIR}/OpenBLAS-${SIMPLE_OPENBLAS_VERSION}.tar.gz")
    file(MAKE_DIRECTORY
        "${SIMPLE_NUMERICS_DOWNLOAD_DIR}"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/lib"
        "${SIMPLE_NUMERICS_INSTALL_PREFIX}/include"
    )

    simple_fetch_pinned_archive(
        "OpenBLAS ${SIMPLE_OPENBLAS_VERSION}"
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${SIMPLE_OPENBLAS_VERSION}/OpenBLAS-${SIMPLE_OPENBLAS_VERSION}.tar.gz"
        "${SIMPLE_OPENBLAS_ARCHIVE}"
        "${SIMPLE_OPENBLAS_SHA256}"
    )
    file(WRITE "${SIMPLE_NUMERICS_RECIPE_STAMP}" "${SIMPLE_OPENBLAS_RECIPE}\n")
    include(ExternalProject)

    set(SIMPLE_OPENBLAS_MAKE_ARGS
        "CC=${CMAKE_C_COMPILER}"
        "FC=${CMAKE_Fortran_COMPILER}"
        ${SIMPLE_OPENBLAS_BUILD_FLAGS}
    )
    # OpenBLAS's compiler probes use mktemp and silently lose HAVE_C11 when
    # an inherited TMPDIR is missing or unwritable. Own the scratch directory
    # inside this build tree and override temp paths only in the child processes.
    set(SIMPLE_OPENBLAS_TMP_DIR "${SIMPLE_OPENBLAS_EP_PREFIX}/scratch")
    set(SIMPLE_OPENBLAS_MAKE_COMMAND
        "${CMAKE_COMMAND}" -E env
        "TMPDIR=${SIMPLE_OPENBLAS_TMP_DIR}"
        "TMP=${SIMPLE_OPENBLAS_TMP_DIR}"
        "TEMP=${SIMPLE_OPENBLAS_TMP_DIR}"
        "${SIMPLE_MAKE_PROGRAM}" ${SIMPLE_OPENBLAS_MAKE_ARGS}
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
        # Targets run sequentially as they update the same static archive,
        # parallel make invocation may corrupt it.
        BUILD_COMMAND "${CMAKE_COMMAND}" -E make_directory "${SIMPLE_OPENBLAS_TMP_DIR}"
            COMMAND ${SIMPLE_OPENBLAS_MAKE_COMMAND} libs
            COMMAND ${SIMPLE_OPENBLAS_MAKE_COMMAND} netlib
        INSTALL_COMMAND "${CMAKE_COMMAND}" -E make_directory "${SIMPLE_OPENBLAS_TMP_DIR}"
            COMMAND ${SIMPLE_OPENBLAS_MAKE_COMMAND} "PREFIX=${SIMPLE_NUMERICS_INSTALL_PREFIX}" install
        BUILD_BYPRODUCTS "${SIMPLE_OPENBLAS_LIBRARY}"
        LOG_DOWNLOAD TRUE
        LOG_BUILD TRUE
        LOG_INSTALL TRUE
        LOG_OUTPUT_ON_FAILURE TRUE
    )
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_openblas_ep)
endif()

# Both flavours stay on the calling thread inside an OpenMP team: code may
# call the BLAS there without a product of its own.
add_compile_definitions(SIMPLE_BLAS_TEAM_SAFE)
set(SIMPLE_BLAS_TEAM_STATUS "ON (${SIMPLE_OPENBLAS_FLAVOUR} OpenBLAS stays on the calling thread)")

list(APPEND SIMPLE_LIBRARIES "${SIMPLE_OPENBLAS_LIBRARY}")
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
# GTK 3 and gtk-fortran (required)
# ------------------------------------------------------------------------------
# GTK 3 comes from the system through pkg-config. gtk-fortran's GTK 3 release is pinned and built
# like OpenBLAS: once per Fortran compiler into a cache outside build/, because
# gfortran module files do not load across compiler versions. SIMPLE_GTK_LIBRARIES
# is for the executables that run GTK code; it is never added to SIMPLE_LIBRARIES,
# so the other executables do not link GTK.
find_package(PkgConfig)
if(NOT PkgConfig_FOUND)
    message(FATAL_ERROR
        "================================================================================\n"
        "  PKG-CONFIG REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE finds GTK 3 through pkg-config.\n"
        "  Install pkg-config via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install pkg-config\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install pkgconf-pkg-config\n"
        "    or: sudo dnf install pkgconf-pkg-config\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S pkgconf\n"
        "  macOS (Homebrew):\n"
        "    brew install pkgconf\n"
        "  macOS (MacPorts):\n"
        "    sudo port install pkgconfig\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
endif()
pkg_check_modules(GTK3 IMPORTED_TARGET gtk+-3.0)
if(NOT GTK3_FOUND)
    message(FATAL_ERROR
        "================================================================================\n"
        "  GTK 3 REQUIRED but NOT FOUND\n"
        "================================================================================\n"
        "  SIMPLE requires the GTK 3 development files to compile;\n"
        "  pkg-config (${PKG_CONFIG_EXECUTABLE}) found no gtk+-3.0 package.\n"
        "  Install them via your system package manager:\n"
        "  Linux (Debian/Ubuntu):\n"
        "    sudo apt-get update\n"
        "    sudo apt-get install libgtk-3-dev\n"
        "  Linux (Fedora/RHEL/CentOS/Rocky):\n"
        "    sudo yum install gtk3-devel\n"
        "    or: sudo dnf install gtk3-devel\n"
        "  Linux (Arch):\n"
        "    sudo pacman -S gtk3\n"
        "  macOS (Homebrew):\n"
        "    brew install gtk+3\n"
        "  macOS (MacPorts):\n"
        "    sudo port install gtk3\n"
        "  GTK 3 in another prefix: add <prefix>/lib/pkgconfig to PKG_CONFIG_PATH.\n"
        "  After installation, reconfigure CMake:\n"
        "    rm -rf build/\n"
        "    cmake -B build\n"
        "================================================================================\n"
    )
endif()

set(SIMPLE_GTK_FORTRAN_VERSION "3.24.49")
set(SIMPLE_GTK_FORTRAN_SHA256 "e8a00fb1b144f56d2d3c84669a1a70256c779ee263206fdc853631823c4f5a91")
# the core interface modules only: no PLplot, high-level API, sketcher or examples
set(SIMPLE_GTK_FORTRAN_BUILD_FLAGS
    "-DCMAKE_BUILD_TYPE=Release"
    "-DCMAKE_Fortran_FLAGS_RELEASE=-O2"
    "-DEXCLUDE_PLPLOT=true"
    "-DNO_BUILD_HL=true"
    "-DNO_BUILD_EXAMPLES=true"
)
# the one target of gtk-fortran's build that SIMPLE links; nothing else of it is built
set(SIMPLE_GTK_FORTRAN_TARGET "gtk-fortran_static")
string(JOIN " " SIMPLE_GTK_FORTRAN_RECIPE
    "gtk-fortran=${SIMPLE_GTK_FORTRAN_VERSION}"
    "target=${SIMPLE_GTK_FORTRAN_TARGET}"
    "system=${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
    "Fortran=${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER}"
    ${SIMPLE_GTK_FORTRAN_BUILD_FLAGS}
)
string(MAKE_C_IDENTIFIER
    "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
    SIMPLE_GTK_CACHE_KEY
)
set(SIMPLE_GTK_CACHE_ROOT "${CMAKE_SOURCE_DIR}/.deps/gtk" CACHE PATH
    "Persistent cache for SIMPLE's pinned gtk-fortran build")
set(SIMPLE_GTK_FORTRAN_PREFIX
    "${SIMPLE_GTK_CACHE_ROOT}/${SIMPLE_GTK_CACHE_KEY}/gtk-fortran-${SIMPLE_GTK_FORTRAN_VERSION}"
)
set(SIMPLE_GTK_FORTRAN_RECIPE_STAMP "${SIMPLE_GTK_FORTRAN_PREFIX}/simple-gtk-fortran-recipe.txt")
set(SIMPLE_GTK_FORTRAN_MODULE_DIR "${SIMPLE_GTK_FORTRAN_PREFIX}/include/gtk-3-fortran")
set(SIMPLE_GTK_FORTRAN_LIBRARY
    "${SIMPLE_GTK_FORTRAN_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtk-3-fortran${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
set(SIMPLE_GTK_FORTRAN_EP_PREFIX "${CMAKE_BINARY_DIR}/_deps/gtk-fortran")

message(STATUS
    "GTK: GTK ${GTK3_VERSION} (pkg-config), pinned gtk-fortran ${SIMPLE_GTK_FORTRAN_VERSION} "
    "(persistent cache: ${SIMPLE_GTK_FORTRAN_PREFIX})"
)

set(SIMPLE_GTK_CACHE_CURRENT FALSE)
if(EXISTS "${SIMPLE_GTK_FORTRAN_LIBRARY}" AND EXISTS "${SIMPLE_GTK_FORTRAN_RECIPE_STAMP}")
    file(READ "${SIMPLE_GTK_FORTRAN_RECIPE_STAMP}" SIMPLE_GTK_CACHED_RECIPE)
    string(STRIP "${SIMPLE_GTK_CACHED_RECIPE}" SIMPLE_GTK_CACHED_RECIPE)
    if(SIMPLE_GTK_CACHED_RECIPE STREQUAL SIMPLE_GTK_FORTRAN_RECIPE)
        set(SIMPLE_GTK_CACHE_CURRENT TRUE)
    endif()
endif()

if(SIMPLE_GTK_CACHE_CURRENT)
    message(STATUS "Reusing the cached gtk-fortran; no gtk-fortran build is scheduled")
else()
    if(EXISTS "${SIMPLE_GTK_FORTRAN_PREFIX}")
        message(STATUS "Cached gtk-fortran recipe changed or is incomplete; rebuilding ${SIMPLE_GTK_FORTRAN_PREFIX}")
        file(REMOVE_RECURSE "${SIMPLE_GTK_FORTRAN_PREFIX}")
    endif()
    # stale ExternalProject stamps in an existing build tree must not suppress the rebuild
    file(REMOVE_RECURSE "${SIMPLE_GTK_FORTRAN_EP_PREFIX}")
    set(SIMPLE_GTK_FORTRAN_ARCHIVE
        "${CMAKE_SOURCE_DIR}/.deps/downloads/gtk-fortran-${SIMPLE_GTK_FORTRAN_VERSION}.tar.gz"
    )
    file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/.deps/downloads" "${SIMPLE_GTK_FORTRAN_PREFIX}")
    simple_fetch_pinned_archive(
        "gtk-fortran ${SIMPLE_GTK_FORTRAN_VERSION}"
        "https://github.com/vmagnin/gtk-fortran/archive/refs/tags/v${SIMPLE_GTK_FORTRAN_VERSION}.tar.gz"
        "${SIMPLE_GTK_FORTRAN_ARCHIVE}"
        "${SIMPLE_GTK_FORTRAN_SHA256}"
    )
    file(WRITE "${SIMPLE_GTK_FORTRAN_RECIPE_STAMP}" "${SIMPLE_GTK_FORTRAN_RECIPE}\n")
    include(ExternalProject)
    # gtk-fortran's own CMake build with SIMPLE's Fortran compiler, limited to the static library:
    # its other products (the shared library, an informational program that needs GLib 2.64, newer
    # than EL8's) are never built. The library and its module files are copied into the cache.
    ExternalProject_Add(simple_gtk_fortran_ep
        URL "${SIMPLE_GTK_FORTRAN_ARCHIVE}"
        URL_HASH "SHA256=${SIMPLE_GTK_FORTRAN_SHA256}"
        PREFIX "${SIMPLE_GTK_FORTRAN_EP_PREFIX}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CMAKE_ARGS
            "-DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}"
            "-DPKG_CONFIG_EXECUTABLE:FILEPATH=${PKG_CONFIG_EXECUTABLE}"
            ${SIMPLE_GTK_FORTRAN_BUILD_FLAGS}
        BUILD_COMMAND "${CMAKE_COMMAND}" --build <BINARY_DIR> --target ${SIMPLE_GTK_FORTRAN_TARGET}
        INSTALL_COMMAND "${CMAKE_COMMAND}" -E copy_directory
                <BINARY_DIR>/src/modules "${SIMPLE_GTK_FORTRAN_MODULE_DIR}"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${SIMPLE_GTK_FORTRAN_PREFIX}/lib"
            COMMAND "${CMAKE_COMMAND}" -E copy
                "<BINARY_DIR>/src/${CMAKE_STATIC_LIBRARY_PREFIX}gtk-3-fortran${CMAKE_STATIC_LIBRARY_SUFFIX}"
                "${SIMPLE_GTK_FORTRAN_LIBRARY}"
        BUILD_BYPRODUCTS "${SIMPLE_GTK_FORTRAN_LIBRARY}"
        LOG_CONFIGURE TRUE
        LOG_BUILD TRUE
        LOG_INSTALL TRUE
        LOG_OUTPUT_ON_FAILURE TRUE
    )
    list(APPEND SIMPLE_DEPENDENCY_TARGETS simple_gtk_fortran_ep)
endif()
set(SIMPLE_GTK_LIBRARIES "${SIMPLE_GTK_FORTRAN_LIBRARY}" PkgConfig::GTK3)

# ------------------------------------------------------------------------------
# Final dependency summary
# ------------------------------------------------------------------------------
message(STATUS "=== SIMPLE Dependency Summary ===")

message(STATUS " BLAS:                YES (${SIMPLE_BLAS_STATUS})")
message(STATUS " BLAS in OMP teams:   ${SIMPLE_BLAS_TEAM_STATUS}")
message(STATUS " FFTW3:               YES")
message(STATUS " GTK 3 / gtk-fortran: YES (GTK ${GTK3_VERSION}, gtk-fortran ${SIMPLE_GTK_FORTRAN_VERSION})")
message(STATUS " JPEG_LIBRARY:        YES")
message(STATUS " LAPACK:              YES (${SIMPLE_BLAS_STATUS})")
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
