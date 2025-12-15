# FindFFTW.cmake
# Finds the FFTW3 library
#
# This will define:
#   FFTW_FOUND - System has FFTW
#   FFTW_INCLUDE_DIRS - The FFTW include directories
#   FFTW_LIBRARIES - The libraries needed to use FFTW
#   FFTW::FFTW - Imported target for FFTW

# Try pkg-config first
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_FFTW QUIET fftw3)
endif()
# Find the include directory
find_path(FFTW_INCLUDE_DIR
    NAMES fftw3.h
    HINTS
        ${PC_FFTW_INCLUDEDIR}
        ${PC_FFTW_INCLUDE_DIRS}
        ${FFTW_ROOT}
        $ENV{FFTW_ROOT}
        $ENV{FFTW_DIR}
        $ENV{FFTWDIR}
    PATH_SUFFIXES include
    PATHS
        /usr/local
        /usr
        /opt/local
        /opt/homebrew
        /opt
)
# Find the library
find_library(FFTW_LIBRARY
    NAMES fftw3 libfftw3
    HINTS
        ${PC_FFTW_LIBDIR}
        ${PC_FFTW_LIBRARY_DIRS}
        ${FFTW_ROOT}
        $ENV{FFTW_ROOT}
        $ENV{FFTW_DIR}
        $ENV{FFTWDIR}
    PATH_SUFFIXES lib lib64
    PATHS
        /usr/local
        /usr
        /opt/local
        /opt/homebrew
        /opt
)
# Handle REQUIRED and QUIET arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_LIBRARY FFTW_INCLUDE_DIR
    VERSION_VAR FFTW_VERSION
)
if(FFTW_FOUND)
    set(FFTW_LIBRARIES ${FFTW_LIBRARY})
    set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
    # Create imported target
    if(NOT TARGET FFTW::FFTW)
        add_library(FFTW::FFTW UNKNOWN IMPORTED)
        set_target_properties(FFTW::FFTW PROPERTIES
            IMPORTED_LOCATION "${FFTW_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
        )
    endif()
    mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)
endif()