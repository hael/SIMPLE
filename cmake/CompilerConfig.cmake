# cmake/CompilerConfig.cmake
#
# Compiler and platform configuration for SIMPLE.
# Requires: CMake >= 3.25, GCC/GFortran >= 14 (Linux), GFortran + AppleClang (macOS).
# ------------------------------------------------------------------------------
# Basic compiler sanity checks
# ------------------------------------------------------------------------------
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "Only GFortran is supported for Fortran. Found: ${CMAKE_Fortran_COMPILER_ID}")
endif()
if(APPLE)
    # On macOS, allow AppleClang for C/C++ if GFortran is used
    if(NOT (CMAKE_C_COMPILER_ID MATCHES "Clang|AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "GNU"))
        message(FATAL_ERROR "C compiler must be GCC or AppleClang on macOS. Found: ${CMAKE_C_COMPILER_ID}")
    endif()
    if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
        message(FATAL_ERROR "C++ compiler must be G++ or AppleClang on macOS. Found: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    # On Linux, enforce GCC/G++
    if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Only GCC is supported for C on Linux. Found: ${CMAKE_C_COMPILER_ID}")
    endif()

    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "Only G++ is supported for C++ on Linux. Found: ${CMAKE_CXX_COMPILER_ID}")
    endif()
endif()

# Require modern GFortran
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14.0")
    message(FATAL_ERROR "GFortran 14.0 or higher required. Found: ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

# ------------------------------------------------------------------------------
# Language standards
# ------------------------------------------------------------------------------
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Platform-specific tuning
# ------------------------------------------------------------------------------
if(APPLE)
    execute_process(
        COMMAND sysctl -n machdep.cpu.brand_string
        OUTPUT_VARIABLE APPLE_CPU
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(APPLE_CPU MATCHES "Apple M")
        set(ARCH_FLAG "-mtune=generic")
    else()
        set(ARCH_FLAG "-mtune=native")
    endif()
else()
    set(ARCH_FLAG "-march=native")
endif()

# ------------------------------------------------------------------------------
# Global preprocessor definitions
# ------------------------------------------------------------------------------
if(APPLE)
    add_compile_definitions(__FreeBSD__ MACOSX)
else()
    add_compile_definitions(__linux__ LINUX)
endif()
add_compile_definitions(
    GNU       # your historic define
    USE_F08=1 # Fortran 2008 support
    $<$<CONFIG:Debug>:_DEBUG>
)

# ------------------------------------------------------------------------------
# Fortran compile flags
#   Use target_compile_options where possible, but keep defaults for all targets.
# ------------------------------------------------------------------------------
# Common Fortran flags (apply to all configurations)
string(APPEND CMAKE_Fortran_FLAGS
       " -cpp -ffree-form -fimplicit-none -ffree-line-length-none"
       " -fno-second-underscore -Wall -Waliasing -Wampersand"
       " -Wsurprising -Wline-truncation"
       " -D__FILENAME__='\"\\$(notdir \\$<)\"'")

# Release flags for Fortran
# for checking vectorization:
#-O3 -march=native -fopenmp -fopt-info-vec-optimized -fopt-info-vec-missed
# inspect vectorization reports: -qopt-report, -fopt-info-vec
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops ${ARCH_FLAG} -fPIC" 
    CACHE STRING "Release flags for Fortran" FORCE)

# Debug flags for Fortran
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fbounds-check -fcheck=all -Wuninitialized -Wunused -fPIC"
    CACHE STRING "Debug flags for Fortran" FORCE)

# ------------------------------------------------------------------------------
# C / C++ compile flags
# ------------------------------------------------------------------------------
set(CMAKE_C_FLAGS_RELEASE "-O3 ${ARCH_FLAG} -fPIC"
    CACHE STRING "Release flags for C" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fPIC"
    CACHE STRING "Debug flags for C" FORCE)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${ARCH_FLAG} -fPIC"
    CACHE STRING "Release flags for C++" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fPIC"
    CACHE STRING "Debug flags for C++" FORCE)
# Suppress Clang deployment version warning on macOS
if(APPLE AND CMAKE_C_COMPILER_ID MATCHES "Clang|AppleClang")
    string(APPEND CMAKE_C_FLAGS " -Wno-overriding-deployment-version")
    string(APPEND CMAKE_CXX_FLAGS " -Wno-overriding-deployment-version")
endif()

# ------------------------------------------------------------------------------
# RPATH for installed shared libraries
# ------------------------------------------------------------------------------
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib;${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
