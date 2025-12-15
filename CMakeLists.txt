cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

# Make sure CMake can see our helper modules early
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Auto-detect GCC/GFortran BEFORE project() if user didn’t set compilers
include(AutoDetectGCC)

# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------
project(SIMPLE
    VERSION 3.0.0
    DESCRIPTION "SIMPLE cryo-EM image processing"
    HOMEPAGE_URL "https://github.com/hael/SIMPLE"
    LANGUAGES Fortran C CXX)

# Export compile_commands.json for tooling
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Build options
# ------------------------------------------------------------------------------
# OpenMP is effectively required (always ON by default)
option(USE_OPENMP   "Use OpenMP for parallelization" ON)

# Future / optional parallel models
option(USE_OPENACC  "Use OpenACC for GPU acceleration (experimental)" OFF)
option(USE_COARRAYS "Use Fortran coarrays for parallelization (experimental)" OFF)
option(USE_MPI      "Use MPI for parallelization" OFF)

option(USE_LIBTIFF  "Use libtiff library" ON)

option(BUILD_TESTS  "Build test executables" ON)
option(BUILD_DOCS   "Build documentation" OFF)
option(INSTALL_GUI  "Install GUI components" ON)

# Always build shared libraries
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared library" FORCE)

# ------------------------------------------------------------------------------
# Install layout & default build type
# ------------------------------------------------------------------------------
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Local install by default (easy testing / packaging)
    set(CMAKE_INSTALL_PREFIX
        "${CMAKE_BINARY_DIR}/install"
        CACHE PATH "Installation prefix" FORCE)
endif()
# Reasonable default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# ------------------------------------------------------------------------------
# Fortran modules & generated headers
# ------------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
file(MAKE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")
# Where various generated files go
set(SIMPLE_BINARY_LIB_DIR "${CMAKE_BINARY_DIR}/lib/simple")
file(MAKE_DIRECTORY "${SIMPLE_BINARY_LIB_DIR}")
include_directories("${SIMPLE_BINARY_LIB_DIR}")
# Git version information (best-effort)
set(SIMPLE_GIT_VERSION "not known")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    find_package(Git QUIET)
    if(Git_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" describe --abbrev=8 --always
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            OUTPUT_VARIABLE SIMPLE_GIT_VERSION
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set(SIMPLE_GIT_VERSION "git${SIMPLE_GIT_VERSION}")
    else()
        set(SIMPLE_GIT_VERSION "-${PROJECT_VERSION}")
    endif()
else()
    set(SIMPLE_GIT_VERSION "-${PROJECT_VERSION}")
endif()
message(STATUS "Git version: ${SIMPLE_GIT_VERSION}")
# Generate SimpleGitVersion.h
configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/GitVersion.h.in"
    "${SIMPLE_BINARY_LIB_DIR}/SimpleGitVersion.h"
    @ONLY
)

# ------------------------------------------------------------------------------
# Compiler configuration (GCC/GFortran 14+, flags, RPATH, etc.)
# ------------------------------------------------------------------------------
include(CompilerConfig)

# ------------------------------------------------------------------------------
# Dependencies (OpenMP, FFTW, TIFF, MPI, etc.)
#   -> populates SIMPLE_LIBRARIES
# ------------------------------------------------------------------------------
include(Dependencies)
# Library target name used in subdirs
set(SIMPLELIB "SIMPLE${PROJECT_VERSION}")

# ------------------------------------------------------------------------------
# Subdirectories
# ------------------------------------------------------------------------------
add_subdirectory(src)         # builds shared library ${SIMPLELIB}
add_subdirectory(production)  # executables + tests
add_subdirectory(scripts)
if(EXISTS "${CMAKE_SOURCE_DIR}/nice" AND INSTALL_GUI)
    add_subdirectory(nice)
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/gui" AND INSTALL_GUI)
    add_subdirectory(gui)
endif()
if(BUILD_DOCS AND EXISTS "${CMAKE_SOURCE_DIR}/doc")
    add_subdirectory(doc)
endif()

# ------------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTS)
    enable_testing()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()

# ------------------------------------------------------------------------------
# Environment setup scripts
# ------------------------------------------------------------------------------
configure_file(
    "${CMAKE_SOURCE_DIR}/scripts/add2.bashrc.sharedtemplate"
    "${CMAKE_BINARY_DIR}/add2.bashrc"
    @ONLY
)
configure_file(
    "${CMAKE_SOURCE_DIR}/scripts/add2.tcshrc.sharedtemplate"
    "${CMAKE_BINARY_DIR}/add2.tcshrc"
    @ONLY
)
install(FILES
    "${CMAKE_BINARY_DIR}/add2.bashrc"
    "${CMAKE_BINARY_DIR}/add2.tcshrc"
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
message(STATUS "========================================")
message(STATUS "SIMPLE Configuration Summary")
message(STATUS "========================================")
message(STATUS "Version:              ${PROJECT_VERSION}")
message(STATUS "Git version:          ${SIMPLE_GIT_VERSION}")
message(STATUS "Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Fortran compiler:     ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
message(STATUS "C compiler:           ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
message(STATUS "C++ compiler:         ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "----------------------------------------")
message(STATUS "OpenMP:               ${USE_OPENMP}")
message(STATUS "OpenACC:              ${USE_OPENACC}")
message(STATUS "Coarrays:             ${USE_COARRAYS}")
message(STATUS "MPI:                  ${USE_MPI}")
message(STATUS "LIBTIFF:              ${USE_LIBTIFF}")
message(STATUS "Build tests:          ${BUILD_TESTS}")
message(STATUS "========================================")

# Add message at end of installation
add_subdirectory(${CMAKE_SOURCE_DIR}/cmake/PostInstall)

