cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(SIMPLE
    VERSION 3.0.0
    DESCRIPTION "SIMPLE cryo-EM image processing"
    HOMEPAGE_URL "https://github.com/hael/SIMPLE"
    LANGUAGES Fortran C CXX)

# CMake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Build options
option(USE_OPENMP "Use OpenMP for parallelization" ON)
option(USE_OPENACC "Use OpenACC for GPU acceleration" OFF)
option(USE_COARRAYS "Use Fortran coarrays for parallelization" OFF)
option(USE_MPI "Use MPI for parallelization" OFF)
option(USE_LIBTIFF "Use libtiff library" ON)
option(BUILD_TESTS "Build test executables" ON)
option(BUILD_DOCS "Build documentation" OFF)
option(INSTALL_GUI "Install GUI components" ON)

# Always build shared library (as requested)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared library" FORCE)

# Installation directories
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation prefix" FORCE)
endif()

# Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Git version information
set(SIMPLE_GIT_VERSION "not known")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=8 --always
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            OUTPUT_VARIABLE SIMPLE_GIT_VERSION
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set(SIMPLE_GIT_VERSION "git${SIMPLE_GIT_VERSION}")
    endif()
else()
    set(SIMPLE_GIT_VERSION "-${PROJECT_VERSION}")
endif()

message(STATUS "Git version: ${SIMPLE_GIT_VERSION}")

# Generate SimpleGitVersion.h
configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/Modules/GitVersion.h.in"
    "${CMAKE_BINARY_DIR}/lib/simple/SimpleGitVersion.h"
    @ONLY
)

# Include the generated header location
include_directories("${CMAKE_BINARY_DIR}/lib/simple")

# Compiler options and flags
include(CompilerOptions)

# Find dependencies
include(Dependencies)

# Build targets
set(SIMPLELIB "SIMPLE${PROJECT_VERSION}")

# Add source directory - builds the library
add_subdirectory(src)

# Add executables
add_subdirectory(production)

# Add scripts
add_subdirectory(scripts)

# Optional subdirectories
if(EXISTS "${CMAKE_SOURCE_DIR}/nice" AND INSTALL_GUI)
    add_subdirectory(nice)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/gui" AND INSTALL_GUI)
    add_subdirectory(gui)
endif()

if(BUILD_DOCS AND EXISTS "${CMAKE_SOURCE_DIR}/doc")
    add_subdirectory(doc)
endif()

# Testing
if(BUILD_TESTS)
    enable_testing()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()

# Generate environment setup scripts
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

# Summary
message(STATUS "========================================")
message(STATUS "SIMPLE Configuration Summary")
message(STATUS "========================================")
message(STATUS "Version:          ${PROJECT_VERSION}")
message(STATUS "Build type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix:   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "C compiler:       ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "C++ compiler:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "----------------------------------------")
message(STATUS "OpenMP:           ${USE_OPENMP}")
message(STATUS "OpenACC:          ${USE_OPENACC}")
message(STATUS "Coarrays:         ${USE_COARRAYS}")
message(STATUS "MPI:              ${USE_MPI}")
message(STATUS "LIBTIFF:          ${USE_LIBTIFF}")
message(STATUS "Build tests:      ${BUILD_TESTS}")
message(STATUS "========================================")