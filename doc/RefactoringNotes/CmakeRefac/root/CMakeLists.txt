cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(SIMPLE
  VERSION 3.0.0
  DESCRIPTION "SIMPLE cryo-EM image processing"
  HOMEPAGE_URL "https://github.com/hael/SIMPLE"
  LANGUAGES Fortran C CXX
)

# ----------------------------------------------------------------------
# Safety
# ----------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not supported. Use cmake -S . -B build")
endif()

# ----------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------
option(USE_OPENMP     "Use OpenMP for parallelization" ON)
option(USE_OPENACC    "Use OpenACC for GPU acceleration" OFF)
option(USE_COARRAYS   "Use Fortran coarrays" OFF)
option(USE_MPI        "Use MPI for parallelization" OFF)
option(USE_LIBTIFF    "Use libtiff library" ON)
option(BUILD_TESTS    "Build tests" ON)
option(BUILD_DOCS     "Build documentation" OFF)
option(INSTALL_GUI    "Install GUI components" ON)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# ----------------------------------------------------------------------
# Install layout
# ----------------------------------------------------------------------
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      "${CMAKE_BINARY_DIR}/install"
      CACHE PATH "Installation prefix" FORCE)
endif()

# ----------------------------------------------------------------------
# CMake modules
# ----------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ----------------------------------------------------------------------
# Compiler + dependencies
# ----------------------------------------------------------------------
include(CompilerOptions)
include(Dependencies)

# ----------------------------------------------------------------------
# Subdirectories (targets live here)
# ----------------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(production)
add_subdirectory(scripts)

if(INSTALL_GUI AND EXISTS "${CMAKE_SOURCE_DIR}/nice")
  add_subdirectory(nice)
endif()

if(INSTALL_GUI AND EXISTS "${CMAKE_SOURCE_DIR}/gui")
  add_subdirectory(gui)
endif()

if(BUILD_DOCS AND EXISTS "${CMAKE_SOURCE_DIR}/doc")
  add_subdirectory(doc)
endif()

# ----------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------
if(BUILD_TESTS)
  enable_testing()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()

# ----------------------------------------------------------------------
# Environment setup scripts
# ----------------------------------------------------------------------
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
  DESTINATION .
)

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
message(STATUS "========================================")
message(STATUS "SIMPLE Configuration Summary")
message(STATUS "Version:        ${PROJECT_VERSION}")
message(STATUS "Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Fortran:        ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "C:              ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "C++:            ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "----------------------------------------")
message(STATUS "OpenMP:         ${USE_OPENMP}")
message(STATUS "OpenACC:        ${USE_OPENACC}")
message(STATUS "Coarrays:       ${USE_COARRAYS}")
message(STATUS "MPI:            ${USE_MPI}")
message(STATUS "LIBTIFF:        ${USE_LIBTIFF}")
message(STATUS "Tests:          ${BUILD_TESTS}")
message(STATUS "========================================")
