cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
include("${CMAKE_SOURCE_DIR}/cmake/Modules/FortranOverride.cmake")
project(SIMPLE Fortran)
#enable_language( CXX C )
include(${CMAKE_SOURCE_DIR}/cmake/Modules/GenericLogging.cmake)

#################################################################
# VERSION 
#################################################################

set(${PROJECT_NAME}_VERSION 2.5)
set(${PROJECT_NAME}_SOVERSION 2.5)

#################################################################
# CONFIGURATION TYPES & BUILD MODE
#################################################################
set(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE RELWITHDEBINFO )
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: NONE DEBUG RELEASE RELWITHDEBINFO"
    FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
endif(NOT CMAKE_BUILD_TYPE)

#################################################################
# OPTIONS
#################################################################
# Set some options the user may choose
option(USE_OPENMP "Use OpenMP for parallelization - debugging may disable OpemMP" ON)
option(USE_OPENACC "Use OpenACC for auto GPU parallelization - debugging may disable OpemMP" ON)
option(USE_CUDA "GPU NVIDIA - in development" OFF)
option(${PROJECT_NAME}_BUILD_TESTS "Build the test examples" ON)
option(${PROJECT_NAME}_BUILD_GUI "Build the simple GUI - not included in release2.5" OFF)
option(BUILD_SHARED_LIBS "Build the simple library as shared (default) or static" ON)
option(ENABLE_LINK_TIME_OPTIMISATION "Enable link-time optimisation" OFF)
OPTION(ENABLE_FAST_MATH_OPTIMISATION "Optimize code for fastmath." OFF)
OPTION(ENABLE_PROFILE_OPTIMISATION "Optimize code for PFO." OFF)
OPTION(ENABLE_AUTO_PARALLELISE "Optimize code for parallel systems." OFF)
#flags for optimization level. You can only turn one of following option to ON, or leave all of them to OFF.
OPTION(ENABLE_PROFILING "enable symbol generation without debug support" OFF)
option(CPACK_GENERATOR "Enable build of distribution packages using CPack" OFF)
#OPTION(ENABLE_OPTIMIZE_MACHINE "Optimize code for this specific computer (typ 5-10% speedup). Binaries may incompatible with some other computers." OFF)
#OPTION(ENABLE_OPTIMIZE_PHI "Optimize code for the Intel compiler on the PHI" OFF)
#OPTION(ENABLE_OPTIMIZE_WINDOWS_VC "Windows visual c++ specific optimization" OFF)
if (APPLE)
IF(CMAKE_SYSTEM MATCHES "Darwin-12*")
        OPTION(ENABLE_OPTIMIZE_MOUNTAIN_LION "enable OSX 10.8/9 support" ON)
ELSE(CMAKE_SYSTEM MATCHES "Darwin-12*")
        OPTION(ENABLE_OPTIMIZE_MOUNTAIN_LION "enable OSX 10.8/9 support" OFF)
ENDIF(CMAKE_SYSTEM MATCHES "Darwin-12*")
endif()

#################################################################
# SET default flags
#################################################################
include(GNUInstallDirs)
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
# Add our local modules to the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(LIB_SRC_PATH ${SRC_PATH}/lib)
set(TEST_SRC_PATH ${SRC_PATH}/test)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
# Have the exec files placed in the bin folder
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# Have the .mod files placed in the lib folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(EXTRA_LIBS "")
set(NUM_JOBS 1)
set(BUILD_NAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE} CACHE STRING "Build name with compiler, type and FFT backend ")

#################################################################
# INCLUDES
#################################################################
## Debugging CMAKE
# include(CMakeAddFortranSubdirectory)
# include(CMakeGraphVizOptions)
include(CheckIncludeFiles)
include(${CMAKE_MODULE_PATH}/SimpleFortranOptions.cmake)


############################################################
# Add definition macros to compilation
############################################################
add_definitions(-DBUILD_NAME=${BUILD_NAME})
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  add_definitions("-D_DEBUG") # use _DEBUG since DEBUG or debug is used as a variable in some files
                              # _DEBUG is defined by Intel Fortran on Windows 
endif()

############################################################
# Define the main targets
############################################################
# Define the executable names
set(SIMPLE_EXE simple_exec)
set(SIMPLE_DISTR_EXE simple_distr_exec)
# Define the library name lib<name>.so
set(SIMPLELIB "SIMPLE${${PROJECT_NAME}_VERSION}")

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(LIB ${CMAKE_BINARY_DIR}/lib)
set(BIN ${CMAKE_BINARY_DIR}/bin)
set(EXESRC ${CMAKE_SOURCE_DIR}/production/simple)
include_directories(${CMAKE_SOURCE_DIR}/simple_utils)
include_directories(${CMAKE_SOURCE_DIR}/src/inc)


set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

#   cmake_add_fortran_subdirectory(
#    ${SRC}                  # name of subdirectory
#    PROJECT libSIMPLE       # project name in subdir top CMakeLists.txt
#    ARCHIVE_DIR ${LIB}      # dir where project places .lib files
#    RUNTIME_DIR ${LIB}      # dir where project places .dll files
#    LIBRARIES <lib>              # names of library targets to import
#    LINK_LIBRARIES          # link interface libraries for LIBRARIES
#     [LINK_LIBS ${FFTW_LIBRARIES} <dep>...]...
#    CMAKE_COMMAND_LINE ...  # extra command line flags to pass to cmake
#    NO_EXTERNAL_INSTALL     # skip installation of external project
#    )


add_subdirectory(${SRC} ${LIB})
set(EXTRA_LIBS
  ${SIMPLELIB}
  ${CMAKE_THREAD_LIBS_INIT}
  ${EXTRA_LIBS})
# link_directories(${LIB})

#################################################################
# The source for the SIMPLE binary and have it placed in the bin folder
#################################################################
add_subdirectory(${EXESRC} ${BIN})

if (${PROJECT_NAME}_BUILD_TESTS)
  add_subdirectory(tests)
endif()

#if (${PROJECT_NAME}_BUILD_GUI)
#  add_subdirectory(gui ${BIN}/gui)
#endif()


#################################################################
# ENABLE TESTING
#################################################################
if(${PROJECT_NAME}_ENABLE_TESTS)
  enable_testing()
  include(CTest)
endif()

#################################################################
# STATIC LIBRARIES
#################################################################
# Try to search first static libraries
if(NOT ${BUILD_SHARED_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
endif()

#################################################################
# Add a distclean target to the Makefile
#################################################################
add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake
)
include(CPack)
include(${CMAKE_MODULE_PATH}/GenericLogging.cmake )
