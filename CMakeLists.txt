cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
include("${CMAKE_SOURCE_DIR}/cmake/Modules/FortranOverride.cmake")
project(SIMPLE LANGUAGES Fortran)

#################################################################
# VERSION and the Shared lib SOVERSION
#################################################################

set(${PROJECT_NAME}_VERSION 2.5)
set(${PROJECT_NAME}_SOVERSION 2.5)

#################################################################
# DEFINE PATHS
#################################################################
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
# Add our local modules to the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#################################################################
# BUILD PATHS
#################################################################
set(SRC_PATH ${CMAKE_SOURCE_DIR}/src)
set(LIB_SRC_PATH ${SRC_PATH}/lib)
set(TEST_SRC_PATH ${SRC_PATH}/test)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

#################################################################
# CONFIGURATION TYPES & BUILD MODE
#################################################################
set(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: NONE DEBUG RELEASE"
    FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
endif(NOT CMAKE_BUILD_TYPE)

#################################################################
# OPTIONS
#################################################################

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
# OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
option(USE_OPENMP "Use OpenMP for parallelization" ON)
option(USE_CUDA "GPU parallelisation" OFF)
option(BUILD_TESTS "Build the test examples" ON)
option(BUILD_GUI "Build the simple GUI" OFF)
option(BUILD_SHARED_LIBS "Build the simple library as shared" ON)

option(LINK_TIME_OPTIMISATION "Enable link-time optimisation" OFF)

#################################################################
# SET default flags
#################################################################
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(EXTRA_LIBS "")
set(NUM_JOBS 1)
include(ProcessorCount)
ProcessorCount(NUM_JOBS)

#################################################################
# INCLUDES
#################################################################
## Debugging CMAKE
# include(CMakeAddFortranSubdirectory)
# include(CMakeGraphVizOptions)
include(CheckIncludeFiles)
include("${CMAKE_SOURCE_DIR}/cmake/Modules/FortranOverride.cmake")
include(${CMAKE_MODULE_PATH}/SimpleFortranOptions.cmake)
#INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
include(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
get_filename_component(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
if(FCNAME STREQUAL "pgf90")
    unset(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
endif(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
set(SIMPLE_EXE simple_exec)
set(SIMPLE_DISTR_EXE simple_distr_exec)

# Define the library name lib<name>.so
set(SIMPLELIB "SIMPLE${${PROJECT_NAME}_VERSION}")

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(LIB ${CMAKE_BINARY_DIR}/lib)
set(BIN ${CMAKE_BINARY_DIR}/bin)
set(EXESRC ${CMAKE_SOURCE_DIR}/production/simple)
include_directories(${CMAKE_SOURCE_DIR}/src/inc)

# Have the .mod files placed in the lib folder
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
if (BUILD_TESTS)
  add_subdirectory(tests)
endif()
#if (BUILD_GUI)
#  add_subdirectory(gui ${BIN}/gui)
#endif()


#################################################################
# ENABLE TESTING
#################################################################
#set(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL} CACHE STRING "Testing Intel MKL" )
#if(${PROJECT_NAME}_ENABLE_TESTS)
#  enable_testing()
#  include(CTest)
#endif()

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

 include(${CMAKE_MODULE_PATH}/GenericLogging.cmake )
