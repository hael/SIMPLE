cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include(FortranOverride)

project(SIMPLE Fortran)
#enable_language( CXX C )

set(PROJECT_BRIEF "SIMPLE is a program package for cryo-EM image processing (https://simplecryoem.com)")

#################################################################
#
# VERSION
#
#################################################################
set(SIMPLE_VERSION_MAJOR 3)
set(SIMPLE_VERSION_MINOR 0)
set(SIMPLE_VERSION_PATCH 0)
set(SIMPLE_VERSION "${SIMPLE_VERSION_MAJOR}.${SIMPLE_VERSION_MINOR}.${SIMPLE_VERSION_PATCH}")

#################################################################
#
# OPTIONS
#
#################################################################
option(USE_OPENMP "Use OpenMP for parallelization - debugging may disable OpemMP" ON )
option(USE_OPENACC "Use OpenACC for auto GPU parallelization - debugging may disable OpenMP" OFF)
#option(USE_CUDA "GPU NVIDIA - in development" OFF)

option(${PROJECT_NAME}_BUILD_TESTS "Build the test examples" ON)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gui)
  option(${PROJECT_NAME}_BUILD_GUI "Build the simple GUI " ON)
endif()
option(BUILD_SHARED_LIBS "Build the simple library as shared or static (default)" OFF)
option(BUILD_DOCS "Enable the generation of PDF manuals and doxygen API documentation" OFF)
option(USE_SYSTEM_JPEG_SQLite "Compile libjpeg9 and sqlite3 from source (default) or use system libraries." OFF)

#################################################################
#
# SET default flags
#
#################################################################
# Add our local modules to the module path
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SCRIPTS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
set(TEST_SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests")
set(PROD_SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/production/simple")
set(EXTRA_LIBS "")
set(NUM_JOBS 1)
include(ProcessorCount)
ProcessorCount(NUM_JOBS)
set(BUILD_NAME "${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}" CACHE STRING "Build name with compiler, type and FFT backend ")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

############################################################
# Define the main targets
############################################################
# Have the .mod files placed in the lib folder
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
set(EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_INSTALL_BINDIR}")
set(LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_LIBDIR}")
set(SCRIPTS_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/scripts")

# Define the executable names
set(SIMPLE_EXE simple_exec)
set(SIMPLE_DISTR_EXE simple_distr_exec)
set(SIMPLE_PRIVATE_EXE simple_private_exec)
# Define the library name lib<name>.so
set(SIMPLELIB "SIMPLE${${PROJECT_NAME}_VERSION}")

# Define some directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/inc)
include_directories(${CMAKE_BINARY_DIR}/lib/simple)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#################################################################
# INCLUDES
#################################################################
# include(CMakeAddFortranSubdirectory)
include(CheckIncludeFiles)
include(SimpleFortranOptions)

############################################################
#
# Add definition macros unique to compilation
#
############################################################
add_definitions("-DBUILD_NAME=\"${BUILD_NAME}\"")
add_definitions(" -D__FILENAME__='\"$(notdir $<)\"' ")
# use _DEBUG since DEBUG or debug is used as a variable in some files
# _DEBUG is defined by Intel Fortran on Windows
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  add_definitions("-D_DEBUG")
endif()


#   cmake_add_fortran_subdirectory(
#    ${SRC}                  # name of subdirectory
#    PROJECT libSIMPLE       # project name in subdir top CMakeLists.txt
#    ARCHIVE_DIR ${LIB}      # dir where project places .lib files
#    RUNTIME_DIR ${LIB}      # dir where project places .dll files
#    LIBRARIES <lib>         # names of library targets to import
#    LINK_LIBRARIES          # link interface libraries for LIBRARIES
#     [LINK_LIBS ${FFTW_LIBRARIES} <dep>...]...
#    CMAKE_COMMAND_LINE ...  # extra command line flags to pass to cmake
#    NO_EXTERNAL_INSTALL     # skip installation of external project
#    )


#################################################################
#
# CREATE libSIMPLE<VERSION>
#
#################################################################

add_subdirectory(${SRC_PATH} lib)
set(EXTRA_LIBS
  ${SIMPLELIB}
  ${CMAKE_THREAD_LIBS_INIT}
  ${EXTRA_LIBS})
link_directories(lib)


#################################################################
#
# ADD SCRIPTS
#
#################################################################
# The script folder contains some handy perl scripts for running SIMPLE
add_subdirectory(${SCRIPTS_PATH})


#################################################################
#
# ENABLE GUI
#
#################################################################
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gui)
  if (${PROJECT_NAME}_BUILD_GUI)
    add_subdirectory(gui  ${CMAKE_INSTALL_BINDIR}/gui)
  endif()
endif()

#################################################################
#
# STATIC LIBRARIES
#
#################################################################
# Try to search first static libraries
if(NOT ${BUILD_SHARED_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES} .a")
endif()

#################################################################
#
# Add a distclean target to the Makefile
#
#################################################################
add_custom_target(distclean
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/distclean.cmake
  )

# here we save a version file
set(VERSION_FILE "${CMAKE_INSTALL_LIBDIR}/simple/CMake/${PROJECT_NAME}ConfigVersion.cmake")
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/ProjectConfigVersion.cmake.in" ${VERSION_FILE} @ONLY )

# modules / depends
FILE(GLOB_RECURSE SIMPLE_CMAKE_FILES
    LIST_DIRECTORIES false
    RELATIVE ${CMAKE_SOURCE_DIR} "*.cmake" )
INSTALL(
  FILES ${SIMPLE_CMAKE_FILES}
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/simple/CMake"
  CONFIGURATIONS DEBUG
  COMPONENT dev)

message( STATUS "Set runtime variables ")
if(NOT BUILD_SHARED_LIBS)
  configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/scripts/add2.bashrc.template"  ${CMAKE_CURRENT_BINARY_DIR}/add2.bashrc @ONLY )
  configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/scripts/add2.tcshrc.template"  ${CMAKE_CURRENT_BINARY_DIR}/add2.tcshrc @ONLY )
else()
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/scripts/add2.bashrc.sharedtemplate"  ${CMAKE_CURRENT_BINARY_DIR}/add2.bashrc @ONLY )
  configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/scripts/add2.tcshrc.sharedtemplate"  ${CMAKE_CURRENT_BINARY_DIR}/add2.tcshrc @ONLY )
  endif()
install(FILES  ${CMAKE_CURRENT_BINARY_DIR}/add2.bashrc ${CMAKE_CURRENT_BINARY_DIR}/add2.tcshrc
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  CONFIGURATIONS DEBUG|RELEASE
  PERMISSIONS   OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ WORLD_WRITE
  )

  if(${BUILD_DOCS} OR VERBOSE)

    #
    # User manual and tutorial
    #
    add_subdirectory(doc)

    #
    # Doxygen
    #
    option(ENABLE_DOXYGEN "Enable a 'make docs' target for Doxygen documentation" ON)
    if (ENABLE_DOXYGEN)
      add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
      include(UseDoxygen OPTIONAL)
      # add_doxygen(
      #     docs
      #     INPUT
      #         "${CMAKE_SOURCE_DIR}/src"
      #     STRIP_FROM_PATH
      #         "${CMAKE_SOURCE_DIR}/src/inc"
      #     STRIP_FROM_INC_PATH
      #         "${CMAKE_SOURCE_DIR}/src/inc"
      #     EXCLUDE_PATTERNS
      #         */gnufor2/*
      #     EXCLUDE_SYMBOLS
      #         *debug*
      #         *verbose*
      #     ALL
      # )
    endif()
  endif()


OPTION(USE_CODE_COVERAGE "Added code coverage flags" OFF)
OPTION(USE_PROFILING "enable symbol generation without debug support" OFF)

option(USE_LINK_TIME_OPTIMISATION "Enable link-time optimisation" OFF)
if(VERBOSE) #only for the fool-hardy
  OPTION(USE_FAST_MATH_OPTIMISATION "Optimize code for fastmath. " OFF)
  OPTION(USE_PROFILE_OPTIMISATION "Optimize code for PFO. (GNU: -fpfo )" OFF)
  OPTION(USE_AUTO_PARALLELISE "Optimize code for parallel systems. [Intel and PGI only]" OFF)
  OPTION(USE_AGGRESSIVE_OPTIMISATION "Flags for aggressive optimization level. Unroll all loops and inline functions."  OFF)
  option(CPACK_GENERATOR "Enable build of distribution packages using CPack" TGZ)

  #
  # CTest
  #
  enable_testing()
  #
  # CPack
  #
  option(CPACK_GENERATOR "Enable build of distribution packages using CPack" OFF)
  include(SimplePackage)
  ## Debugging CMAKE
  include(CMakeGraphVizOptions)
  include(GenericLogging)
endif()



#################################################################
#
# ADD PRODUCTION EXECUTABLES
#
#################################################################
# The source for the SIMPLE binary and have it placed in the bin folder
add_subdirectory(production bin)


# Add message at end of installation
add_subdirectory(${CMAKE_SOURCE_DIR}/cmake/PostInstall)
