# - Finds OpenMP support
# This module can be used to detect OpenMP support in a compiler.
# If the compiler supports OpenMP, the flags required to compile with
# openmp support are set.
#
# This module was modified from the standard FindOpenMP module to find Fortran
# flags.
#
# The following variables are set:
#   OpenMP_Fortran_FLAGS - flags to add to the Fortran compiler for OpenMP
#                          support.  In general, you must use these at both
#                          compile- and link-time.
#   OMP_NUM_PROCS - the max number of processors available to OpenMP

#=============================================================================
# Copyright 2009 Kitware, Inc.
# Copyright 2008-2009 Andr\`e Rigland Brodtkorb <Andre.Brodtkorb@ifi.uio.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Modified by Michael Eager 2018 Monash University
#
#

INCLUDE (${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
UNSET (OPENMP_FOUND CACHE)
SET (OpenMP_Fortran_FLAG_CANDIDATES
  #Gnu
  "-fopenmp"
  #Intel
  "-qopenmp"
  #IBM XL C/c++
  "-qsmp"
  #Portland Group
  "-mp"
  #Clang
  "-fopenmp=libomp"
  #Sun
  "-xopenmp"
  #HP
  "+Oopenmp"
  #PathScale, Intel
  "-openmp"
  #Microsoft Visual Studio
  "/openmp"
  #Intel windows
  "/Qopenmp"
  #Empty, if compiler automatically accepts openmp
  " "
  )
set(OMP_FLAG_GNU "-fopenmp")

set (IFM3 $ENV{SLURM_CLUSTER_NAME})
if("${IFM3} " STREQUAL "m3 ")
  message(STATUS "Testing OpenMP on MASSIVE ")
  message(STATUS "Warning: if cuda module loaded after gcc or openmpi, LDFLAGS will override libstdc++ location. ")

  string(FIND "$ENV{LDFLAGS}" "-L/usr/lib64" OMP_ERROR )
  if(NOT ${OMP_ERROR} EQUAL -1)
    message(STATUS " Warning: /usr/lib64 in LDFLAGS ")
	  message(STATUS " Warning: LDFLAGS=  $ENV{LDFLAGS}" )
  endif()
  #  message(STATUS " ${CMAKE_C_COMPILER} ${CMAKE_Fortran_COMPILER}")
  #    find_package(OpenMP REQUIRED)
  #    if(OPENMP_FOUND)
  #    endif()
  add_definitions("-DMPI_NOF08_MODULE=1")
endif() # building on MASSIVE


if(NOT OPENMP_FOUND)
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "Intel")
    add_definitions("-DMPI_NOF08_MODULE=1")
    if("${CMAKE_${LANG}_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
      set(OMP_FLAG_Intel "-openmp")
    else()
      set(OMP_FLAG_Intel "-qopenmp")
    endif()
  endif()
  set(OMP_FLAG_PGI "-mp")
  IF (DEFINED OpenMP_Fortran_FLAGS)
    SET (OpenMP_Fortran_FLAG_CANDIDATES)
  ENDIF (DEFINED OpenMP_Fortran_FLAGS)
  set(OpenMP_Fortran_TEST_SOURCE_ORIGINAL
    "
      program test
      use omp_lib
      integer :: n
      n = omp_get_num_threads()
      end program test
  ")
  #
  # Use the !$ to force the compiler to use preprocessor directive
  set(OpenMP_Fortran_TEST_SOURCE
    "
      program TestOpenMP
!$        use omp_lib
          write(*,'(I0)',ADVANCE='NO') omp_get_num_procs()
      end
  ")

  FILE (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90
    "${OpenMP_Fortran_TEST_SOURCE}")

  # check fortran compiler. also determine number of processors
  FOREACH (FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})
    SET (SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    SET (CMAKE_REQUIRED_FLAGS "${FLAG}")
    UNSET (OpenMP_FLAG_DETECTED CACHE)
    MESSAGE (STATUS "Try OpenMP Fortran flag = [${FLAG}]")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS "-DOpenMP_FLAG_DETECTED ${CMAKE_REQUIRED_FLAGS}")
    if(APPLE)
      try_compile(OpenMP_FLAG_DETECTED  ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} -DCMAKE_LINKER_STATIC_FLAGS=""
        COMPILE_OUTPUT_VARIABLE OUTPUT
        RUN_OUTPUT_VARIABLE OMP_NUM_PROCS_INTERNAL)
    else()
	    message(STATUS "OpenMP testing : ${CMAKE_REQUIRED_FLAGS}")
      TRY_RUN (OpenMP_RUN_FAILED OpenMP_FLAG_DETECTED ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90
        COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} -DCMAKE_LINKER_STATIC_FLAGS=""
        COMPILE_OUTPUT_VARIABLE OUTPUT
        RUN_OUTPUT_VARIABLE OMP_NUM_PROCS_INTERNAL)
    endif(APPLE)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(STATUS "OpenMP compilation output : ${OUTPUT}")
      message(STATUS "OpenMP execution output   : ${OMP_NUM_PROCS_INTERNAL}")
    endif()
    IF (OpenMP_FLAG_DETECTED)
      FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Determining if the Fortran compiler supports OpenMP passed with "
        "the following output:\n${OMP_NUM_PROCS_INTERNAL}\n\n")
      SET (OpenMP_FLAG_DETECTED 1)
      # IF (OpenMP_RUN_FAILED)
      #     MESSAGE (FATAL_ERROR "OpenMP found, but test code did not run")
      # ENDIF (OpenMP_RUN_FAILED)
      STRING(REGEX MATCH "^[^0-9]*$" OMP_ERROR "${OMP_NUM_PROCS_INTERNAL}")

      if ("${OMP_ERROR} " STREQUAL " ")
        SET (OMP_NUM_PROCS ${OMP_NUM_PROCS_INTERNAL} CACHE
          STRING "Number of processors OpenMP may use" FORCE)
        SET (OpenMP_Fortran_FLAGS_INTERNAL "${FLAG}")
      else()
        message(STATUS "OMP NUM PROCS ERROR ${OMP_ERROR}" )
        STRING(REGEX REPLACE ".* \([0-9]+\)$" "\\1" OMP_NUM_PROCS_INTERNAL "${OMP_NUM_PROCS_INTERNAL}")
        message(STATUS "OMP NUM PROCS OUTPUT set to ${OMP_NUM_PROCS_INTERNAL}" )
        SET (OMP_NUM_PROCS ${OMP_NUM_PROCS_INTERNAL} CACHE
          STRING "Number of processors OpenMP may use" FORCE)
        SET (OpenMP_Fortran_FLAGS_INTERNAL "${FLAG}")
      endif()
      BREAK ()
    ELSE ()
      FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining if the Fortran compiler supports OpenMP failed with "
        "the following output:\n${OMP_NUM_PROCS_INTERNAL}\n\n")
      SET (OpenMP_FLAG_DETECTED 0)
    ENDIF (OpenMP_FLAG_DETECTED)
  ENDFOREACH (FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})
  SET (OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS_INTERNAL}"
    CACHE STRING "Fortran compiler flags for OpenMP parallelism")
endif(NOT OPENMP_FOUND)

UNSET (OpenMP_Version_DETECTED CACHE)
MESSAGE (STATUS "Try OpenMP version")


# if(APPLE)
# FILE (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMPVersion.f90
#   "
#       program TestOpenMPVersion
# !$ use omp_lib    !    only include when _OPENMP defined
#           write(*,'(I6)',ADVANCE='NO') openmp_version
#       end program TestOpenMPVersion
#   ")
#     try_compile(OpenMP_Version_DETECTED ${CMAKE_BINARY_DIR}
#       ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMPVersion.f90
#       COMPILE_DEFINITIONS ${OpenMP_Fortran_FLAGS}
#       CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} -DCMAKE_LINKER_STATIC_FLAGS="" -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
#       COMPILE_OUTPUT_VARIABLE OUTPUT
#       RUN_OUTPUT_VARIABLE OMP_VERSION_INTERNAL)
#       message (STATUS " OpenMP version compilation output: ${OUTPUT}")
#       message (STATUS " OpenMP version runtime output: ${OMP_VERSION_INTERNAL}")

#  else()
FILE (WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMPVersion.f90
  "
      program TestOpenMPVersion
          include 'omp_lib.h' !    only include when _OPENMP defined
          write(*,'(I6)',ADVANCE='NO') openmp_version
      end program TestOpenMPVersion
  ")

#	message(STATUS " TESTING  ")
TRY_RUN (OpenMP_RUN_FAILED OpenMP_Version_DETECTED ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMPVersion.f90
  COMPILE_DEFINITIONS ${OpenMP_Fortran_FLAGS}
  CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} -DCMAKE_LINKER_STATIC_FLAGS="" -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  COMPILE_OUTPUT_VARIABLE OUTPUT
  RUN_OUTPUT_VARIABLE OMP_VERSION_INTERNAL)
       message (STATUS " OpenMP version compilation flag: ${OpenMP_RUN_FAILED}")
       message (STATUS " OpenMP version runtime flag: ${OpenMP_Version_DETECTED}")

#       message (STATUS " OpenMP version compilation output: ${OUTPUT}")
       message (STATUS " OpenMP version runtime output: ${OMP_VERSION_INTERNAL}")

#endif()
IF (OpenMP_Version_DETECTED)
  FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining the Fortran compiler version of OpenMP passed with "
    "the following output:\n${OMP_VERSION_INTERNAL}\n\n")
  SET (OpenMP_Version_DETECTED 1)
  # IF (OpenMP_RUN_FAILED)
  #     MESSAGE (FATAL_ERROR "OpenMP found, but test code did not run")
  # ENDIF (OpenMP_RUN_FAILED)
message(STATUS " openmp version test output: ${OMP_VERSION_INTERNAL}")
  if(OMP_VERSION_INTERNAL STREQUAL "")
    message(FATAL_ERROR "OpenMP version unexpected error. Output empty")
  endif()
  # Check for other error mesages
  string(REGEX REPLACE "[\n\t]" "" OMP_VERSION_INTERNAL "${OMP_VERSION_INTERNAL}")
message(STATUS " openmp version test output: ${OMP_VERSION_INTERNAL}")
  STRING(REGEX MATCH "^[^0-9]" OMP_ERROR "${OMP_VERSION_INTERNAL}")
  message(STATUS " openmp version test numerical output: ${OMP_VERSION_INTERNAL}: regex ${OMP_ERROR}")
  if ("${OMP_ERROR} " STREQUAL " ")
     message(STATUS "Setting OpenMP_Fortran_Version to ${OMP_VERSION_INTERNAL}")
    SET (OpenMP_Fortran_VERSION ${OMP_VERSION_INTERNAL}  CACHE
      STRING " OpenMP version " FORCE)
  else()
    message(STATUS " OMP VERSION ERROR ${OMP_ERROR}" )
    STRING(REGEX REPLACE ".*\(2[0-9]*\)$" "\\1" OMP_VERSION_INTERNAL "${OMP_VERSION_INTERNAL}")
    message(STATUS " OMP VERSION OUTPUT set to ${OMP_VERSION_INTERNAL}" )
    SET (OpenMP_Fortran_VERSION ${OMP_VERSION_INTERNAL} CACHE
      STRING " OpenMP version in Fortran header " FORCE)
  endif()

ELSE ()
  FILE (APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining the Fortran compiler version of OpenMP failed with "
    "the following output:\n${OMP_VERSION_INTERNAL}\n\n")
  SET (OpenMP_Version_DETECTED 0)
ENDIF (OpenMP_Version_DETECTED)
UNSET (OpenMP_Version_DETECTED CACHE)

if(OpenMP_Fortran_VERSION)
  message (STATUS "OpenMP Version string cleanup")
  string(REGEX REPLACE ".*\(2[0-9]*\)$" "\\1" OpenMP_Fortran_VERSION "${OpenMP_Fortran_VERSION}")
  message (STATUS "OpenMP Version string = ${OpenMP_Fortran_VERSION}")
else()
  message(FATAL_ERROR " OpenMP version string not defined : ${OMP_VERSION_INTERNAL}")
endif(OpenMP_Fortran_VERSION)


mark_as_advanced(OpenMP_Fortran_FLAGS)
# mark_as_advanced(OpenMP_Fortran_VERSION)


# handle the standard arguments for FIND_PACKAGE
find_package_handle_standard_args (OpenMP_Fortran DEFAULT_MSG
  OpenMP_Fortran_FLAGS)
