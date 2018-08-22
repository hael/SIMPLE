# This overrides the default CMake Debug and Release compiler options.
# The user can still specify different options by setting the
# CMAKE_Fortran_FLAGS_[RELEASE,DEBUG] variables (on the command line or in the
# CMakeList.txt). This files serves as better CMake defaults and should only be
# modified if the default values are to be changed. Project specific compiler
# flags should be set in the CMakeList.txt by setting the CMAKE_Fortran_FLAGS_*
# variables.

# Override CMakeDetermineFortranCompiler default fortran of f95
if(NOT $ENV{FC} STREQUAL "")
  set(CMAKE_Fortran_COMPILER_NAMES $ENV{FC})
else()
  set(CMAKE_Fortran_COMPILER_NAMES gfortran mpif90 mpifort)
  # set(ENV{FC} "gfortran")
endif()
# Override preprocessor in CMakeDetermineCompiler default cpp

if(NOT $ENV{CPP} STREQUAL "")
  set(CMAKE_CPP_COMPILER_NAMES $ENV{CPP})
else()
  find_file (CMAKE_CPP_COMPILER_NAMES
    NAMES cpp- cpp-6 cpp6 cpp-5 cpp5 cpp-4.9 cpp
    PATHS /sw/bin /usr/local/bin /opt/local/bin /usr/bin
    #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
    DOC "GNU cpp preprocessor "
    #  [NO_DEFAULT_PATH]
    #  [NO_CMAKE_ENVIRONMENT_PATH]
    #  [NO_CMAKE_PATH]
    # NO_SYSTEM_ENVIRONMENT_PATH
    #  [NO_CMAKE_SYSTEM_PATH]
    #  [CMAKE_FIND_ROOT_PATH_BOTH |
    #   ONLY_CMAKE_FIND_ROOT_PATH |
    #   NO_CMAKE_FIND_ROOT_PATH]
    )
  if(NOT EXISTS ${CMAKE_CPP_COMPILER_NAMES})
    set(CMAKE_CPP_COMPILER_NAMES cpp-5)
  endif()
  set(ENV{CPP} ${CMAKE_CPP_COMPILER_NAMES})
endif()

enable_language(Fortran C CXX)
include(CMakeDetermineFortranCompiler)
include(CMakeDetermineCompiler)
# include(CMakeFortranInformation)

# Try to identify the ABI and configure it into CMakeFortranCompiler.cmake
include(${CMAKE_ROOT}/Modules/CMakeDetermineCompilerABI.cmake)
CMAKE_DETERMINE_COMPILER_ABI(Fortran ${CMAKE_ROOT}/Modules/CMakeFortranCompilerABI.F)

# Test Fortran 2008 support by using a F2008 specific construct.
if(NOT DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_F08)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008")
  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerF08.f90 "
! BLOCK should be rejected without F2008.
PROGRAM main
  IMPLICIT NONE
  BLOCK
    INTEGER :: i
  END BLOCK
END PROGRAM
")
  try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_F08 ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerF08.f90
    OUTPUT_VARIABLE OUTPUT)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F08)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 -- yes")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the Fortran compiler supports Fortran 2008 passed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_F08 1)
  else()
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 -- no")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determining if the Fortran compiler supports Fortran 2008 failed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_F08 0)
  endif()
  unset(CMAKE_Fortran_COMPILER_SUPPORTS_F08 CACHE)
endif()
# Test Fortran 2008 support by using a F2008 specific construct.
if(NOT DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_USC4)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 Unicode 10646")
  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranUnicodeSupport.f90 "
    program test_iso_10646_support
        use iso_fortran_env ,only: output_unit, error_unit
        implicit none
        integer, parameter :: UCS4_K = selected_char_kind('ISO_10646')
        if ( UCS4_K == -1 ) then
            write(error_unit,'(A)') 'Your compiler does not support ISO 10646/UCS4 characters!'
            write(error_unit,'(A)') 'JSON-Fortran must/will be configured to use the DEFAULT'
            write(error_unit,'(A)') 'character set. (Should be ASCII on a reasonable system.)'
            stop 2
        else
            write(error_unit,'(A)') 'Congratulations! Your compiler supports ISO 10646/UCS4!'
            write(error_unit,'(A)') 'JSON-Fortran may be configured to enable UCS4 support.'
            write(output_unit,'(A)') 'UCS4_SUPPORTED'
        end if
    end program test_iso_10646_support
")
  try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_USC4 ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranUnicodeSupport.f90
    OUTPUT_VARIABLE OUTPUT)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_USC4)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Unicode 10646 -- yes")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the Fortran compiler supports Fortran Unicode 10646 passed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_USC4 1)
  else()
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Unicode 10646 -- no")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determining if the Fortran compiler supports Fortran 2008 failed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_USC4 0)
  endif()
  unset(CMAKE_Fortran_COMPILER_SUPPORTS_USC4 CACHE)
endif()
# Test Fortran 2008 support for compiler_version and compiler_options.
if(NOT DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran 2008 iso_fortran_env compiler_version and compiler_options")
  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortran_F08_ISOENVSupport.f90 "
    program test_iso_support
        use iso_fortran_env
        implicit none
    character(len=:),allocatable :: str
!character(len=:),external :: compiler_version, compiler_options
    str=compiler_version()
        if ( .not. allocated(str) ) then
            write(error_unit,'(A)') 'Your compiler does not support ISO_FORTRAN_ENV compiler_version()'
            stop 2
        endif
    if (len (str) < 1) then
           write(error_unit,'(A)') 'ISO_FORTRAN_ENV compiler_version() returned string has length zero'
        stop 2
    endif
    write(error_unit,'(A)') 'Your compiler supports ISO_FORTRAN_ENV compiler_version() '
    write(error_unit,'(A)') ' Compiler version: '// trim(str)
    deallocate(str)
    str=compiler_options()
        if ( .not. allocated(str) ) then
            write(error_unit,'(A)') 'Your compiler does not support ISO_FORTRAN_ENV compiler_options()'
            stop 2
        endif
    if (len (str) < 1) then
            write(error_unit,'(A)') 'ISO_FORTRAN_ENV compiler_options() returned string has length zero'
        stop 2
    endif
    write(error_unit,'(A)') 'Your compiler supports ISO_FORTRAN_ENV compiler_options() '
    write(error_unit,'(A)') ' Compiler options: '// trim(str)
end program test_iso_support
")
  try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortran_F08_ISOENVSupport.f90
    OUTPUT_VARIABLE OUTPUT)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports ISO_FORTRAN_ENV compiler_version/options -- yes")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the Fortran compiler supports Fortran2008 ISO_FORTRAN_ENV compiler_version passed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV 1)
  else()
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports ISO_FORTRAN_ENV compiler_version/options -- no")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determining if the Fortran compiler supports Fortran2008 ISO_FORTRAN_ENV compiler_version failed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV 0)
  endif()
  unset(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV CACHE)
endif()



# Test Fortran preprocessor support for variadic macros
if(NOT DEFINED CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports variadic macros")
  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCPPCompilerVariadic.f90 "
#define c99_count(...)    _c99_count1 ( , ##__VA_ARGS__)/* */
#define _c99_count1(...)  _c99_count2 (__VA_ARGS__,10,9,8,7,6,5,4,3,2,1,0)
#define _c99_count2(_,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,n,...) n
program dummyprog
 integer i
 integer,parameter :: nv=c99_count (__VA_ARGS__);
 character(255)::p_tokens= #__VA_ARGS__ ;
 i = 5
end program dummyprog
")
  try_compile(CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCPPCompilerVariadic.f90
    OUTPUT_VARIABLE OUTPUT)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports variadic macros -- yes")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the Fortran compiler supports variadic macros passed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC 1)
  else()
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports Fortran variadic macros -- no")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determining if the Fortran compiler supports variadic macros failed with "
      "the following output:\n${OUTPUT}\n\n")
    set(CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC 0)
  endif()
  unset(CMAKE_Fortran_COMPILER_SUPPORTS_VARIADIC CACHE)
endif()





# If user specifies the build type, use theirs, otherwise use release
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "")
endif()

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create separate directory for build files.
")
endif()


# Look at system to see what if any options are available because
# of build environment
#include(SystemDefines)

# Turn on all compiler warnings
#include(EnableAllWarnings)

# Bring in helper functions for dealing with CACHE INTERNAL variables
include(CacheInternalHelpers)

# We want to create static libraries
set(BUILD_SHARED_LIBS OFF)
include(GNUInstallDirs)
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
endif()
# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
  SET(GNUNATIVE "-mtune=native")
ELSE()
  SET(GNUNATIVE "-march=native")
ENDIF(APPLE)

# include(CheckCXXCompilerFlag)

# check_cxx_compiler_flag(-std=gnu++11  HAS_FLAG_STD_GNUCXX11)

# if(NOT HAS_FLAG_STD_GNUCXX11)
#     check_cxx_compiler_flag(-std=c++11    HAS_FLAG_STD_CXX11)
# endif()

# if(HAS_FLAG_STD_GNUCXX11)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
# endif()

# if(HAS_FLAG_STD_CXX11)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# endif()


option(USE_GNU_EXTENSIONS "Enable GNU extensions in C " OFF)
if(USE_GNU_EXTENSIONS)
  set(C_DIALECT gnu)
else()
  set(C_DIALECT c)
endif()

###########  SETTING UP PROCESSING FLAGS ################
#include(PlatformDefines)

message( STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gfortran
  set(preproc  "-cpp  -Wp,C,CC,no-endif-labels")                                                # preprocessor flags
  set(dialect  "-ffree-form  -fimplicit-none  -ffree-line-length-none  -fno-second-underscore")  # language style
  set(warn     "-Waliasing -Wampersand  -Wsurprising -Wline-truncation -Wreal-q-constant ")
  #-Wconversion, -Wc-binding-type, -Wintrinsics-std, -Wtabs, -Wintrinsic-shadow -Wtarget-lifetime
  #-Wcompare-reals, -Wunused-parameter and -Wdo-subscript
  set(checks   "-frange-check -fstack-protector -fstack-check -fbounds-check ")                # checks
  set(forspeed "-O3 -ffrontend-optimize -fno-stack-protector -fno-lifetime-dse")                # optimisation
  set(forpar   "" )# -fopenmp  -Wp,-fopenmp")                                                   # parallel flags
  set(target   "${GNUNATIVE} -fPIC ")                                                           # target platform
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F08 EQUAL 1)
    set(target "${target} -std=f2008 -fall-intrinsics -Wintrinsics-std")
  else()
    set(target "${target} -std=f2003 -fall-intrinsics -Wintrinsics-std")
  endif()

  set(common   "${preproc} ${dialect} ${target} ${warn}")

  set(warnDebug "-Wall -Wcompare-reals -Wtype-limits -Wuninitialized -Wunused-but-set-parameter -Wimplicit-interface -Wno-unused-dummy-argument -Wno-unused-value -Wno-surprising")     # extra warning flag
  set(fordebug "-Og -g -pedantic -fno-inline -fno-f2c -Og -ggdb -fbacktrace  ${warnDebug} ${checks}")    # debug flags
  # -O0 -g3 -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wimplicit-interface
  # -Wimplicit-procedure -Wunderflow -Wuninitialized -fcheck=all -fmodule-private -fbacktrace -dump-core -finit-real=nan -ffpe-trap=invalid,zero,overflow
  #
  set(cstd   "-std=${C_DIALECT}11" )
  set(cppstd "-std=${C_DIALECT}++14" )

  option(GFORTRAN_EXTRA_CHECKING "Use extra checks in commandline " OFF)

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  # pgfortran
  message(STATUS " PGI Compiler settings: default USE_CUDA=ON")
  set(USE_CUDA ON)
  set(preproc  "-Mpreprocess ")
  set(dialect  "-Mfreeform  -Mextend -Mnosecond_underscore -Mlarge_arrays ") #-Mstandard -Mallocatable=03
  set(checks   "-Mdclchk -Mchkptr -Mchkstk -Mdepchk -Munixlogical -Mflushz -Mdaz -Mfpmisalign")
  set(warn     "-Minform=warn -Minfo=all,ftn ") # ${checks}")
  # bounds checking cannot be done in CUDA fortran or OpenACC GPU
  set(fordebug "-g ${warn}  -traceback -gopt -Mcuda=debug -Mneginfo=all,ftn -Mpgicoff -traceback -Mprof  ")
  set(forspeed "-O3 -fast " ) # -Munroll -O4  -Mipa=fast -fast -Mcuda=fastmath,unroll -Mvect=nosizelimit,short,simd,sse  ")
  set(forpar   "-mp -acc ")   # -Mconcur=bind,allcores -Mcuda=cuda8.0,cc60,flushz,fma
  set(target   "-m64 ")  #
  set(common   "${preproc} ${dialect} ${target}")
  # further PGI options
  option(PGI_EXTRACT_ALL  "PGI --Extract subprograms for inlining (-Mextract)" OFF)
  option(PGI_LARGE_FILE_SUPPORT  "PGI -- Link with library directory for large file support (-Mlfs)" OFF)
  option(PGI_CUDA_MANAGED_MEMORY "Use CUDA Managed Memory" OFF)
  option(PGI_CUDA_IOMUTEX "Use mutex for IO calls (-Miomutex)" ON)
  option(PGI_CHECKING "Use extra checks in commandline " OFF)
  option(PGI_EXTRA_FAST "Use extra compile options to speed up code e.g. -Munroll -Mvect" OFF)
  #
  option(USE_OPENACC_ONLY "Enable OpenACC without OpenMP (OpenMP on by default)" OFF)
  message(STATUS "In PGI: FFTW should be set with one of the following environment variables: FFTWDIR,FFTW_DIR, or FFTW_ROOT ")
  if(NOT "$ENV{FFTW_DIR}" STREQUAL "")
    set(FFTWDIR "$ENV{FFTW_DIR}")
  elseif (NOT "$ENV{FFTWDIR}" STREQUAL "")
    set(FFTWDIR "$ENV{FFTWDIR}")
  elseif (NOT "$ENV{FFTW_ROOT}" STREQUAL "")
    set(FFTWDIR "$ENV{FFTW_ROOT}")
  else()
    set(FFTWDIR "/usr/local/pgi/src/fftw/")
  endif()

  set(cstd   "-c1x" )
  set(cppstd "--c++14")

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # ifort
  # set(FC "ifort" CACHE PATH "Intel Fortran compiler")
  set(preproc  "-fpp")

  set(dialect  "-free -implicitnone -list-line-len=264 -diag-disable 6477  -diag-disable 406 -gen-interfaces -assume no2underscore -assume buffered_io -assume realloc_lhs")
  set(checks   "-check bounds -check uninit")

  set(warn     "-warn all")
  set(fordebug "-g -debug -O0 -ftrapuv -debug all -check all ${warn} -assume byterecl -align sequence -traceback")
  set(forspeed "-O3 -fp-model fast=2 -inline all -unroll-aggressive -no-fp-port")
  set(forpar   " ")
  set(target   "-no-prec-div -fPIC -xHost -traceback ")
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F08 EQUAL 1)
    set(target "${target} -std08")
  else()
    set(target "${target} -std03")
  endif()

  set(common   "${preproc} ${dialect} ${checks} ${target}")
  # else()
  #   message(" Fortran compiler not supported. Set FC environment variable")
  if(NOT "$ENV{MKLROOT}" STREQUAL "")
    set(MKLROOT $ENV{MKLROOT})
  else()
    message( "MKLROOT must be set using INTEL mklvars ")
  endif()
  if(NOT "$ENV{INTEL_DIR}" STREQUAL "")
    set(INTEL_DIR $ENV{INTEL_DIR})
    set(INTEL_TARGET_ARCH $ENV{INTEL_TARGET_ARCH})
    set(INTEL_TARGET_PLATFORM $ENV{INTEL_TARGET_PLATFORM})
  else()
    message( "INTEL_DIR must be set using INTEL compilervars ")
  endif()
  if (CMAKE_Fortran_COMPILER STREQUAL "mpiif*")
    if(NOT "$ENV{I_MPI_ROOT}" STREQUAL "")
      set(I_MPI_ROOT $ENV{I_MPI_ROOT})
    else()
      message( "I_MPI_ROOT must be set using INTEL mpivars ")
    endif()
  endif()
  option(INTEL_OMP_OVERRIDE "Allow intel compiler to override limits when compiling OpenMP" OFF)
  set(cstd "-std=c11" )
  set(cppstd "-std=c++14" )
endif ()


option(LARGE_FILE_SUPPORT  "GNU -- Link with library directory for large file support (-mcmodel=large)" ON)

string(TOUPPER "${CMAKE_Fortran_COMPILER_ID}" ID_STRING)
message( STATUS "CMAKE_C_FLAGS_RELEASE_INIT: ${CMAKE_C_FLAGS_RELEASE_INIT} (Old version)")
message( STATUS "CMAKE_CXX_FLAGS_RELEASE_INIT: ${CMAKE_CXX_FLAGS_RELEASE_INIT} (Old version)")
set(CMAKE_C_FLAGS_RELEASE_INIT " ${CMAKE_C_FLAGS_RELEASE_INIT} ${cstd} -D${ID_STRING}"  CACHE STRING "Default C release flags -- this cannot be edited " FORCE)
set(CMAKE_CXX_FLAGS_RELEASE_INIT "${CMAKE_CXX_FLAGS_RELEASE_INIT} ${cppstd} -D${ID_STRING}"  CACHE STRING "Default CXX release flags -- this cannot be edited" FORCE)
message( STATUS "CMAKE_C_FLAGS_RELEASE_INIT: ${CMAKE_C_FLAGS_RELEASE_INIT}")
message( STATUS "CMAKE_CXX_FLAGS_RELEASE_INIT: ${CMAKE_CXX_FLAGS_RELEASE_INIT}")
message( STATUS "CMAKE_Fortran_FLAGS_RELEASE_INIT: ${CMAKE_Fortran_FLAGS_RELEASE_INIT} (Old version)")
message( STATUS "CMAKE_Fortran_FLAGS_DEBUG_INIT: ${CMAKE_Fortran_FLAGS_DEBUG_INIT}  (Old version)")
set(CMAKE_Fortran_FLAGS_RELEASE_INIT "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${common} ${forspeed} ${forpar} -D${ID_STRING}"  CACHE STRING "Default release flags -- this cannot be edited (see cmake/FortranOverride.cmake) -- Use CMAKE_Fortran_FLAGS_RELEASE instead" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG_INIT  " ${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${common} ${fordebug} ${forpar} -D${ID_STRING}"  CACHE STRING "Default debug flags -- this cannot be edited (see cmake/FortranOverride.cmake -- Use CMAKE_Fortran_FLAGS_DEBUG instead)" FORCE)
message( STATUS "CMAKE_Fortran_FLAGS_RELEASE_INIT: ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
message( STATUS "CMAKE_Fortran_FLAGS_DEBUG_INIT: ${CMAKE_Fortran_FLAGS_DEBUG_INIT}")
#
# Make recent cmake not spam about stuff
if(POLICY CMP0063)
  cmake_policy(SET CMP0063 OLD)
endif()
if(POLICY CMP0004)
  cmake_policy(SET CMP0004 OLD)
endif()

#message(STATUS Executing: ENV{PATH} $ENV{PATH})
## Remove EMAN bin and possible library paths from PATH
# message(STATUS  "ENV{PATH} -- $ENV{PATH}")
file(TO_CMAKE_PATH "$ENV{PATH}" TMPPATH)
# message(STATUS " FortranOverride Before REGEX: TMPPATH ${TMPPATH} ")
string(REGEX REPLACE  "[^;]\+(EMAN|eman)[2]*;" "" TMPPATH "${TMPPATH}")
# cmake > 3.6   list(FILTER  TMPPATH EXCLUDE  REGEX "[^;]\+(EMAN|eman)[2]*[^;]\+")
#message(STATUS " FortranOverride After REGEX: TMPPATH ${TMPPATH} ")
string (REPLACE ";" ":" TMPPATH "${TMPPATH}")
#file(TO_NATIVE_PATH "${TMPPATH}" TMPPATH)
#message(STATUS " FortranOverride After TO_CMAKE_PATH: TMPPATH ${TMPPATH} ")
set(ENV{PATH} "${TMPPATH}")
#message(STATUS " FortranOverride After TO_CMAKE_PATH: $ENV{PATH} ")



if (NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
  file(TO_CMAKE_PATH   "$ENV{LD_LIBRARY_PATH}" TMPPATH)
  string(REGEX REPLACE  "[^:]\+(eman|EMAN)[2]*[^:]\+" "" TMPPATH "${TMPPATH}")
  file(TO_CMAKE_PATH "${TMPPATH}" TMPPATH)
  set(ENV{LD_LIBRARY_PATH} ${TMPPATH})
elseif(NOT "$ENV{DYLD_LIBRARY_PATH}" STREQUAL "")
  file(TO_CMAKE_PATH "$ENV{DYLD_LIBRARY_PATH}" TMPPATH)
  string(REGEX REPLACE  "[^:]\+(eman|EMAN)[2]*[^:]\+" "" TMPPATH "${TMPPATH}")
  file(TO_CMAKE_PATH "${TMPPATH}" TMPPATH)
  set(ENV{DYLD_LIBRARY_PATH} ${TMPPATH})
endif()
if (NOT "$ENV{FFTW_ROOT}" STREQUAL "")
  file(TO_CMAKE_PATH  "$ENV{FFTW_ROOT}" TMPPATH)
  string(REGEX REPLACE  "[^:]\+(eman|EMAN)[2]*[^:]\+" "" TMPPATH "${TMPPATH}")
  file(TO_CMAKE_PATH "${TMPPATH}" TMPPATH)
  set(ENV{FFTW_ROOT} ${TMPPATH})
elseif (NOT "$ENV{FFTW_DIR}" STREQUAL "")
  file(TO_CMAKE_PATH  "$ENV{FFTW_DIR}" TMPPATH)
  string(REGEX REPLACE  "[^:]\+(eman|EMAN)[2]*[^:]\+" "" TMPPATH "${TMPPATH}")
  file(TO_CMAKE_PATH "${TMPPATH}" TMPPATH)
  set(ENV{FFTW_DIR} ${TMPPATH})
elseif (NOT "$ENV{FFTW_ROOT}" STREQUAL "")
  file(TO_CMAKE_PATH "$ENV{FFTWDIR}" TMPPATH)
  string(REGEX REPLACE  "[^:]\+(eman|EMAN)[2]*[^:]\+" "" TMPPATH "${TMPPATH}")
  file(TO_CMAKE_PATH "${TMPPATH}" TMPPATH)
  set(ENV{FFTWDIR} ${TMPPATH})
endif()

# message(STATUS  "XXXXX ENV{PATH} -- $ENV{PATH}")

# execute_process( COMMAND cmake -E env which ${CMAKE_Fortran_COMPILER}
#   OUTPUT_VARIABLE GFC_VERSION
#   ERROR_VARIABLE GFC_VERSION_ERROR
#   RESULT_VARIABLE GFC_VERSION_RESULT)
# message(STATUS " whichFO ${CMAKE_Fortran_COMPILER}  :res ${GFC_VERSION_RESULT} :err ${GFC_VERSION_ERROR} :out ${GFC_VERSION}")


# Make sure the build type is uppercase
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

#################################################################
# CONFIGURATION TYPES & BUILD MODE
#################################################################
set(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE RELWITHDEBINFO )
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: NONE DEBUG RELEASE RELWITHDEBINFO"
    FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE RELWITHDEBINFO TESTING)
endif(NOT CMAKE_BUILD_TYPE)

if(BT STREQUAL "RELEASE")
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are DEBUG, RELEASE, RELWITHDEBINFO or TESTING."
    FORCE)
elseif(BT STREQUAL "DEBUG")
  set (CMAKE_BUILD_TYPE DEBUG CACHE STRING
    "Choose the type of build, options are DEBUG, RELEASE, RELWITHDEBINFO or TESTING."
    FORCE)
elseif(BT STREQUAL "RELWITHDEBINFO")
  set (CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING
    "Choose the type of build, options are DEBUG, RELEASE, RELWITHDEBINFO  or TESTING."
    FORCE)
elseif(BT STREQUAL "TESTING")
  set (CMAKE_BUILD_TYPE TESTING CACHE STRING
    "Choose the type of build, options are DEBUG, RELEASE, RELWITHDEBINFO or TESTING."
    FORCE)
elseif(NOT BT)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are DEBUG, RELEASE, RELWITHDEBINFO or TESTING."
    FORCE)
  message(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
else()
  message(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, RELWITHDEBINFO or TESTING")
endif(BT STREQUAL "RELEASE")


IF (LINUX)
  option(MAP_TEXT_HUGE_PAGES "Remap hot static code onto huge pages" ON)
ENDIF()
