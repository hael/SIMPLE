# This overrides the default CMake Debug and Release compiler options.
# The user can still specify different options by setting the
# CMAKE_Fortran_FLAGS_[RELEASE,DEBUG] variables (on the command line or in the
# CMakeList.txt). This files serves as better CMake defaults and should only be
# modified if the default values are to be changed. Project specific compiler
# flags should be set in the CMakeList.txt by setting the CMAKE_Fortran_FLAGS_*
# variables.
if(NOT $ENV{FC} STREQUAL "")
  set(CMAKE_Fortran_COMPILER_NAMES $ENV{FC})
else()
  set(CMAKE_Fortran_COMPILER_NAMES gfortran)
  set(ENV{FC} "gfortran")
  set(ENV{CC} "gcc")
  set(ENV{CXX} "g++")
endif()

set (CLANG_FAIL_MSG  "FATAL ERROR: SIMPLE cannot support Clang.
Set FC,CC,CXX environment variables to GNU compilers, e.g. in bash:
export FC=/sw/bin/gfortran
export CC=/sw/bin/gcc
export CXX=/sw/bin/g++

Clang has overridden Gfortran links /usr/bin/gfortran GRRR!
In PATH prepend  /usr/local/bin (Homebrew) or /opt/local/bin (MacPorts) or /sw/bin (FINK)
In LD_LIBRARY_PATH prepend the appropriate lib path.
OR set FC, CC CXX variables with absolute paths.
    ")


if(APPLE)
  # Try setting the GNU compiler
  __darwin_compiler_gnu()

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
    message( FATAL_ERROR "${CLANG_FATAL_MSG}" )
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
    message(STATUS "Making sure your Mac OS X GNU compiler points to the correct binary")
    if(Fortran_COMPILER_NAME MATCHES "gfortran*")
      execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
        OUTPUT_VARIABLE ACTUAL_FC_TARGET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
        message( FATAL_ERROR  "${CLANG_FATAL_MSG}")
      endif()
      if(NOT $ENV{CPP} STREQUAL "")
        set(CMAKE_CPP_COMPILER $ENV{CPP})
      else()
        set(CMAKE_CPP_COMPILER cpp)
      endif()
      execute_process(COMMAND ${CMAKE_CPP_COMPILER} --version
        OUTPUT_VARIABLE ACTUAL_FC_TARGET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
        message( FATAL_ERROR  "${CLANG_FATAL_MSG}  -- CPP compiler ${CMAKE_CPP_COMPILER} links to Clang")
      endif()
    endif()
  endif()
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_FC_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
    message( FATAL_ERROR "${CLANG_FATAL_MSG} -- C++ compiler")
  endif()
  execute_process(COMMAND ${CMAKE_C_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_FC_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
    message( FATAL_ERROR "${CLANG_FATAL_MSG} -- C compiler")
  endif()

endif(APPLE)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create separate directory for build files.
")
endif()


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # gfortran
    set(dialect  "-ffree-form -cpp -fimplicit-none  -ffree-line-length-none")                 # language style
    set(checks   "-fcheck-array-temporaries  -frange-check -ffpe-trap=invalid,zero,overflow -fstack-protector -fstack-check") # checks
    set(warn     "-Wall -Wextra -Wimplicit-interface  -Wline-truncation")                     # warning flags
    set(fordebug "-pedantic -fno-inline -fno-f2c -Og -ggdb -fbacktrace -fbounds-check")       # debug flags
    set(forspeed "-O3 -ffast-math -finline-functions -funroll-all-loops -fno-f2c ")           # optimisation
    set(forpar   "-fopenmp -pthread ")                                                         # parallel flags
    set(target   "-march=native -fPIC")                                                       # target platform
    set(common   "${dialect} ${checks} ${target} ${warn} ")
#
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    # pgfortran
    set(dialect  "-Mpreprocess -Mfreeform  -Mstandard -Mallocatable=03")
    set(checks   "-Mdclchk  -Mchkptr -Mchkstk  -Munixlogical -Mlarge_arrays -Mflushz -Mdaz -Mfpmisalign")
    set(warn     "-Minform=warn")
    # bounds checking cannot be done in CUDA fortran or OpenACC GPU
    set(fordebug "-Minfo=all,ftn  -traceback -gopt -Mneginfo=all,ftn -Mnodwarf -Mpgicoff -traceback -Mprof -Mbound -C")
    set(forspeed "-Munroll -O4  -Mipa=fast -fast -Mcuda=fastmath,unroll -Mvect=nosizelimit,short,simd,sse -mp -acc ")
    set(forpar   "-Mconcur -Mconcur=bind,allcores -Mcuda=cuda8.0,cc60,flushz,fma ")
    set(target   " -m64 -fPIC ")
    set(common   " ${dialect} ${checks} ${target} ${warn}  -DPGI")
#
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # ifort
    # set(FC "ifort" CACHE PATH "Intel Fortran compiler")
    set(dialect  "-fpp -free -implicitnone -std08  -80")
    set(checks   "-check bounds -check uninit -assume buffered_io -assume byterecl -align sequence  -diag-disable 6477  -gen-interfaces ") # -mcmodel=medium -shared-intel
    set(warn     "-warn all")
    set(fordebug "-debug -O0 -ftrapuv -debug all -check all")
    set(forspeed "-O3 -fp-model fast=2 -inline all -unroll-aggressive ")
    set(forpar   "-qopenmp")
    set(target   "-xHOST -no-prec-div -static -fPIC")
    set(common   "${dialect} ${checks} ${target} ${warn} -DINTEL")
# else()
#   message(" Fortran compiler not supported. Set FC environment variable")
  endif ()
  #
  
   set(CMAKE_Fortran_FLAGS_RELEASE_INIT "${common} ${forspeed} ${forpar} ")
   set(CMAKE_Fortran_FLAGS_DEBUG_INIT   "${common} ${fordebug} ${forpar} -g ")
#
# Make recent cmake not spam about stuff
if(POLICY CMP0063)
    cmake_policy(SET CMP0063 OLD)
endif()
if(POLICY CMP0004)
    cmake_policy(SET CMP0004 OLD)
endif()
#
set(ENV_PATH ENV{PATH})

