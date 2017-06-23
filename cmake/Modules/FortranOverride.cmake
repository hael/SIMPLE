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

