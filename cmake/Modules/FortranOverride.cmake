# This overrides the default CMake Debug and Release compiler options.
# The user can still specify different options by setting the
# CMAKE_Fortran_FLAGS_[RELEASE,DEBUG] variables (on the command line or in the
# CMakeList.txt). This files serves as better CMake defaults and should only be
# modified if the default values are to be changed. Project specific compiler
# flags should be set in the CMakeList.txt by setting the CMAKE_Fortran_FLAGS_*
# variables.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # gfortran
    set(dialect  "-ffree-form -std=f2008 -fimplicit-none")                     #language style
    set(checks   "-fbounds-check -fcheck-array-temporaries -fmax-errors=1 ")   #checks
    set(warn     "-Wall -Wextra -Wimplicit-interface ")                        #warning flags
    set(fordebug "-DTRACE -fno-inline -fno-f2c -Og -fbacktrace")               #debug flags
    set(forspeed "-ffast-math -funroll-all-loops -fno-f2c -O3")                #optimisation
    set(forpar   "-fopenmp -pthread")                                          #parallel flags
    set(target   "-march=native -fPIC")                                        #platform
    set(common   "${dialect} ${checks} ${target} ${warn} -DGNU")               #general

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    # pgfortran
    set(dialect  "-Mpreprocess -Mfreeform  -Mstandard -Mallocatable=03")
    set(checks   "-Mdclchk  -Mchkptr -Mchkstk  -Munixlogical -Mlarge_arrays -Mflushz -Mdaz -Mfpmisalign")
    set(warn     "-Minform=warn")
    # bounds checking cannot be done in CUDA fortran or OpenACC GPU
    set(fordebug "-Minfo=all,ftn  -traceback -gopt -Mneginfo=all,ftn -Mnodwarf -Mpgicoff -traceback -Mprof -Mbound -C")
    set(forspeed "-Munroll -O4  -Mipa=fast -fast -Mcuda=fastmath,unroll -Mvect=nosizelimit,short,simd,sse -mp -acc ")
    set(forpar   "-Mconcur -Mconcur=bind,allcores -Mcuda=cuda8.0,cc60,flushz,fma ")
    set(target   "-tp=p7-64 -m64 -fPIC ")
    set(common   " ${dialect} ${checks} ${target} ${warn}  -DPGI")

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # ifort
    set(dialect  "-cpp -free -implicitnone -stand f08  -80")
    set(checks   "-check bounds -check uninit")
    set(warn     "-warn all")
    set(fordebug "-debug -O0 -ftrapuv -debug all -check all")
    set(forspeed "-O3 -fp-model fast=2")
    set(forpar   "-openmp")
    set(target   "-xHOST -no-prec-div -static -fPIC")
    set(common   "${dialect} ${checks} ${target} ${warn} -DINTEL")
  else()
    message(" Fortran compiler not supported. Set FC environment variable")
  endif ()
    
    set(CMAKE_Fortran_FLAGS_RELEASE_INIT "${common} ${forspeed} ${forpar} ")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT   "${common} ${fordebug} ${forpar} -g ")


# Make recent cmake not spam about stuff
if(POLICY CMP0063)
    cmake_policy(SET CMP0063 OLD)
endif()
