
####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
include(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)
enable_language(Fortran C)
#include(${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
#include(${CMAKE_ROOT}/Modules/CMakeDetermineCCompiler.cmake)

if(NOT $ENV{CPP} STREQUAL "")
  set(TMP_CPP_COMPILER $ENV{CPP})
else()
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set(TMP_CPP_COMPILER fpp)
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
    set(TMP_CPP_COMPILER "pgcc -E")
  endif()
endif()

set (CLANG_FAIL_MSG  "FATAL ERROR: SIMPLE cannot support Clang Preprocessor.
Set FC,CC, and CPP environment variables to GNU compilers, e.g. in bash:
export FC=/sw/bin/gfortran
export CC=/sw/bin/gcc
export CPP=/sw/bin/cpp-5

Clang has overridden Gfortran links /usr/bin/gfortran GRRR!
Make sure you prepend  /usr/local/bin (Homebrew) or /opt/local/bin (MacPorts) or /sw/bin (FINK) to PATH.
In LD_LIBRARY_PATH prepend the appropriate lib path.
    ")

message(STATUS "Making sure your preprocessor points to the correct binary")
execute_process(COMMAND "${TMP_CPP_COMPILER}" --version
  OUTPUT_VARIABLE ACTUAL_FC_TARGET
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
  message(STATUS "Found cpp is actually Clang -- Trying other paths")
  find_file (
    TMP_CPP_COMPILER
    NAMES cpp- cpp-4.9 cpp-5 cpp-6 cpp5 cpp6
    PATHS /usr/local/bin /opt/local/bin /sw/bin /usr/bin
    #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
    DOC "Searing for GNU cpp preprocessor "
    #  [NO_DEFAULT_PATH]
    #  [NO_CMAKE_ENVIRONMENT_PATH]
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    #  [NO_CMAKE_SYSTEM_PATH]
    #  [CMAKE_FIND_ROOT_PATH_BOTH |
    #   ONLY_CMAKE_FIND_ROOT_PATH |
    #   NO_CMAKE_FIND_ROOT_PATH]
    )
  if(NOT EXISTS "${TMP_CPP_COMPILER}")
    message( FATAL_ERROR  "Cannot find GNU cpp compiler -- 
${CLANG_FATAL_MSG}")
  endif()
endif()
set(CMAKE_CPP_COMPILER ${TMP_CPP_COMPILER})



#if(APPLE)
#  message(STATUS  "Try setting the GNU compiler")
#    include(${CMAKE_ROOT}/Modules/Compiler/GNU-Fortran.cmake)
#    include(${CMAKE_ROOT}/Modules/Platform/Darwin-GNU-Fortran.cmake)
#    __darwin_compiler_gnu(Fortran)

# if (CMAKE_Fortran_COMPILER_ID STREQUAL "Clang" OR "${CMAKE_GENERATOR}" MATCHES "XCode")

#   message( FATAL_ERROR "${CLANG_FATAL_MSG}" )

# elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
message(STATUS "Making sure your Fortran compiler points to the correct binary")
if(Fortran_COMPILER_NAME MATCHES "gfortran*")
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_FC_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
    message(STATUS "gfortran points to Clang -- Trying other paths")
    find_file (
      CMAKE_Fortran_COMPILER
      NAMES gfortran- gfortran-4.9 gfortran-5 gfortran-6 gfortran5 gfortran6
      PATHS /usr/local/bin /opt/local/bin /sw/bin /usr/bin
      #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
      DOC "Searing for GNU gfortran preprocessor "
      #  [NO_DEFAULT_PATH]
      #  [NO_CMAKE_ENVIRONMENT_PATH]
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      #  [NO_CMAKE_SYSTEM_PATH]
      #  [CMAKE_FIND_ROOT_PATH_BOTH |
      #   ONLY_CMAKE_FIND_ROOT_PATH |
      #   NO_CMAKE_FIND_ROOT_PATH]
      )
    if(NOT EXIST "${CMAKE_Fortran_COMPILER}")
      message( FATAL_ERROR  "Cannot find ${CMAKE_Fortran_COMPILER} -- 
${CLANG_FATAL_MSG}")
    endif()
  endif()
endif()


set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS} "f03;F03;f08;F08")
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
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE RELWITHDEBINFO)
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


#################################################################
# Setting up options
#################################################################
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message(STATUS "Fortran compiler: ${Fortran_COMPILER_NAME}")

#################################################################
# FFLAGS depend on the compiler and the build type
#################################################################
# set(EXTRA_FLAGS "${EXTRA_FLAGS} -fPIC")
include(ProcessorCount)
ProcessorCount(NUM_JOBS)

# if (${ENABLE_HARD_OPTIMIZE})
#   if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" ) #AND Fortran_COMPILER_NAME MATCHES "gfortran*")
#     set(EXTRA_FLAGS "${EXTRA_FLAGS}  -O3 -ffast-math -funroll-all-loops")
#   elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR Fortran_COMPILER_NAME MATCHES "pgf*")
#     set(EXTRA_FLAGS "${EXTRA_FLAGS}  -O3 -fast ")
#   elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
#     set(EXTRA_FLAGS "${EXTRA_FLAGS}  -O3 -inline all -unroll-aggressive -fast ")
#   endif ()
# endif ()

# #if (${${PROJECT_NAME}_ENABLE_OPENMP})
# #  find_package( OpenMP )
# #  if (${OPENMP_FOUND})
# #    set(EXTRA_FLAGS ${OpenMP_Fortran_FLAGS})
# #    add_definitions("-DOPENMP")
# #  else ()
# #    option(${PROJECT_NAME}_ENABLE_OPENMP "To activate the OpenMP extensions for Fortran" OFF)
# #  endif ()
# #endif ()



##############################################
# Linker            (FROM FACEBOOK HHVM)
#############################################
set(GOLD_FOUND FALSE)
mark_as_advanced(GOLD_FOUND)

#############################################
## DEBUG is used as a variable in some files so use _ as prefix
#############################################
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  add_definitions("-D_DEBUG")
endif()
#############################################

add_definitions(-D${CMAKE_Fortran_COMPILER_ID})
message(STATUS "Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gfortran
  set(dialect  "-ffree-form -cpp -fimplicit-none  -ffree-line-length-none")                 # language style
  set(checks   "-fcheck-array-temporaries  -frange-check -ffpe-trap=invalid,zero,overflow -fstack-protector -fstack-check") # checks
  set(warn     "-Wall -Wextra -Wimplicit-interface  -Wline-truncation")                     # warning flags
  set(fordebug "-pedantic -fno-inline -fno-f2c -Og -ggdb -fbacktrace -fbounds-check")       # debug flags
  set(forspeed "-O3 -ffast-math -finline-functions -funroll-all-loops -fno-f2c ")           # optimisation
  set(forpar   "-fopenmp -pthread ")                                                         # parallel flags
  set(target   "-march=native -fPIC")                                                       # target platform
  set(common   "${dialect} ${checks} ${target} ")
  #
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  # pgfortran
  set(dialect  "-Mpreprocess -Mfreeform  -Mstandard -Mallocatable=03 -Mextend")
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
  #-fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -std03
  set(warn     "-warn all")
  set(fordebug "-debug -O0 -ftrapuv -debug all -check all")
  set(forspeed "-O3 -fp-model fast=2 -inline all -unroll-aggressive ")
  set(forpar   "-qopenmp")
  set(target   "-xHOST -no-prec-div -static -fPIC")
  set(common   "${dialect} ${checks} ${target} ${warn} -DINTEL")
  # else()
  #   message(" Fortran compiler not supported. Set FC environment variable")
endif ()
set(CMAKE_Fortran_FLAGS_RELEASE_INIT "${common} ${forspeed} ${forpar} " )
set(CMAKE_Fortran_FLAGS_DEBUG_INIT   "${common} ${warn} ${fordebug} ${forpar} -g ")

#############################################
## COMPLER SPECIFIC SETTINGS
#############################################
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" ) #AND Fortran_COMPILER_NAME MATCHES "gfortran*")
  #############################################
  #
  ## GNU fortran
  #
  #############################################
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion OUTPUT_VARIABLE GFC_VERSION)
  if (NOT (GFC_VERSION VERSION_GREATER 4.8 OR GFC_VERSION VERSION_EQUAL 4.8))
    message(FATAL_ERROR "${PROJECT_NAME} requires gfortran 4.8 or greater.")
  endif ()
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -cpp -fimplicit-none -fall-intrinsics -ffree-line-length-none ")
  set(CMAKE_CPP_COMPILER_FLAGS           "-E -w -C -CC -P")

  set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${EXTRA_FLAGS} ")
  # set(CMAKE_Fortran_FLAGS_DEBUG         "${CMAKE_Fortran_FLAGS_DEBUG_INIT} -O0 -g3 -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wuninitialized -fcheck=all -fmodule-private -fbacktrace -dump-core -finit-real=nan " CACHE STRING "" FORCE)
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL     "-Os ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_DEBUG_INIT}")
  
  # #  CMAKE_EXE_LINKER_FLAGS
  # if (LINK_TIME_OPTIMISATION)
  #   set(CMAKE_EXE_LINKER_FLAGS             "${CMAKE_EXE_LINKER_FLAGS_INIT} -flto ")
  #   set(CMAKE_SHARED_LINKER_FLAGS           "${CMAKE_EXE_LINKER_FLAGS_INIT} -flto -flto=${NUM_JOBS}")
  # endif(LINK_TIME_OPTIMISATION)

  # ## use gold as linker (from HVVM)

  # find_program(GOLD_EXECUTABLE NAMES gold ld.gold DOC "path to gold")
  # mark_as_advanced(GOLD_EXECUTABLE)
  # if(GOLD_EXECUTABLE)
  #   set(GOLD_FOUND TRUE)
  #   execute_process(COMMAND ${GOLD_EXECUTABLE} --version
  #     OUTPUT_VARIABLE GOLD_VERSION
  #     OUTPUT_STRIP_TRAILING_WHITESPACE)
  #   message(STATUS "Found gold: ${GOLD_EXECUTABLE}")
  #   add_definitions(" -fuse-ld=gold -Wl,--threads")
  # else()
  #   message(STATUS "Could not find gold linker. Using the default")
  # endif()

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  #############################################
  #
  ## INTEL fortran
  #
  #############################################
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -fpp -I${MKLROOT}/include  -assume realloc_lhs -assume source_include")
  #
  set(CMAKE_AR                           "xiar")
  set(CMAKE_CPP_COMPILER                 "fpp")
  set(CMAKE_CPP_COMPILER_FLAGS           "  -noJ -B -C -P")
  set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS_INIT} ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG_INIT} -O0 -debug all -check all -warn all -extend-source 132 -traceback -gen-interfaces")
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL     "-Os ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "-O3 ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_DEBUG_INIT}")
  # set(CMAKE_EXE_LINKER_FLAGS             "${CMAKE_EXE_LINKER_FLAGS_INIT}")
  # set(CMAKE_SHARED_LINKER_FLAGS          "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${EXTRA_LIBS}")
  # set(CMAKE_STATIC_LINKER_FLAGS          "${CMAKE_STATIC_LINKER_FLAGS_INIT} ${EXTRA_LIBS}")
  # if (LINK_TIME_OPTIMISATION)
  #   set(CMAKE_EXE_LINKER_FLAGS           "${CMAKE_EXE_LINKER_FLAGS} -ipo-separate -ipo-jobs=${NUM_JOBS}")
  #   set(CMAKE_SHARED_LINKER_FLAGS        "${CMAKE_SHARED_LINKER_FLAGS} -ip -ipo-separate -ipo-jobs=${NUM_JOBS}")
  #   set(CMAKE_STATIC_LINKER_FLAGS        "${CMAKE_STATIC_LINKER_FLAGS} -ip -ipo")
  # endif(LINK_TIME_OPTIMISATION)
  add_definitions("-DOPENMP") 

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR Fortran_COMPILER_NAME MATCHES "pgfortran.*")
  #############################################
  #
  ## Portland Group fortran
  #
  #############################################
  message(STATUS "NVIDIA PGI Linux compiler")
  set(PGICOMPILER ON)
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -Mpreprocess -module ${CMAKE_Fortran_MODULE_DIRECTORY} -I${CMAKE_Fortran_MODULE_DIRECTORY}")
  # NVIDIA PGI Linux compiler
  set(CMAKE_CPP_COMPILER                "pgcc -E ")
  set(CMAKE_CPP_COMPILER_FLAGS                " ")
  set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS_INIT} ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG_INIT} ")
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-gopt ${CMAKE_Fortran_FLAGS_RELEASE_INIT} -Mneginfo=all")
  # set(CMAKE_EXE_LINKER_FLAGS             "${CMAKE_EXE_LINKER_FLAGS_INIT} -acclibs -cudalibs -Mcudalib=cufft,curand")
  # #  CMAKE_SHARED_LINKER_FLAG

  ################################################################
  # CUDA PGI Default options
  ################################################################
  # default PGI library
  set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)
  set(CUDA_rt_LIBRARY  /usr/lib/x86_64-linux-gnu/librt.so)
  ## PGI does not have the OpenMP proc_bind option
  add_definitions("-Dproc_bind\\(close\\)=\"\"")  # disable proc_bind in OMP 
  # #  add_definitions("-module ${CMAKE_Fortran_MODULE_DIRECTORY}") # pgc++ doesn't have -module
  #   set(CMAKE_SHARED_LINKER_FLAGS          "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${EXTRA_LIBS}  -module ${CMAKE_Fortran_MODULE_DIRECTORY}")
  #   set(CMAKE_STATIC_LINKER_FLAGS        "${CMAKE_STATIC_LINKER_FLAGS_INIT} ${EXTRA_LIBS}  -module ${CMAKE_Fortran_MODULE_DIRECTORY}")

  #   #  CMAKE_EXE_LINKER_FLAGS
  #   if (LINK_TIME_OPTIMISATION)
  #     set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS} -Mipa=fast ")
  #     set(CMAKE_EXE_LINKER_FLAGS           "${CMAKE_EXE_LINKER_FLAGS} -Mipa=fast")
  #     set(CMAKE_SHARED_LINKER_FLAGS        "${CMAKE_SHARED_LINKER_FLAGS} -Mipa=fast")
  #     set(CMAKE_STATIC_LINKER_FLAGS        "${CMAKE_STATIC_LINKER_FLAGS} -Mipa=fast")
  #   endif(LINK_TIME_OPTIMISATION)


elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Clang")

  #############################################
  ## APPLE Clang
  #############################################
  message ("Clang is not supported.  Please use GNU toolchain with either Homebrew, MacPorts or Fink.")
  # find_package(LLVM REQUIRED CONFIG)

  # set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS_RELEASE_INIT} -Wno-mismatched-tags -Qunused-arguments")
  # if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  #   # In OSX, clang requires "-stdlib=libc++" to support C++11
  #   set(CMAKE_Fortran_FLAGS              "${CMAKE_Fortran_FLAGS} -stdlib=f2003")
  #   set(CMAKE_EXE_LINKER_FLAGS           "-stdlib=libc++")
  # endif()

else ()
  #############################################
  ## UNKNOWN fortran
  #############################################
  message (STATUS " Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID}")
  message (STATUS " Set environment variable FC to fortran compiler and rebuild cache.")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif () # COMPILER_ID



if (IMAGE_TYPE_DOUBLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIMAGETYPEDOUBLE" )
endif()


################################################################
# Compiler-specific C++11/Modern fortran activation.
################################################################
# IF(UNIX)
# if ()
# elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Clang")
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-mismatched-tags -Qunused-arguments")
#   if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     # In OSX, clang requires "-stdlib=libc++" to support C++11
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stdlib=f95")
#     set(SIMPLE_EXTRA_LINKER_FLAGS "-stdlib=libc++")
#   endif ()
# else ()
#   message(FATAL_ERROR "Your C++ compiler does not support C++11.")
# endif ()

# ELSE(UNIX)
#   IF(WIN32)
#     SET(GUI "Win32")
#   ELSE(WIN32)
#     SET(GUI "Unknown")
#   ENDIF(WIN32)
# ENDIF(UNIX)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

if (APPLE)
  message(STATUS "Applied CUDA OpenMP macOS workaround")
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS_BACKUP "${CMAKE_SHARED_LIBRARY_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -Wno-unused-function")
  string(REGEX REPLACE "-fopenmp[^ ]*" "" CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS}")
endif()


################################################################
# FFTW  -- MKL core already inlcuded in Intel config
################################################################
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # Append MKL FFTW interface libs
  if(BUILD_SHARED_LIBS)
    set(EXTRA_LIBS   ${EXTRA_LIBS} -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl )
  else()
    set(EXTRA_LIBS    ${EXTRA_LIBS} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
  endif()
  set(BUILD_NAME "${BUILD_NAME}_MKL" )

else()
  if (FFTW_DIR)
    set(FFTW_ROOT ${FFTW_DIR})
  endif()
  find_package(FFTW QUIET)
  if (NOT FFTW_FOUND)
    message(FATAL_ERROR "Unable to find FFTW")
  else()
    message(STATUS "fftw3 found")
    message(STATUS "lib: ${FFTW_LIBRARIES}")
    include_directories(" ${FFTW_INCLUDE_DIRS} ")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES})
  endif()
  set(BUILD_NAME "${BUILD_NAME}_FFTW" )
endif()


# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
  SET(GNUNATIVE "-mtune=native")
ELSE()
  SET(GNUNATIVE "-march=native")
ENDIF(APPLE)

################################################################
# Generic Flags 
################################################################

# # Don't add underscores in symbols for C-compatability
# message(STATUS "Testing flag no underscore")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#   Fortran
#   "-fnosecond-underscore"     # GNU
#   "-Mnosecond_underscore" # PGI
#   "-assume nounderscore"  # Intel
#   )


# Optimize for the host's architecture
# message(STATUS "Testing flag host arch")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#   Fortran 
#   "-xHost"        # Intel
#   "/QxHost"       # Intel Windows
#   ${GNUNATIVE}    # GNU
#   "-tp=x64"      # Portland Group - generic 64 bit platform
#   )



# # Preprocessing
# # message(STATUS "Testing flag cpp")
# # SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
# #   Fortran REQUIRED
# #   "-fpp"                # Intel
# #   "-cpp"                # GNU
# #   "-Mpreprocess"        # PGI
# #   )
# # free form
# message(STATUS "Testing flag free-form")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#   Fortran 
#   "-ffree-form"       # GNU
#   "-free"             # Intel
#   "-Mfreeform"        # PGI
#   )
# # line length
# message(STATUS "Testing flag no line length")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#   Fortran
#   "-ffree-line-length-none"    # GNU
#   "-Mextend"                   # PGI
#   "-extend-source"            # Intel
#   "-list-line-len=264"        # Intel
#   )

# ###################
# ### DEBUG FLAGS ###
# ###################
# ## Disable optimizations
# message(STATUS "Testing debug flags")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#   Fortran 
#   "-O0" # All compilers not on Windows
#   "/Od" # Intel Windows
#   )

# # Turn on all warnings
# message(STATUS "Testing warn all flags")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#   Fortran
#   "/warn:all" # Intel Windows
#   "-Wall"     # GNU
#   "-warn all" # Intel
#   )

# # Traceback
# message(STATUS "Testing flags traceback ")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#   Fortran
#   "-traceback"   # Intel/Portland Group
#   "/traceback"   # Intel Windows
#   "-fbacktrace"  # GNU (gfortran)
#   "-ftrace=full" # GNU (g95)
#   )

# # Check array bounds
# message(STATUS "Testing flags array bounds check")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#   Fortran 
#   "-fbounds-check" # GNU (Old style)
#   "-fcheck=bounds" # GNU (New style)
#   "-Mbounds"       # Portland Group
#   "/check:bounds"  # Intel Windows
#   "-check bounds"  # Intel
#   )


#####################
### TESTING FLAGS ###
#####################

# Optimizations
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
#  Fortran REQUIRED
#  "-O2" # All compilers not on Windows
#  "/O2" # Intel Windows
#  )

#####################
### RELEASE FLAGS ###
#####################

# Unroll loops
#message(STATUS "Testing flags unroll")
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#  Fortran
#  "-unroll-aggressive"        # Intel
#  "-funroll-loops"            # GNU, Intel, Clang
#  "/unroll"                   # Intel Windows
#  "-Munroll"                  # Portland Group
#  )


if(ENABLE_LINK_TIME_OPTIMISATION)
  # Interprocedural (link-time) optimizations
  # See IPO performance issues https://software.intel.com/en-us/node/694501
  #  may need to use -ipo-separate in case of multi-file problems),  using the -ffat-lto-objects compiler option is provided for GCC compatibility
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-ipo-separate -ffat-lto-object"   # Intel (Linux OS X)
    "/Qipo-separate"                   # Intel Windows
    "-flto "                           # GNU
    "-Mipa=fast,libinline"    # Portland Group
    )

  # Single-file optimizations
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-ip-dir=.profiling"           # Intel
    "/Qip-dir=_profiling"          # Intel Windows
    "-fprofile-dir=.profiling" # GNU
    # PGI
    )
endif(ENABLE_LINK_TIME_OPTIMISATION)
if(ENABLE_AUTOMATIC_VECTORIZATION)
  # Interprocedural (link-time) optimizations
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-m"              # Intel ( may need to use -ipo-separate in case of multi-file problems)
    "/Qipo"             # Intel Windows
    "-flto "            # GNU
    "-Mvect"    # Portland Group
    )
endif()

# # Fast math code
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#   Fortran
#   "-fastmath"        # Intel
#   "-ffast-math"      # GNU
#   "-Mcuda=fastmath"  # Portland Group
#   )

# # Vectorize code
if (ENABLE_FAST_MATH_OPTIMISATION)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-fast"             # Intel, PGI
    "/Ofast"     # Intel Windows
    "-Ofast"            # GNU
    )
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-fastmath"        # Intel
    "-ffast-math"      # GNU
    "-Mcuda=fastmath"  # Portland Group
    )
endif(ENABLE_FAST_MATH_OPTIMISATION)

# Auto parallelize
if (ENABLE_AUTO_PARALLELISE)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-parallel -opt-report-phase=par -opt-report:5"            # Intel (Linux)
    "/Qparallel /Qopt-report-phase=par /Qopt-report:5"         # Intel (Windows)
    "-Mconcur"                                                 # PGI
    )
endif()

# Auto parallelize with OpenACC
if (USE_OPENACC)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-acc"                 # PGI
    "-fopenacc"            # GNU
    "/acc"
    )
endif()

# # Instrumentation
# if (USE_INSTRUMENTATION)
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#   Fortran
#   "-Minstrument -Mpfi"      # PGI
#   "-finstrument "           # GNU
#   )
# endif()

# Profile-feedback optimisation
if(ENABLE_PROFILE_OPTIMISATION)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-Mpfo"                 # PGI
    "-fpfo "                # GNU
    "-prof-gen"             # Intel (Linux, OS X)
    "/Qprof-gen"            # Intel (Windows)
    )
endif()

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  ")
#set(CMAKE_FCPP_FLAGS " -C -P ") # Retain comments due to fortran slash-slash
#set(CMAKE_Fortran_CREATE_PREPROCESSED_SOURCE "${CMAKE_FCPP_COMPILER} <DEFINES> <INCLUDES> <FLAGS> -E <SOURCE> > <PREPROCESSED_SOURCE>")

add_definitions(" -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# Override Fortran preprocessor
set(CMAKE_Fortran_COMPILE_OBJECT "grep --silent '#include' <SOURCE> && ( ${CMAKE_CPP_COMPILER} ${CMAKE_CPP_COMPILER_FLAGS} -DOPENMP <DEFINES> <INCLUDES> <SOURCE> > <OBJECT>.f90 &&  <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> \${FFLAGS} -c <OBJECT>.f90 -o <OBJECT> ) || <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> \${FFLAGS} -c <SOURCE> -o <OBJECT>")


# Option for code coverage
option(CODE_COVERAGE "Build code coverage results, requires GCC compiler (forces Debug build)" OFF)
if(CODE_COVERAGE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG
      "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)

    # Ensure that CDash targets are always enabled if coverage is enabled.
    if (NOT CDASH_SUPPORT)
      get_property(HELP_STRING CACHE CDASH_SUPPORT PROPERTY HELPSTRING)
      set(CDASH_SUPPORT ON CACHE BOOL "${HELP_STRING}" FORCE)
      message(STATUS "Enabling CDash targets as coverage has been enabled.")
    endif()
  endif(CMAKE_COMPILER_IS_GNUCXX)
endif()

