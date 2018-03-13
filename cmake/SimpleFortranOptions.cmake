#------------------------------------------------------------------------------!
# SIMPLE              Elmlund & Elmlund Lab          simplecryoem.com          !
#------------------------------------------------------------------------------!

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
include(SetCompileFlag)

#include(${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
#include(${CMAKE_ROOT}/Modules/CMakeDetermineCCompiler.cmake)
## Double check cpp is not Clang
if(NOT CMAKE_CPP_COMPILER)
  if(NOT $ENV{CPP} STREQUAL "")
    set(TMP_CPP_COMPILER $ENV{CPP})
  else()
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(TMP_CPP_COMPILER fpp)
    elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
      set(TMP_CPP_COMPILER "pgcc -E")
    endif()
  endif()
else()
  set(TMP_CPP_COMPILER ${CMAKE_CPP_COMPILER})
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
      DOC "Searching for GNU gfortran preprocessor "
      )
    if(NOT EXIST "${CMAKE_Fortran_COMPILER}")
      message( FATAL_ERROR  "Cannot find ${CMAKE_Fortran_COMPILER} --
${CLANG_FATAL_MSG}")
    endif()
  endif()
endif()

get_filename_component(FORTRAN_PARENT_DIR ${CMAKE_Fortran_COMPILER} PATH)
message(STATUS "FORTRAN_PARENT_DIR ${FORTRAN_PARENT_DIR}")

message(STATUS "Making sure your preprocessor points to the correct binary")
if(TMP_CPP_COMPILER MATCHES "cpp*")
execute_process(COMMAND "${TMP_CPP_COMPILER}" --version
  OUTPUT_VARIABLE ACTUAL_FC_TARGET
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
  message(STATUS "Found cpp is actually Clang -- Trying other paths")
  find_file (
    TMP_CPP_COMPILER cpp-
    NAMES  cpp-6 cpp6 cpp-5 cpp5 cpp-4.9 cpp
    PATHS ${FORTRAN_PARENT_DIR} /usr/local/bin /opt/local/bin /sw/bin /usr/bin
    #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
    DOC "Searching for GNU cpp preprocessor "
    )
  if(NOT EXISTS "${TMP_CPP_COMPILER}")
    message( FATAL_ERROR  "Cannot find GNU cpp compiler --
${CLANG_FATAL_MSG}")
  endif()
endif()
endif()
set(CMAKE_CPP_COMPILER ${TMP_CPP_COMPILER})

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS} "f03;F03;f08;F08")

if(CMAKE_INSTALL_LIBDIR MATCHES "lib64")
  set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "" FORCE)
endif()



#figure out our git version
option(UPDATE_GIT_VERSION_INFO "update git version info in source tree" ON)
mark_as_advanced(UPDATE_GIT_VERSION_INFO)
if(UPDATE_GIT_VERSION_INFO)
	include(GitInfo)
endif()


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

string(TOUPPER "${CMAKE_Fortran_COMPILER_ID}" ID_STRING)
add_definitions(-D${ID_STRING})
message(STATUS "Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")



################################################################
# Generic Flags
################################################################


# Optimize for the host's architecture
# message(STATUS "Testing flag host arch")
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#   Fortran
#   "-xHost"        # Intel
#   "/QxHost"       # Intel Windows
#   ${GNUNATIVE}    # GNU
#   "-tp=x64"      # Portland Group - generic 64 bit platform
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

###################
### DEBUG FLAGS ###
###################
## Disable optimizations
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
  if (NOT (GFC_VERSION VERSION_GREATER 4.9 OR GFC_VERSION VERSION_EQUAL 4.9))
    message(FATAL_ERROR "${PROJECT_NAME} requires gfortran version 4.9 or above")
  endif ()
  set(EXTRA_FLAGS "${EXTRA_FLAGS}")
  # -fopenmp-simd -mfma -mfma4 -faggressive-function-elimination")
  # -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
  # -Wunused-parameter -Wunused-variable -Wuninitialized ")

  set(CMAKE_CPP_COMPILER_FLAGS           "-E -C -CC -w -Wno-endif-labels -fopenmp") # only seen by preprocessor if #include.*timer is present
  set(CMAKE_Fortran_FLAGS                " ${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG          " ${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${EXTRA_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}" )
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL     "-Os ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        " ${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_FLAGS} ")

  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE_INIT} \
 ${EXTRA_FLAGS} \
-O0 -g -pedantic  -Wextra -Wvector-operation-performance \
-Wopenmp-simd -Wstack-protector -Wpedantic -Wunsafe-loop-optimizations \
-Wshadow -Wsystem-headers -Warray-bounds -Wsuggest-attribute=pure \
-Wsuggest-final-methods  -Wsurprising -Wuse-without-only \
-Wintrinsic-shadow -Wno-unused-dummy-argument")
##  -fcheck-array-temporaries -frange-check -fstack-protector -fstack-check -fbounds-check \
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
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include  -assume source_include -diag-enable=openmp,vec,par,error  -diag-disable=warn -sox -qoverride-limits")
  # -diag-file-append=diagnostics.txt
  set(CMAKE_AR                           "xiar")
  set(CMAKE_CPP_COMPILER                 "fpp")
  set(CMAKE_CPP_COMPILER_FLAGS           " -noJ -B -C ")  #Recognize C,C++,F90 style comments.
  set(CMAKE_Fortran_FLAGS                " ${CMAKE_Fortran_FLAGS_RELEASE_INIT}  ${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG          " ${EXTRA_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${CMAKE_Fortran_FLAGS}")
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL    "-Os ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_RELEASE}  ${EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${EXTRA_FLAGS}")
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
  ## NVIDIA PGI Linux compiler
  #
  #############################################
  message(STATUS "NVIDIA PGI Linux compiler")
  set(PGICOMPILER ON)
  set(USE_LINK_TIME_OPTIMISATION ON)
  if(PGI_CHECKING)
     set (EXTRA_FLAGS "${EXTRA_FLAGS} -Mdclchk -Mchkptr -Mchkstk -Mdepchk -Munixlogical -Mflushz -Mdaz -Mfpmisalign  -Minfo=all,ftn -Mneginfo=all")
  endif()
  if (USE_OPENACC_ONLY)
   set(EXTRA_FLAGS "${EXTRA_FLAGS}  -acc")
    add_definitions(" -DUSE_OPENACC ")
   else()
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -mp")
  endif()
 if(PGI_EXTRA_FAST)
     set (EXTRA_FLAGS "${EXTRA_FLAGS} -Munroll -O4  -fast -Mcuda=fastmath,unroll -Mvect=nosizelimit,short,simd,sse  ")
  endif()

  if(PGI_CUDA_IOMUTEX)
     set (EXTRA_FLAGS "${EXTRA_FLAGS} -Miomutex")
  endif()
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -module ${CMAKE_Fortran_MODULE_DIRECTORY} -I${CMAKE_Fortran_MODULE_DIRECTORY}")
  # NVIDIA PGI Linux compiler
#  set(CMAKE_AR                           "pgfortran")
  set(CMAKE_CPP_COMPILER                 "pgcc -E ")
  set(CMAKE_CPP_COMPILER_FLAGS           "  ")
  set(CMAKE_Fortran_FLAGS                " ${EXTRA_FLAGS} ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG          " ${EXTRA_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG_INIT}  ${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS}")
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "${EXTRA_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-gopt ${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_DEBUG_INIT} -Mneginfo=all")
  set(CMAKE_EXE_LINKER_FLAGS             "${CMAKE_Fortran_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_INIT} -acclibs -cudalibs -Mcudalib=cufft")
  set(CMAKE_SHARED_LINKER_FLAGS          "${CMAKE_SHARED_LINKER_FLAGS_INIT} -acclibs -cudalibs -Mcudalib=cufft")
  set(CMAKE_STATIC_LINKER_FLAGS          "${CMAKE_SHARED_LINKER_FLAGS_INIT} ")

  ################################################################
  # CUDA PGI Default options
  ################################################################
  # default PGI library
  set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)
  set(CUDA_rt_LIBRARY  /usr/lib/x86_64-linux-gnu/librt.so)
  ## PGI 16.10 does not have the OpenMP proc_bind option
  # add_definitions("-Dproc_bind\\(close\\)=\"\"")  # disable proc_bind in PGI  OMP

  if (PGI_LARGE_FILE_SUPPORT)
    set(CMAKE_EXE_LINKER_FLAGS          "-Mlfs ${CMAKE_EXE_LINKER_FLAGS}")
  endif()


#  if (PGI_EXTRACT_ALL)
#    set(CMAKE_Fortran_FLAGS           " -Minline ${CMAKE_Fortran_FLAGS}")
#  endif()
   if (PGI_CUDA_MANAGED_MEMORY)
     set(CMAKE_Fortran_FLAGS           "-ta=nvidia:managed ${CMAKE_Fortran_FLAGS}")
   endif()
   if (PGI_CUDA_PINNED_MEMORY)
     set(CMAKE_Fortran_FLAGS           "-ta=nvidia:pinned ${CMAKE_Fortran_FLAGS}")
   endif()
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


#####################
### RELEASE FLAGS ###
#####################

# Unroll loops
#message(STATUS "Testing flags unroll")
if(USE_AGGRESSIVE_OPTIMISATION)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  Fortran
  "-unroll-aggressive"        # Intel
  "-funroll-all-loops"        # GNU, Intel, Clang
  "/unroll"                   # Intel Windows
  "-Munroll "                 # Portland Group
  )
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  Fortran
  "-inline-level=2"           # Intel
  "-finline-functions"        # GNU, Intel, Clang
  "/unroll"                   # Intel Windows
  "-Minline=maxsize=100,reshape,smallsize=10"      # Portland Group
  )
endif()
if(USE_LINK_TIME_OPTIMISATION)
  # Interprocedural (link-time) optimizations
  # See IPO performance issues https://software.intel.com/en-us/node/694501
  #  may need to use -ipo-separate in case of multi-file problems),  using the -ffat-lto-objects compiler option is provided for GCC compatibility
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-ipo-separate -ffat-lto-object"   # Intel (Linux OS X)
    "/Qipo-separate"                   # Intel Windows
    "-flto "                           # GNU
    "-Mipa=fast,libinline,vestigial,reaggregation"    # Portland Group
    )

  #Profile optimizations
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-ip-dir=.profiling"           # Intel
    "/Qip-dir=_profiling"          # Intel Windows
    "-fprofile-dir=.profiling" # GNU
    "-Mpfo "# PGI
    )
endif(USE_LINK_TIME_OPTIMISATION)
#if(USE_AUTOMATIC_VECTORIZATION)
#   # Interprocedural (link-time) optimizations
#   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#     Fortran
#     "-m"              # Intel ( may need to use -ipo-separate in case of multi-file problems)
#     "/Qipo"             # Intel Windows
#     "-flto "            # GNU
#     "-Mvect"    # Portland Group
#     )
# endif()


# # Vectorize code
if (USE_FAST_MATH_OPTIMISATION)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-fast"             # Intel, PGI
    "/Ofast"     # Intel Windows
    "-Ofast"            # GNU
    )
  # Fast math code
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-fastmath"        # Intel
    "-ffast-math"      # GNU
    "-Mcuda=fastmath"  # Portland Group
    )
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-ffp-contract=fast"              # GNU
    "-Mfma -Mvect=assoc,tile,fuse,gather,simd,partial,prefetch"        # Portland Group
    )   # PGI
endif(USE_FAST_MATH_OPTIMISATION)

if  (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
# Auto parallelize
 if (USE_AUTO_PARALLELISE)
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
     Fortran
     "-parallel -opt-report-phase=par -opt-report:5"            # Intel (Linux)
     "/Qparallel /Qopt-report-phase=par /Qopt-report:5"         # Intel (Windows)
     "-Mconcur=bind,allcores,cncall"                            # PGI
     "-floop-parallelize-all -ftree-parallelize-loops=4"        # GNU
     )
 endif()
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
if(USE_PROFILE_OPTIMISATION)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran
    "-Mpfo"                 # PGI
    "-fpfo "                # GNU
    "-prof-gen"             # Intel (Linux, OS X)
    "/Qprof-gen"            # Intel (Windows)
    )
endif()




if (IMAGE_TYPE_DOUBLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIMAGETYPEDOUBLE" )
endif()

if(APPLE)
  add_definitions("-DMACOSX")
else()
  add_definitions("-DLINUX")
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

if (APPLE AND USE_CUDA)
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
  # https://software.intel.com/en-us/articles/intel-mkl-main-libraries-contain-fftw3-interfaces
  # https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
  set(MKLROOT $ENV{MKLROOT})
  if( NOT "${MKLROOT}" STREQUAL "")
    if(SIMPLE_INTEL_8byte_INTERFACE)
      set(INTEL_INTERFACE "ilp64")  # 8 byte integer interface
    else()
    # default interface
      set(INTEL_INTERFACE "lp64") # 4 byte integer interface
    endif()

    if(BUILD_SHARED_LIBS)
      set(EXTRA_LIBS ${EXTRA_LIBS} -L${MKLROOT}/lib/intel64 -limf -lmkl_intel_${INTEL_INTERFACE} -lmkl_intel_thread -lmkl_core -lmkl_rt -liomp5 -lpthread -lm -ldl )
    else()
     set(EXTRA_LIBS ${EXTRA_LIBS} -Wl,--start-group  ${MKLROOT}/lib/intel64/libmkl_cdft_core.a  ${MKLROOT}/lib/intel64/libmkl_blas95_${INTEL_INTERFACE}.a ${MKLROOT}/lib/intel64/libmkl_intel_${INTEL_INTERFACE}.a ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a  ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L${MKLROOT}/lib/intel64/ -liomp5 -lpthread -lm -ldl)
    endif()
    set(BUILD_NAME "${BUILD_NAME}_MKL" )
  else()
    message( FATAL_ERROR "MKLROOT undefined. Set Intel MKL environment variables by sourcing the relevant shell script (typically in /opt/intel/bin), for example: source /opt/intel/bin/compilervars.sh intel64 lp64; source /opt/intel/bin/compilervars.sh intel64 lp64")
  endif()
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
    include_directories(" ${FFTW_INCLUDE_DIRS}")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES})
  endif()
  set(BUILD_NAME "${BUILD_NAME}_FFTW" )
endif()

#
#  JPEG, TIFF, PNG img support through LibGD
#
if(BUILD_WITH_LIBGD)
  find_package(GD QUIET)
  if (NOT GD_FOUND)
    message(STATUS "Unable to find LibGD imaging library. Please install libgd-dev (https://libgd.github.io) which also depends on libjpeg, libpng and zlib.
DEB:  sudo apt-get install libgd-dev
FINK: fink install gd2-shlibs gd2
BREW: brew install gd
COMPILE FROM SOURCE: git clone https://github.com/libgd/libgd
Warning -- libgd will pull latest stable libjpeg that may conflict with libjpeg9.
")
    set(BUILD_WITH_LIBGD OFF)
  else()
    message(STATUS "LibGD found")
    message(STATUS "lib: ${GD_LIBRARIES}")
    add_definitions(" -D_LIBGD ")
    include_directories(" ${GD_INCLUDE_DIRS}")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${GD_LIBRARIES} ${ZLIB_LIBRARY_RELEASE})
    set(EXTRA_LIBS ${EXTRA_LIBS} -ljpeg)
  endif()
endif()


#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  ")\
#set(CMAKE_FCPP_FLAGS " -C -P ") # Retain comments due to fortran slash-slash
#set(CMAKE_Fortran_CREATE_PREPROCESSED_SOURCE "${CMAKE_FCPP_COMPILER} <DEFINES> <INCLUDES> <FLAGS> -E <SOURCE> > <PREPROCESSED_SOURCE>")

# add_definitions(" -D__FILENAME__='\"$(notdir $<)\"' ")
# add_definitions(" -DHALT\\(X\\)='call simple_stop(X, __FILENAME__, __LINE__)'")
# add_definitions(" -Dallocchk\\(X\\)='call alloc_errchk(X,alloc_stat,\"$(notdir $<)\",__LINE__);' ")


# Override Fortran preprocessor
# block constructs (F2008), unlimited polymorphism and variadic macros (not included in F2003 -- but is part of C99 )
 if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")

 set(CMAKE_Fortran_COMPILE_OBJECT "grep --silent -E '#include.*timer' <SOURCE> && ( ${CMAKE_CPP_COMPILER} ${CMAKE_CPP_COMPILER_FLAGS} -DOPENMP <DEFINES> <INCLUDES> <SOURCE> > <OBJECT>.f90 &&  <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <OBJECT>.f90 -o <OBJECT> ) ||  <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <SOURCE> -o <OBJECT>")
 else()
# #elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
 set(CMAKE_Fortran_COMPILE_OBJECT "grep --silent -E '#include.*timer' <SOURCE> && ( ${CMAKE_CPP_COMPILER} ${CMAKE_CPP_COMPILER_FLAGS} -DOPENMP <DEFINES> <INCLUDES> <SOURCE> > <OBJECT>.f90 &&  <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <OBJECT>.f90 -o <OBJECT> ) || <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <SOURCE> -o <OBJECT>")
 endif()


if(USE_PROFILING)

  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg -fno-omit-frame-pointer")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -fno-omit-frame-pointer")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -fno-omit-frame-pointer")
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
      string(REGEX REPLACE "-f[pP][iI][cC]" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  endif()
endif()


# Option for code coverage
#if(VERBOSE OR ${BUILD_WITH_COVERAGE})
#  option(USE_CODE_COVERAGE "Build code coverage results, requires GCC compiler (forces Debug build)" OFF)
  if(USE_CODE_COVERAGE)
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" )
      set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -g -Og -pg -fprofile-arcs -ftest-coverage")
      set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS}       -g -Og -pg -fprofile-arcs -ftest-coverage")
      set(CMAKE_CXX_FLAGS       "${CMAKE_CXX_FLAGS}     -g -Og -pg -fprofile-arcs -ftest-coverage")
     # set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "" FORCE)
      SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -g -Og -pg -fprofile-arcs -ftest-coverage")
      SET(CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -Og -g -pg -fprofile-arcs -ftest-coverage")
      # Ensure that CDash targets are always enabled if coverage is enabled.
      if (NOT CDASH_SUPPORT)
        get_property(HELP_STRING CACHE CDASH_SUPPORT PROPERTY HELPSTRING)
        set(CDASH_SUPPORT ON CACHE BOOL "${HELP_STRING}" FORCE)
        message(STATUS "Enabling CDash targets as coverage has been enabled.")
      endif()
    endif()
  endif()
#endif()


set(CMAKE_INSTALL_DO_STRIP FALSE)
if(APPLE)
  #https://cmake.org/Wiki/CMake_RPATH_handling
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF("${isSystemDir}" STREQUAL "-1")

elseif(UNIX)
  set(CMAKE_INSTALL_DO_STRIP FALSE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

endif()
