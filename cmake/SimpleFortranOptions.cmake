#------------------------------------------------------------------------------!
# SIMPLE              Elmlund & Elmlund Lab          simplecryoem.com          !
#------------------------------------------------------------------------------!

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
include(SetCompileFlag)
include(ZSetParallelLibrary)

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
Set FC,CC, and CPP environment variables to GNU compilers, e.g. using Fink in bash :
export FC=/sw/bin/gfortran
export CC=/sw/bin/gcc
export CXX=/sw/bin/g++
export CPP=/sw/bin/cpp-5

Clang has overridden Gfortran links /usr/bin/gfortran GRRR!
Make sure you prepend  /usr/local/bin (Homebrew) or /opt/local/bin (MacPorts) or /sw/bin (FINK) to PATH.
In LD_LIBRARY_PATH prepend the appropriate lib path.
    ")

option(USE_GNU_EXTENSIONS "Enable GNU extensions in C " OFF)
if(USE_GNU_EXTENSIONS)
  set(C_DIALECT gnu)
else()
  set(C_DIALECT c)
endif()

message(STATUS "Making sure your MPI Fortran compiler enables the USE_MPI flag")
if(Fortran_COMPILER_NAME MATCHES "mpi*")
  set(USE_MPI ON)
endif()

message(STATUS "Making sure your Fortran compiler points to the correct binary")
if(Fortran_COMPILER_NAME MATCHES "gfortran*")
message(STATUS "Making sure your Fortran compiler points to the correct binary")
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_FC_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "GFORTRAN version: ${ACTUAL_FC_TARGET}")
  if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
    message(STATUS "WARNING gfortran points to Clang -- Trying other paths")
    find_file (
      CMAKE_Fortran_COMPILER
      NAMES gfortran- gfortran-13 gfortran-12 gfortran-11 gfotran-10 gfortran-9 gfortran-8 gfortran-7 gfortran-6 gfortran-5 gfortran-4.9 gfortran13 gfortran12 gfortran11 gfortran10 gfortran9 gfortran8 gfortran7 gfortran6 gfortran5
      PATHS /usr/local/bin /opt/local/bin /sw/bin /opt/homebrew/bin /usr/bin
      #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
      DOC "Searching for GNU gfortran preprocessor "
      )
    if(NOT EXISTS "${CMAKE_Fortran_COMPILER}")
      message( FATAL_ERROR  "Cannot find ${CMAKE_Fortran_COMPILER} --
${CLANG_FATAL_MSG}")
    endif()
    message(STATUS "GFORTRAN replaced with ${CMAKE_Fortran_COMPILER")

  endif()
endif()
get_filename_component(GFORTRAN_ABSPATH ${CMAKE_Fortran_COMPILER} REALPATH)
message(STATUS "FORTRAN ABSPATH ${GFORTRAN_ABSPATH}")
get_filename_component(FORTRAN_PARENT_DIR ${GFORTRAN_ABSPATH} DIRECTORY)
message(STATUS "FORTRAN_PARENT_DIR ${FORTRAN_PARENT_DIR}")

message(STATUS "Making sure your C compiler points to the correct binary")

  execute_process(COMMAND ${CMAKE_C_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_C_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GCC version: ${ACTUAL_C_TARGET}")

  if(ACTUAL_C_TARGET MATCHES "Clang|clang")
    message(STATUS "WARNING gcc points to Clang -- Attempting other paths, starting with ${FORTRAN_PARENT_DIR}")
    find_file (
      CMAKE_C_COMPILER_NEW
      NAMES gcc- gcc-13 gcc-12 gcc-11 gcc-10 gcc-9 gcc-8 gcc-7 gcc-6 gcc-5 gcc-4.9 gcc-fsf-6 gcc-fsf-5 gcc13 gcc12 gcc11 gcc10 gcc9 gcc8 gcc7 gcc6 gcc5 gcc4.9
      HINTS ${FORTRAN_PARENT_DIR}
      PATHS  /sw/bin /usr/local/bin /opt/local/bin /opt/homebrew/bin /usr/bin
      #  [PATH_SUFFIXES suffix1 [suffix2 ...]]
      DOC "Searching for GNU gcc preprocessor, starting with ${FORTRAN_PARENT_DIR} "
      NO_CMAKE_PATH NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH
)
    message (STATUS " Found ${CMAKE_C_COMPILER_NEW}")
    if(EXISTS "${CMAKE_C_COMPILER_NEW}")
    set(CMAKE_C_COMPILER "${CMAKE_C_COMPILER_NEW}" CACHE FILEPATH "GNU gcc compiler " FORCE)
    message(STATUS "C compiler points to ${CMAKE_C_COMPILER}")
else()
      message( FATAL_ERROR  "Cannot find GNU ${CMAKE_C_COMPILER_NEW} --
${CLANG_FATAL_MSG}")
    endif()
  endif()

message(STATUS "Making sure your C++ compiler points to the correct binary")

  execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version
    OUTPUT_VARIABLE ACTUAL_CXX_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "G++ version: ${ACTUAL_CXX_TARGET}")
  if(ACTUAL_CXX_TARGET MATCHES "Clang|clang")
    message(STATUS "WARNING g++ points to Clang -- Trying other paths")
    find_file (
      CMAKE_CXX_COMPILER_NEW
      NAMES g++- g++-13 g++-12 g++-11 g++-10 g++-9 g++-8 g++-7 g++-6 g++-5 g++-4.9 g++-fsf-6 g++-fsf-5 g++13 g++12 g++11 g++10 g++9 g++8 g++7 g++6 g++5 g++4.9 g++
      PATHS ${FORTRAN_PARENT_DIR} /sw/bin /usr/local/bin /opt/local/bin /opt/homebrew/bin /usr/bin
      DOC "Searching for GNU g++ preprocessor "
NO_DEFAULT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH
      )
    if(EXISTS "${CMAKE_CXX_COMPILER_NEW}")
    set(CMAKE_CXX_COMPILER "${CMAKE_CXX_COMPILER_NEW}" CACHE FILEPATH "GNU C++ COMPILER" FORCE)
    message(STATUS "C++ compiler points to ${CMAKE_CXX_COMPILER}")
else()
      message( FATAL_ERROR  "Cannot find GNU ${CMAKE_CXX_COMPILER} --
${CLANG_FATAL_MSG}")
    endif()
  endif()


message(STATUS "Making sure your preprocessor points to the correct binary")
if(TMP_CPP_COMPILER MATCHES "cpp*")
execute_process(COMMAND "${TMP_CPP_COMPILER}" --version
  OUTPUT_VARIABLE ACTUAL_FC_TARGET
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(ACTUAL_FC_TARGET MATCHES "Clang|clang")
  message(STATUS "Found cpp is actually Clang -- Trying other paths")
  find_file (
    TMP_CPP_COMPILER cpp-
    NAMES  cpp-10 cpp-9 cpp-8 cpp-7 cpp-6 cpp-5 cpp-4.9 cpp10 cpp9 cpp8 cpp7 cpp6 cpp5 cpp4.9 cpp
    PATHS ${FORTRAN_PARENT_DIR} /sw/bin /usr/local/bin /opt/local/bin /usr/bin
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
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE_INIT: ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")

################################################################
# Json-fortran definitions
################################################################
#add_definitions("-DJF_REAL32=4")
#add_definitions("-DJF_INT32=4")
if(CMAKE_Fortran_COMPILER_SUPPORTS_USC4 EQUAL 1)
  add_definitions("-DUSE_USC4=1")
endif()

################################################################
# fortran version
################################################################

if(CMAKE_Fortran_COMPILER_SUPPORTS_F08_ISOENV EQUAL 1)
  add_definitions("-DUSE_F08_ISOENV=1")
endif()
if(CMAKE_Fortran_COMPILER_SUPPORTS_F08 EQUAL 1)
  add_definitions("-DUSE_F08=1")
endif()
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
if(NOT BUILD_SHARED_LIBS)
  message(STATUS "Build static executables.")
  message(STATUS "On glibc-based systems, OpenMP enabled applications cannot be statically linked due to limitations of the underlying pthreads-implementation. It might be possible to get a working solution if -Wl,--whole-archive -lpthread -Wl,--no-whole-archive is added to the command line. However, this is not supported by gcc and thus not recommended.")
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  add_definitions("-fPIE -fPIC")
endif()


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
#if(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
#else()
if (USE_OPENACC)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
    Fortran REQUIRED
    "-fopenacc"            # GNU
    "-acc"
    "/acc"
    )
  add_definitions("-DOPENACC")  ## FIXME above
endif()
#endif()
#include(FindOpenMP_Fortran)


if (USE_OPENMP)
 # message(STATUS "in USE_OPENMP : ${CMAKE_Fortran_FLAGS_RELEASE}")
  # SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  #      Fortran REQUIRED
  #      "-fopenmp"            # GNU
  #      "-qopenmp"            # Intel
  #      "-openmp"             # depreciated Intel Open MP flag
  #      "/mp"
  #      )
  #    message(STATUS "leaving USE_OPENMP : ${CMAKE_Fortran_FLAGS_RELEASE}")
#   endif()
# endif()

#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DOPENMP ")
  add_definitions("-DOPENMP")  ## FIXME above
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




if(USE_MPI)
  find_package(MPI REQUIRED)
 # if((${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
 #     if(MPI_Fortran_FOUND)
        # do something different for Intel MPI

 #     else()

  if(MPI_Fortran_FOUND)
    message(STATUS "MPI found -- ${MPI_COMPILE_FLAGS}")
    message(STATUS "MPI_Fortran_COMPILER ${MPI_Fortran_COMPILER}")
    message(STATUS "MPI_Fortran_COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}")
    message(STATUS "MPI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES}")
    message(STATUS "MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH}")
    message(STATUS "MPI_Fortran_LINK_FLAGS ${MPI_Fortran_LINK_FLAGS}")
    message(STATUS "MPI_Fortran mpi_f08_mod_path  ${mpi_f08_mod_path}")
    # message(STATUS "MPI_C_COMPILER ${MPI_C_COMPILER}")
    # message(STATUS "MPI_C_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS}")
    # message(STATUS "MPI_C_LIBRARIES ${MPI_C_LIBRARIES}")
    # message(STATUS "MPI_C_INCLUDE_PATH ${MPI_C_INCLUDE_PATH}")
    # message(STATUS "MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS}")

    # message(STATUS "MPI_CXX_COMPILER ${MPI_CXX_COMPILER}")
    # message(STATUS "MPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS}")
    # message(STATUS "MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES}")
    # message(STATUS "MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH}")
    # message(STATUS "MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS}")

    # For backward compatibility with older versions of FindMPI, these
    # variables are set, but deprecated:
    message(STATUS "MPI_LIBRARY ${MPI_LIBRARY}")
    message(STATUS "MPI_COMPILE_FLAGS ${MPI_COMPILE_FLAGS}")
    message(STATUS "MPI_EXTRA_LIBRARY ${MPI_EXTRA_LIBRARY}")
    message(STATUS "MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH}")
    message(STATUS "MPI_CXX_LINK_FLAGS ${MPI_LINK_FLAGS}")

   # get_filename_component(mpi_f08_mod_path ${MPI_EXTRA_LIBRARY} DIRECTORY)
   # include_directories(SYSTEM ${mpi_f08_mod_path})

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    #-I${libdir} -I${includedir} -I${includedir}/openmpi/opal/mca/event/libevent2021/libevent -I${includedir}/openmpi/opal/mca/event/libevent2021/libevent/include   -pthread
    foreach(mpi_inc_paths ${MPI_Fortran_INCLUDE_PATH})
      include_directories(SYSTEM ${mpi_inc_paths})
      #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -isystem ${mpi_inc_paths}")
    endforeach()
    foreach(mpi_mod_paths ${mpi_f08_mod_path})
      link_directories(${mpi_mod_paths})
      include_directories(SYSTEM ${mpi_mod_paths})
    endforeach()
    foreach(mpi_exe_flags ${MPI_COMPILE_FLAGS}  ${MPI_LINK_FLAGS} ${MPI_LIBRARIES})
      #/ME  strip udev from linker list, this allows for full static compilation
      if( NOT ${mpi_exe_flags} STREQUAL "-ludev" )
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${mpi_exe_flags}")
      endif()
    endforeach()
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  ${MPI_LINK_FLAGS}")
   # set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
    add_definitions("-DUSING_MPI")
  else(MPI_Fortran_FOUND)
    message(FATAL_ERROR "Unable to find MPI -- mpif90 and libraries were not found in system paths")
  endif(MPI_Fortran_FOUND)

endif(USE_MPI)


if(USE_CUDA)
  if (APPLE)
    message(STATUS "Apple MacOSX CUDA unsupported")
    set(USE_CUDA OFF)
    #     set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    #     set(CMAKE_SHARED_LIBRARY_CXX_FLAGS_BACKUP "${CMAKE_SHARED_LIBRARY_CXX_FLAGS}")
    #     set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -Wno-unused-function")
    #     string(REGEX REPLACE "-fopenmp[^ ]*" "" CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS}")
  else()
    find_package(Threads)
    find_package(CUDA)

    if(CUDA_FOUND)
      message(STATUS "CUDA Found ")
      set(EXTRA_LIBS ${EXTRA_LIBS}  ${CUDA_LIBRARIES})
      include_directories( ${CUDA_INCLUDE_DIRS} )
      add_definitions("-DUSING_CUDA")
      if(NOT BUILD_SHARED_LIBS)
        set(CUDA_USE_STATIC_CUDA_RUNTIME ON)
      endif()
      set(CUDA_PROPAGATE_HOST_FLAGS OFF)
      CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/cuda/kernels)
      set(CUDA_HOST_COMPILER /usr/bin/gcc )
      CUDA_INCLUDE_DIRECTORIES(/usr/include)

      if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        set(CUDA_VERBOSE_BUILD ON)
        set(CUDA_NVCC_FLAGS -g ${CUDA_NVCC_FLAGS})
      endif()

      if( CUDA_TOOLKIT_ROOT_DIR )
        CUDA_INCLUDE_DIRECTORIES( ${CUDA_TOOLKIT_ROOT_DIR} )
      endif()

    else()
      set(USE_CUDA OFF)
      unset (CUDA_INCLUDE_DIRS CACHE)
      unset (CUDA_TOOLKIT_ROOT_DIR CACHE)
      unset (CUDA_LIBRARIES CACHE)
      #   message(FATAL_ERROR "Unable to find CUDA -- set CUDA_TOOLKIT_ROOT_DIR in shell environment")
    endif()
  endif(APPLE)
endif(USE_CUDA)

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
    else() #STATIC
     set(EXTRA_LIBS ${EXTRA_LIBS} -Wl,--start-group  ${MKLROOT}/lib/intel64/libmkl_cdft_core.a  ${MKLROOT}/lib/intel64/libmkl_blas95_${INTEL_INTERFACE}.a ${MKLROOT}/lib/intel64/libmkl_intel_${INTEL_INTERFACE}.a ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a  ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group ${INTEL_DIR}/lib/intel64/libiomp5.a  -lm -ldl)
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

#############################################
# TIFF library & dependencies
#############################################
if(USE_LIBTIFF)
  set(TIFFLIBS "0")
  set(EXTRA_TIFF_DEPENDENCIES "0")
  find_package(TIFF 4.3.0 QUIET)
  if( TIFF_FOUND )
    add_definitions("-DTIFF430=1")
  endif()
  find_package(TIFF 4.0.10)
  if( TIFF_FOUND )
    set(EXTRA_TIFF_DEPENDENCIES "1")
  else()
    find_package(TIFF)
  endif()
  if(TIFF_FOUND)
    set(TIFFLIBS "1")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${TIFF_LIBRARY} )
    # JBIG
    find_library(JBIG_LIBRARY jbig)
    message("-- JBIG_LIBRARY: ${JBIG_LIBRARY}")
    if( JBIG_LIBRARY STREQUAL "JBIG_LIBRARY-NOTFOUND" )
      set(TIFFLIBS "0")
    else()
      set(EXTRA_LIBS ${EXTRA_LIBS} ${JBIG_LIBRARY} )
    endif()
    # JPEG
    find_library(JPEG_LIBRARY jpeg)
    message("-- JPEG_LIBRARY: ${JPEG_LIBRARY}")
    if( JPEG_LIBRARY STREQUAL "JPEG_LIBRARY-NOTFOUND" )
      set(TIFFLIBS "0")
    else()
      set(EXTRA_LIBS ${EXTRA_LIBS} ${JPEG_LIBRARY} )
    endif()
    # LZMA
    find_library(LZMA_LIBRARY lzma)
    message("-- LZMA_LIBRARY: ${LZMA_LIBRARY}")
    if( LZMA_LIBRARY STREQUAL "LZMA_LIBRARY-NOTFOUND" )
      set(TIFFLIBS "0")
    else()
      set(EXTRA_LIBS ${EXTRA_LIBS} ${LZMA_LIBRARY} )
    endif()
    # ZLIB
    find_library(Z_LIBRARY z)
    message("-- Z_LIBRARY: ${Z_LIBRARY}")
    if( Z_LIBRARY STREQUAL "Z_LIBRARY-NOTFOUND" )
      set(TIFFLIBS "0")
    else()
      set(EXTRA_LIBS ${EXTRA_LIBS} ${Z_LIBRARY} )
    endif()
    if( EXTRA_TIFF_DEPENDENCIES STREQUAL "1")
      if(APPLE)
        # Does not seem required for macos 10+
      else()
        # ZSTDLIB
        find_library(ZSTD_LIBRARY zstd)
        message("-- ZSTD_LIBRARY: ${ZSTD_LIBRARY}")
        if( ZSTD_LIBRARY STREQUAL "ZSTD_LIBRARY-NOTFOUND" )
          set(TIFFLIBS "0")
        else()
          set(EXTRA_LIBS ${EXTRA_LIBS} ${ZSTD_LIBRARY} )
        endif()
        # WEBP
        find_library(WEBP_LIBRARY webp)
        message("-- WEBP_LIBRARY: ${WEBP_LIBRARY}")
        if( WEBP_LIBRARY STREQUAL "WEBP_LIBRARY-NOTFOUND" )
          set(TIFFLIBS "0")
        else()
          set(EXTRA_LIBS ${EXTRA_LIBS} ${WEBP_LIBRARY} )
        endif()
        # Deflate
        find_package(TIFF 4.2.0 QUIET)
        if( TIFF_FOUND )
          find_library(Deflate_LIBRARY deflate)
          message("-- Deflate_LIBRARY: ${Deflate_LIBRARY}")
          if( Deflate_LIBRARY STREQUAL "Deflate_LIBRARY-NOTFOUND" )
            set(TIFFLIBS "0")
          else()
            set(EXTRA_LIBS ${EXTRA_LIBS} ${Deflate_LIBRARY} )
          endif()
        endif()
      endif()
    endif()
  endif()
  if( TIFFLIBS STREQUAL "1")
    add_definitions("-DUSING_TIFF=1")
    include_directories(${TIFF_INCLUDE_DIR})
    set(EXTRA_LIBS ${EXTRA_LIBS} -lm )
  else()
    message(STATUS "ERROR: Dependencies for TIFF support were not all found (see above)")
    message(STATUS "ERROR: Install dependencies or check for missing soft links and INCLUDE directories")
    message(FATAL_ERROR "CMake configuration cannot proceed")
  endif()
endif()

#############################################
## COMPLER SPECIFIC SETTINGS
#############################################
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*" ) #OR Fortran_COMPILER_NAME MATCHES "mpif(90|ort)")
  #############################################
  #
  ## GNU fortran
  #
  #############################################

  message(STATUS "Executing: ${CMAKE_Fortran_COMPILER} -dumpversion")
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
    OUTPUT_VARIABLE GFC_VERSION
    ERROR_VARIABLE GFC_VERSION_STDERR
    RESULT_VARIABLE GFC_VERSION_RESULT)
 # message(STATUS "${CMAKE_Fortran_COMPILER} -dumpversion :  ${GFC_VERSION} : ${GFC_VERSION_RESULT} : ${GFC_VERSION_STDERR}")
  if( RESULT_VARIABLE EQUAL 1 )
    message(FATAL_ERROR "${PROJECT_NAME} requires ${CMAKE_Fortran_COMPILER} to compile")
  endif()
  if (NOT (GFC_VERSION VERSION_GREATER 4.9 OR GFC_VERSION VERSION_EQUAL 4.9))
    message(FATAL_ERROR "${PROJECT_NAME} requires gfortran version 4.9 or above")
  endif ()
  set(EXTRA_FLAGS "${EXTRA_FLAGS}")
  # -fopenmp-simd -mfma -mfma4 -faggressive-function-elimination")
  # -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
  # -Wunused-parameter -Wunused-variable -Wuninitialized ")

  set(CMAKE_CPP_COMPILER_FLAGS           "-E -C -CC -w -Wno-endif-labels -fopenmp") # only seen by preprocessor if #include.*timer is present
  set(CMAKE_Fortran_FLAGS                "${EXTRA_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_FLAGS} " )
  set(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_FLAGS} ")

  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE_INIT} \
 ${EXTRA_FLAGS} \
-O0 -g -pedantic  -Wcompare-reals -Wvector-operation-performance \
-Wopenmp-simd -Wstack-protector -Wpedantic -Wunsafe-loop-optimizations \
-Wshadow -Wsystem-headers -Warray-bounds -Wsuggest-attribute=pure \
-Wsuggest-final-methods  -Wsurprising -Wno-use-without-only \
-Wintrinsic-shadow -Wno-unused-dummy-argument -Wno-unused-variable")
##  -fcheck-array-temporaries -frange-check -fstack-protector -fstack-check -fbounds-check \
  # #  CMAKE_EXE_LINKER_FLAGS
  # if (LINK_TIME_OPTIMISATION)
  #   set(CMAKE_EXE_LINKER_FLAGS             "${CMAKE_EXE_LINKER_FLAGS_INIT} -flto ")
  #   set(CMAKE_SHARED_LINKER_FLAGS           "${CMAKE_EXE_LINKER_FLAGS_INIT} -flto -flto=${NUM_JOBS}")
  # endif(LINK_TIME_OPTIMISATION)

  if(GFORTRAN_EXTRA_CHECKING)
    set(CMAKE_Fortran_FLAGS "-frange-check -fstack-protector -fstack-check -fbounds-check -fcheck-array-temporaries")
  endif()


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

  set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE_INIT}")
  set(CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG_INIT}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG_INIT}")
  #if(BUILD_SHARED_LIBS STREQUAL "OFF")
  #  set(CMAKE_EXE_LINKER_FLAGS   "${CMAKE_EXE_LINKER_FLAGS} -static ")
  #  set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -Wl,static-libgfortran")
 #   set(LARGE_FILE_SUPPORT OFF)
 # endif()
  if(LARGE_FILE_SUPPORT)
      if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64") 
          set(CMAKE_Fortran_FLAGS             "${CMAKE_Fortran_FLAGS} -mcmodel=small")
     else()
          set(CMAKE_Fortran_FLAGS             "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
     endif()
  endif()

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*" OR Fortran_COMPILER_NAME MATCHES "mpiifort*")
  #############################################
  #
  ## INTEL fortran
  #
  #############################################
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include/  -assume source_include -diag-enable=openmp,vec,par,error   -I${INTEL_DIR}/include/${INTEL_TARGET_ARCH}")
  if(INTEL_OMP_OVERRIDE)
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -sox -qoverride-limits  -diag-disable=warn")
endif()
  # -diag-file-append=diagnostics.txt
  set(CMAKE_AR                           "xiar")
  set(CMAKE_CPP_COMPILER                 "fpp")
  set(CMAKE_CPP_COMPILER_FLAGS           " -noJ -B -C ")  #Recognize C,C++,F90 style comments.
  set(CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${CMAKE_Fortran_FLAGS_DEBUG}")
  # set(CMAKE_Fortran_FLAGS_MINSIZEREL    "-Os ${CMAKE_Fortran_FLAGS_RELEASE_INIT}")
  set(CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_RELEASE} ")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE_INIT} ${CMAKE_Fortran_FLAGS_DEBUG_INIT} ${EXTRA_FLAGS}")
  # set(CMAKE_EXE_LINKER_FLAGS       "${CMAKE_EXE_LINKER_FLAGS_INIT}")
  # set(CMAKE_SHARED_LINKER_FLAGS    "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${EXTRA_LIBS}")
  # set(CMAKE_STATIC_LINKER_FLAGS    "${CMAKE_STATIC_LINKER_FLAGS_INIT} ${EXTRA_LIBS}")
  # if (LINK_TIME_OPTIMISATION)
  #   set(CMAKE_EXE_LINKER_FLAGS     "${CMAKE_EXE_LINKER_FLAGS} -ipo-separate -ipo-jobs=${NUM_JOBS}")
  #   set(CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -ip -ipo-separate -ipo-jobs=${NUM_JOBS}")
  #   set(CMAKE_STATIC_LINKER_FLAGS  "${CMAKE_STATIC_LINKER_FLAGS} -ip -ipo")
  # endif(LINK_TIME_OPTIMISATION)

  set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE_INIT}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
  if(NOT BUILD_SHARED_LIBS)
    set(CMAKE_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} ")
    if(Fortran_COMPILER_NAME MATCHES "mpiif*")
      set(CMAKE_EXE_LINKER_FLAGS       "${CMAKE_EXE_LINKER_FLAGS} -static_mpi")
    endif()
    if(USING_OPENMP)
      set(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} ${INTEL_DIR}/lib/${INTEL_TARGET_ARCH}/libiomp5.a")
    endif()
    set(LARGE_FILE_SUPPORT OFF) # disable mcmodel in static mode
  else()
    set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -traceback -shared-intel ")
  endif()
  if(LARGE_FILE_SUPPORT)
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64") 
        set(CMAKE_Fortran_FLAGS             "${CMAKE_Fortran_FLAGS} -mcmodel=small")
    else()
        set(CMAKE_Fortran_FLAGS             "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
    endif()
  endif()
  message(STATUS "Intel CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ")
  ## end Intel compilation section

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
     set (EXTRA_FLAGS "${EXTRA_FLAGS}")
  endif()
  if (USE_OPENACC_ONLY)
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -acc")
    add_definitions(" -DUSE_OPENACC ")
   else()
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -mp")
  endif()
 if(PGI_EXTRA_FAST)
     set (EXTRA_FLAGS "${EXTRA_FLAGS} -Munroll -O4  -fast -Mcuda=fastmath,unroll -Mvect=nosizelimit,short,simd,sse  ")
  endif()

  if(PGI_CUDA_IOMUTEX)
     set(CMAKE_Fortran_FLAGS  " ${CMAKE_Fortran_FLAGS} -Miomutex")
  endif()
  set(CMAKE_Fortran_FLAGS  " ${CMAKE_Fortran_FLAGS} -module ${CMAKE_Fortran_MODULE_DIRECTORY} -I${CMAKE_Fortran_MODULE_DIRECTORY}")
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


if (IMAGE_TYPE_DOUBLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIMAGETYPEDOUBLE" )
endif()

if(APPLE)
  add_definitions("-DMACOSX")
else()
  add_definitions("-DLINUX")
endif()


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
    set(CMAKE_Fortran_COMPILE_OBJECT "grep --silent -E '#include.*timer' <SOURCE> && ( ${CMAKE_CPP_COMPILER} ${CMAKE_CPP_COMPILER_FLAGS} $ENV{FFLAGS} -DOPENMP <DEFINES> <INCLUDES> <SOURCE> > <OBJECT>.f90 &&  <CMAKE_Fortran_COMPILER> $ENV{FFLAGS} <DEFINES> <INCLUDES> <FLAGS> -c <OBJECT>.f90 -o <OBJECT> ) || <CMAKE_Fortran_COMPILER> $ENV{FFLAGS} <DEFINES> <INCLUDES> <FLAGS> -c <SOURCE> -o <OBJECT> ")
    # uncomment below to show compilation times
    #set(CMAKE_Fortran_COMPILE_OBJECT "grep --silent -E '#include.*timer' <SOURCE> && ( ${CMAKE_CPP_COMPILER} ${CMAKE_CPP_COMPILER_FLAGS} -DOPENMP <DEFINES> <INCLUDES> <SOURCE> > <OBJECT>.f90 &&  <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <OBJECT>.f90 -o <OBJECT> ) ||(TIME='PTIME %C %E' time <CMAKE_Fortran_COMPILER> <DEFINES> <INCLUDES> <FLAGS> -c <SOURCE> -o <OBJECT> 2>&1 | awk '/PTIME/ {cmd=\"basename -- \" \$\$\(NF-1\)\;cmd|getline out\;print \$\$\(NF-1\),\$\$NF\;close(cmd)\;}')")
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
if (COVERALLS)
    include(Coveralls)
    set(USE_CODE_COVERAGE ON)
endif()
if(USE_CODE_COVERAGE)
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" )
    set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -g -Og -pg -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS         "${CMAKE_C_FLAGS}       -g -Og -pg -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS       "${CMAKE_CXX_FLAGS}     -g -Og -pg -fprofile-arcs -ftest-coverage")
    # set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "" FORCE)
    SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -g -Og -pg -fprofile-arcs -ftest-coverage -lgcov --coverage")

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

  # when building, don't use the install RPATH
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


set (GRAPHVIZ_GRAPH_TYPE "digraph")
set (GRAPHVIZ_GRAPH_NAME "GG")
set (GRAPHVIZ_GRAPH_HEADER "node [n fontsize = \"8.0\"];")
set (GRAPHVIZ_NODE_PREFIX  "n")
set (GRAPHVIZ_EXECUTABLES TRUE)
set (GRAPHVIZ_STATIC_LIBS TRUE)
set (GRAPHVIZ_SHARED_LIBS TRUE)
set (GRAPHVIZ_MODULE_LIBS TRUE)
set (GRAPHVIZ_EXTERNAL_LIBS FALSE )
set (GRAPHVIZ_IGNORE_TARGETS ".*(.build|.required|.proxy)")
set (GRAPHVIZ_GENERATE_PER_TARGET FALSE)
set (GRAPHVIZ_GENERATE_DEPENDERS FALSE)
