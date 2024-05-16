#Look for MKL_FFT libraries ( amd FFTE for Fortran in the future)
find_package(MKL)
${MKLROOT}
MKLLIB /lib/intel64
module load comp-intel/2018.0.128

# Make sure that the MKL wrapper is compatible with FFTW version
#LIB_DIRS += -L/opt/intel/mkl/$(MKL_VER)/lib/$(MKL_ARCH)
#INC_DIRS += -I/opt/intel/mkl/$(MKL_VER)/include
#INC_DIRS += -I/opt/intel/mkl/$(MKL_VER)/include/fftw
#LIB_FILES += -lpthread -lfftw3xc_intel -lmkl -lm

#ifort -O3                                                 \
#      -o fftw_example.exe fftw_example.f                  \
#      -I/nasa/intel/Compiler/2018.0.128/mkl/include/fftw  \
#      -L/nasa/intel/Compiler/2018.0.128/mkl/lib/intel64   \
#      -lfftw2xf_single_intel                              \
#      -mkl
# https://www.nas.nasa.gov/hecc/support/kb/mkl-fftw-interface_204.html
# - Find MKL_FFT
# Find the Intel MKL_FFT includes and library
# It uses the FFTW wrappers that translater the FFTW3 calls to their corresponding MKL_FFT implementations

find_path(MKL_FFT_INCLUDE_DIRS
MKL_FFT_LIBRARIES MKL_FFT_INCLUDE_DIRS )
    NAMES fftw3.f03
    HINTS
        ${MKLROOT}/include/fftw
    PATHS
        ${MRLROOT}/include/fftw
)
mark_as_advanced(MKL_FFT_INCLUDE_DIRS)

find_library(MKL_FFT_SINGLE_PRECISION_LIBRARIES
    NAMES fftw3f libfftw3f libfftw3f-3
    HINTS
        ${MKLROOT}share/mkl/interfaces/fftw2xf/wrappers/
# Interfaces
/mkl/2024.1/share/mkl/interfaces
        ${FFTW_ROOT}/lib
        ${FFTW_ROOT}/.libs
        ${FFTW_ROOT}
        ${FFTWDIR}/lib
        $ENV{FFTW_ROOT}/lib
        $ENV{FFTW_ROOT}/.libs
        $ENV{FFTWLIB}
        $ENV{FFTWDIR}/lib
        ENV FFTW_ROOT
        ENV FFTWLIB
    PATHS
        ${FFTWDIR}/lib
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib               # Homebrew
        /opt/local/lib               # Macports
        /usr/opt/local/lib
        /sw/lib                      # Fink
    DOC "FFTW dynamic library -- single precision, serial"
)
mark_as_advanced(FFTW_SINGLE_PRECISION_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( MKL_FFT DEFAULT_MSG MKL_FFT_LIBRARIES MKL_FFT_INCLUDE_DIRS )

if(NOT MKL_FFT_FOUND)
    message( STATUS "Error MKL_FFT not found")
    message( STATUS "FindMKL_FFT looked for single precision libraries -- :  ${FFTW_SINGLE_PRECISION_LIBRARIES}" )
    message( STATUS "FindMKL_FFT looked for double precision libraries -- :  ${FFTW_DOUBLE_PRECISION_LIBRARIES}" )
    message( STATUS "FindMKL_FFT looked for single precision threaded libraries -- :  ${FFTW_SINGLE_PRECISION_THREADED_LIBRARIES}" )
    message( STATUS "FindMKL_FFT looked for double precision threaded libraries -- :  ${FFTW_DOUBLE_PRECISION_THREADED_LIBRARIES}" )
endif()
