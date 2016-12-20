/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 09th June 2015
 *
 *      May 2015
 *
 *   Header file for the FFT library under CUDA for cpp code to obtain
 *   information from the GPU cards
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

/* GPU header files */
#if defined (CUDA)
#include "cufft.h"
//#include "cufftw.h"
//#include "cufftXt.h"
#endif

#define FFTW_FORWARD        -1
#define FFTW_BACKWARD       +1

#ifndef _GET_FFT123D_GPU_H_
#define _GET_FFT123D_GPU_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- routines used to interface back to fortran
   *
   *  FORTRAN API -
   * 
   */

#if defined (CUDA) /*preprossing for the CUDA environment */

/* gatherer methods for cpu */
  /* 1D */
  int gather_fft_1D_Z2Z_gpu_cpp(int *nx,
				double complex *in, double complex *out,
				int *sign);
  /* 2D */
  int gather_fft_2D_Z2Z_gpu_cpp(int *nx, int *ny,
				cufftDoubleComplex *in, cufftDoubleComplex *out,
				int *sign);
  int gather_fft_2D_C2C_gpu_cpp(int *nx, int *ny,
				cufftComplex *in, cufftComplex *out,
				int *sign);
  int gather_fft_2D_D2Z_gpu_cpp(int *nx, int *ny,
				const cufftDoubleReal *in,
                                cufftDoubleComplex *out);
  int gather_fft_2D_S2C_gpu_cpp(int *nx, int *ny,
				const cufftReal *in,
                                cufftComplex *out);
  int gather_fft_2D_Z2D_gpu_cpp(int *nx, int *ny,
				const cufftDoubleComplex *in,
				cufftDoubleReal *out);
  int gather_fft_2D_C2S_gpu_cpp(int *nx, int *ny,
				const cufftComplex *in,
                                cufftReal *out);
  /* 3D */
  int gather_fft_3D_Z2Z_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleComplex *in,
				cufftDoubleComplex *out,
				int *sign);
  int gather_fft_3D_C2C_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftComplex *in,
                                cufftComplex *out,
				int *sign);
  int gather_fft_3D_D2Z_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleReal *in,
				cufftDoubleComplex *out);
  int gather_fft_3D_S2C_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftReal *in,
                                cufftComplex *out);
  int gather_fft_3D_Z2D_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleComplex *in,
				cufftDoubleReal *out);
  int gather_fft_3D_C2S_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftComplex *in,
                                cufftReal *out);
/* helper methods from CUDA */
/*error and warning handlers methods */
  int get_warning_message_fft123D_gpu();
  int get_error_id_fft123D_gpu(cudaError_t error_id);
  int get_cuda_cores_error_fft123D_gpu(int ncores);
/* fft handlers methods */

  static const char *_cuFFTGetErrorEnum(cufftResult error) {
    switch (error)
      {
      case CUFFT_SUCCESS        : return "CUFFT_SUCCESS";
      case CUFFT_INVALID_PLAN   : return "CUFFT_INVALID_PLAN";
      case CUFFT_ALLOC_FAILED   : return "CUFFT_ALLOC_FAILED";
      case CUFFT_INVALID_TYPE   : return "CUFFT_INVALID_TYPE";
      case CUFFT_INVALID_VALUE  : return "CUFFT_INVALID_VALUE";
      case CUFFT_INTERNAL_ERROR : return "CUFFT_INTERNAL_ERROR";
      case CUFFT_EXEC_FAILED    : return "CUFFT_EXEC_FAILED";
      case CUFFT_SETUP_FAILED   : return "CUFFT_SETUP_FAILED";
      case CUFFT_INVALID_SIZE   : return "CUFFT_INVALID_SIZE";
      case CUFFT_UNALIGNED_DATA : return "CUFFT_UNALIGNED_DATA";
      case CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
      case CUFFT_INVALID_DEVICE : return "CUFFT_INVALID_DEVICE";
      case CUFFT_PARSE_ERROR    : return "CUFFT_PARSE_ERROR";
      case CUFFT_NO_WORKSPACE   : return "CUFFT_NO_WORKSPACE";
      case CUFFT_NOT_IMPLEMENTED: return "CUFFT_NOT_IMPLEMENTED";
      case CUFFT_LICENSE_ERROR  : return "CUFFT_LICENSE_ERROR";
      }
    return "<unknown>";
  }

#endif /* CUDA */

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _GET_FFT123D_GPU_H_ */
