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
/* CPU header files */
#include <fftw3.h>

#ifndef _GET_FFT123D_CPU_H_
#define _GET_FFT123D_CPU_H_
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

/* gatherer methods for cpu */
  /* 1D */
  int gather_fft_1D_Z2Z_cpu_cpp(int *nx,
				double complex *in, double complex *out,
				int *sign);
  /* 2D */
  int gather_fft_2D_Z2Z_cpu_cpp(int *nx, int *ny,
				double complex *in, double complex *out,
				int *sign);
  int gather_fft_2D_C2C_cpu_cpp(int *nx, int *ny,
                                fftwf_complex *in, fftwf_complex *out,
                                int *sign);
  int gather_fft_2D_D2Z_cpu_cpp(int *nx, int *ny,
				double *in,
                                double complex *out);
  int gather_fft_2D_S2C_cpu_cpp(int *nx, int *ny,
				float *in,
                                fftwf_complex *out);
  int gather_fft_2D_Z2D_cpu_cpp(int *nx, int *ny,
				double complex *in,
                                double *out);
  int gather_fft_2D_C2S_cpu_cpp(int *nx, int *ny,
				fftwf_complex *in,
                                float *out);
  /* 3D */
  int gather_fft_3D_Z2Z_cpu_cpp(int *nx, int *ny, int *nz,
				double complex *in,
                                double complex *out,
				int *sign);
  int gather_fft_3D_C2C_cpu_cpp(int *nx, int *ny, int *nz,
				fftwf_complex *in,
                                fftwf_complex *out,
				int *sign);
  int gather_fft_3D_D2Z_cpu_cpp(int *nx, int *ny, int *nz,
				double *in,
                                double complex *out);
  int gather_fft_3D_S2C_cpu_cpp(int *nx, int *ny, int *nz,
				float *in,
                                fftwf_complex *out);
  int gather_fft_3D_Z2D_cpu_cpp(int *nx, int *ny, int *nz,
				double complex *in,
                                double *out);
  int gather_fft_3D_C2S_cpu_cpp(int *nx, int *ny, int *nz,
				fftwf_complex *in,
                                float *out);
/* helper methods */
/*error and warning handlers methods */
  int get_warning_message_fft123D_cpu();
/* fft handlers methods */
/*testers methods */
  int test_fft_1D_cpu_cpp(double *xi, double *xf, int *nstep,
			  double *x, double *fx, double complex *fh);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _GET_FFT123D_CPU_H_ */
