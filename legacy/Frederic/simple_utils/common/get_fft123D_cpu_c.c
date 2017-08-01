/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 19th June 2015
 *
 *      June 2015
 *
 *   code to obtain FFT from the cpu using the FFTW library0
 *
 * @precisions normal z -> s d c
 */

#include "simple.h"
/* CPU header files */
#include <fftw3.h>

/* 3D cufft */
/* Z2Z */
int gather_fft_3d_z2z_cpu_c_(int *nx, int *ny, int *nz,
			     double complex *h_in, double complex *h_out,
			     int *sign) {
  int rc = 0;

  double complex *in = (double complex*) h_in; 
  double complex *out = (double complex*) h_out;

  rc = gather_fft_3D_Z2Z_cpu_cpp(nx, ny, nz, in, out,sign);

  return rc;
}
void GATHER_FFT_3D_Z2Z_CPU_C_() __attribute__((weak,alias("gather_fft_3d_z2z_cpu_c_")));
void gather_fft_3d_z2z_cpu_c__() __attribute__((weak,alias("gather_fft_3d_z2z_cpu_c_")));
void GATHER_FFT_3D_Z2Z_CPU_C__() __attribute__((weak,alias("gather_fft_3d_z2z_cpu_c_")));
/* C2C */
int gather_fft_3d_c2c_cpu_c_(int *nx, int *ny, int *nz,
			     float complex *h_in, float complex *h_out,
			     int *sign) {
  int rc = 0;

  fftwf_complex *in = (float complex*) h_in;
  fftwf_complex *out = (float complex*) h_out;

  rc = gather_fft_3D_C2C_cpu_cpp(nx, ny, nz, in, out,sign);

  return rc;
}
void GATHER_FFT_3D_C2C_CPU_C_() __attribute__((weak,alias("gather_fft_3d_c2c_cpu_c_")));
void gather_fft_3d_c2c_cpu_c__() __attribute__((weak,alias("gather_fft_3d_c2c_cpu_c_")));
void GATHER_FFT_3D_C2C_CPU_C__() __attribute__((weak,alias("gather_fft_3d_c2c_cpu_c_")));
/* D2Z */
int gather_fft_3d_d2z_cpu_c_(int *nx, int *ny, int *nz,
			     double *h_in, double complex *h_out) {
  int rc = 0;

  double *in = (double*) h_in; 
  double complex *out = (double complex*) h_out;

  rc = gather_fft_3D_D2Z_cpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_D2Z_CPU_C_() __attribute__((weak,alias("gather_fft_3d_d2z_cpu_c_")));
void gather_fft_3d_d2z_cpu_c__() __attribute__((weak,alias("gather_fft_3d_d2z_cpu_c_")));
void GATHER_FFT_3D_D2Z_CPU_C__() __attribute__((weak,alias("gather_fft_3d_d2z_cpu_c_")));
/* S2C */
int gather_fft_3d_s2c_cpu_c_(int *nx, int *ny, int *nz,
			     float *h_in, float complex *h_out) {
  int rc = 0;

  float *in = (float*) h_in;
  fftwf_complex *out = (float complex*) h_out;

  rc = gather_fft_3D_S2C_cpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_S2C_CPU_C_() __attribute__((weak,alias("gather_fft_3d_s2c_cpu_c_")));
void gather_fft_3d_s2c_cpu_c__() __attribute__((weak,alias("gather_fft_3d_s2c_cpu_c_")));
void GATHER_FFT_3D_S2C_CPU_C__() __attribute__((weak,alias("gather_fft_3d_s2c_cpu_c_")));
/* Z2D */
int gather_fft_3d_z2d_cpu_c_(int *nx, int *ny, int *nz,
			     double complex *h_in, double *h_out) {
  int rc = 0;

  double complex *in = (double complex*) h_in;
  double *out = (double*) h_out; 

  rc = gather_fft_3D_Z2D_cpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_Z2D_CPU_C_() __attribute__((weak,alias("gather_fft_3d_z2d_cpu_c_")));
void gather_fft_3d_z2d_cpu_c__() __attribute__((weak,alias("gather_fft_3d_z2d_cpu_c_")));
void GATHER_FFT_3D_Z2D_CPU_C__() __attribute__((weak,alias("gather_fft_3d_z2d_cpu_c_")));
/* C2S */
int gather_fft_3d_c2s_cpu_c_(int *nx, int *ny, int *nz,
			     float complex *h_in, float *h_out) {
  int rc = 0;

  fftwf_complex *in = (float complex*) h_in;
  float *out = (float*) h_out;

  rc = gather_fft_3D_C2S_cpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_C2S_CPU_C_() __attribute__((weak,alias("gather_fft_3d_c2s_cpu_c_")));
void gather_fft_3d_c2s_cpu_c__() __attribute__((weak,alias("gather_fft_3d_c2s_cpu_c_")));
void GATHER_FFT_3D_C2S_CPU_C__() __attribute__((weak,alias("gather_fft_3d_c2s_cpu_c_")));
/* 2D cufft */
/* C2C */
int gather_fft_2d_c2c_cpu_c_(int *nx, int *ny,
			     float complex *h_in, float complex *h_out,
			     int *sign) {
  int rc = 0;

  fftwf_complex *in = (float complex*) h_in;
  fftwf_complex *out = (float complex*) h_out;

  rc = gather_fft_2D_C2C_cpu_cpp(nx, ny, in, out, sign);

  return rc;
}
void GATHER_FFT_2D_C2C_CPU_C_() __attribute__((weak,alias("gather_fft_2d_c2c_cpu_c_")));
void gather_fft_2d_c2c_cpu_c__() __attribute__((weak,alias("gather_fft_2d_c2c_cpu_c_")));
void GATHER_FFT_2D_C2C_CPU_C__() __attribute__((weak,alias("gather_fft_2d_c2c_cpu_c_")));
/* Z2Z */
int gather_fft_2d_z2z_cpu_c_(int *nx, int *ny,
			     double complex *h_in, double complex *h_out,
			     int *sign) {
  int rc = 0;

  double complex *in = (double complex*) h_in; 
  double complex *out = (double complex*) h_out;

  rc = gather_fft_2D_Z2Z_cpu_cpp(nx, ny, in, out, sign);

  return rc;
}
void GATHER_FFT_2D_Z2Z_CPU_C_() __attribute__((weak,alias("gather_fft_2d_z2z_cpu_c_")));
void gather_fft_2d_z2z_cpu_c__() __attribute__((weak,alias("gather_fft_2d_z2z_cpu_c_")));
void GATHER_FFT_2D_Z2Z_CPU_C__() __attribute__((weak,alias("gather_fft_2d_z2z_cpu_c_")));
/* D2Z */
int gather_fft_2d_d2z_cpu_c_(int *nx, int *ny,
			     double *h_in, double complex *h_out) {
  int rc = 0;

  double *in = (double*) h_in;
  double complex *out = (double complex*) h_out;

  rc = gather_fft_2D_D2Z_cpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_D2Z_CPU_C_() __attribute__((weak,alias("gather_fft_2d_d2z_cpu_c_")));
void gather_fft_2d_d2z_cpu_c__() __attribute__((weak,alias("gather_fft_2d_d2z_cpu_c_")));
void GATHER_FFT_2D_D2Z_CPU_C__() __attribute__((weak,alias("gather_fft_2d_d2z_cpu_c_")));
/* S2C */
int gather_fft_2d_s2c_cpu_c_(int *nx, int *ny,
			     float *h_in, float complex *h_out) {
  int rc = 0;

  float *in = (float*) h_in;
  fftwf_complex *out = (float complex*) h_out;

  rc = gather_fft_2D_S2C_cpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_S2C_CPU_C_() __attribute__((weak,alias("gather_fft_2d_s2c_cpu_c_")));
void gather_fft_2d_s2c_cpu_c__() __attribute__((weak,alias("gather_fft_2d_s2c_cpu_c_")));
void GATHER_FFT_2D_S2C_CPU_C__() __attribute__((weak,alias("gather_fft_2d_s2c_cpu_c_")));
/* Z2D */
int gather_fft_2d_z2d_cpu_c_(int *nx, int *ny,
			     double complex *h_in, double *h_out) {
  int rc = 0;

  double complex *in = (double complex*) h_in;
  double *out = (double*) h_out;

  rc = gather_fft_2D_Z2D_cpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_Z2D_CPU_C_() __attribute__((weak,alias("gather_fft_2d_z2d_cpu_c_")));
void gather_fft_2d_z2d_cpu_c__() __attribute__((weak,alias("gather_fft_2d_z2d_cpu_c_")));
void GATHER_FFT_2D_Z2D_CPU_C__() __attribute__((weak,alias("gather_fft_2d_z2d_cpu_c_")));
/* C2S */
int gather_fft_2d_c2s_cpu_c_(int *nx, int *ny,
			     float complex *h_in, float *h_out) {
  int rc = 0;

  fftwf_complex *in = (float complex*) h_in;
  float *out = (float*) h_out;

  rc = gather_fft_2D_C2S_cpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_C2S_CPU_C_() __attribute__((weak,alias("gather_fft_2d_c2s_cpu_c_")));
void gather_fft_2d_c2s_cpu_c__() __attribute__((weak,alias("gather_fft_2d_c2s_cpu_c_")));
void GATHER_FFT_2D_C2S_CPU_C__() __attribute__((weak,alias("gather_fft_2d_c2s_cpu_c_")));

/* 1D fftw */
/* Z2Z */
int gather_fft_1d_z2z_cpu_c_(int *nx,
			     double complex *in, double complex *out,
			     int *sign) {
  int rc = 0;

  rc = gather_fft_1D_Z2Z_cpu_cpp(nx,in,out,sign);

  return rc;
}
void GATHER_FFT_1D_Z2Z_CPU_C_() __attribute__((weak,alias("gather_fft_1d_z2z_cpu_c_")));
void gather_fft_1d_z2z_cpu_c__() __attribute__((weak,alias("gather_fft_1d_z2z_cpu_c_")));
void GATHER_FFT_1D_Z2Z_CPU_C__() __attribute__((weak,alias("gather_fft_1d_z2z_cpu_c_")));
