/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 21st May 2015
 *
 *      May 2015
 *
 *   code to obtain information from the GPU cards
 *
 * @precisions normal z -> s d c
 */

#include "simple.h"

#if defined (CUDA)
#include "cufft.h"
//#include "cufftw.h"
//#include "cufftXt.h"
#endif

#if defined (CUDA)
/* 3D cufft */
/* Z2Z */
int gather_fft_3d_z2z_gpu_c_(int *nx, int *ny, int *nz,
			     const double complex *h_in, double complex *h_out,
			     int *sign) {
  int rc = 0;

  cufftDoubleComplex *in = (cufftDoubleComplex*) h_in; 
  cufftDoubleComplex *out = (cufftDoubleComplex*) h_out;

  rc = gather_fft_3D_Z2Z_gpu_cpp(nx, ny, nz, in, out,sign);

  return rc;
}
void GATHER_FFT_3D_Z2Z_GPU_C_() __attribute__((weak,alias("gather_fft_3d_z2z_gpu_c_")));
void gather_fft_3d_z2z_gpu_c__() __attribute__((weak,alias("gather_fft_3d_z2z_gpu_c_")));
void GATHER_FFT_3D_Z2Z_GPU_C__() __attribute__((weak,alias("gather_fft_3d_z2z_gpu_c_")));
/* C2C */
int gather_fft_3d_c2c_gpu_c_(int *nx, int *ny, int *nz,
			     const float complex *h_in, float complex *h_out,
			     int *sign) {
  int rc = 0;

  cufftComplex *in = (cufftComplex*) h_in;
  cufftComplex *out = (cufftComplex*) h_out;

  rc = gather_fft_3D_C2C_gpu_cpp(nx, ny, nz, in, out,sign);

  return rc;
}
void GATHER_FFT_3D_C2C_GPU_C_() __attribute__((weak,alias("gather_fft_3d_c2c_gpu_c_")));
void gather_fft_3d_c2c_gpu_c__() __attribute__((weak,alias("gather_fft_3d_c2c_gpu_c_")));
void GATHER_FFT_3D_C2C_GPU_C__() __attribute__((weak,alias("gather_fft_3d_c2c_gpu_c_")));
/* D2Z */
int gather_fft_3d_d2z_gpu_c_(int *nx, int *ny, int *nz,
			     const double *h_in, double complex *h_out) {
  int rc = 0;

  cufftDoubleReal *in = (cufftDoubleReal*) h_in; 
  cufftDoubleComplex *out = (cufftDoubleComplex*) h_out;

  rc = gather_fft_3D_D2Z_gpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_D2Z_GPU_C_() __attribute__((weak,alias("gather_fft_3d_d2z_gpu_c_")));
void gather_fft_3d_d2z_gpu_c__() __attribute__((weak,alias("gather_fft_3d_d2z_gpu_c_")));
void GATHER_FFT_3D_D2Z_GPU_C__() __attribute__((weak,alias("gather_fft_3d_d2z_gpu_c_")));
/* S2C */
int gather_fft_3d_s2c_gpu_c_(int *nx, int *ny, int *nz,
			     const float *h_in, float complex *h_out) {
  int rc = 0;

  cufftReal *in = (cufftReal*) h_in;
  cufftComplex *out = (cufftComplex*) h_out;

  rc = gather_fft_3D_S2C_gpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_S2C_GPU_C_() __attribute__((weak,alias("gather_fft_3d_s2c_gpu_c_")));
void gather_fft_3d_s2c_gpu_c__() __attribute__((weak,alias("gather_fft_3d_s2c_gpu_c_")));
void GATHER_FFT_3D_S2C_GPU_C__() __attribute__((weak,alias("gather_fft_3d_s2c_gpu_c_")));
/* Z2D */
int gather_fft_3d_z2d_gpu_c_(int *nx, int *ny, int *nz,
			     const double complex *h_in, double *h_out) {
  int rc = 0;

  cufftDoubleComplex *in = (cufftDoubleComplex*) h_in;
  cufftDoubleReal *out = (cufftDoubleReal*) h_out; 

  rc = gather_fft_3D_Z2D_gpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_Z2D_GPU_C_() __attribute__((weak,alias("gather_fft_3d_z2d_gpu_c_")));
void gather_fft_3d_z2d_gpu_c__() __attribute__((weak,alias("gather_fft_3d_z2d_gpu_c_")));
void GATHER_FFT_3D_Z2D_GPU_C__() __attribute__((weak,alias("gather_fft_3d_z2d_gpu_c_")));
/* C2S */
int gather_fft_3d_c2s_gpu_c_(int *nx, int *ny, int *nz,
			     const float complex *h_in, float *h_out) {
  int rc = 0;

  cufftComplex *in = (cufftComplex*) h_in;
  cufftReal *out = (cufftReal*) h_out;

  rc = gather_fft_3D_C2S_gpu_cpp(nx, ny, nz, in, out);

  return rc;
}
void GATHER_FFT_3D_C2S_GPU_C_() __attribute__((weak,alias("gather_fft_3d_c2s_gpu_c_")));
void gather_fft_3d_c2s_gpu_c__() __attribute__((weak,alias("gather_fft_3d_c2s_gpu_c_")));
void GATHER_FFT_3D_C2S_GPU_C__() __attribute__((weak,alias("gather_fft_3d_c2s_gpu_c_")));

/* 2D cufft */
/* Z2Z */
int gather_fft_2d_z2z_gpu_c_(int *nx, int *ny,
			     double complex *h_in, double complex *h_out,
			     int *sign) {
  int rc = 0;

  cufftDoubleComplex *in = (cufftDoubleComplex*) h_in; 
  cufftDoubleComplex *out = (cufftDoubleComplex*) h_out;

  rc = gather_fft_2D_Z2Z_gpu_cpp(nx, ny, in, out, sign);

  return rc;
}
void GATHER_FFT_2D_Z2Z_GPU_C_() __attribute__((weak,alias("gather_fft_2d_z2z_gpu_c_")));
void gather_fft_2d_z2z_gpu_c__() __attribute__((weak,alias("gather_fft_2d_z2z_gpu_c_")));
void GATHER_FFT_2D_Z2Z_GPU_C__() __attribute__((weak,alias("gather_fft_2d_z2z_gpu_c_")));
/* C2C */
int gather_fft_2d_c2c_gpu_c_(int *nx, int *ny,
			     float complex *h_in, float complex *h_out,
			     int *sign) {
  int rc = 0;

  cufftComplex *in = (cufftComplex*) h_in;
  cufftComplex *out = (cufftComplex*) h_out;

  rc = gather_fft_2D_C2C_gpu_cpp(nx, ny, in, out, sign);

  return rc;
}
void GATHER_FFT_2D_C2C_GPU_C_() __attribute__((weak,alias("gather_fft_2d_c2c_gpu_c_")));
void gather_fft_2d_c2c_gpu_c__() __attribute__((weak,alias("gather_fft_2d_c2c_gpu_c_")));
void GATHER_FFT_2D_C2C_GPU_C__() __attribute__((weak,alias("gather_fft_2d_c2c_gpu_c_")));
/* D2Z */
int gather_fft_2d_d2z_gpu_c_(int *nx, int *ny,
			     const double *h_in, double complex *h_out) {
  int rc = 0;

  const cufftDoubleReal *in = (cufftDoubleReal*) h_in;
  cufftDoubleComplex *out = (cufftDoubleComplex*) h_out;

  rc = gather_fft_2D_D2Z_gpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_D2Z_GPU_C_() __attribute__((weak,alias("gather_fft_2d_d2z_gpu_c_")));
void gather_fft_2d_d2z_gpu_c__() __attribute__((weak,alias("gather_fft_2d_d2z_gpu_c_")));
void GATHER_FFT_2D_D2Z_GPU_C__() __attribute__((weak,alias("gather_fft_2d_d2z_gpu_c_")));
/* S2C */
int gather_fft_2d_s2c_gpu_c_(int *nx, int *ny,
			     const float *h_in, float complex *h_out) {
  int rc = 0;

  const cufftReal *in = (cufftReal*) h_in;
  cufftComplex *out = (cufftComplex*) h_out;

  rc = gather_fft_2D_S2C_gpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_S2C_GPU_C_() __attribute__((weak,alias("gather_fft_2d_s2c_gpu_c_")));
void gather_fft_2d_s2c_gpu_c__() __attribute__((weak,alias("gather_fft_2d_s2c_gpu_c_")));
void GATHER_FFT_2D_S2C_GPU_C__() __attribute__((weak,alias("gather_fft_2d_s2c_gpu_c_")));
/* Z2D */
int gather_fft_2d_z2d_gpu_c_(int *nx, int *ny,
			     const double complex *h_in, double *h_out) {
  int rc = 0;

  const cufftDoubleComplex *in = (cufftDoubleComplex*) h_in;
  cufftDoubleReal *out = (cufftDoubleReal*) h_out;

  rc = gather_fft_2D_Z2D_gpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_Z2D_GPU_C_() __attribute__((weak,alias("gather_fft_2d_z2d_gpu_c_")));
void gather_fft_2d_z2d_gpu_c__() __attribute__((weak,alias("gather_fft_2d_z2d_gpu_c_")));
void GATHER_FFT_2D_Z2D_GPU_C__() __attribute__((weak,alias("gather_fft_2d_z2d_gpu_c_")));
/* C2S */
int gather_fft_2d_c2s_gpu_c_(int *nx, int *ny,
			     const float complex *h_in, float *h_out) {
  int rc = 0;

  const cufftComplex *in = (cufftComplex*) h_in;
  cufftReal *out = (cufftReal*) h_out;

  rc = gather_fft_2D_C2S_gpu_cpp(nx, ny, in, out);

  return rc;
}
void GATHER_FFT_2D_C2S_GPU_C_() __attribute__((weak,alias("gather_fft_2d_c2s_gpu_c_")));
void gather_fft_2d_c2s_gpu_c__() __attribute__((weak,alias("gather_fft_2d_c2s_gpu_c_")));
void GATHER_FFT_2D_C2S_GPU_C__() __attribute__((weak,alias("gather_fft_2d_c2s_gpu_c_")));

/* 1D cufft */
/* Z2Z */
int gather_fft_1d_z2z_gpu_c_(int *nx, double complex *in, double complex *out,
                             int *sign) {
  int rc = 0;

  rc = gather_fft_1D_Z2Z_gpu_cpp(nx,in,out,sign);

  return rc;
}
void GATHER_FFT_1D_Z2Z_GPU_C_() __attribute__((weak,alias("gather_fft_1d_z2z_gpu_c_")));
void gather_fft_1d_z2z_gpu_c__() __attribute__((weak,alias("gather_fft_1d_z2z_gpu_c_")));
void GATHER_FFT_1D_Z2Z_GPU_C__() __attribute__((weak,alias("gather_fft_1d_z2z_gpu_c_")));

#endif
