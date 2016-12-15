/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 28th April 2015
 *      Monash University
 *      April 2015
 *
 *   Header files for the kernels on GPU.
 *
 * @precisions normal z -> s d c
 */
#include "simple.h"

#ifndef _COMMONBLAS_ZZ2D_H_
#define _COMMONBLAS_ZZ2D_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined (CUDA) /*preprossing for the CUDA environment */

#define ZZ2DGEMM( name ) int zz2dgemm_kernel_##name(double *C, const cuDoubleComplex *A, const cuDoubleComplex *B, int m, int n, int k, int lda, int ldb, int ldc, double alpha, double beta)
  ZZ2DGEMM( N_N                       );
  ZZ2DGEMM( T_N                       );
  ZZ2DGEMM( N_T                       );
  ZZ2DGEMM( T_T                       );

#define ZZ2DGEMM_SUMSQ( name ) int zz2dgemm_kernel_sumsq_##name(double *C, const cuDoubleComplex *A, int m, int n, int k, int lda, int ldb, int ldc, double alpha, double beta)
  ZZ2DGEMM_SUMSQ( N_N                 );
  ZZ2DGEMM_SUMSQ( T_T                 );

#define POLARFT( name ) int polarft_corr_##name(deviceDetails_t * s_devD, polar_corr_calc_t *s_polar, float *r, const cuFloatComplex *A, const cuFloatComplex *B, int npart, int nrot, int nk, float alpha, bench_t *s_bench, debug_gpu_t *s_debug_gpu)
  POLARFT( N_N                        );
  POLARFT( F_N                        );
  POLARFT( P_N                        );
  POLARFT( X_N                        );

#define GENCORRALL( name )int polarft_gencorrAll_##name(deviceDetails_t * s_devD, polar_corr_calc_t *s_polar, float *r, float *cormat3d, const cuFloatComplex *A, const cuFloatComplex *B, const float *sqsums_A, const float *sqsums_B, int npart, int nrot, int nk, float alpha, bench_t *s_bench, debug_gpu_t *s_debug_gpu)
  GENCORRALL( Z_N                        );

#define MULTIGPUS( name )int polarft_multi_GPUs_##name(deviceDetails_t * s_devD, polar_corr_calc_t *s_polar, float *r, float *cormat3d, const cuFloatComplex *A, const cuFloatComplex *B, const float *sqsums_A, const float *sqsums_B, int npart, int nrot, int nk, float alpha, bench_t *s_bench, debug_gpu_t *s_debug_gpu)
  MULTIGPUS( Z_M                        );

#define KRNLOPTI( name )int polarft_krnl_Opti_##name(deviceDetails_t * s_devD, polar_corr_calc_t *s_polar, float *r, float *cormat3d, const cuFloatComplex *A, const cuFloatComplex *B, const float *sqsums_A, const float *sqsums_B, int npart, int nrot, int nk, float alpha, bench_t *s_bench, debug_gpu_t *s_debug_gpu)
  KRNLOPTI( O_N                        );

#define CARTE2DFTEXT( name )int carte2d_ftExt_corr_##name(deviceDetails_t * s_devD, polar_corr_calc_t *s_carte, float *r, cuFloatComplex *shmat, const cuFloatComplex *A, const cuFloatComplex *B, int vx, int vy, int vz, float alpha, bench_t *s_bench, debug_gpu_t *s_debug_gpu)
  CARTE2DFTEXT( C_N                        );
  CARTE2DFTEXT( C_F                        );

#endif /* CUDA */

#ifdef __cplusplus
}
#endif

#endif /* _COMMONBLAS_ZZ2D_H_ */

