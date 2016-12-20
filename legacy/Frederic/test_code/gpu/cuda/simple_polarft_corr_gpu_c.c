/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th Sep 2015
 *
 *      September 2015
 *
 *      c code for the MACOSX code this wrapps around the the cpp for LINUX
 *      files methods for the cpp files to extract the information from the GPU
 *      devices.
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

#if defined (CUDA) /*preprossing for the CUDA environment */

/* wrapper method for the polarft corr calculator */
int get_polarft_corr_gpu_c_(deviceDetails_t * s_devD,
                            polar_corr_calc_t *s_polar,
                            char *TRANSA, char *TRANSB,
                            float *r,
                            const float complex *A, 
                            const float complex *B, 
                            int *npart, int *nrot, int *nk,
                            int *lda, int *ldb, int *ldc,
                            float *alpha,
                            bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
  int rc = RC_SUCCESS;

  rc = get_polarft_corr_gpu(s_devD, s_polar, TRANSA, TRANSB,
                            r,
                            A, B, 
                            npart, nrot, nk,
                            lda, ldb, ldc,
                            alpha,s_bench, s_debug_gpu);

  return rc;
}
void GET_POLARFT_CORR_GPU_C_() __attribute__((weak,alias("get_polarft_corr_gpu_c_")));
void get_polarft_corr_gpu_c__() __attribute__((weak,alias("get_polarft_corr_gpu_c_")));
void GET_POLARFT_CORR_GPU_C__() __attribute__((weak,alias("get_polarft_corr_gpu_c_")));

/* wrapper method for the polarft corrAll calculator */
int get_polarft_gencorrall_gpu_c_(deviceDetails_t * s_devD,
                                  polar_corr_calc_t *s_polar,
                                  char *TRANSA, char *TRANSB,
                                  float *r, float *cormat3d,
                                  const float complex *A,
                                  const float complex *B,
                                  const float *sqsums_A,
                                  const float *sqsums_B,
                                  int *npart, int *nrot, int *nk,
                                  int *lda, int *ldb, int *ldc,
                                  float *alpha,
                                  bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
  int rc = RC_SUCCESS;

  rc = get_polarft_gencorrall_gpu(s_devD, s_polar, TRANSA, TRANSB,
                                  r, cormat3d,
                                  A, B,
                                  sqsums_A, sqsums_B,
                                  npart, nrot, nk,
                                  lda, ldb, ldc,
                                  alpha, s_bench, s_debug_gpu);
  
  return rc;
}
void GET_POLARFT_GENCORRALL_GPU_C_() __attribute__((weak,alias("get_polarft_gencorrall_gpu_c_")));
void get_polarft_gencorrall_gpu_c__() __attribute__((weak,alias("get_polarft_gencorrall_gpu_c_")));
void GET_POLARFT_GENCORRALL_GPU_C__() __attribute__((weak,alias("get_polarft_gencorrall_gpu_c_")));

/* wrapper method for the polarft corrAll calculator */
int get_polarft_multi_gpus_gpu_c_(deviceDetails_t * s_devD,
                                  polar_corr_calc_t *s_polar,
                                  char *TRANSA, char *TRANSB,
                                  float *r, float *cormat3d,
                                  const float complex *A,
                                  const float complex *B,
                                  const float *sqsums_A,
                                  const float *sqsums_B,
                                  int *npart, int *nrot, int *nk,
                                  int *lda, int *ldb, int *ldc,
                                  float *alpha,
                                  bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = RC_SUCCESS;
  
  rc = get_polarft_multi_GPUs_gpu(s_devD, s_polar, TRANSA, TRANSB,
                                  r, cormat3d,
                                  A, B,
                                  sqsums_A, sqsums_B,
                                  npart, nrot, nk,
                                  lda, ldb, ldc,
                                  alpha, s_bench, s_debug_gpu);

  return rc;
}
void GET_POLARFT_MULTI_GPUS_GPU_C_() __attribute__((weak,alias("get_polarft_multi_gpus_gpu_c_")));
void get_polarft_multi_gpus_gpu_c__() __attribute__((weak,alias("get_polarft_multi_gpus_gpu_c_")));
void GET_POLARFT_MULTI_GPUS_GPU_C__() __attribute__((weak,alias("get_polarft_multi_gpus_gpu_c_")));

/* wrapper method for the polarft corrAll calculator */
int get_polarft_krnl_opti_gpu_c_(deviceDetails_t * s_devD,
                                  polar_corr_calc_t *s_polar,
                                  char *TRANSA, char *TRANSB,
                                  float *r, float *cormat3d,
                                  const float complex *A,
                                  const float complex *B,
                                  const float *sqsums_A,
                                  const float *sqsums_B,
                                  int *npart, int *nrot, int *nk,
                                  int *lda, int *ldb, int *ldc,
                                  float *alpha,
                                  bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = RC_SUCCESS;
  
  rc = get_polarft_krnl_Opti_gpu(s_devD, s_polar, TRANSA, TRANSB,
                                 r, cormat3d,
                                 A, B,
                                 sqsums_A, sqsums_B,
                                 npart, nrot, nk,
                                 lda, ldb, ldc,
                                 alpha, s_bench, s_debug_gpu);

  return rc;
}
void GET_POLARFT_KRNL_OPTI_GPU_C_() __attribute__((weak,alias("get_polarft_krnl_opti_gpu_c_")));
void get_polarft_krnl_opti_gpu_c__() __attribute__((weak,alias("get_polarft_krnl_opti_gpu_c_")));
void GET_POLARFT_KRNL_OPTI_GPU_C__() __attribute__((weak,alias("get_polarft_krnl_opti_gpu_c_")));

/* wrapper method for the polarft corrAll calculator */
int get_carte2d_ftext_corr_gpu_c_(deviceDetails_t * s_devD,
                                  polar_corr_calc_t *s_carte,
                                  char *TRANSA, char *TRANSB,
                                  float *r,
                                  float complex *shmat,
                                  const float complex *A,
                                  const float complex *B,
                                  int *vx, int *vy, int *vz,
                                  int *lda, int *ldb, int *ldc,
                                  float *alpha,
                                  bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = RC_SUCCESS;

  rc = get_carte2d_ftExt_corr_gpu(s_devD, s_carte, TRANSA, TRANSB,
                                  r,
                                  shmat,
                                  A, B,
                                  vx, vy, vz,
                                  lda, ldb, ldc,
                                  alpha, s_bench, s_debug_gpu);

  return rc;
}
void GET_CARTE2D_FTEXT_CORR_GPU_C_() __attribute__((weak,alias("get_carte2d_ftext_corr_gpu_c_")));
void get_carte2d_ftext_corr_gpu_c__() __attribute__((weak,alias("get_carte2d_ftext_corr_gpu_c_")));
void GET_CARTE2D_FTEXT_CORR_GPU_C__() __attribute__((weak,alias("get_carte2d_ftext_corr_gpu_c_")));



#endif /* CUDA */
