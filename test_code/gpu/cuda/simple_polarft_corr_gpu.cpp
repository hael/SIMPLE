/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th Sep 2015
 *
 *      May 2015
 *
 *   code to obtain the correlator between pFT on GPU cards
 *
 * @precisions normal z -> s d c
 */

//MAGMA include files from {path}/magma-1.6.1/include
#if defined (CUDA) && defined (MAGMA) 
#include "magma.h"

#include "polarft_gpu.h"

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************/
/**
 *  FORTRAN API - math functions (simple interface)
 **/

#if defined (CUDA) /*preprossing for the CUDA environment */

  /* wrapper method for the polarft corr calculator */
  int get_polarft_corr_gpu(deviceDetails_t * s_devD,
                           polar_corr_calc_t *s_polar,
                           char *TRANSA, char *TRANSB,
                           float *r,
                           const float complex *A_in, 
                           const float complex *B_in, 
                           int *npart, int *nrot, int *nk,
                           int *lda, int *ldb, int *ldc,
                           float *alpha,
                           bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
    int rc = RC_SUCCESS;

    cuFloatComplex *A = (cuFloatComplex*) A_in;
    cuFloatComplex *B = (cuFloatComplex*) B_in;

    rc = polarft_corr_Hadmr_gpu_(s_devD, s_polar, *TRANSA, *TRANSB,
                                       r,
                                       A, B,
                                       *npart, *nrot, *nk,
                                       *lda, *ldb, *ldc,
                                       *alpha,s_bench, s_debug_gpu);
    return rc;
  }

  /* wrapper method for the polarft corrAll calculator */
  int get_polarft_gencorrall_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_polar,
                                 char *TRANSA, char *TRANSB,
                                 float *r, float *cormat3d_in,
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 const float *sqsums_A_in,
                                 const float *sqsums_B_in,
                                 int *npart, int *nrot, int *nk,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
    int rc = RC_SUCCESS;

    cuFloatComplex *A = (cuFloatComplex*) A_in;
    cuFloatComplex *B = (cuFloatComplex*) B_in;
    float *sqsums_A = (float*) sqsums_A_in;
    float *sqsums_B = (float*) sqsums_B_in;
    float *cormat3d = (float*) cormat3d_in;

    rc = polarft_gencorrAll_gpu_(s_devD, s_polar, *TRANSA, *TRANSB,
                                       r, cormat3d,
                                       A, B,
                                       sqsums_A,sqsums_B,
                                       *npart, *nrot, *nk,
                                       *lda, *ldb, *ldc,
                                       *alpha,
                                       s_bench, s_debug_gpu);

    return rc;
  }

  /* wrapper method for the polarft corrAll calculator */
  int get_polarft_multi_GPUs_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_polar,
                                 char *TRANSA, char *TRANSB,
                                 float *r, float *cormat3d_in,
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 const float *sqsums_A_in,
                                 const float *sqsums_B_in,
                                 int *npart, int *nrot, int *nk,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
    int rc = RC_SUCCESS;

    cuFloatComplex *A = (cuFloatComplex*) A_in;
    cuFloatComplex *B = (cuFloatComplex*) B_in;
    float *sqsums_A = (float*) sqsums_A_in;
    float *sqsums_B = (float*) sqsums_B_in;
    float *cormat3d = (float*) cormat3d_in;

    rc = polarft_multi_GPUs_gpu_(s_devD, s_polar, *TRANSA, *TRANSB,
                                       r, cormat3d,
                                       A, B,
                                       sqsums_A,sqsums_B,
                                       *npart, *nrot, *nk,
                                       *lda, *ldb, *ldc,
                                       *alpha,
                                       s_bench, s_debug_gpu);

    return rc;
  }
  /* wrapper method for the polarft corrAll calculator */
  int get_polarft_krnl_Opti_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_polar,
                                 char *TRANSA, char *TRANSB,
                                 float *r, float *cormat3d_in,
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 const float *sqsums_A_in,
                                 const float *sqsums_B_in,
                                 int *npart, int *nrot, int *nk,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
    int rc = RC_SUCCESS;

    cuFloatComplex *A = (cuFloatComplex*) A_in;
    cuFloatComplex *B = (cuFloatComplex*) B_in;
    float *sqsums_A = (float*) sqsums_A_in;
    float *sqsums_B = (float*) sqsums_B_in;
    float *cormat3d = (float*) cormat3d_in;

    rc = polarft_krnl_Opti_gpu_(s_devD, s_polar, *TRANSA, *TRANSB,
                                r, cormat3d,
                                A, B,
                                sqsums_A,sqsums_B,
                                *npart, *nrot, *nk,
                                *lda, *ldb, *ldc,
                                *alpha,
                                s_bench, s_debug_gpu);

    return rc;
  }
  /* wrapper method for the polarft corrAll calculator */
  int get_carte2d_ftExt_corr_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_carte,
                                 char *TRANSA, char *TRANSB,
                                 float *r,
                                 float complex *shmat_in, 
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 int *vx, int *vy, int *vz,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu) {
    int rc = RC_SUCCESS;

    cuFloatComplex     *A = (cuFloatComplex*) A_in;
    cuFloatComplex     *B = (cuFloatComplex*) B_in;
    cuFloatComplex *shmat = (cuFloatComplex*) shmat_in;

    rc = carte2d_ftExt_corr_gpu_(s_devD, s_carte, *TRANSA, *TRANSB,
                                 r,
                                 shmat,
                                 A, B,
                                 *vx, *vy, *vz,
                                 *lda, *ldb, *ldc,
                                 *alpha,
                                 s_bench, s_debug_gpu);

    return rc;
  }

  /* the aliases for external access */

  extern "C" int get_polarft_corr_gpu_()       __attribute__((weak,alias("get_polarft_corr_gpu")));
  extern "C" int get_polarft_gencorrall_gpu_() __attribute__((weak,alias("get_polarft_gencorrall_gpu")));
  extern "C" int get_polarft_multi_GPUs_gpu_() __attribute__((weak,alias("get_polarft_multi_GPUs_gpu")));
  extern "C" int get_polarft_krnl_Opti_gpu_()  __attribute__((weak,alias("get_polarft_krnl_Opti_gpu")));
  extern "C" int get_carte2d_ftExt_corr_gpu_() __attribute__((weak,alias("get_carte2d_ftExt_corr_gpu")));
  
#endif /* CUDA */

#ifdef __cplusplus
}
#endif

#endif /* defined (CUDA) && defined (MAGMA) */

