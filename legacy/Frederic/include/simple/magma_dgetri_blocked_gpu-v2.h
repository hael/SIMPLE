/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 1st Apr 2015
 *      Monash University
 *      Apr 2015
 *
 * @precisions normal z -> s d c
 */

#include "simple_magma_common.h"

#ifndef _MAGMA_DGETRI_BLOCKED_GPU_V2_
#define _MAGMA_DGETRI_BLOCKED_GPU_V2_

//#define PRECISION_d

#define lapackf77_dtrti2  dtrti2_

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Lapack routines used for the block code
*/
  void lapackf77_dtrti2( char *uplo, char *diag, 
			 magma_int_t *n, double *a,
			 magma_int_t *lda, magma_int_t *info);

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/

#if defined (CUDA) && defined (MAGMA) /*preprossing for the CUDA and MAGMA environment */

  magma_int_t magma_dgetri_block_gpu_v2( magma_int_t n, double *dA,
					 magma_int_t ldda,
					 magma_int_t *ipiv, double *work,
					 magma_int_t *lwork,
					 magma_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Some getters
*/

#endif /* CUDA and MAGMA */

#ifdef __cplusplus
}
#endif

//#undef PRECISION_d
#endif /* _MAGMA_DGETRI_BLOCKED_GPU-V2_ */

