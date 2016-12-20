/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 1st Apr 2015
 *      Monash University
 *      Apr 2015
 *
 *   code to get the identity in CUDA.
 *
 * @precisions normal z -> s d c
 */

/*preprossing for the CUDA and MAGMA environment */
#if defined (CUDA) && defined (MAGMA) 

//MAGMA include files from {path}/magma-1.6.1/include
#include "magma.h"

//#include "magma_zgetri_blocked_gpu-v2.h"
//#include "magma_dgetri_blocked_gpu-v2.h"
#include "simple_math_gpu.h"

//Magama definition of device pointer for a complex pointer
#define magma_zdevptr(ptr_) ((magmaDoubleComplex*)(uintptr_t)(*(ptr_)))
#define magma_ddevptr(ptr_) ((double*)(uintptr_t)(*(ptr_)))

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *  FORTRAN API - math functions (simple interface)
 *
*/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA function definitions / Data on CPU
*/

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/
  void MAGMA_DGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *dA,
			 magma_int_t *ldda,
			 magma_int_t *ipiv, magma_int_t *info ) {
    magma_dgetrf_gpu(*m, *n, magma_ddevptr(dA), *ldda, ipiv, info );
  }

  void MAGMA_ZGETRF_GPU( magma_int_t *m, magma_int_t *n, devptr_t *dA,
			 magma_int_t *ldda,
			 magma_int_t *ipiv, magma_int_t *info ) {
    magma_zgetrf_gpu(*m, *n, magma_zdevptr(dA), *ldda, ipiv, info );
  }

  void MAGMA_DGETRI_BLOCK_GPU_V2( magma_int_t *n, devptr_t *A,
				  magma_int_t *lda, 
				  magma_int_t *ipiv, double *work, 
				  magma_int_t *lwork, 
				  magma_int_t *info) {
    double *d_a = (double *)(*A);
    magma_dgetri_block_gpu_v2( *n, d_a, *lda, ipiv, work, lwork, info);
  }

  void MAGMA_ZGETRI_BLOCK_GPU_V2( magma_int_t *n, devptr_t *A,
				  magma_int_t *lda, 
				  magma_int_t *ipiv, cuDoubleComplex *work,
				  magma_int_t *lwork, 
				  magma_int_t *info) {
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    magma_zgetri_block_gpu_v2( *n, d_a, *lda, ipiv, work, lwork, info);
  }
  
  extern "C" void magma_dgetrf_gpu_(magma_int_t, magma_int_t, devptr_t, magma_int_t, magma_int_t, magma_int_t) __attribute__((weak,alias("MAGMA_DGETRF_GPU")));
  extern "C" void magma_dgetri_block_gpu_v2_(magma_int_t, devptr_t, magma_int_t, magma_int_t, double, magma_int_t, magma_int_t) __attribute__((weak,alias("MAGMA_DGETRI_BLOCK_GPU_V2")));
  extern "C" void magma_zgetrf_gpu_(magma_int_t, magma_int_t, devptr_t, magma_int_t, magma_int_t, magma_int_t) __attribute__((weak,alias("MAGMA_ZGETRF_GPU")));
  extern "C" void magma_zgetri_block_gpu_v2_(magma_int_t, devptr_t, magma_int_t, magma_int_t, cuDoubleComplex, magma_int_t, magma_int_t) __attribute__((weak,alias("MAGMA_ZGETRI_BLOCK_GPU_V2")));


#ifdef __cplusplus
}
#endif

#endif /* (CUDA) && (MAGMA) */
