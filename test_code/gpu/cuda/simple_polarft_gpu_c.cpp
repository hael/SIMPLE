/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 22nd Apr 2015
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

#include "polarft_gpu.h"

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
/*
  void CORR_CUDA_GPU(int *rot1, int *rot2,
                          self1%nradial/2,           &
                          self1%khp,                 &
                          self1%klp,                 &
                          self1%pft(i1,1:self1%khp), &
                          self2%pft(i2,j),           &
                          r,                         &
                          sumasq,                    &
                          sumbsq                     )
*/
  void ZZ2DGEMM_ELMTWS_TESLA_GPU( char *TRANSA, char *TRANSB,
				  int *m, int *n, int *k,
				  double *alpha,
				  const devptr_t *A, int *lda,
				  const devptr_t *B, int *ldb,
				  double *beta,
				  devptr_t *C, int *ldc)
  {
    const cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    const cuDoubleComplex *d_b = (cuDoubleComplex *)(*B);
    double *d_c          = (double *)(*C);
    zz2dgemm_ElmtWs_tesla_gpu_(*TRANSA, *TRANSB, *m, *n, *k, 
			       *alpha, 
			       d_a, *lda, 
			       d_b, *ldb, *beta, 
			       d_c, *ldc);
  }
  void ZZ2DGEMM_ELMTWS_TESLA_SUMSQ_GPU( char *TRANSA, char *TRANSB,
					int *m, int *n, int *k,
					double *alpha,
					const devptr_t *A, int *lda,
					const devptr_t *B, int *ldb,
					double *beta,
					devptr_t *C, int *ldc)
  {
    const cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    const cuDoubleComplex *d_b = (cuDoubleComplex *)(*B);
    double *d_c          = (double *)(*C);
    zz2dgemm_ElmtWs_tesla_sumsq_gpu_(*TRANSA, *TRANSB, *m, *n, *k, 
				     *alpha, 
				     d_a, *lda, 
				     d_b, *ldb, *beta, 
				     d_c, *ldc);
  }


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA function definitions / Data on GPU
*/

  extern "C" void zz2dgemm_elmtws_tesla_gpu_( char, char, int, int, int, double , const devptr_t, int, const devptr_t, int, double, devptr_t, int) __attribute__((weak,alias("ZZ2DGEMM_ELMTWS_TESLA_GPU")));
  extern "C" void zz2dgemm_elmtws_tesla_sumsq_gpu_( char, char, int, int, int, double , const devptr_t, int, const devptr_t, int, double, devptr_t, int) __attribute__((weak,alias("ZZ2DGEMM_ELMTWS_TESLA_SUMSQ_GPU")));

#ifdef __cplusplus
}
#endif

#endif /* (CUDA) && (MAGMA) */
