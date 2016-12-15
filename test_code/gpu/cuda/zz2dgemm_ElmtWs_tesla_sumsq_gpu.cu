/**********************************************************************
 * SIMPLE routine for CUDA BLAS, including a CUDA kernel and wappers  *
 **********************************************************************
 *
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 27th Apr 2015
 *      Monash University
 *      April 2015
 *
 *      Routine which calculates the product of two complex matrices 
 *      elementwise and returns only the real part of the product
 *      ie
 *      if A = x_1 + i y_1 and B = x_2 + i y_2
 *      then it returns a double product matrix element wise C  
 *      Re(A B^{*}) = x_1 * x_2 + y_1 * y_2 = (double*) C
 *      it also calculates the 
 *
 * @precisions normal z -> s d c
 */
#include "common_magma.h"
#include "commonblas_zz2d.h"
#include "polarft_gpu.h"
#include "simple.h"

//#define zz2dgemm_ElmtWs_tesla_gpu dzgemm_gpu

extern "C" void
zz2dgemm_ElmtWs_tesla_sumsq_gpu_( char TRANSA, char TRANSB,
				  int m, int n, int k,
				  double alpha, 
				  const cuDoubleComplex *A, int lda, 
				  const cuDoubleComplex *B, int ldb,
				  double beta,
				  double *C, int ldc)
{
  int rc = 0; // return code
  double *dC;

#if defined (CUDA) /*preprossing for the CUDA environment */

  dC = C;

  if(m==0 || n==0  || ( ( ( alpha == 0.0 ) ) || k==0 ) && ( beta == 0.0 ) ) return;
    
  TRANSA = toupper( TRANSA ); 	
  TRANSB = toupper( TRANSB ); 	

  if(ldc < m ) {rc = -1; return ;}
  if(TRANSA=='N')
    {
      if(TRANSB=='N')
	{ 
	  if(lda < m ) {rc = -1; return ;}
	  if(ldb < k ) {rc = -1; return ;}
	  /*====================================================================
	    =================C = alpha * A * B + beta * C ======================
	    ==================================================================*/
	  rc = zz2dgemm_kernel_sumsq_N_N(dC,A, m, n, k, lda, ldb, ldc, alpha,
                                         beta);
	}
      else 
	{
	  if(lda < m ) {rc = -1; return ;}
	  if(ldb < n ) {rc = -1; return ;}
	  /*====================================================================
	    ===================C = alpha * A * B^T + beta * C===================
	    ==================================================================*/
	  //TODO: if needed, not yet implemented 
	}
    } 
  else
    {
      if(TRANSB=='N') {
	if(lda < k ) {rc = -1; return ;}
	if(ldb < k ) {rc = -1; return ;}
	/*======================================================================
	  ===================C = alpha * A^T * B + beta * C=====================
	  ====================================================================*/
	//TODO: if needed, not yet implemented 
      }
      else
	{
	  if(lda < k) {rc = -1; return ;}
	  if(ldb < n ) {rc = -1; return ;}
	  /*====================================================================
	    ===================C = alpha * A^T* B^T + beta * C==================
	    ==================================================================*/
	  rc = zz2dgemm_kernel_sumsq_T_T(dC,A, m, n, k, lda, ldb, ldc, alpha,
                                         beta);
	}
    }

  cuCtxSynchronize();
#else

  rc = get_warning_message_corr_Hadmr_gpu();

#endif /* CUDA */

}/* zz2dgemm_ElmtWs_tesla_gpu */
