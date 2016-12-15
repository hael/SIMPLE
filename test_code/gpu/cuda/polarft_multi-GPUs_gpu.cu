/**********************************************************************
 * SIMPLE routine for CUDA BLAS, including a CUDA kernel and wappers  *
 **********************************************************************
 *
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 8th Feb 2016
 *      Monash University
 *      March 2016
 *
 *      Wrapper routine which calculates the product of two complex matrices 
 *      elementwise and returns only the real part of the product
 *      ie
 *      if A = x_1 + i y_1 and B = x_2 + i y_2
 *      then it returns a double product matrix element wise C  
 *      Re(A B^{*}) = x_1 * x_2 + y_1 * y_2 = (double*) C
 *      it also calculates the sum and the normailsation 
 *
 * @precisions normal z -> s d c
 */
#include "common_magma.h"
#include "commonblas_zz2d.h"
#include "polarft_gpu.h"
#include "simple.h"
#include "get_deviceQuery_gpu.h"

//#define debug false
//#define debug_high false

extern "C" int
polarft_multi_GPUs_gpu_(deviceDetails_t * s_devD,
                        polar_corr_calc_t *s_polar,
                        char TRANSA, char TRANSB,
                        float *r, float *cormat3d,
                        const cuFloatComplex *A,
                        const cuFloatComplex *B,
                        const float *sqsums_A,
                        const float *sqsums_B,
                        int npart, int nrot, int nk,
                        int lda, int ldb, int ldc,
                        float alpha,
                        bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = RC_SUCCESS; // return code
#if defined (CUDA) /*preprossing for the CUDA environment */
  
  if ( get_bool(s_debug_gpu->debug_i) == "true" ) {
    rc = print_s_debug_struct(s_debug_gpu);
    rc = print_s_bench_struct(s_bench);
    rc = print_s_devD_struct(s_devD);
  }

  if(npart==0 || nrot==0  || ( ( ( alpha == 0.0 ) ) || nk==0 ) )
    {return rc = RC_FAIL;}

  TRANSA = toupper( TRANSA ); 	
  TRANSB = toupper( TRANSB ); 	

  if(ldc < npart ) {return rc = RC_FAIL;}
  if(TRANSA=='Z') {
    if(TRANSB=='M')
      {
        if(lda < npart ) {return rc = RC_FAIL;}
        if(ldb < nrot ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/
        
        rc = polarft_multi_GPUs_Z_M(s_devD,
                                    s_polar,
                                    r, cormat3d,
                                    A, B,
                                    sqsums_A,
                                    sqsums_B,
                                    npart, nrot, nk,
                                    alpha,
                                    s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i)=="true"){
          rc = print_s_polar_struct(s_polar);}

      }
  }
  
  cuCtxSynchronize();

#else

  rc = get_warning_message_corr_Hadmr_gpu();

#endif /* CUDA */

  return rc;

}/* polarft_multi_GPUs_tesla_gpu */
