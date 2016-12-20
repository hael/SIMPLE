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

//#define debug false
//#define debug_high false

extern "C" int
carte2d_ftExt_corr_gpu_(deviceDetails_t * s_devD,
                        polar_corr_calc_t *s_carte,
                        char TRANSA, char TRANSB,
                        float *r,
                        cuFloatComplex *shmat,
                        const cuFloatComplex *A,
                        const cuFloatComplex *B,
                        int vx, int vy, int vz,
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

  if(vx==0 || vy==0  || ( ( ( alpha == 0.0 ) ) || vz==0 ) )
    {return rc = RC_FAIL;}

  TRANSA = toupper( TRANSA ); 	
  TRANSB = toupper( TRANSB ); 	

  if(ldc < vx ) {return rc = RC_FAIL;}
  if(TRANSA=='C') {
    if(TRANSB=='N')
      {
        if(lda < vx ) {return rc = RC_FAIL;}
        if(ldb < vy ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/

        rc = carte2d_ftExt_corr_C_N(s_devD,
                                    s_carte,
                                    r,
                                    shmat,
                                    A, B,
                                    vx, vy, vz,
                                    alpha,
                                    s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i)=="true"){
          rc = print_s_polar_struct(s_carte);}

      }
  }
  if(TRANSA=='C') {
    if(TRANSB=='F')
      {
        if(lda < vx ) {return rc = RC_FAIL;}
        if(ldb < vy ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/

        rc = carte2d_ftExt_corr_C_F(s_devD,
                                    s_carte,
                                    r,
                                    shmat,
                                    A, B,
                                    vx, vy, vz,
                                    alpha,
                                    s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i)=="true"){
          rc = print_s_polar_struct(s_carte);}

      }
  }
  
  cuCtxSynchronize();

#else

  rc = get_warning_message_corr_Hadmr_gpu();

#endif /* CUDA */

  return rc;

}/* carte2d_ftExt_corr_gpu */
