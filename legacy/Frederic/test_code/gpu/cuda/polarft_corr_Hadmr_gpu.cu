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

//#define debug false
//#define debug_high false

extern "C" int
polarft_corr_Hadmr_gpu_(deviceDetails_t * s_devD,
                        polar_corr_calc_t *s_polar,
                        char TRANSA, char TRANSB,
                        float *r,
                        const cuFloatComplex *A,
                        const cuFloatComplex *B,
                        int npart, int nrot, int nk,
                        int lda, int ldb, int ldc,
                        float alpha,
                        bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = RC_SUCCESS; // return cide

  if (get_bool(s_debug_gpu->debug_i) == "true" ) {
    rc = print_s_debug_struct(s_debug_gpu);
    rc = print_s_bench_struct(s_bench);
    rc = print_s_devD_struct(s_devD);
  }

#if defined (CUDA) /*preprossing for the CUDA environment */

  if(npart==0 || nrot==0  || ( ( ( alpha == 0.0 ) ) || nk==0 ) )
    {return rc = RC_FAIL;}

  TRANSA = toupper( TRANSA ); 	
  TRANSB = toupper( TRANSB ); 	

  if(ldc < npart ) {return rc = RC_FAIL;}
  if(TRANSA=='N') {
    if(TRANSB=='N') { 
      if(lda < npart ) {return rc = RC_FAIL;}
      if(ldb < nrot ) {return rc = RC_FAIL;}
      /*======================================================================
        ===================C = alpha * A * B            ======================
        ====================================================================*/
      rc = polarft_corr_N_N(s_devD, s_polar,
                            r,
                            A, B,
                            npart, nrot, nk,
                            alpha,
                            s_bench, s_debug_gpu);

      if (get_bool(s_debug_gpu->debug_i) == "true"){
        rc = print_s_polar_struct(s_polar);}

    }
  }
  else if(TRANSA=='F') {
    if(TRANSB=='N')
      {
        if(lda < npart ) {return rc = RC_FAIL;}
        if(ldb < nrot ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/

        rc = polarft_corr_F_N(s_devD, s_polar,
                              r,
                              A, B,
                              npart, nrot, nk,
                              alpha,
                              s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i) == "true"){
          rc = print_s_polar_struct(s_polar);}

      }
  }
  else if(TRANSA=='P') {
    if(TRANSB=='N')
      {
        if(lda < npart ) {return rc = RC_FAIL;}
        if(ldb < nrot ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/

        rc = polarft_corr_P_N(s_devD, s_polar,
                              r,
                              A, B,
                              npart, nrot, nk,
                              alpha,
                              s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i) == "true"){
          rc = print_s_polar_struct(s_polar);}

      }
  }
  else if(TRANSA=='X') {
    if(TRANSB=='N')
      {
        if(lda < npart ) {return rc = RC_FAIL;}
        if(ldb < nrot ) {return rc = RC_FAIL;}
        /*====================================================================
          ===================C = alpha * A * B             ===================
          ==================================================================*/

        rc = polarft_corr_X_N(s_devD, s_polar,
                              r,
                              A, B,
                              npart, nrot, nk,
                              alpha,
                              s_bench, s_debug_gpu);

        if(get_bool(s_debug_gpu->debug_i) == "true"){
          rc = print_s_polar_struct(s_polar);}

      }
  }
  
  cuCtxSynchronize();

#else

  rc = get_warning_message_corr_Hadmr_gpu();
  
#endif /* CUDA */

  return rc;

}/* polarft_corr_Hadmr_tesla_gpu */
