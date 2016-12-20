/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2016
 *                           Modified: 25th of September 2016
 *
 *      September 2016
 *
 *   code to perform timing tests.
 *
 * @precisions normal z -> s d c
 */
#include "simple.h"
//#include "waveletsD.h"

/***************************************************************************/
/**
 *  FORTRAN API - wavelet functions (simple interface)
 **/

/* //////////////////////////////////////////////////////////////////////////// 
 -- C function wrapper function that gets the wavelet for a given type
 integer, float and double
*/
int test_transform_sig_2d_sym_c_(char *wname, int *J_in) {
  int rc = RC_SUCCESS;
  int J = *J_in;
  rc =  test_transform_sig_2d_sym(wname, J);
  return rc;
}
void TEST_TRANSFORM_SIG_2D_SYM_C_() __attribute__((weak,alias("test_transform_sig_2d_sym_c_")));
void test_transform_sig_2d_sym_c__() __attribute__((weak,alias("test_transform_sig_2d_sym_c_")));
void TEST_TRANSFORM_SIG_2D_SYM_C__() __attribute__((weak,alias("test_transform_sig_2d_sym_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function wrapper function that gets the wavelet for a given type
 integer, float and double
*/
int test_transform_sig_1d_sym_c_(char *wname, int *J_in) {
  int rc = RC_SUCCESS;
  int J = *J_in;
  rc =  test_transform_sig_1d_sym(wname, J);
  return rc;
}
void TEST_TRANSFORM_SIG_1D_SYM_C_() __attribute__((weak,alias("test_transform_sig_1d_sym_c_")));
void test_transform_sig_1d_sym_c__() __attribute__((weak,alias("test_transform_sig_1d_sym_c_")));
void TEST_TRANSFORM_SIG_1D_SYM_C__() __attribute__((weak,alias("test_transform_sig_1d_sym_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function wrapper function that gets the wavelet for a given type
 integer, float and double
*/
int get_wavelet_coef_c_(char *name) {
  int rc = RC_SUCCESS;
  rc = get_wavelet_coef(name);
  return rc;
}
void GET_WAVELET_COEF_C_() __attribute__((weak,alias("get_wavelet_coef_c_")));
void get_wavelet_coef_c__() __attribute__((weak,alias("get_wavelet_coef_c_")));
void GET_WAVELET_COEF_C__() __attribute__((weak,alias("get_wavelet_coef_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function wrapper function to get the coefficients of the wavelets
 integer, float and double
*/
int filtcoef_c_(char *name, double *lp1, double *hp1,
                double *lp2, double *hp2) {
  int rc = RC_SUCCESS;

  //TODO:need recast of double type array to a vector<double> type array
  
  return rc;
}
void FILTCOEF_C_() __attribute__((weak,alias("filtcoef_c_")));
void filtcoef_c__() __attribute__((weak,alias("filtcoef_c_")));
void FILTCOEF_C__() __attribute__((weak,alias("filtcoef_c_")));
