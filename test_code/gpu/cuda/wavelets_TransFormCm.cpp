/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2016
 *                           Modified: 25th of September 2016
 *
 *      September 2016
 *
 *   common code to perform wavelet transform in 1 and 2D.
 *
 * @precisions normal z -> s d c
 */

#include "waveletsD.h"
//#include "get_fft123D_cpu.h"
//#include "get_fft123D_gpu.h" TODO: need to insert the cuFFT later if needed
// in that case fftw3.h mist be taken care of.
/* CPU header files */
#include <fftw3.h>

using namespace std;

/*Header for the MACOSX clock environment */
#if defined (MACOSX)
#endif

#ifdef __cplusplus
extern "C" {
#endif
/******************************************************************************/
/**
 *  FORTRAN API - timming functions (simple interface)
 **/
/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculate the 1D inverse wavelet transform based on a given 
   wavelet:
   sig = input signal
   name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
   db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
   bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
   bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5.
*/
  /* function to return the maxvalue fo the 2D array*/
  int maxval(vector<vector<double> > &arr, double &max) {
    int rc = RC_SUCCESS;
    max = 0;
    for (unsigned int i =0; i < arr.size(); i++) {
      for (unsigned int j =0; j < arr[0].size(); j++) {
	if (max <= arr[i][j]) {max = arr[i][j];}
      }
    }
    return rc;
  }
  /* function to return the maxvalue fo the 1D array*/
  int maxval1(vector<double> &arr, double &max) {
    int rc = RC_SUCCESS;
    max = 0;
    for (unsigned int i =0; i < arr.size(); i++) {
      if (max <= arr[i]) {max = arr[i];}
    }
    return rc;
  }

  int downsamp(vector<double> &sig, int M, vector<double> &sig_d){
    int rc = RC_SUCCESS;
    int len = sig.size();
    double len_n = ceil( (double) len / (double) M);
    for (int i = 0; i < (int) len_n; i++) {
      double temp = sig[i*M];
      sig_d.push_back(temp);
    }
    return rc;
  }
  /* symmetric extension*/
  int symm_1D_ext(vector<double> &sig, int a) {
    int rc = RC_SUCCESS;
    unsigned int len = sig.size();
    for (int i =0; i < a; i++) {
      double temp1= sig[i * 2];
      double temp2= sig[len - 1];
      sig.insert(sig.begin(),temp1);
      sig.insert(sig.end(),temp2);
    }
    return rc;
  }
  /**/
  int upsamp(vector<double> &sig, int M, vector<double> &sig_u) {
    int rc = RC_SUCCESS;
    int len = sig.size();
    double len_n = ceil( (double) len * (double) M);
    
    for (int i = 0; i < (int) len_n; i++) {
      if ( i % M == 0) {
        double temp = sig[i / M];
        sig_u.push_back(temp);
      }
      else
        {
          sig_u.push_back(0);
        }
    }
    
    return rc;
  }
  /* vector sum */
  int vecsum(vector<double> &a, vector<double> &b, vector<double> &c){
    int rc = RC_SUCCESS;
    c.resize(a.size());
    transform (a.begin(), a.end(), b.begin(), b.begin(), op_sum);
    c = b;
    return rc;
  }
  
  double op_sum(double i, double j) {
    return (i+j);
  }
  /* outputs the dimensions of dwt*/
  int dwt_output_dim_sym(vector<int> &length,vector<int> &length2, int J) {
    int rc = RC_SUCCESS;
    unsigned int sz=length.size();
    int rows = length[sz-2];
    int cols = length[sz-1];
    for (int i =0; i < J; i++) {
      rows =(int) ceil((double) rows/ 2.0);
      cols =(int) ceil((double) cols/ 2.0);
    }
    for (int i =0; i < J + 1; i++) {
      length2.push_back(rows);
      length2.push_back(cols);
      rows = rows * 2;
      cols = cols*2;
    }
    return rc;
  }
  /*dispersion of the DWT*/
  int dispDWT(vector<double> &output,vector<vector<double> > &dwtdisp,
	      vector<int> &length , vector<int> &length2, int J) {
    int rc = RC_SUCCESS;
    int sum = 0;

    for (int iter =0; iter < J; iter++) {
      int d_rows=length[2*iter]-length2[2*iter];
      int d_cols=length[2*iter+1]-length2[2*iter + 1];

      int rows_n =length[2 * iter];
      int cols_n = length[2 * iter + 1];
      vector< vector<double> > dwt_output(2*rows_n, vector<double>(2*cols_n));
      if (iter == 0) {
	for(int i =0; i < rows_n; i++){
	  for (int j =0; j < cols_n; j++){
	    dwt_output[i][j]=output[i*cols_n + j];
	  }
	}
	
	for(int i =0; i < rows_n; i++){
	  for (int j = cols_n; j < cols_n * 2; j++){
	    dwt_output[i][j]= output[rows_n * cols_n + i * cols_n + (j - cols_n)];
	  }
	}

	for(int i = rows_n; i < rows_n * 2; i++){
	  for (int j =0; j < cols_n; j++){
	    dwt_output[i][j]=output[2 * rows_n * cols_n+ (i - rows_n) * cols_n + j];
	  }
	}

	
	for(int i = rows_n; i < rows_n * 2; i++){
	  for (int j = cols_n; j < cols_n * 2; j++){
	    dwt_output[i][j]=output[3 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
	  }
	}
      } else {
	for(int i =0; i < rows_n; i++){
	  for (int j = cols_n; j < cols_n * 2; j++){
	    dwt_output[i][j]= output[sum + i * cols_n + (j - cols_n)];
	  }
	}
	
	for(int i = rows_n; i < rows_n * 2; i++){
	  for (int j =0; j < cols_n; j++){
	    dwt_output[i][j]=output[sum + rows_n * cols_n+ (i - rows_n) * cols_n + j];
	  }
	}

	
	for(int i = rows_n; i < rows_n * 2; i++){
	  for (int j = cols_n; j < cols_n * 2; j++){
	    dwt_output[i][j]=output[sum + 2 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
	  }
	}
	
      }

      int rows_x = length2[2*iter];
      int cols_x =length2[2*iter +1];
      
      int d_cols2 = (int) ceil( (double) (d_cols - 1) / 2.0);
      int d_rows2 = (int) ceil( (double) (d_rows - 1) / 2.0);
      if (iter ==0) {
        for(int i =0; i < rows_x; i++){
	  for (int j =0; j < cols_x; j++){
	    if (i + d_rows -1 < 0){
	      dwtdisp[i][j]=0;
	    }
	    else if (j + d_cols -1 < 0){
	      dwtdisp[i][j]=0;
	    } else {
	      dwtdisp[i][j]=dwt_output[i+d_rows -1][j+d_cols -1];
	    }
	  }
	}
      }
      for(int i =0; i < rows_x; i++){
	for (int j = cols_x; j < cols_x * 2; j++){
	  if (i + d_rows2 < 0) {
	    dwtdisp[i][j]=0;
	  }
	  else if (j + 2* (d_cols -1) +1 > (signed) dwt_output[0].size() - 1){
	    dwtdisp[i][j]=0;
	  } else {
	    dwtdisp[i][j]= dwt_output[i+d_rows2 ][j + 2* (d_cols -1)+1 ];
	  }
	}
      }
      
      for(int i = rows_x; i < rows_x * 2; i++){
	for (int j =0; j < cols_x; j++){
	  if (i + 2* (d_rows -1) + 1 > (signed) dwt_output.size() - 1){
	    dwtdisp[i][j]=0;
	  }
	  else if (j + d_cols2 < 0){
	    dwtdisp[i][j]=0;
	  } else {	    
	    dwtdisp[i][j]=dwt_output[i+2 * (d_rows - 1) + 1 ][j+d_cols2 ];
	  }
	}
      }

      for(int i = rows_x; i < rows_x * 2; i++){
	for (int j = cols_x; j < cols_x * 2; j++){
	  
	  if (i +  (d_rows -1) + 1 + d_rows2 > (signed) dwt_output.size() - 1){
	    dwtdisp[i][j]=0;
	  }
	  else if (j + (d_cols -1) + 1 + d_cols2  > (signed) dwt_output[0].size() - 1){
	    dwtdisp[i][j]=0;
	  } else {
	    dwtdisp[i][j]=dwt_output[i +  (d_rows -1) + 1 + d_rows2 ][j + (d_cols -1) + 1 + d_cols2 ];
	  }
	}
      }
      if (iter == 0) {
	sum+= 4*rows_n*cols_n;
      } else {
	sum+= 3*rows_n * cols_n;
      }
      
    } 
    return rc;
  }

  /*Aliases*/
  extern "C" void maxval_() __attribute__((weak,alias("maxval")));
  extern "C" void maxval1_() __attribute__((weak,alias("maxval1")));
  extern "C" void downsamp_() __attribute__((weak,alias("downsamp")));
  extern "C" void symm_1d_ext_() __attribute__((weak,alias("symm_1D_ext")));
  extern "C" void upsamp_() __attribute__((weak,alias("upsamp")));
  extern "C" void vecsum_() __attribute__((weak,alias("vecsum")));
  extern "C" void op_sum_() __attribute__((weak,alias("op_sum")));
  extern "C" void dwt_output_dim_sym_() __attribute__((weak,alias("dwt_output_dim_sym")));

#ifdef __cplusplus
}
#endif
