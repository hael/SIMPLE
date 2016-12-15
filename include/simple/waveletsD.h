/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 10th of August 2016
 *                           Modified: 10th of August 2016
 *
 *      August 2016
 *
 *   code to perform timing tests.
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
/*image processing library OpenCV*/
//#include "cv.h"
//#include "highgui.h"
//#include "cxcore.h"

using namespace std;

#ifndef _WAVELETSD_H_
#define _WAVELETSD_H_
#define PRECISION_z

typedef struct {
  string name;
  vector<double> lp1; vector<double> hp1;
  vector<double> lp2; vector<double> hp2;
} wavelets_coef_t;

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- routines used to interface back to fortran
   *
   *  FORTRAN API -
   * 
   */
  /*testers function calls*/
  int test_transform_sig_1d_sym(char *wname_in, int J);
  int test_transform_sig_2d_sym(char *wname_in, int J);
  /*getter for the filter coefficients*/
  int get_wavelet_coef(char *wname_in);
  /*filter coefficients*/
  int filtcoef(string name, vector<double> &lp1, vector<double> &hp1,
               vector<double> &lp2, vector<double> &hp2);
  /* common functions for the testers*/
  int maxval(vector<vector<double> > &arr, double &max);
  int maxval1(vector<double> &arr, double &max);
  /*2D symmetric transform*/
  int dwt_2D_sym(vector<vector<double> > &origsig, int J,
		 string name, vector<double> &dwt_output,
		 vector<double> &flag , vector<int> &length);
  int dwt2_2D_sym(string name,vector<vector<double> > &signal,
		  vector<vector<double> >  &cLL,
		  vector<vector<double> >  &cLH,
		  vector<vector<double> >  &cHL,
		  vector<vector<double> > &cHH);
  int inv_dwt_2D_sym(vector<double>  &dwtop,vector<double> &flag, string name,
		     vector<vector<double> > &idwt_output, vector<int> &length);
  /* 1D symmetric transform*/
  int dwt_1D_sym(vector<double> &sig, int J, string name,
                 vector<double> &dwt_output,
                 vector<double> &flag, vector<int> &length );
  int dwt1_1D_sym(string wname, vector<double> &signal,
                  vector<double> &cA, vector<double> &cD);
  int dwt1_1D_sym_m(string wname, vector<double> &signal,
		    vector<double> &cA, vector<double> &cD);
  int inv_dwt_1D_sym(vector<double> &dwtop,vector<double> &flag, string name,
                     vector<double> &idwt_output, vector<int> &length);
  int inv_dwt_1D_sym_m(string wname, vector<double> &idwt_output,
		       vector<double> &app, vector<double> &detail);
  /* convoluters*/
  int convolate_fft(vector<double> &a, vector<double> &b, vector<double> &c);
  int convolate_fft_m(vector<double> &a, vector<double> &b, vector<double> &c);
  int myconvolate_fft(vector<double> &a, vector<double> &b, vector<double> &c);
  /* common methods */
  int symm_1D_ext(vector<double> &sig, int a);
  int upsamp(vector<double> &sig, int M, vector<double> &sig_u);
  int vecsum(vector<double> &a, vector<double> &b, vector<double> &c);
  double op_sum(double i, double j);
  int downsamp(vector<double> &sig, int M, vector<double> &sig_d);
  int dwt_output_dim_sym(vector<int> &length,vector<int> &length2, int J);
  int dispDWT(vector<double> &output,vector<vector<double> > &dwtdisp,
	      vector<int> &length , vector<int> &length2, int J);
  //int transform_sig_1D(char *name_in );
  /* 1D transform */
  int dwt_1D(vector<double> &sig, int J, string name,
             vector<double> &dwt_output, vector<double> &flag,
             vector<int> &length );

  /* class to return the coefficients for the wavelets*/
  class wavelets_coef {
  private:
  public:
    /*global variables*/
    string name_in;
    vector<double> lp1_in; vector<double> hp1_in;
    vector<double> lp2_in; vector<double> hp2_in;
    /*setters*/
    void set_wavelets_name(string name);
    void set_wavelets_lp1(vector<double> lp1);
    void set_wavelets_hp1(vector<double> hp1);
    void set_wavelets_lp2(vector<double> lp2);
    void set_wavelets_hp2(vector<double> hp2);
    /*getters*/
    string get_wavelets_name();
    vector<double> get_wavelets_lp1();
    vector<double> get_wavelets_hp1();
    vector<double> get_wavelets_lp2();
    vector<double> get_wavelets_hp2();
    /*constructor*/
    wavelets_coef(string name_in,
                  vector<double> lp1_in, vector<double> hp1_in,
                  vector<double> lp2_in, vector<double> hp2_in);    
    /*destructor*/
    ~wavelets_coef();
  };

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _WAVELETSD_H_ */
