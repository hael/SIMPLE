/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2016
 *                           Modified: 25th of september 2016
 *
 *      September 2016
 *
 *   code to perform timing tests.
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
 *  FORTRAN API - Discrete wavelet transform 1D and 2D case (simple interface)
 **/

/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculate the 1D inverse wavelet transform based on a given
   wavelet for multi diemnsional case.
    sig = input signal
    name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
    db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
    bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
    bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5. 
*/

  int inv_dwt_1D_sym_m(string wname, vector<double> &idwt_output,
		       vector<double> &app, vector<double> &detail) {
    int rc = RC_SUCCESS;
    int U = 2; /* Upsampling Factor */
    vector<double> lpd1,hpd1, lpr1, hpr1;

    rc = filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
    int lf = lpr1.size();

    /* Operations in the Low Frequency branch of the Synthesis Filter Bank*/
    vector<double> X_lp;
    vector<double> cA_up;
    rc = upsamp(app, U,cA_up );
    cA_up.pop_back();
    rc = convolate_fft_m(cA_up, lpr1, X_lp);

    /* Operations in the High Frequency branch of the Synthesis Filter Bank*/

    vector<double> X_hp;
    vector<double> cD_up;
    rc = upsamp(detail, U, cD_up);
    cD_up.pop_back();
    rc = convolate_fft_m(cD_up, hpr1, X_hp);

    rc = vecsum(X_lp,X_hp,idwt_output);

    idwt_output.erase(idwt_output.begin(),idwt_output.begin()+lf-2);
    idwt_output.erase(idwt_output.end()-(lf - 2),idwt_output.end());

   return rc;
}
  
/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculate the 1D inverse wavelet transform based on a given
   wavelet:
    sig = input signal
    name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
    db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
    bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
    bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5. 
*/
  int inv_dwt_1D_sym(vector<double> &dwtop,vector<double> &flag, string name,
                     vector<double> &idwt_output, vector<int> &length) {
    int rc = RC_SUCCESS;
    int J =(int) flag[1];
    unsigned int lf;

    vector<double> app;
    vector<double> detail;
    unsigned int app_len = length[0];
    unsigned int det_len = length[1];

    vector<double>::iterator dwt;
    dwt = dwtop.begin();
    app.assign(dwt,dwtop.begin()+app_len);
    detail.assign(dwtop.begin()+app_len, dwtop.begin()+ 2* app_len);

    for (int i = 0; i < J; i++) {

      int U = 2; // Upsampling Factor
      vector<double> lpd1,hpd1, lpr1, hpr1;

      rc = filtcoef(name,lpd1,hpd1,lpr1,hpr1);
      lf = lpr1.size();

      // Operations in the Low Frequency branch of the Synthesis Filter Bank
      vector<double> X_lp;
      vector<double> cA_up;
      rc = upsamp(app, U,cA_up );
      cA_up.pop_back();
      rc = convolate_fft(cA_up, lpr1, X_lp);

      // Operations in the High Frequency branch of the Synthesis Filter Bank

      vector<double> X_hp;
      vector<double> cD_up;
      rc = upsamp(detail, U, cD_up);
      cD_up.pop_back();
      rc = convolate_fft(cD_up, hpr1, X_hp);

      app_len += det_len;
      vecsum(X_lp,X_hp,idwt_output);

      idwt_output.erase(idwt_output.begin(),idwt_output.begin()+lf-2);
      idwt_output.erase(idwt_output.end()-(lf - 2),idwt_output.end());

      app.clear();
      detail.clear();
      if ( i < J - 1 ) {
        det_len = length[i+2];
        for (unsigned int l = 0; l < det_len;l++) {
          double temp = dwtop[app_len + l];
          detail.push_back(temp);
        }
      }
      app = idwt_output;

      for (int iter1 = 0; iter1 < (int) (app.size() - det_len);iter1++) {
        app.pop_back();
      }
    }

    /* Remove ZeroPadding */

    int zerop =(int) flag[0];
    idwt_output.erase(idwt_output.end()- zerop,idwt_output.end());
    return 0;
  }

/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculatet the 1D wavelet transform based on a given wavelet:
    sig = input signal
    name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
    db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
    bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
    bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5.
    
*/
  int dwt_1D_sym(vector<double> &signal, int J, string name,
                 vector<double> &dwt_output,
                 vector<double> &flag, vector<int> &length ) {
    int rc = RC_SUCCESS;
    
    unsigned int temp_len = signal.size();
    if ( (temp_len % 2) != 0) {
      double temp =signal[temp_len - 1];
      signal.push_back(temp);
      flag.push_back(1);
      temp_len++;
    } else {
      flag.push_back(0);
    }
    length.push_back(temp_len);
    flag.push_back(double(J));

    vector<double> original_copy, appx_sig, det_sig;
    original_copy = signal;

    //  Storing Filter Values for GnuPlot
    vector<double> lp1,hp1,lp2,hp2;
    filtcoef(name,lp1,hp1,lp2,hp2);
    for (int iter = 0; iter < J; iter++) {
      rc = dwt1_1D_sym(name,signal, appx_sig, det_sig);
      dwt_output.insert(dwt_output.begin(),det_sig.begin(),det_sig.end());
      int l_temp = det_sig.size();
      length.insert(length.begin(),l_temp);
      
      if (iter == J-1 ) {
        dwt_output.insert(dwt_output.begin(),appx_sig.begin(),appx_sig.end());
        int l_temp = appx_sig.size();
        length.insert(length.begin(),l_temp);
      }
      signal.clear();
      signal = appx_sig;
      appx_sig.clear();
      det_sig.clear();
    }
    signal = original_copy;
    
    return rc;
  }/*end of dwt_1D_sym*/
  
  /* submethod for the 1D symmetric wavelet*/
  int dwt1_1D_sym(string wname, vector<double> &signal,
                  vector<double> &cA, vector<double> &cD) {
    int rc = RC_SUCCESS;
    vector<double> lp1, hp1, lp2, hp2;

    rc = filtcoef(wname,lp1,hp1,lp2,hp2);
    int D = 2; // Downsampling Factor is 2
    int lf = lp1.size();
    rc = symm_1D_ext(signal,lf-1);

    vector<double> cA_undec;
    //sig value
    rc = convolate_fft(signal,lp1,cA_undec);
    cA_undec.erase(cA_undec.begin(),cA_undec.begin()+lf);
    cA_undec.erase(cA_undec.end()-lf+1,cA_undec.end());
    rc = downsamp(cA_undec, D, cA);
    //High Pass Branch Computation
    vector<double> cD_undec;
    rc = convolate_fft(signal,hp1,cD_undec);
    cD_undec.erase(cD_undec.begin(),cD_undec.begin()+lf);
    cD_undec.erase(cD_undec.end()-lf+1,cD_undec.end());
    rc = downsamp(cD_undec,D,cD);

    rc = filtcoef(wname,lp1,hp1,lp2,hp2);

    return rc;
  }/* dwt1_1D_sym */

  /* 1D wavelet transform for multi dimensional called from the 2D and higher*/
  int dwt1_1D_sym_m(string wname, vector<double> &signal,
		    vector<double> &cA, vector<double> &cD) {
    int rc = RC_SUCCESS;
    vector<double> lp1, hp1, lp2, hp2;

    rc = filtcoef(wname,lp1,hp1,lp2,hp2);
    int D = 2; // Downsampling Factor is 2
    int lf = lp1.size();
    symm_1D_ext(signal,lf-1);

    vector<double> cA_undec;
    //sig value
    rc = convolate_fft_m(signal,lp1,cA_undec);
    cA_undec.erase(cA_undec.begin(),cA_undec.begin()+lf);
    cA_undec.erase(cA_undec.end()-lf+1,cA_undec.end());
    downsamp(cA_undec, D, cA);

    //High Pass Branch Computation

    vector<double> cD_undec;
    rc = convolate_fft_m(signal,hp1,cD_undec);
    cD_undec.erase(cD_undec.begin(),cD_undec.begin()+lf);
    cD_undec.erase(cD_undec.end()-lf+1,cD_undec.end());
    downsamp(cD_undec,D,cD);

    filtcoef(wname,lp1,hp1,lp2,hp2);

    return rc;
  }/* end of dwt1_1D_sym_m method*/

  /* 1D disctrete wavelet transform the non symmetric case */
  int dwt_1D(vector<double> &sig, int J, string name, vector<double> &dwt_output
             , vector<double> &flag, vector<int> &length ) {
    int rc = RC_SUCCESS;
    //TODO: need to implement the 1D discrete wavelet transform
    return rc;
  }

  /*Aliases*/
  extern "C" void dwt_1d_() __attribute__((weak,alias("dwt_1D")));
  extern "C" void dwt_1d_sym_() __attribute__((weak,alias("dwt_1D_sym")));
  extern "C" void dwt1_1d_sym_m_() __attribute__((weak,alias("dwt1_1D_sym_m")));
  extern "C" void inv_dwt_1d_sym_() __attribute__((weak,alias("inv_dwt_1D_sym")));
  extern "C" void inv_dwt_1d_sym_m_() __attribute__((weak,alias("inv_dwt_1D_sym_m")));
  
#ifdef __cplusplus
}
#endif
