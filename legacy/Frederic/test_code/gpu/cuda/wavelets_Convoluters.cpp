/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2016
 *                           Modified: 25th of September 2016
 *
 *      september 2016
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
/* Fourier transform plan declaratiomns*/
fftw_plan plan_forward_inp,plan_forward_filt, plan_backward;
static unsigned int transient_size_of_fft = 0;
/*Header for the MACOSX clock environment */
#if defined (MACOSX)
#endif


//convfftm
//convolate_fft

#ifdef __cplusplus
extern "C" {
#endif
/******************************************************************************/
/**
 *  FORTRAN API - timming functions (simple interface)
 **/

/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculate the 1D inverse wavelet transform based on a given wavelet:
    sig = input signal
    name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
    db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
    bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
    bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5.
    
*/
  
  /* method to convulate the fourier transform to the signal*/
  int convolate_fft(vector<double> &a, vector<double> &b, vector<double> &c) {
    int rc = RC_SUCCESS;

    fftw_complex *inp_data, *filt_data, *inp_fft, *filt_fft;
    fftw_complex *temp_data, *temp_ifft;
    fftw_plan plan_forward_inp,plan_forward_filt, plan_backward;

    int sz = a.size() + b.size() - 1;
    inp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    inp_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    temp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    temp_ifft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    plan_forward_inp  = fftw_plan_dft_1d( sz, inp_data, inp_fft, FFTW_FORWARD, FFTW_ESTIMATE );
    plan_forward_filt  = fftw_plan_dft_1d( sz, filt_data, filt_fft, FFTW_FORWARD, FFTW_ESTIMATE );
    plan_backward = fftw_plan_dft_1d( sz, temp_data, temp_ifft, FFTW_BACKWARD, FFTW_ESTIMATE );

    for (unsigned int i =0; i < sz; i++) {
      if (i < a.size()) {
        inp_data[i][0] = a[i];
      } else {
        inp_data[i][0] = 0.0;
        
      }
      inp_data[i][1] = 0.0;
      if (i < b.size()) {
        filt_data[i][0] = b[i];
      } else {
        filt_data[i][0] = 0.0;
        
      }
      filt_data[i][1] = 0.0;

    }

    fftw_execute(plan_forward_inp);
    fftw_execute(plan_forward_filt);

    for (unsigned int i =0; i < sz; i++){
      temp_data[i][0] = inp_fft[i][0]*filt_fft[i][0] -
                        inp_fft[i][1]*filt_fft[i][1] ;
      temp_data[i][1] = inp_fft[i][0]*filt_fft[i][1] +
                        inp_fft[i][1]*filt_fft[i][0];
    }

    fftw_execute(plan_backward);

    for (unsigned int i = 0; i < sz; i++) {
        double temp1;
        temp1 = temp_ifft[i][0] / (double) sz;
        c.push_back(temp1);
    }
    
    fftw_free(inp_data);
    fftw_free(filt_data);
    fftw_free(inp_fft);
    fftw_free(filt_fft);
    fftw_free(temp_data);
    fftw_free(temp_ifft);
    fftw_destroy_plan(plan_forward_inp);
    fftw_destroy_plan(plan_forward_filt);
    fftw_destroy_plan(plan_backward);

    return rc;
  }
  /*convoluter for multi dimensional wavelets */
  int convolate_fft_m(vector<double> &a, vector<double> &b, vector<double> &c) {
    int rc = RC_SUCCESS;
    fftw_complex *inp_data, *filt_data, *inp_fft, *filt_fft;
    fftw_complex *temp_data, *temp_ifft;

    unsigned int sz = a.size() + b.size() - 1;
    inp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    inp_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    filt_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    temp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
    temp_ifft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );

    if (sz != transient_size_of_fft) {

      if (transient_size_of_fft != 0) {
	fftw_destroy_plan(plan_forward_inp);
	fftw_destroy_plan(plan_forward_filt);
	fftw_destroy_plan(plan_backward);
      }
      plan_forward_inp  =
	fftw_plan_dft_1d( sz, inp_data, inp_fft, FFTW_FORWARD, FFTW_MEASURE );
      plan_forward_filt  =
	fftw_plan_dft_1d( sz, filt_data, filt_fft, FFTW_FORWARD, FFTW_MEASURE );
      plan_backward =
	fftw_plan_dft_1d( sz, temp_data, temp_ifft, FFTW_BACKWARD, FFTW_MEASURE );
      transient_size_of_fft = sz;
    }

    for (unsigned int i =0; i < sz; i++) {
      if (i < a.size()) {
        inp_data[i][0] = a[i];
      } else {
	inp_data[i][0] = 0.0;
      }
      inp_data[i][1] = 0.0;
      if (i < b.size()) {
	filt_data[i][0] = b[i];
      } else {
	filt_data[i][0] = 0.0;
      }
      filt_data[i][1] = 0.0;
    }

    fftw_execute_dft( plan_forward_inp,inp_data, inp_fft);
    fftw_execute_dft( plan_forward_filt,filt_data, filt_fft);

    for (unsigned int i =0; i < sz; i++){
      temp_data[i][0] =
	inp_fft[i][0]*filt_fft[i][0] - inp_fft[i][1]*filt_fft[i][1];
      temp_data[i][1] =
	inp_fft[i][0]*filt_fft[i][1] + inp_fft[i][1]*filt_fft[i][0];
    }

    fftw_execute_dft( plan_backward, temp_data, temp_ifft);
    for (unsigned int i = 0; i < sz; i++) {
      double temp1;
      temp1 = temp_ifft[i][0] / (double) sz;
      c.push_back(temp1);
    }
    fftw_free(inp_data);
    fftw_free(filt_data);
    fftw_free(inp_fft);
    fftw_free(filt_fft);
    fftw_free(temp_data);
    fftw_free(temp_ifft);

    return rc;
  } /* end convolate_fft_m method */ 

  /*Aliases*/
  extern "C" void convolate_fft_() __attribute__((weak,alias("convolate_fft")));
  extern "C" void convolate_fft_m_() __attribute__((weak,alias("convolate_fft_m")));


#ifdef __cplusplus
}
#endif
