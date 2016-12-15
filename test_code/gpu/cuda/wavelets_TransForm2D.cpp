/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2016
 *                           Modified: 25th of September 2016
 *
 *      September 2016
 *
 *   code to perform 2D wavelet transforms.
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
   -- function to calculate the 2D inverse wavelet transform based on a given
   wavelet:
   sig = input signal
   name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
   db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
   bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
   bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5. 
*/

  int inv_dwt_2D_sym(vector<double>  &dwtop,vector<double> &flag, string name,
		     vector<vector<double> > &idwt_output, vector<int> &length){
    int rc = RC_SUCCESS;
    int J =(int) flag[0];
    int rows =length[0];
    int cols =length[1];

    int sum_coef =0;
    vector<double> lp1,hp1,lp2,hp2;
    rc = filtcoef(name,lp1,hp1,lp2,hp2);
    unsigned int lf = lp1.size();
    vector<vector<double> >  cLL(rows, vector<double>(cols));


    for (int iter=0; iter < J; iter++) {

      int rows_n = length[2*iter];
      int cols_n = length[2*iter + 1];

      vector<vector<double> >  cLH(rows_n, vector<double>(cols_n));
      vector<vector<double> >  cHL(rows_n, vector<double>(cols_n));
      vector<vector<double> >  cHH(rows_n, vector<double>(cols_n));

      for (int i = 0 ; i < rows_n; i++ ){
	for (int j = 0; j < cols_n; j++){
	  if (iter == 0) {
	    cLL[i][j] = dwtop[sum_coef+ i * cols_n + j];
	    
	    cLH[i][j] = dwtop[sum_coef+ rows_n * cols_n+ i * cols_n + j];
	    
	    cHL[i][j] = dwtop[sum_coef+ 2 * rows_n * cols_n + i * cols_n + j];
	    
	    cHH[i][j] = dwtop[sum_coef+ 3* rows_n * cols_n + i * cols_n + j];
	  } else {
	    
	    cLH[i][j] = dwtop[sum_coef+  i * cols_n + j];
	    
	    cHL[i][j] = dwtop[sum_coef+ rows_n * cols_n + i * cols_n + j];
	    
	    cHH[i][j] = dwtop[sum_coef+ 2* rows_n * cols_n + i * cols_n + j];
	    
	  }
	}
      }
      //      temp_A = cLL;
      //  	idwt2_sym(nm,idwt_output2, cA, cH,cV,cD);
      
      unsigned int len_x = cLH.size();
      unsigned int len_y = cLH[0].size();
      
      // Row Upsampling and Column Filtering at the first LP Stage
      vector<vector<double> > cL(2 *len_x - lf + 2,vector<double>(len_y ));
      vector<vector<double> > cH(2 * len_x - lf +2,vector<double>(len_y ));

      if (iter ==0) {
        for (unsigned int j =0; j < len_y; j++) {
	  
	  vector<double> sigLL,sigLH,oup;

	  for (unsigned int i=0;i <  len_x;i++) {
	    
	    double temp1 = cLL[i][j];
	    double temp2 = cLH[i][j];
	    sigLL.push_back(temp1);
	    sigLH.push_back(temp2);
	  }
	  inv_dwt_1D_sym_m(name,oup,sigLL,sigLH);
	  
	  for (int i=0;i < (int) oup.size();i++) {
	    cL[i][j] = oup[i];
	  }
	  
	}
      } else{
	unsigned int rows1 =cLH.size();
	unsigned int cols1 =cLH[0].size();
	
	for (unsigned int j =0; j < cols1;j++){
	  vector<double> temp_L1,temp_L2,oup;
	  for (unsigned int i =0; i < rows1; i++){
	    double temp = cLL[i][j];
	    temp_L1.push_back(temp);
	    
	    double temp2 = cLH[i][j];
	    temp_L2.push_back(temp2);
	  }
	  inv_dwt_1D_sym_m(name,oup,temp_L1,temp_L2);
	  
	  for (unsigned int i =0; i < oup.size(); i++){
	    cL[i][j]=oup[i];
	  }
	  
	}
      }
      
      for (unsigned int j =0; j < len_y; j++) {
	vector<double> sigHL,sigHH,oup2;
	for (unsigned int i=0;i <  len_x;i++) {
	  double temp3 = cHL[i][j];
	  double temp4 = cHH[i][j];
	  sigHL.push_back(temp3);
	  sigHH.push_back(temp4);
	}

	inv_dwt_1D_sym_m(name,oup2,sigHL,sigHH);
	
	for (int i=0;i < (int) oup2.size();i++) {
	  cH[i][j] = oup2[i];
	}
      }

      vector<vector<double> >signal(2*len_x-lf +2,vector<double>(2*len_y-lf+2));
      for (unsigned int i =0; i < 2 * len_x - lf +2; i++) {
	vector<double> sigL,sigH,oup;
	for (unsigned int j=0;j <  len_y;j++) {
	  double temp5 = cL[i][j];
	  double temp6 = cH[i][j];
	  sigL.push_back(temp5);
	  sigH.push_back(temp6);
	}

	inv_dwt_1D_sym_m(name,oup,sigL,sigH);
	
	for (int j=0;j < (int) oup.size();j++) {
	  signal[i][j] = oup[j];
	}
      }
      idwt_output = signal;
      if (iter ==0) {
	sum_coef+= 4 *rows_n * cols_n;
      } else {
	sum_coef+= 3 *rows_n * cols_n;
      }
      cLL = signal;
    }
    return rc;
  }/*end of  inv_dwt_2D_sym method*/
 
/* /////////////////////////////////////////////////////////////////////////////
   -- function to calculate the 2D wavelet transform based on a given wavelet:
    sig = input signal
    name =haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10, db11, db12,
    db13, db14, db15, bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,
    bior2.6,bior2.8, bior3.1, bior3.3, bior3.5, bior3.7, bior3.9,
    bior4.4, bior5.5, bior6.8, coif1, coif2, coif3, coif4, coif5.
    
*/
  int dwt_2D_sym(vector<vector<double> > &origsig, int J,
		 string name, vector<double> &dwt_output,
		 vector<double> &flag , vector<int> &length) {
    int rc = RC_SUCCESS;
    vector<vector<double> >  sig = origsig;
    int rows_n = sig.size();
    int cols_n = sig[0].size();
    vector<vector<double> > original_copy(rows_n,vector<double>(cols_n));

    original_copy = sig;
    int Max_Iter;
    Max_Iter = min((int) ceil(log( double(sig.size()))/log (2.0)),(int) ceil(log( double(sig[0].size()))/log (2.0)));
    if ( Max_Iter < J) {
      cout << J << " Iterations are not possible with signals of this dimension "  << endl;
      exit(1);
    }
    vector<double> lp1,hp1,lp2,hp2;

    flag.push_back(double(J));
    length.insert(length.begin(),cols_n);
    length.insert(length.begin(),rows_n);
    // Flag Values
    /*
      double temp = (double) (sig2.size() - sig.size()); // Number of zeropad rows
      flag.push_back(temp);
      double temp2 = (double) (sig2[0].size() - sig[0].size());// Number of zpad cols
      flag.push_back(temp2);
      flag.push_back((double) J); // Number of Iterations
    */
    int sum_coef = 0;
    for (int iter = 0; iter < J; iter++) {
      filtcoef(name,lp1,hp1,lp2,hp2);
      unsigned int lf = lp1.size();

      rows_n =(int) floor((double)(rows_n + lf -1)/2);
      cols_n =(int) floor((double) (cols_n + lf -1)/2);
      length.insert(length.begin(),cols_n);
      length.insert(length.begin(),rows_n);

      vector<vector<double> >  cA(rows_n, vector<double>(cols_n));
      vector<vector<double> >  cH(rows_n, vector<double>(cols_n));
      vector<vector<double> >  cV(rows_n, vector<double>(cols_n));
      vector<vector<double> >  cD(rows_n, vector<double>(cols_n));

      if (iter == 0) {
	dwt2_2D_sym(name,original_copy,cA,cH,cV,cD);
      } else {
	dwt2_2D_sym(name,original_copy,cA,cH,cV,cD);
      }
      vector<double>   temp_sig2;

      original_copy = cA;
      if (iter == J-1) {
	for(int i =0; i < rows_n; i++){
	  for (int j =0; j < cols_n; j++){
	    double temp=cA[i][j];
	    temp_sig2.push_back(temp);
	  }
	}
      }
      for(int i =0; i < rows_n; i++){
	for (int j = cols_n; j < cols_n * 2; j++){
	  double temp =cH[i][j - cols_n];
	  temp_sig2.push_back(temp);
	}
      }
      for(int i = rows_n; i < rows_n * 2; i++){
	for (int j =0; j < cols_n; j++){
	  double temp=cV[i - rows_n][j];
	  temp_sig2.push_back(temp);
	}
      }
      for(int i = rows_n; i < rows_n * 2; i++){
	for (int j = cols_n; j < cols_n * 2; j++){
	  double temp =cD[i- rows_n][j - cols_n];
	  temp_sig2.push_back(temp);
	}
      }
      dwt_output.insert(dwt_output.begin(),temp_sig2.begin(),temp_sig2.end());
      sum_coef += 4 * rows_n * cols_n;
    }
    /*
      ofstream dwt2out("dwt2out.dat");
      for (unsigned int i= 0; i < dwt_output.size(); i++){
      dwt2out << dwt_output[i] <<endl;
      }
    */
    return rc;
  } /*end of the method  dwt_2D_sym */

  /* submethod for the 2D symmetric wavelet*/
  int dwt2_2D_sym(string name,vector<vector<double> > &signal,
		  vector<vector<double> >  &cLL,
		  vector<vector<double> >  &cLH,
		  vector<vector<double> >  &cHL,
		  vector<vector<double> > &cHH) {
    //Analysis
    int rc = RC_SUCCESS;
    int rows = signal.size();
    int cols = signal[0].size();
    int cols_lp1 = cLL[0].size();
    int cols_hp1 = cLL[0].size();
    vector<double> lp1,hp1,lp2,hp2;
    rc = filtcoef(name, lp1,hp1,lp2,hp2);
    vector<vector<double> > lp_dn1(rows, vector<double>( cols_lp1));
    vector<vector<double> > hp_dn1(rows, vector<double>( cols_hp1));

    // Implementing row filtering and column downsampling in each branch.
    for (int i =0; i < rows; i++) {
      vector<double> temp_row,oup_lp,oup_hp;
      for (int j=0;j <  cols;j++) {
	double temp = signal[i][j];
	temp_row.push_back(temp);
      }
      rc = dwt1_1D_sym_m(name,temp_row,oup_lp,oup_hp);
      for (int j=0;j < (int) oup_lp.size();j++) {
	lp_dn1[i][j] = oup_lp[j];
	hp_dn1[i][j] = oup_hp[j];
      }  
    }

    cols =cols_lp1;
    // Implementing column filtering and row downsampling in Low Pass branch.
    for (int j =0; j < cols; j++) {
      vector<double> temp_row3,oup_lp,oup_hp;
      for (int i=0;i <  rows;i++) {
	double temp = lp_dn1[i][j];
	temp_row3.push_back(temp);
      }
      rc = dwt1_1D_sym_m(name,temp_row3,oup_lp,oup_hp);
      for (int i=0;i < (int) oup_lp.size();i++) {
	cLL[i][j] = oup_lp[i];
	cLH[i][j] = oup_hp[i];
      }
    }

    // Implementing column filtering and row downsampling in High Pass branch.
    for (int j =0; j < cols; j++) {
      vector<double> temp_row5,oup_lp,oup_hp;
      for (int i=0;i <  rows;i++) {
	double temp = hp_dn1[i][j];
	temp_row5.push_back(temp);
      }
      rc = dwt1_1D_sym_m(name,temp_row5,oup_lp,oup_hp);
      for (int i=0;i < (int) oup_lp.size();i++) {
	cHL[i][j] = oup_lp[i];
	cHH[i][j] = oup_hp[i];
      }
    }
    return rc;
  } /* end of the dwt2_2D_sym */
 
  /*Aliases*/
  extern "C" void dwt_2d_sym_() __attribute__((weak,alias("dwt_2D_sym")));
  extern "C" void dwt2_2d_sym_() __attribute__((weak,alias("dwt2_2D_sym")));
  extern "C" void inv_dwt_2d_sym_() __attribute__((weak,alias("inv_dwt_2D_sym")));


  
#ifdef __cplusplus
}
#endif
