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

#include "waveletsD.h"
using namespace std;
//using namespace cv;

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
  int test_transform_sig_2d_sym(char *wname_in, int J) {
    /* in the argument list add: file name or data coming in, J, */
    int rc = RC_SUCCESS;
    string name = string(wname_in);
    printf(ANSI_COLOR_BRIGHT_YELLOW);
    std::cout<<"Wavelet 2D transform for wavelet: "
      ANSI_COLOR_BRIGHT_RED<<name<<std::endl;
    printf(ANSI_COLOR_RESET);
    /*
    IplImage* img = cvLoadImage("Lena.png");

    //checking if the image exits or not if not throw error and stopper
    if (!img){
      cout << " Can't read Image. Try Different Format." << endl;
      exit(1);
    }
    // getting the diemensions of the image
    int height, width;
    height = img->height;
    width = img->width;
    int nc = img->nChannels;
    int pix_depth = img->depth;
    CvSize size;
    size.width =width;
    size.height=height;
    // print image info
    cout << "depth: " << pix_depth <<  " Channels: " << nc << endl;
    // tag original image and pop it to screen
    cvNamedWindow("Original Image", CV_WINDOW_AUTOSIZE);
    cvShowImage("Original Image", img);
    cvWaitKey();
    cvDestroyWindow("Original Image");
    cvSaveImage("orig.bmp",img);

    int rows =(int) height;
    int cols =(int) width;
    //getting the Mat class from the cv namespace
    Mat matimg(img);
    //creating a matrix image of type double
    vector<vector<double> > vec1(rows, vector<double>(cols));
    int k =1;
    for (int i=0; i < rows; i++) {
      for (int j =0; j < cols; j++){
	unsigned char temp;
	temp = ((uchar*) matimg.data +
		i * matimg.step)[j  * matimg.elemSize() + k ];
	vec1[i][j] = (double) temp;
      }  
    }
    // get the filter coefficients from the fixed filtcoef method 
    vector<double>l1,h1,l2,h2;
    rc = filtcoef(name, l1, h1, l2, h2);
    //getting the 2D wavelet transform for symmetric
    vector<int> length;
    vector<double> output,flag;
    //int J = 3; // TODO: need to insert the J factor from the argument list 
    dwt_2D_sym(vec1,J,name,output,flag,length);

    double max;
    vector<int> length2;
    //This algorithm computes DWT of image of any given size. Together with
    //convolution and subsampling operations it is clear that subsampled
    //images are of different length than dyadic length images. In order to
    //compute the "effective" size of DWT we do additional calculations.
    
    rc = dwt_output_dim_sym(length,length2,J); 
    //length2 is gives the integer vector that contains the size of subimages
    //that will combine to form the displayed output image. The last two 
    //entries of length2 gives the size of DWT ( rows_n by cols_n) 
    

    int siz = length2.size();
    int rows_n=length2[siz-2];
    int cols_n = length2[siz-1];

    vector<vector< double> > dwtdisp(rows_n, vector<double>(cols_n));
    rc = dispDWT(output,dwtdisp, length ,length2, J);
    //dispDWT returns the 2D object dwtdisp which will be displayed using 
    //OPENCV's image handling functions
    

    vector<vector<double> >  dwt_output= dwtdisp;
    rc = maxval(dwt_output,max);
    //max value is needed to take care of overflow which happens because
    //of convolution operations performed on unsigned 8 bit images
    //isplaying Scaled Image. Creating Image in OPENCV.
    
    IplImage *cvImg; // image used for output 
    CvSize imgSize;  // size of output image  

    imgSize.width = cols_n;
    imgSize.height = rows_n;

    //cvImg = cvCreateImage( imgSize, 8, 1 );
    cvImg = cvCreateImage( imgSize, pix_depth, 1 );//pixel_depth should be 8
    //dwt_hold is created to hold the dwt output as further operations need to
    //be carried out on dwt_output in order to display scaled images.
    
    vector<vector<double> > dwt_hold(rows_n, vector<double>( cols_n));
    dwt_hold = dwt_output;
    //Setting coefficients of created image to the scaled DWT output values 
    for (int i = 0; i < imgSize.height; i++ ) {
      for (int j = 0; j < imgSize.width; j++ ){
	if ( dwt_output[i][j] <= 0.0){
	  dwt_output[i][j] = 0.0;
	}
	if ( i <= (length2[0]) && j <= (length2[1]) ) {
	  ((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
	    (char) ( (dwt_output[i][j] / max) * 255.0);
	} else {
	  ((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
	    (char) (dwt_output[i][j]) ;
	}
      }
    }
    // saving image of the DWT to disk as bmp image 
    cvNamedWindow( "DWT Image", 1 ); // creation of a visualisation window
    cvShowImage( "DWT Image", cvImg ); // image visualisation
    cvWaitKey();
    cvDestroyWindow("DWT Image");
    cvSaveImage("dwt.bmp",cvImg);

    // now getting the inverse DWT for the 2D image 
    vector<vector<double> > idwt_output(rows, vector<double>(cols));
    inv_dwt_2D_sym(output,flag, name, idwt_output,length);

    //Displaying Reconstructed Image
    IplImage *dvImg;
    CvSize dvSize; // size of output image

    dvSize.width = idwt_output[0].size();
    dvSize.height = idwt_output.size();

    cout << idwt_output.size() << idwt_output[0].size() << endl;
    dvImg = cvCreateImage( dvSize, pix_depth, 1 );

    for (int i = 0; i < dvSize.height; i++ )
      for (int j = 0; j < dvSize.width; j++ )
	((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
	  (char) (idwt_output[i][j])  ;

    cvNamedWindow( "Reconstructed Image", 1 );
    cvShowImage( "Reconstructed Image", dvImg );
    cvWaitKey();
    cvDestroyWindow("Reconstructed Image");
    cvSaveImage("recon.bmp",dvImg);
    */
    return rc;
  }  /* end of test_transform_sig_2d_sym */
  
/* /////////////////////////////////////////////////////////////////////////////
   -- function to test 1D symmteric trasnformer for a given wavelet in real spc
*/
  int test_transform_sig_1d_sym(char *wname_in, int J) {
    /* in the argument list add: file name or data coming in, J, */ 
    int rc = RC_SUCCESS;
    string name = string(wname_in);
    printf(ANSI_COLOR_BRIGHT_YELLOW);
    std::cout<<"Wavelet 1D transform for wavelet: "
      ANSI_COLOR_BRIGHT_RED<<name<<std::endl;
    printf(ANSI_COLOR_RESET);

    cout<<"Enter a signal file, fu_real102.asc: "<<endl;
    char inp[50];
    cin >> inp;
    vector<double> sig;
    ifstream sig_inp(inp);
    if (!sig_inp.good() ) {
      cout<<"The file does not exits"<<endl;
    } else {
      cout<<"The file does exits"<<endl;
    }
    while (sig_inp) {
      double temp;
      sig_inp>>temp;
      sig.push_back(temp);
    }
    sig.pop_back();
    //getting the original signal
    vector<double> original;
    original = sig;

    //int J = 6;

    vector<double> dwt_1D_sym_output, flag;
    //getting the decomposition J-level DWT
    vector<int> length;

    rc = dwt_1D_sym(sig, J, name, dwt_1D_sym_output, flag, length );
    ofstream dwtout("fu_dwtout.asc");
    for (unsigned int i = 0; i < dwt_1D_sym_output.size(); i++){
      dwtout << dwt_1D_sym_output[i] << endl;
    }

    //now getting the inverse wavelet transform

    //Perform J-Level IDWT
    vector<double> output;
    inv_dwt_1D_sym(dwt_1D_sym_output, flag,name,output,length);
    
    ofstream sig1("fu_recon.asc");
    ofstream diff("diff.asc");

    cout <<" Reconstructed signal size" << output.size() << endl;
    for (unsigned int i = 0; i < output.size(); i++){
      sig1 << output[i] << endl;
      diff << output[i] - original[i] << endl;
      
    }
    
    return rc;
  }
/* /////////////////////////////////////////////////////////////////////////////
   -- function definitions for the Daubechies wavelets
*/
  int get_wavelet_coef(char *name_in) {
    int rc = RC_SUCCESS;
    string name = string(name_in);
    printf(ANSI_COLOR_BRIGHT_YELLOW);
    std::cout<<"Wavelet cofficients for wavelet: "
      ANSI_COLOR_BRIGHT_RED<<name<<std::endl;
    vector<double> lp1,hp1;
    vector<double> lp2,hp2;
    /*Data structure for the coeficients*/
    wavelets_coef_t s_wcoef;
    wavelets_coef *p_wcoef;

    rc = filtcoef(name,lp1,hp1,lp2,hp2);

    /*Filling in the data structure and instantiating the coefficient class*/
    s_wcoef.name = name;
    s_wcoef.lp1 = lp1; s_wcoef.hp1 = hp1;
    s_wcoef.lp2 = lp2; s_wcoef.hp2 = hp2;
    p_wcoef = new wavelets_coef(name,lp1,hp1,lp2,hp2); 
    /*Printing out the coefficients initialising the iterators for the loop*/
    std:: vector<double>::iterator ilp1;
    std:: vector<double>::iterator ihp1;
    std:: vector<double>::iterator ilp2;
    std:: vector<double>::iterator ihp2;
    printf(ANSI_COLOR_BRIGHT_GREEN);
    std::cout<<std::setw(10)<<"lp1"<<std::setw(12)<<"hp1"
             <<std::setw(12)<<"lp2"<<std::setw(12)<<"hp2"<<'\n';
    printf(ANSI_COLOR_BRIGHT_BLUE);
    for (ilp1=lp1.begin(),ihp1=hp1.begin(),ilp2=lp2.begin(),ihp2=hp2.begin() ;
         ilp1 != lp1.end() ; ++ilp1,++ihp1,++ilp2,++ihp2) {
      std::cout<<std::setw(10)<<*ilp1<<"  "<<std::setw(10)<<*ihp1<<"  "
               <<std::setw(10)<<*ilp2<<"  "<<std::setw(10)<<*ihp2<<'\n';
    }
    printf(ANSI_COLOR_RESET);
    //Destroy the object at later stage p_wcoef->~wavelets_coef();
    return rc;
  } /*end get_wavelet_coef method*/

  /*Aliases*/
  extern "C" void get_wavelet_coef_() __attribute__((weak,alias("get_wavelet_coef")));
  extern "C" void test_transform_sig_1d_sym_() __attribute__((weak,alias("test_transform_sig_1d_sym")));
  extern "C" void test_transform_sig_2d_sym_() __attribute__((weak,alias("test_transform_sig_2d_sym")));
  
#ifdef __cplusplus
}
#endif
