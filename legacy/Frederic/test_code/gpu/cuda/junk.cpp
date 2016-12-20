vector<double> c_lp1 = p_wcoef->get_wavelets_lp1();
vector<double> c_hp1 = p_wcoef->get_wavelets_hp1();
vector<double> c_lp2 = p_wcoef->get_wavelets_lp2();
vector<double> c_hp2 = p_wcoef->get_wavelets_hp2();

std:: vector<double>::iterator jlp1;
std:: vector<double>::iterator jhp1;
std:: vector<double>::iterator jlp2;
std:: vector<double>::iterator jhp2;

std::cout<<std::setw(10)<<"lp1"<<'\n';
for (jlp1=c_lp1.begin(); jlp1 != c_lp1.end() ; ++jlp1) {
  std::cout<<std::setw(10)<<*jlp1<<'\n';
 }
std::cout<<std::setw(10)<<"hp1"<<'\n';
for (jhp1=c_hp1.begin(); jhp1 != c_hp1.end() ; ++jhp1) {
  std::cout<<std::setw(10)<<*jhp1<<'\n';
 }
std::cout<<std::setw(10)<<"lp2"<<'\n';
for (jlp2=c_lp2.begin(); jlp2 != c_lp2.end() ; ++jlp2) {
  std::cout<<std::setw(10)<<*jlp2<<'\n';
 }
std::cout<<std::setw(10)<<"hp2"<<'\n';
for (jhp2=c_hp2.begin(); jhp2 != c_hp2.end() ; ++jhp2) {
  std::cout<<std::setw(10)<<*jhp2<<'\n';
 }
    
for (ilp1=c_lp1.begin(), ihp1=c_hp1.begin(),
       ilp2=c_lp2.begin(), ihp2=c_hp2.begin() ;
     ilp1 != c_lp1.end() ;
     ++ilp1,++ihp1,++ilp2,++ihp2) {
  std::cout<<std::setw(10)<<*ilp1<<"  "<<std::setw(10)<<*ihp1<<"  "
           <<std::setw(10)<<*ilp2<<"  "<<std::setw(10)<<*ihp2<<'\n';
 }

for(vector<double>::iterator i=lp1.begin() ; i != lp1.end(); i++ ) {
  std::cout << *i << '\n';
 }

int num = 0;
int nelm = 1;
int length_num = 1;

printf("*iter: %i, *iint: %i, num: %i, length_num: %i, indexed_int: %i ",
       *iter, *iint, num, length_num, indexed_int);

//recomuting the length of the indexed int.
//length = get_length_of_string(indexed_int);


//printf("length+1: %i ,-----,",length+1);

//printf(ANSI_COLOR_BRIGHT_GREEN"The length of indexed_int: %i\n",indexed_int);


int get_int2char_pos(char *char_out, int *int_in) {
  int rc = RC_SUCCESS;
  int length = get_length_of_string(int_in);
  int size_int_in = (length) * sizeof(char);
  char *buff = (char*)malloc(size_int_in);
  sprintf(buff,"%i",*int_in);
  for (int il=0 ; il < length ; il++) {
    char_out[il] = buff[il];
  }
  free(buff);
  return rc;
}

int convert_int2char_neg(char *char_out, int int_in) {
  int rc = RC_SUCCESS; //return code

  int length = get_lenght_of_string(int_in);
  int size_int_in = (length) * sizeof(char);
  char *buff = (char*)malloc(size_int_in);
  int iil;

  sprintf(buff,"%i",int_in);
  //strcpy(char_out,buff);

  free(buff);

  return rc;
}

  if ( my_int < 0 ) size_int_in = (length) * sizeof(char);
  if ( my_int >= 0 ) size_int_in = (length) * sizeof(char);

if ( my_int+i < 0 ) rc = convert_int2char_neg(char_out,my_int+i);
rc = convert_int2char_pos(char_out,my_int+i);


  if ( int_in >= 0 ) {


  } else if (int_in < 0 ) {
    //printf("input integer is < 0 skip\n");
    //exit;
    ///*
    sprintf(buff,"%i",int_in);
    //printf("buff: %s\n",buff);
    //copying over the buffer into the char_out
    if (int_in < 0 ) iil = 0;
    if (int_in >= 0 ) iil = 1;
    //for (int il=iil ; il < length ; il++) {
      //char_out[il+1] = buff[il];
       strcpy(char_out,buff);
       //}
    //*/
  }




    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6);


    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6);


    err = cudaMalloc((void**)&d_in, npts_in*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts_in*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_Z2D);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2D(plan_fwd, d_in, d_out);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}

      err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_fwd);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}

      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      // rc = get_warning_message_fft123D_gpu(); }


/*-----2D----------------1st of June 2015-------------------------------------*/

      for (int i=0; i<nptx-(nptx-3) ; i++) {
	for (int j=0; j<npty-(npty-3) ; j++) {
	  printf(ANSI_COLOR_YELLOW"direction == CUFFT_FORWARD(-1)= %i,"
		 ANSI_COLOR_GREEN"in[%i][%i]=(%15.8f,%15.8f)"
		 ANSI_COLOR_BLUE" out[%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
		 direction,i,j, in[j*nptx+i].x, in[j*nptx+i].y,
		 i,j,out[j*nptx+i].x,out[j*nptx+i].y );
	}
      }

/*-----3D----------------1st of June 2015-------------------------------------*/
      for (int i=0; i<nptx-(nptx-2) ; i++) {
	for (int j=0; j<npty-(npty-2) ; j++) {
	  for (int k=0; k<nptz-(nptz-2) ; k++) {
	    printf(ANSI_COLOR_GREEN"in[%i][%i][%i]=(%15.8f,%15.8f)"
		   ANSI_COLOR_BLUE" out[%i][%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
		   i,j,k, in[(j+npty*k)*nptx+i].x, in[(j+npty*k)*nptx+i].y,
		   i,j,k,out[(j+npty*k)*nptx+i].x,out[(j+npty*k)*nptx+i].y );
	  }
	}
      }

/*-----------------------30th of June 2015------------------------------------*/


/* c */

/* 2D cufft Z2Z */
int gather_fft_2d_gpu_c_(int *nx, int *ny,
			 double complex *h_in, double complex *h_out,
			 int *sign) {
  int rc = 0;

  cufftDoubleComplex *in = (cufftDoubleComplex*) h_in; 
  cufftDoubleComplex *out = (cufftDoubleComplex*) h_out;

  rc = gather_fft_2D_gpu_cpp(nx, ny,in, out,sign);

  return rc;
}

/* h */

int gather_fft_2D_gpu_cpp(int *nx, int *ny,
			  cufftDoubleComplex *in, cufftDoubleComplex *out,
			  int *sign);

  int gather_fft_2D_gpu_cpp(int *nx, int *ny,
			    cufftDoubleComplex *d_in, cufftDoubleComplex *out,
			    int *sign);

/*cpp*/

  /* 2D fourier transform */
  int gather_fft_2D_gpu_cpp(int *nx, int *ny,
			    cufftDoubleComplex *in, cufftDoubleComplex *out,
			    int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;
    
    printf("nptx=%i, npty=%i, npts=%i, nptx * npty=%i\n",nptx, npty, npts, nptx*npty);
    printf("nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu\n",nptx*npty*sizeof(cufftDoubleComplex));

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleComplex));
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleComplex));
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);

    if ( direction == CUFFT_FORWARD ) {
 
      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_fwd, d_in, d_out, direction);
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

      for (int i=0; i<nptx-(nptx-3) ; i++) {
	for (int j=0; j<npty-(npty-3) ; j++) {
	  printf(ANSI_COLOR_GREEN"in[%i][%i]=(%15.8f,%15.8f)"
		 ANSI_COLOR_BLUE" out[%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
		 i,j, in[j*nptx+i].x, in[j*nptx+i].y,
		 i,j,out[j*nptx+i].x,out[j*nptx+i].y );
	}
      }

      cufft_err = cufftDestroy(plan_fwd);
      err = cudaFree(d_in);
      err = cudaFree(d_out);

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan2d(&plan_bwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_bwd, d_in, d_out, direction);
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

      cufft_err = cufftDestroy(plan_bwd);

      err = cudaFree(d_in);
      err = cudaFree(d_out);

    } else { rc = get_warning_message_fft123D_gpu(); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }

  /* 2D fourier transform */
  int gather_fft_2D_gpu_cpp(int *nx, int *ny,
			    cufftDoubleComplex *d_in, cufftDoubleComplex *out,
			    int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    cufftDoubleComplex *dIn;
    cufftDoubleComplex *d_out;
    
    printf("nptx=%i, npty=%i, npts=%i, nptx * npty=%i\n",nptx, npty, npts, nptx*npty);
    printf("nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu\n",nptx*npty*sizeof(cufftDoubleComplex));

    dIn = d_in; /*mapping the pointer to avoid data corruption*/

    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleComplex));

    if ( direction == CUFFT_FORWARD ) {
 
      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_fwd, dIn, d_out, CUFFT_FORWARD);

      cublasGetMatrix (nptx, npty, sizeof(cufftDoubleComplex), d_out, nptx, out, npty);

      cufft_err = cufftDestroy(plan_fwd);

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan2d(&plan_bwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_bwd, dIn, d_out, CUFFT_INVERSE);

      cublasGetMatrix (nptx, npty, sizeof(cufftDoubleComplex), d_out, nptx, out, npty);

      cufft_err = cufftDestroy(plan_bwd);

    } else { rc = get_warning_message_fft123D_gpu(); }

    
    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }




/*----------------------------------------------------------------------------*/


      for (int i=0; i<nptx-(nptx-2) ; i++) {
	for (int j=0; j<npty-(npty-2) ; j++) {
	  for (int k=0; k<nptz-(nptz-2) ; k++) {
	    printf("in[%i][%i][%i]=(%15.8f,%15.8f), out[%i][%i][%i]=(%15.8f,%15.8f)\n",
		   i,j,k, in[(j+npty*k)*nptx+i].x, in[(j+npty*k)*nptx+i].y,
		   i,j,k,out[(j+npty*k)*nptx+i].x,out[(j+npty*k)*nptx+i].y );
	  }
	}
      }


      for (int i=0; i<nptx-(nptx-2) ; i++) {
	for (int j=0; j<npty-(npty-2) ; j++) {
	  for (int k=0; k<npty-(nptz-2) ; k++) {
	    printf("in[%i][%i][%i]=(%15.8f,%15.8f), out[%i][%i][%i]=(%15.8f,%15.8f)\n",
		   i,j,k,creal(in[(j+npty*k)*nptx+i]),cimag(in[(j+npty*k)*nptx+i]),
		   i,j,k,creal(out[(j+npty*k)*nptx+i]),cimag(out[(j+npty*k)*nptx+i]) );
	  }
	}
      }

      for (int i=0; i<nptx-(nptx-3) ; i++) {
	for (int j=0; j<npty-(npty-3) ; j++) {
	  printf("in[%i][%i]=(%15.8f,%15.8f), out[%i][%i]=(%15.8f,%15.8f)\n",
		 i,j,creal(in[j*nptx+i]),cimag(in[j*nptx+i]),
		 i,j,creal(out[j*nptx+i]),cimag(out[j*nptx+i]) );
	}
      }


    for (int i=0; i<nptx-(nptx-3) ; i++) {
      for (int j=0; j<npty-(npty-3) ; j++) {
	printf("in[%i][%i]=(%15.8f,%15.8f)",i,j,creal(in[i*nptx+j]),cimag(in[i*nptx+j]));
      }
    }

//------------------------Faulty method-----------------------------------------

  int gather_fft_2D_gpu_cpp(int *nx, int *ny,
			    cufftDoubleComplex *in, cufftDoubleComplex *out,
			    int *sign) {
    int rc = 0; /* the return code from the function */
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    //int i,j;

    printf("nptx=%i, npty=%i, npts=%i, nptx * npty=%i\n",nptx, npty, npts, nptx * npty);
    printf("nptx * npty*sizeof(cufftDoubleComplex)=%lu\n",nptx * npty*sizeof(cufftDoubleComplex));

    cufftDoubleComplex *h_in  = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*nptx*npty);
    cufftDoubleComplex *h_out = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*nptx*npty);
    
    for (int i=0 ; i < nptx ; i++) {
      for (int j=0 ; j < npty ; j++) {
	h_in[i + nptx * j].x = in[i + nptx * j].x;
	h_in[i + nptx * j].y = in[i + nptx * j].y;
      }
    }

    for (int i=0; i<nptx-(nptx-3) ; i++) {
      for (int j=0; j<npty-(npty-3) ; j++) {
	printf(ANSI_COLOR_GREEN"in[%i][%i]=(%15.8f,%15.8f)" 
	       ANSI_COLOR_GREEN" h_in[%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
	       i,j,in[j*nptx+i].x,in[j*nptx+i].y,
	       i,j,h_in[j*nptx+i].x,h_in[j*nptx+i].y );
      }
    }

    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;

    err = cudaMalloc((void**)&d_in, nptx * npty*sizeof(cufftDoubleComplex));
    err = cudaMalloc((void**)&d_out, nptx * npty*sizeof(cufftDoubleComplex));
    err = cudaMemcpy(d_in, h_in, nptx * npty * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);

    if ( direction == CUFFT_FORWARD ) {

      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_fwd, d_in, d_out, direction);
      err = cudaMemcpy(h_out, d_out, nptx * npty*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

      for (int i = 0; i < nptx; i++) {
	for (int j = 0; j < npty; j++)
	  {
	    out[i+nptx*j].x = h_out[i+nptx*j].x;
	    out[i+nptx*j].y = h_out[i+nptx*j].y;
	  }
      }

      cufft_err = cufftDestroy(plan_fwd);
      err = cudaFree(d_in);
      err = cudaFree(d_out);

      for (int i=0; i<nptx-(nptx-3) ; i++) {
	for (int j=0; j<npty-(npty-3) ; j++) {
	  printf(ANSI_COLOR_GREEN"in[%i][%i]=(%15.8f,%15.8f)" 
		 ANSI_COLOR_BLUE" out[%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
		 i,j,in[j*nptx+i].x,in[j*nptx+i].y,
		 i,j,out[j*nptx+i].x,out[j*nptx+i].y );
	}
      }

    } else if ( direction == CUFFT_INVERSE) {
      
      cufft_err = cufftPlan2d(&plan_bwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_bwd, d_in, d_out, direction);
      err = cudaMemcpy(h_out, d_out, nptx * npty*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

      cufft_err = cufftDestroy(plan_bwd);
      err = cudaFree(d_in);
      err = cudaFree(d_out);

      for (int i = 0; i < nptx; ++i) {
	for (int j = 0; j < npty; ++j)
	  {
	    out[i+nptx*j].x = h_out[i+nptx*j].x;
	    out[i+nptx*j].y = h_out[i+nptx*j].y;
	  }
      }

      for (int i=0; i<nptx-(nptx-3) ; i++) {
	for (int j=0; j<npty-(npty-3) ; j++) {
	  printf(ANSI_COLOR_GREEN"in[%i][%i]=(%15.8f,%15.8f)" 
		 ANSI_COLOR_BLUE" out[%i][%i]=(%15.8f,%15.8f)\n" ANSI_COLOR_RESET,
		 i,j,in[j*nptx+i].x,in[j*nptx+i].y,
		 i,j,out[j*nptx+i].x,out[j*nptx+i].y );
	}
      }

    } else { rc = get_warning_message_fft123D_gpu(); }


    free(h_in);
    free(h_out);

    cudaDeviceReset();

    return rc;
  }
//------------------------------------------------------------------------------


  swicth(sizeof(T)) {
  case sizeof(cufftDoubleComplex):
    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;
  default:
  }

  if (sizeof(T) == sizeof(cufftDoubleComplex) ){
    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;
  }





    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;

    int npts = *nstep;
    int direction = *sign;
    /* fourier plan */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    //Temporary test data
    int SIGNAL_SIZE = 32;
    cufftComplex *h_in;
    h_in = (cufftComplex*)malloc(sizeof(cufftComplex)*SIGNAL_SIZE);
    for (unsigned int i=0 ; i < SIGNAL_SIZE ; i++) {
      h_in[i].x = rand() / (float)RAND_MAX;
      h_in[i].y = sin(i*4.0*atan(1.0)*2.0/SIGNAL_SIZE);
    }

    cufftComplex *d_in;
    cudaMalloc((void**)&d_in, SIGNAL_SIZE*sizeof(cufftComplex));
    cudaMemcpy(d_in, h_in, SIGNAL_SIZE * sizeof(cufftComplex), cudaMemcpyHostToDevice);
    cufftPlan1d(&plan_fwd, SIGNAL_SIZE, CUFFT_C2C, 1);

    cudaDeviceReset();

    /* Calculating the fourier transform */ 
    if ( direction == FFTW_FORWARD ) {
      //cufft_err = cufftPlan1d(&plan_fwd, npts, CUFFT_Z2Z, 1);
      //if ( cufft_err != CUFFT_SUCCESS) {rc = get_error_fft123D_gpu(cufft_err, __FUNCTION__, __FILE__, __LINE__);}
      //cufft_err = cufftExecZ2Z(plan_fwd, (cufftDoubleComplex*)d_in,(cufftDoubleComplex*)d_out,CUFFT_FORWARD);
      //if ( cufft_err != CUFFT_SUCCESS) {rc = get_error_fft123D_gpu(cufft_err, __FUNCTION__, __FILE__, __LINE__);}
      //cufftDestroy(plan_fwd);
      //if ( cufft_err != CUFFT_SUCCESS) {rc = get_error_fft123D_gpu(cufft_err, __FUNCTION__, __FILE__, __LINE__);}
      //plan_fwd = fftw_plan_dft_1d(nstepi, fz, fh, direction, FFTW_ESTIMATE);
      //fftw_execute(plan_fwd);
      //fftw_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      //plan_bwd = fftw_plan_dft_1d(nstepi, fz, fh, FFTW_BACKWARD, FFTW_ESTIMATE);
      //fftw_execute(plan_bwd);
      //fftw_destroy_plan ( plan_bwd );
    }






    for (unsigned int i = 0; i < npts-(npts-3); ++i) {
      printf("%i %15.8f %15.8f \n",i,h_in[i].x,h_in[i].y);
    }
    for (unsigned int i = (npts-3); i < npts; ++i) {
      printf("%i %15.8f %15.8f \n",i,h_in[i].x,h_in[i].y);
    }

    /* device variables */
    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;

    /* first allocating the memory on device for both in and out */ 
    err = cudaMalloc((void**)&d_in, npts * sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts * sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    /* Transferring the data to gpu */
    err = cudaMemcpy(d_in, in, npts * sizeof(cufftDoubleComplex),
		     cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_out, out, npts * sizeof(cufftDoubleComplex), 
		     cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }




    /* transfering back to the host from device */
    err = cudaMemcpy(out, d_out, npts * sizeof(cuDoubleComplex), 
		     cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    /* freeing the cuda allocated resources */
    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }




printf("sizeof complex=%lu\n",sizeof(complex));
printf("sizeof double complex=%lu\n",sizeof(double complex));
printf("sizeof cuDoubleComplex=%lu\n",sizeof(cuDoubleComplex));
printf("sizeof cufftDoubleComplex=%lu\n",sizeof(cufftDoubleComplex));

for (unsigned int i = (npts-3); i < npts; ++i)
  {
    printf("%i %15.8f %15.8f \n",i,creal(in[i]),cimag(in[i]));
  }


fftw_complex *in_c;
fftw_complex *out_c;


in_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);
out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);

fprintf(testfile,"%15.5f  %15.5f %15.5f\n", x[istep], creal(in_c[istep]), creal(out_c[istep]));
fprintf(testfile,"%15.5f\n", creal(out_c[istep]));

fftw_free(in_c);
fftw_free(out_c);

FFTW_FORWARD

in_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);
out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);

for (istep = 0 ; istep < nstepi ; istep++ ) { 
  in_c[istep] = fz[istep]; 
 }

printf("direction= %d, FFTW_FORWARD= %d, FFTW_BACKWARD= %d\n",direction, FFTW_FORWARD, FFTW_BACKWARD);

printf("direction= %d, FFTW_FORWARD= %d, FFTW_BACKWARD= %d\n",direction, FFTW_FORWARD, FFTW_BACKWARD);

fftw_complex *in_back_c;

in_back_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);


plan_bwd = fftw_plan_dft_1d(nstepi, fh, in_back_c, FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_execute(plan_bwd);

FILE * testfile;
testfile = fopen("Testfile.log","a");

for (istep = 0 ; istep < nstepi ; istep++ ) {
  //fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
  //fprintf(testfile,"%15.5f  %15.5f %15.5f\n", x[istep], creal(in_c[istep]), creal(out_c[istep]));
  //fprintf(testfile,"%15.5f %15.5f\n", x[istep], creal(in_c[istep]));
  //fprintf(testfile,"%15.5f\n", creal(out_c[istep]));
  //fprintf(testfile,"%d %15.5f\n", direction, creal(fh[istep]));
  //fprintf(testfile,"%15.5f\n", creal(in_back_c[istep])/nstepi);
 }

fclose(testfile);

fftw_free(in_back_c);




/* the method */

  int gather_fft_1D_cpu_cpp(double *xi, double *xf, int *nstep,
			    double *x, double *fx, double complex *fh) {
    int rc = 0;

    int istep;
    int nstepi = *nstep;

    fftw_complex *in_c;
    fftw_complex *out_c;

    fftw_complex *in_back_c;

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    in_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);
    //out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);
    in_back_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);

    for (istep = 0 ; istep < nstepi ; istep++ ) { 
      in_c[istep] = fx[istep]; 
    }

    plan_fwd = fftw_plan_dft_1d(nstepi, in_c, fh, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);

    plan_bwd = fftw_plan_dft_1d(nstepi, fh, in_back_c, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_bwd);

    FILE * testfile;
    testfile = fopen("Testfile.log","a");

    for (istep = 0 ; istep < nstepi ; istep++ ) {
      //fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
      //fprintf(testfile,"%15.5f  %15.5f %15.5f\n", x[istep], creal(in_c[istep]), creal(out_c[istep]));
      //fprintf(testfile,"%15.5f %15.5f\n", x[istep], creal(in_c[istep]));
      //fprintf(testfile,"%15.5f\n", creal(out_c[istep]));
      fprintf(testfile,"%15.5f %15.5f\n", x[istep], creal(in_back_c[istep])/nstepi);
    }

    fclose(testfile);

    fftw_destroy_plan ( plan_fwd );
    fftw_destroy_plan ( plan_bwd );
    fftw_free(in_c);
    fftw_free(out_c);
    fftw_free(in_back_c);

    return rc;
  }

  int test_fft_1D_cpu_cpp(double *xi, double *xf, int *nstep,
			  double *x, double *fx, double complex *fh) {
    int rc = 0;

    int nc;
    int istep;
    int nstepi = *nstep;

    fftw_complex in_c[nstepi];
    fftw_complex out_c[nstepi];
    fftw_complex in_back_c[nstepi];

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    for (istep = 0 ; istep < nstepi ; istep++ ) { 
      in_c[istep] = fx[istep]; 
    }

    plan_fwd = fftw_plan_dft_1d(nstepi, in_c, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);

    plan_bwd = fftw_plan_dft_1d(nstepi, out_c, in_back_c, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_bwd);

    FILE * testfile;
    testfile = fopen("Testfile.log","a");

    for (istep = 0 ; istep < nstepi ; istep++ ) {
      //fprintf(testfile,"  %15.5f  %15.5f %15.5f\n", x[istep], creal(in_c[istep]), creal(out_c[istep]));
      //fprintf(testfile,"  %15.5f %15.5f\n", x[istep], creal(in_c[istep]));
      //fprintf(testfile,"  %15.5f\n", creal(out_c[istep]));
      fprintf(testfile,"  %15.5f %15.5f\n", x[istep], creal(in_back_c[istep])/nstepi);
      //      fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
    }

    fclose(testfile);

    fftw_destroy_plan ( plan_fwd );
    fftw_destroy_plan ( plan_bwd );

    return rc;
  }


  int gather_fft_1D_cpu_cpp(double *xi, double *xf, int *nstep,
			    double *x, double *fx, double complex *fh) {
    int rc = 0;

    /*
    double xii = *xi;
    double xfi = *xf;
    int nstepi = *nstep;

    double *in;
    fftw_complex *out;

    int nc = (nstepi / 2) + 1; 
    int istep;

    
    in = (double*)fftw_malloc(sizeof(double)*nstepi);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nc);

    for (istep = 0 ; istep < nstepi ; istep++ ) { in[istep] = fx[istep]; }

    plan_fwd = fftw_plan_dft_r2c_1d( nstepi, in, out, FFTW_ESTIMATE );

    fftw_execute ( plan_fwd );
    */
    int nc;
    int istep;
    int nstepi = *nstep;

    fftw_complex in_c[nstepi];
    fftw_complex out_c[nstepi];

    fftw_complex in_back_c[nstepi];

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    for (istep = 0 ; istep < nstepi ; istep++ ) { 
      in_c[istep] = fx[istep]; 
      //      in_c[istep][1] = 0.0;
    }

    plan_fwd = fftw_plan_dft_1d(nstepi, in_c, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);

    plan_bwd = fftw_plan_dft_1d(nstepi, out_c, in_back_c, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_bwd);

    FILE * testfile;
    testfile = fopen("Testfile.log","a");

    for (istep = 0 ; istep < nstepi ; istep++ ) {
      //fprintf(testfile,"  %12.5f  %12.5f %12.5f\n", x[istep], creal(in_c[istep]), creal(out_c[istep]));
      //fprintf(testfile,"  %12.5f %12.5f\n", x[istep], creal(in_c[istep]));
      //fprintf(testfile,"  %12.5f\n", creal(out_c[istep]));
      fprintf(testfile,"  %12.5f %12.5f\n", x[istep], creal(in_back_c[istep])/nstepi);
      //      fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
    }

    fclose(testfile);

    fftw_destroy_plan ( plan_fwd );
    fftw_destroy_plan ( plan_bwd );
    //fftw_free(in_c);
    //fftw_free(out_c);

    return rc;
  }









    FILE * testfile;
    int istep;
    printf("xii=%f xfi=%f nstepi=%i,Exit at Line %i in file %s %s\n",xii,xfi,nstepi,__LINE__,__FILE__,__FUNCTION__);
   
    testfile = fopen("Testfile.log","a");

    for (istep = 0 ; istep < nstepi ; istep++ ) {
      fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
    }

    fclose(testfile); 


//**************************from simple_polarft_gpu_c.cpp *******************/

  void MAGMABLAS_DZGEMM_TESLA_GPU( char *TRANSA, char *TRANSB,
				   int *m, int *n, int *k,
				   cuDoubleComplex *alpha,
				   const devptr_t *A, int *lda,
				   const devptr_t *B, int *ldb,
				   cuDoubleComplex *beta,
				   devptr_t *C, int *ldc)
  {
    double *d_a = (double *)(*A);
    cuDoubleComplex *d_b = (cuDoubleComplex *)(*B);
    cuDoubleComplex *d_c = (cuDoubleComplex *)(*C);
    magmablas_dzgemm_tesla_gpu_(*TRANSA, *TRANSB, *m, *n, *k, *alpha, d_a, *lda, d_b, *ldb, *beta, d_c, *ldc);
  }
  extern "C" void dzgemm_tesla_gpu_( char, char, int, int, int, cuDoubleComplex , const devptr_t, int, const devptr_t, int, cuDoubleComplex, devptr_t, int) __attribute__((weak,alias("MAGMABLAS_DZGEMM_TESLA_GPU")));


  void ZZ2DGEMM_ELMTWS_TESLA_SUMSQ_GPU( char *TRANSA, char *TRANSB,
					int *m, int *n, int *k,
					double *alpha,
					const devptr_t *A, int *lda,
					const devptr_t *B, int *ldb,
					double *beta,
					devptr_t *C, int *ldc,
					devptr_t *sumasq, devptr_t *sumbsq)
  {
    cuDoubleComplex *d_a = (cuDoubleComplex *)(*A);
    cuDoubleComplex *d_b = (cuDoubleComplex *)(*B);
    double *d_c          = (double *)(*C);
    double *d_sumasq     = (double *)(*sumasq);
    double *d_sumbsq     = (double *)(*sumbsq);
    zz2dgemm_ElmtWs_tesla_gpu_(*TRANSA, *TRANSB, *m, *n, *k, 
			       *alpha, 
			       d_a, *lda, 
			       d_b, *ldb, *beta, 
			       d_c, *ldc, 
			       d_sumasq, d_sumbsq);

//***************************************************************************/


//From the magma_zgetri_blocked_gpu-v2.cpp

typedef struct {
  double re;
  double im;
}ComplexDouble;

/* ////////////////////////////////////////////////////////////////////////////
   -- Defines a set of function for the error call.
*/
#define CUDA_CALL(x) do { if((x) != cudaSuccess){ \
      printf("Error %i:, failed with message %s, at Line %i in file %s %s\n",x,cudaGetErrorString(x),__LINE__,__FILE__,__FUNCTION__); \
      break;/*return EXIT_FAILURE;*/}} while(0)


#define  inAT(i,j)  ( dAT +  (i) * n  * ldda  + (j) * n  )    //blocks of size n first 


// Flops formula
#define PRECISION_z
#if defined(PRECISION_z) || defined(PRECISION_c)
#define FLOPS(m, n) ( 6. * FMULS_GETRI(m, n) + 2. * FADDS_GETRI(m, n) )
#else
#define FLOPS(m, n) (      FMULS_GETRI(m, n) +      FADDS_GETRI(m, n) )
#endif

//from magma_z_GPU-v2.cpp

//#include <ctype.h>
//#include <stdio.h>

#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

#include "timming_c.h"

//global objects declared for project

//magma_int_t kb_val;
//kbZlasyfGPU *the_kb = new kbZlasyfGPU(kb_val);

//int gettimeofday(struct timeval *restrict tp, void *restrict tzp);

void GETTIMEOFDAY_C( double Time[2] )
{
  struct timeval tp;
  struct timezone tz;

  gettimeofday(&tp, &tz);
  
  Time[0] = tp.tv_sec;
  Time[1] = tp.tv_usec;
}

void ELAPSED_TIME_C( double start_Time[2], double end_Time[2], double *elapsed_time)
{
  *elapsed_time = (end_Time[0] - start_Time[0]) + (end_Time[1] - start_Time[1])/1000000.0;
}

extern "C" void gettimeofday_c_(double) __attribute__((alias("GETTIMEOFDAY_C")));
extern "C" void elapsed_time_c_(double, double, double ) __attribute__((alias("ELAPSED_TIME_C")));



//from magma_z_GPU-v2.h file

//#include "f2c.h"
//#include "blaswrap.h"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

#include <objectDec_GPU_funcdec.h>

#define lapackf77_zgetrf  zgetrf_
#define lapackf77_zgetri  zgetri_
#define lapackf77_ztrti2  FORTRAN_NAME( ztrti2, ZTRTI2 )

//global objects declared for the project.

extern kbZlasyfGPU *the_kb;

void lapackf77_zgetri( magma_int_t *n, cuDoubleComplex *a, 
		       magma_int_t *lda, magma_int_t *ipiv, 
		       cuDoubleComplex *work, magma_int_t *lwork,
		       magma_int_t *info);


