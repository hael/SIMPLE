// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//include the cufft library
#include <cufft.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

//declarations
void runtest(int argc, char**argv);

int main(int argc, char **argv) {
  runtest(argc,argv);
}
void runtest(int argc, char**argv) {

  int SIGNAL_SIZE = 32;
  cufftHandle plan_fwd;
  cufftHandle plan_bwd;
  cufftComplex *h_in;
  cufftComplex *h_out = (cufftComplex*)malloc(sizeof(cufftComplex) * SIGNAL_SIZE);
  cufftComplex *h_in_rev = (cufftComplex*)malloc(sizeof(cufftComplex) * SIGNAL_SIZE);
  
  h_in = (cufftComplex*)malloc(sizeof(cufftComplex)*SIGNAL_SIZE);
  for (unsigned int i=0 ; i < SIGNAL_SIZE ; i++) {
    h_in[i].x = rand() / (float)RAND_MAX;
    h_in[i].y = sin(i*4.0*atan(1.0)*2.0/SIGNAL_SIZE);
  }

  printf("pi: %f\n",4.0*atan(1.0));
  
  cufftComplex *d_in;
  cufftComplex *d_in_rev;
  cufftComplex *d_out;
  cudaMalloc((void**)&d_in, SIGNAL_SIZE*sizeof(cufftComplex));
  cudaMalloc((void**)&d_out, SIGNAL_SIZE*sizeof(cufftComplex));

  cudaMemcpy(d_in, h_in, SIGNAL_SIZE * sizeof(cufftComplex), cudaMemcpyHostToDevice);
  //tranform data
  cufftPlan1d(&plan_fwd, SIGNAL_SIZE, CUFFT_C2C, 1);
  cufftExecC2C(plan_fwd, (cufftComplex *)d_in, (cufftComplex *)d_out, CUFFT_FORWARD);
  //copy trans into h_out from device 
  cudaMemcpy(h_out, d_out, SIGNAL_SIZE*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
  //transform back
  cudaMalloc((void**)&d_in_rev, SIGNAL_SIZE*sizeof(cufftComplex));
  cufftPlan1d(&plan_bwd, SIGNAL_SIZE, CUFFT_C2C, 1);
  cufftExecC2C(plan_fwd, (cufftComplex *)d_out, (cufftComplex *)d_in_rev, CUFFT_INVERSE);
  cudaMemcpy(h_in_rev, d_in_rev, SIGNAL_SIZE*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

  // check result
  for (unsigned int i = 0; i < SIGNAL_SIZE; ++i)
    {
      h_out[i].x = h_out[i].x / (float)SIGNAL_SIZE;
      h_out[i].y /= (float)SIGNAL_SIZE;

      h_in_rev[i].x = h_in_rev[i].x / (float)SIGNAL_SIZE;
      h_in_rev[i].y /= (float)SIGNAL_SIZE;

      printf( ANSI_COLOR_GREEN "data: %15.8f %15.8f"
	      ANSI_COLOR_BLUE" Fourier %15.8f %15.8f"
	      ANSI_COLOR_RED" Inverse %15.8f %15.8f\n",
	     h_in[i].x, h_in[i].y, 
	     h_out[i].x, h_out[i].y, 
	     h_in_rev[i].x, h_in_rev[i].y);
      //printf("1 Error %g %g \n", fabs(h_signal[i].x - h_reversed_signal[i].x), fabs(h_signal[i].y - h_reversed_signal[i].y));
    }

  cudaDeviceReset();
}
