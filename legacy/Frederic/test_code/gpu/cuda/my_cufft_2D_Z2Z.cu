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

  int nx = 12000;
  int ny = 12000;

  int npts = nx * ny;
  /* the error handlers from the cuda library */
  cudaError_t err;
  cufftResult cufft_err;
  /* the plan for the cuFFT */
  cufftHandle plan_fwd;
  cufftHandle plan_bwd;

  printf("nx=%i, ny=%i, npts=%i, nx * ny=%i\n",nx, ny, npts, nx * ny);
  printf("nx*ny*sizeof(cufftDoubleComplex)=%lu\n",nx * ny*sizeof(cufftDoubleComplex));

  cufftDoubleComplex *h_in;
  cufftDoubleComplex *h_out = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * npts);
  cufftDoubleComplex *h_in_rev = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex) * npts);
  
  h_in = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)* npts);
  for (unsigned int i=0 ; i < nx ; i++) {
    for (unsigned int j=0 ; j < ny ; j++) {
      h_in[i+nx*j].x = rand() / (float)RAND_MAX;
      h_in[i+nx*j].y = sin(i*4.0*atan(1.0)*2.0/npts);
    }
  }

  printf("pi: %f\n",4.0*atan(1.0));
  
  cufftDoubleComplex *d_in;
  cufftDoubleComplex *d_in_rev;
  cufftDoubleComplex *d_out;

  cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleComplex));
  cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleComplex));
  cudaMemcpy(d_in, h_in, npts * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);

  //tranform data
  cufftPlan2d(&plan_fwd, nx, ny, CUFFT_Z2Z);
  cufftExecZ2Z(plan_fwd, (cufftDoubleComplex *)d_in, (cufftDoubleComplex *)d_out, CUFFT_FORWARD);
  //copy trans into h_out from device 
  cudaMemcpy(h_out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

  cufft_err = cufftDestroy(plan_fwd);
  err = cudaFree(d_in);

  //transform back
  cudaMalloc((void**)&d_in_rev, npts*sizeof(cufftDoubleComplex));
  cufftPlan2d(&plan_bwd, nx, ny, CUFFT_Z2Z);
  cufftExecZ2Z(plan_bwd, (cufftDoubleComplex *)d_out, (cufftDoubleComplex *)d_in_rev, CUFFT_INVERSE);
  cudaMemcpy(h_in_rev, d_in_rev, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);

  cufft_err = cufftDestroy(plan_bwd);
  err = cudaFree(d_in_rev);
  err = cudaFree(d_out);

  // check result
  for (unsigned int i = 0; i < nx-(nx-3); ++i) {
    for (unsigned int j = 0; j < ny-(ny-3); ++j)
      {
	h_out[i+nx*j].x = h_out[i+nx*j].x / (float)npts;
	h_out[i+nx*j].y /= (float)npts;

	h_in_rev[i+nx*j].x = h_in_rev[i+nx*j].x / (float)npts;
	h_in_rev[i+nx*j].y /= (float)npts;

	printf( ANSI_COLOR_GREEN "data: %15.8f %15.8f"
		ANSI_COLOR_BLUE" Fourier %15.8f %15.8f"
		ANSI_COLOR_RED" Inverse %15.8f %15.8f\n",
		h_in[i+nx*j].x, h_in[i+nx*j].y, 
		h_out[i+nx*j].x, h_out[i+nx*j].y, 
		h_in_rev[i+nx*j].x, h_in_rev[i+nx*j].y);

      }
  }

  free(h_in);
  free(h_out);
  free(h_in_rev);

  cudaDeviceReset();
}
