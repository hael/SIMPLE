/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 17th June 2015
 *
 *      June 2015
 *
 *   code to obtain fftw on gpu
 *
 * @precisions normal z -> s d c
 */

#include "get_fft123D_gpu.h"

//#include <map>

#if defined (CUDA)
cudaDeviceProp deviceProp_fft123D;

/* getting the error message from the cuFFT api using a template typename */
template <typename T>
int get_error_fft123D_gpu(T result, int const line, const char *const file, 
			  char const *const func) {
  int rc = 0;
  if (result) {
    fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n",
	    file, line, static_cast<unsigned int>(result),
	    _cuFFTGetErrorEnum(result), func);
    cudaDeviceReset();
  }
  rc = static_cast<unsigned int>(result);
}
/* template metode to alloc and set different types */
template <typename T>
int t_allocset_fft_gpu(int npts_in, T *h_in, T *h_out,
    		       T *d_in, T *d_out) {
  int rc = 0;
  int npts = npts_in;
  /* the error handlers from the cuda library */
  cudaError_t err;

  err = cudaMalloc((void**)&d_in, npts*sizeof(T));
  err = cudaMalloc((void**)&d_out, npts*sizeof(T));
  err = cudaMemcpy(d_in, h_in, npts * sizeof(T), cudaMemcpyHostToDevice);

  return rc;
}
/* template metode to get and dealloc different types */
template <typename T>
int t_dealcset_fft_gpu(cufftHandle plan_in,
		       int *npts_in, T *h_in, T *h_out,
		       T *d_in, T *d_out,
		       int *sign) {
  int rc = 0;

  int direction = *sign;
  int npts = *npts_in;
  /* the error handlers from the cuda library */
  cudaError_t err;
  cufftResult cufft_err;

  err = cudaMemcpy(h_out, d_out, npts*sizeof(T), cudaMemcpyDeviceToHost);
  cufft_err = cufftDestroy(plan_in);
  err = cudaFree(d_in);
  err = cudaFree(d_out);

  return rc;
}

#endif

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************/
/**
 *  FORTRAN API - math functions (simple interface)
 **/

#if defined (CUDA) /*preprossing for the CUDA environment */

  /* initialising the data structure */
  /* freeing the data structure from memory */
  /* determine if peer to peer is allowed */
  /* getting the combinations for p2p support */
  /* getting the error correcting code memory device true or flase*/
  /* getting the threads details */
  /* getting the number of registers */
  /* getting # of Multiprocessors, CUDA cores/MP and total # of CUDA cores */ 
  // General check for CUDA GPU SM Capabilities
  /* getting the driver version */
  /* getting the runtime version */
  /* getting the name of the device */
  /* get the number of devices avaible on the system */
  /* getting the warning message from the cuda */

  /* 3D fourier transform */
  /* Z2Z */
  int gather_fft_3D_Z2Z_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleComplex *in, cufftDoubleComplex *out,
				int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;
    
    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    if ( direction == CUFFT_FORWARD ) {

      cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_Z2Z);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2Z(plan_fwd, d_in, d_out, direction);
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

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan3d(&plan_bwd, nptx, npty, nptz, CUFFT_Z2Z);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2Z(plan_bwd, d_in, d_out, direction);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_bwd);
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      
    } else { rc = get_warning_message_fft123D_gpu(); }

    return rc;
  }
  /* C2C */
  int gather_fft_3D_C2C_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftComplex *in, cufftComplex *out,
				int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    cufftComplex *d_in;
    cufftComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftComplex), nptx*npty*nptz*sizeof(cufftComplex)/1.e6,
	   nptx*npty*nptz*sizeof(cufftComplex), nptx*npty*nptz*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    if ( direction == CUFFT_FORWARD ) {

      cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_C2C);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecC2C(plan_fwd, d_in, d_out, direction);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_fwd);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan3d(&plan_bwd, nptx, npty, nptz, CUFFT_C2C);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecC2C(plan_bwd, d_in, d_out, direction);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_bwd);
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else { rc = get_warning_message_fft123D_gpu(); }

    return rc;
  }
  /* D2Z */
  int gather_fft_3D_D2Z_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleReal *in, cufftDoubleComplex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_out = nptx * npty * (nptz / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftDoubleReal *d_in;
    cufftDoubleComplex *d_out;
    
    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_out, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex), nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts_out*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_D2Z);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecD2Z(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts_out*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftDestroy(plan_fwd);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* S2C */
  int gather_fft_3D_S2C_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftReal *in, cufftComplex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_out = nptx * npty * (nptz / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftReal *d_in;
    cufftComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_out, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*(nptz/2+1)*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(cufftComplex), nptx*npty*nptz*sizeof(cufftComplex)/1.e6,
	   nptx*npty*(nptz/2+1)*sizeof(cufftComplex), nptx*npty*(nptz/2+1)*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts_out*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftReal), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_R2C);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecR2C(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts_out*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftDestroy(plan_fwd);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* Z2D */
  int gather_fft_3D_Z2D_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftDoubleComplex *in, cufftDoubleReal *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_in = nptx * npty * (nptz / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftDoubleComplex *d_in;
    cufftDoubleReal *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_in, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex), nptx*npty*(nptz/2+1)*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*npty*nptz*sizeof(cufftDoubleComplex), nptx*npty*nptz*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts_in*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts_in*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_Z2D);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecZ2D(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleReal), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* C2S */
  int gather_fft_3D_C2S_gpu_cpp(int *nx, int *ny, int *nz,
				const cufftComplex *in, cufftReal *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_in = nptx * npty * (nptz / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftComplex *d_in;
    cufftReal *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_in, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*(nptz/2+1)*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*(nptz/2+1)*sizeof(cufftComplex), nptx*npty*(nptz/2+1)*sizeof(cufftComplex)/1.e6,
	   nptx*npty*nptz*sizeof(cufftComplex), nptx*npty*nptz*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts_in*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts_in*sizeof(cufftComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan3d(&plan_fwd, nptx, npty, nptz, CUFFT_C2R);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecC2R(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts*sizeof(cufftReal), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }

  /* 2D fourier transform */
  /* Z2Z */
  int gather_fft_2D_Z2Z_gpu_cpp(int *nx, int *ny,
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

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty=%i\n",
	   nptx, npty, npts, nptx*npty, npts, nptx*npty);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(cufftDoubleComplex), nptx*npty*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*npty*sizeof(cufftDoubleComplex), nptx*npty*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    if ( direction == CUFFT_FORWARD ) {
 
      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_Z2Z);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2Z(plan_fwd, d_in, d_out, direction);
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

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan2d(&plan_bwd, nptx, npty, CUFFT_Z2Z);
      cufft_err = cufftExecZ2Z(plan_bwd, d_in, d_out, direction);
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_bwd);
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else { rc = get_warning_message_fft123D_gpu(); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* C2C */
  int gather_fft_2D_C2C_gpu_cpp(int *nx, int *ny,
				cufftComplex *in, cufftComplex *out,
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

    cufftComplex *d_in;
    cufftComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty=%i\n",
	   nptx, npty, npts, nptx*npty, npts, nptx*npty);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(cufftComplex), nptx*npty*sizeof(cufftComplex)/1.e6,
	   nptx*npty*sizeof(cufftComplex), nptx*npty*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    if ( direction == CUFFT_FORWARD ) {

      cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_C2C);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecC2C(plan_fwd, d_in, d_out, direction);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_fwd);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan2d(&plan_bwd, nptx, npty, CUFFT_C2C);
      cufft_err = cufftExecC2C(plan_bwd, d_in, d_out, direction);
      err = cudaMemcpy(out, d_out, npts*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

      cufft_err = cufftDestroy(plan_bwd);
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else { rc = get_warning_message_fft123D_gpu(); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* D2Z */
  int gather_fft_2D_D2Z_gpu_cpp(int *nx, int *ny,
				const cufftDoubleReal *in, cufftDoubleComplex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_out = nptx * (npty / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftDoubleReal *d_in;
    cufftDoubleComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*(npty/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_out, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*(npty/2+1)*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(cufftDoubleComplex), nptx*npty*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*(npty/2+1)*sizeof(cufftDoubleComplex), nptx*(npty/2+1)*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftDoubleReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts_out*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_D2Z);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecD2Z(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts_out*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftDestroy(plan_fwd);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* S2C */
  int gather_fft_2D_S2C_gpu_cpp(int *nx, int *ny,
				const cufftReal *in, cufftComplex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_out = nptx * (npty / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftReal *d_in;
    cufftComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*(npty/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_out, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*(npty/2+1)*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(cufftComplex), nptx*npty*sizeof(cufftComplex)/1.e6,
	   nptx*(npty/2+1)*sizeof(cufftComplex), nptx*(npty/2+1)*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts*sizeof(cufftReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts_out*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts*sizeof(cufftReal), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_R2C);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecR2C(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts_out*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftDestroy(plan_fwd);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* Z2D */
  int gather_fft_2D_Z2D_gpu_cpp(int *nx, int *ny,
				const cufftDoubleComplex *in, cufftDoubleReal *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_in = nptx * (npty / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftDoubleComplex *d_in;
    cufftDoubleReal *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_in, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*(npty/2+1)*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*(npty/2+1)*sizeof(cufftDoubleComplex), nptx*(npty/2+1)*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*npty*sizeof(cufftDoubleComplex), nptx*npty*sizeof(cufftDoubleComplex)/1.e6);

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
    err = cudaMemcpy(out, d_out, npts*sizeof(cufftDoubleReal), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* C2S */
  int gather_fft_2D_C2S_gpu_cpp(int *nx, int *ny,
				const cufftComplex *in, cufftReal *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_in = nptx * (npty / 2 + 1);

    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;

    cufftComplex *d_in;
    cufftReal *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_in, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*(npty/2+1)*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(cufftComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*(npty/2+1)*sizeof(cufftComplex), nptx*(npty/2+1)*sizeof(cufftComplex)/1.e6,
	   nptx*npty*sizeof(cufftComplex), nptx*npty*sizeof(cufftComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, npts_in*sizeof(cufftComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, npts*sizeof(cufftReal));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMemcpy(d_in, in, npts_in*sizeof(cufftComplex), cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cufft_err = cufftPlan2d(&plan_fwd, nptx, npty, CUFFT_C2R);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    cufft_err = cufftExecC2R(plan_fwd, d_in, d_out);
    if ( cufft_err != CUFFT_SUCCESS) {
      rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
    err = cudaMemcpy(out, d_out, npts*sizeof(cufftReal), cudaMemcpyDeviceToHost);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    err = cudaFree(d_in);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaFree(d_out);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }
  /* 1D fourier transform */
  /* Z2Z */
  int gather_fft_1D_Z2Z_gpu_cpp(int *nx,
				double complex *in, double complex *out,
				int *sign) {
    int rc = 0; /* the return code from the function */
    int nptx = *nx;

    int npts = nptx;

    int direction = *sign;
    /* the error handlers from the cuda library */
    cudaError_t err;
    cufftResult cufft_err;
    /* the plan for the cuFFT */
    cufftHandle plan_fwd;
    cufftHandle plan_bwd;

    cufftDoubleComplex *d_in;
    cufftDoubleComplex *d_out;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i"
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i\n",nptx, npts, npts);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*sizeof(cufftDoubleComplex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*sizeof(cufftDoubleComplex), nptx*sizeof(cufftDoubleComplex)/1.e6,
	   nptx*sizeof(cufftDoubleComplex), nptx*sizeof(cufftDoubleComplex)/1.e6);

    err = cudaMalloc((void**)&d_in, nptx*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
    err = cudaMalloc((void**)&d_out, nptx*sizeof(cufftDoubleComplex));
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    err = cudaMemcpy(d_in, in, nptx * sizeof(cufftDoubleComplex),
		     cudaMemcpyHostToDevice);
    if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    if ( direction == CUFFT_FORWARD ) {

      cufft_err = cufftPlan1d(&plan_fwd, nptx, CUFFT_Z2Z, 1);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2Z(plan_fwd, d_in, d_out, CUFFT_FORWARD);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, nptx*sizeof(cufftDoubleComplex),
		       cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      cufft_err = cufftDestroy(plan_fwd);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else if ( direction == CUFFT_INVERSE) {

      cufft_err = cufftPlan1d(&plan_bwd, nptx, CUFFT_Z2Z, 1);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      cufft_err = cufftExecZ2Z(plan_bwd, d_in, d_out, CUFFT_INVERSE);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaMemcpy(out, d_out, nptx*sizeof(cufftDoubleComplex),
		       cudaMemcpyDeviceToHost);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      cufft_err = cufftDestroy(plan_bwd);
      if ( cufft_err != CUFFT_SUCCESS) {
	rc = get_error_fft123D_gpu(cufft_err,__LINE__,__FILE__,__FUNCTION__);}
      err = cudaFree(d_in);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }
      err = cudaFree(d_out);
      if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_fft123D_gpu(err); }

    } else { rc = get_warning_message_fft123D_gpu(); }

    cuCtxSynchronize();
    cudaDeviceReset();

    return rc;
  }

  /* getting the warning message from the cuda */
  int get_warning_message_fft123D_gpu() {
    int rc = 0;

    printf("***************************WARNING*****************************\n");
    printf("You need to compile with -DCUDA to acces the CUDA environment  \n");
    printf("computation using GPU                                          \n");
    printf("***************************************************************\n");
    printf("\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    printf("\n");
    rc = -1;

    return rc;
  }
  /* getting the error id from the cuda */
  int get_error_id_fft123D_gpu(cudaError_t error_id) {
    int rc = 0;
    printf("cudaDriverGetVersion returned %d\n-> %s\n", 
	     (int)error_id, cudaGetErrorString(error_id));
    printf("Result = FAIL\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    rc = (int)error_id;
    exit(EXIT_FAILURE);
    return rc;
  }
  /* getting the cuda cores error */
  int get_cuda_cores_error_fft123D_gpu(int ncores) {
    int rc = 0;
    printf("There are no CUDA cores available on system %d\n",ncores);
    printf("Result = FAIL\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    rc = -1;
    exit(EXIT_FAILURE);
    return rc;
  }/* --end of get_cuda_cores_error--*/
  /* Beginning of GPU Architecture definitions */

  /* the aliases for external access */

  /* 1D */
  extern "C" int gather_fft_1d_z2z_gpu_cpp_() __attribute__((weak,alias("gather_fft_1D_Z2Z_gpu_cpp")));
  /* 2D */
  extern "C" int gather_fft_2d_c2c_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_C2C_gpu_cpp")));
  extern "C" int gather_fft_2d_s2c_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_S2C_gpu_cpp")));
  extern "C" int gather_fft_2d_c2s_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_C2S_gpu_cpp")));
  //double precision
  extern "C" int gather_fft_2d_z2z_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_Z2Z_gpu_cpp")));
  extern "C" int gather_fft_2d_d2z_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_D2Z_gpu_cpp")));
  extern "C" int gather_fft_2d_z2d_gpu_cpp_() __attribute__((weak,alias("gather_fft_2D_Z2D_gpu_cpp")));
  /* 3D */
  extern "C" int gather_fft_3d_c2c_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_C2C_gpu_cpp")));
  extern "C" int gather_fft_3d_s2c_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_S2C_gpu_cpp")));
  extern "C" int gather_fft_3d_c2s_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_C2S_gpu_cpp")));
  //double precision
  extern "C" int gather_fft_3d_z2z_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_Z2Z_gpu_cpp")));
  extern "C" int gather_fft_3d_d2z_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_D2Z_gpu_cpp")));
  extern "C" int gather_fft_3d_z2d_gpu_cpp_() __attribute__((weak,alias("gather_fft_3D_Z2D_gpu_cpp")));

#endif /* CUDA */

#ifdef __cplusplus
}
#endif
