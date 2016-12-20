/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 25th of September 2015
 *      Monash University
 *      Spetember 2015
 *
 *      Routine which calculates the product (element wise) of 
 *      cuDoubleComplex matrix A and cuDoubleComplex matrix B. and takes the
 *      real part of the product and puts it into a matrix of type double
 *      C = cuCreal(x) * cuCreal(y) + cuCimag(x) * cuCimag(y)
 *      C = Re( A * conjg(B) ) element wise
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/
#include "common_magma.h"
#include "commonblas_zz2d.h"
#include "simple_cuDoubleComplex.h"
#include "polarft_gpu.h"
#include "simple.h"

#define imin(a,b) (a<b?a:b)
//#define debug false
//#define debug_high false
//#define debug_write false
#if defined (CUDA) /*preprossing for the OPENCL environment */

/* doing the r product Re(A*conjg(B)) */
extern "C" __global__ void
pprod_3D_mat( float *C,
              const cuFloatComplex *A,
              const cuFloatComplex *B,
              int npart, int nrot, int nk,
              float alpha)
{

  int ipart = blockIdx.x * blockDim.x + threadIdx.x;
  int  irot = blockIdx.y * blockDim.y + threadIdx.y;
  int    ik = blockIdx.z * blockDim.z + threadIdx.z;

  if ( ipart < npart) {
    if (irot < nrot ){
      if (ik < nk ){
        C[(irot+nrot*ik)*npart+ipart ] =
          cuReCCstarmulf( A[(irot+nrot*ik)*npart+ipart], 
                          B[(irot+nrot*ik)*npart+ipart] );
      }
    }
  }
  __syncthreads();
}

/* main kernel entry */
extern "C" int
polarft_corr_P_N(deviceDetails_t * s_devD,
                 polar_corr_calc_t *s_polar,
                 float *r,
                 const cuFloatComplex *A,
                 const cuFloatComplex *B,
                 int npart, int nrot, int nk,
                 float alpha,
                 bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = 0;
  /*numbers of points to be considered */
  int npts = npart * nrot * nk;
  /* grid and threads block definition */
  int nx, ny, nz;          //Dimension of the threads
  int gridx, gridy, gridz; //Dimension of the 3D Grid
  /* setting the values into the object */
  pfts_Sizes *p_pfts_Sizes = new pfts_Sizes(npart,nrot,nk);
  mesh_3D *p_mesh_3D = new mesh_3D(s_polar);
  /* the error handlers from the cuda library */
  cudaError_t err;
  /*device allocation */
  cuFloatComplex *d_A = NULL;
  cuFloatComplex *d_B = NULL;
  float *d_C = NULL;   //device pointer for the r product
  /* size of the element in consideration */
  int size_m;
  int size_reC;
  /* indexer */
  int i,j;
  /*start of the execution commands */
  cuFloatComplex *C = (cuFloatComplex*)malloc(npts);

  float *reC = (float*)malloc(npts);

  /* nx=16; ny=16; nz=4; dim3 threads(16,16,4); */
  nx=s_polar->nx; ny=s_polar->ny; nz=s_polar->nz;
  dim3 threads(nx,ny,nz);
  gridx = npart/(float)nx+(npart%nx!=0);
  gridy =  nrot/(float)ny+( nrot%ny!=0);
  gridz =    nk/(float)nz+(   nk%nz!=0);
  dim3 grid(gridx,gridy,gridz);
  /*dim3 grid(npart/16.0+(npart%16!=0),
    nrot/16.0+( nrot%16!=0),
    nk/4.0+(     nk%4!=0) );
  */
  
  if (s_debug_gpu->debug_i == true ) {
    rc = print_3D_mesh(0,p_mesh_3D,s_polar,p_pfts_Sizes,nx,ny,nz);

    /* printing input variables on screen */
    if ( s_debug_gpu->debug_high_i == true ) {
      if ( s_debug_gpu->debug_write_i == true ) {
        rc = print_function_header_P_N_info(s_polar,
                                            r,
                                            C, A, B, 
                                            npart, nrot, nk,
                                            alpha);
      }
    }
    
  }
  //allocating the memory on GPU
  size_m = npts*sizeof(cuFloatComplex);
  err = cudaMalloc((void**)&d_A, size_m);
  err = cudaMalloc((void**)&d_B, size_m);
  size_reC = npts*sizeof(float);
  err = cudaMalloc((void**)&d_C, size_reC);
  //uploading the matrices A, B to GPU
  err = cudaMemcpy(d_A, A, size_m, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_B, B, size_m, cudaMemcpyHostToDevice);
  //computing the r product
  pprod_3D_mat<<<grid,threads>>>( d_C, d_A, d_B, npart, nrot, nk, alpha);
  //synchronizing CPU with GPU
  cuCtxSynchronize();

  float *h_C = (float*)malloc(size_reC);
  err = cudaMemcpy(h_C, d_C, size_reC, cudaMemcpyDeviceToHost);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }

  //To be OpenMP if cannot get the double sum to work on GPU
  //#pragma omp parallel default(shared) private (i,j)
  //#pragma omp for schedule(auto)
  for (i=0; i<npart ; i++) {
    r[i] = 0.0;
    for (j=0; j<nrot*nk ; j++) {
      r[i] += h_C[j*npart+i];
    }
  }

  if (s_debug_gpu->debug_i == true ) {
    for (i=0; i<npart/*-(npart-5)*/ ; i++) {
      //      printf(ANSI_COLOR_BRIGHT_GREEN "vector float r[%i]=%f \n" ANSI_COLOR_RESET,
      //     i,r[i]);
    }
  }

  //freeing the temporary resources on CPU
  free(h_C);
  //free(sumb_vec);

  //freeing resources on device
  err = cudaFree(d_A);
  err = cudaFree(d_B);
  err = cudaFree(d_C);

  return rc;
 
} /* End of polarft_corr_P_N */

#endif /* CUDA */
