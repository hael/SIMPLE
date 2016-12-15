/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 27th of January 2016
 *      Monash University
 *      January 2016
 *
 *      Routine which calculates the product (element wise) of 
 *      cuDoubleComplex matrix A and cuDoubleComplex matrix B. and takes the
 *      real part of the product and puts it into a matrix of type double
 *      C = cuCreal(x) * cuCreal(y) + cuCimag(x) * cuCimag(y)
 *      C = Re( A * conjg(B) ) element wise
 *      The product is then summed over the 2nd and 3rd dimensions
 *      to produce a vector of correlators.
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
#if defined (CUDA) /*preprossing for the CUDA environment */

/* doing the r product Re(A*conjg(B)) */
extern "C" __global__ void
xprod_3D_mat( float *C,
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
/* summing the 3D matrix treating it a 1D vec */
extern "C" __global__ void
mysum_23D(float *sumA, float *A, int chunk, int npart){
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int sharedIdx = threadIdx.x;
  __shared__ float shared_A[256];
  int offset = -1;
  
  float temp = 0.0;
  while ( i < chunk) {
    offset++;
    temp += A[i*npart + offset];//A[i];
    i += blockDim.x * gridDim.x;
  }
  shared_A[sharedIdx] = temp;
  __syncthreads();
  int j = blockDim.x / 2;
  while ( j != 0 ) {
    if ( sharedIdx < j ) shared_A[sharedIdx] += shared_A[sharedIdx + j];
    __syncthreads();
    j /= 2;
  }
  if ( sharedIdx == 0 ) sumA[blockIdx.x] = shared_A[0];
  //atomicAdd(temp,2.0);//atomicAdd(&shared_A[0],2.0);

}

/* main kernel entry */
extern "C" int
polarft_corr_X_N(deviceDetails_t * s_devD,
                 polar_corr_calc_t *s_polar,
                 float *r,
                 const cuFloatComplex *A,
                 const cuFloatComplex *B,
                 int npart, int nrot, int nk,
                 float alpha,
                 bench_t *s_bench, debug_gpu_t *s_debug_gpu)
{
  int rc = 0;
  /* numbers of points to be considered */
  int npts = npart * nrot * nk; //Total # of pts
  /* grid and threads block definition */
  int nx, ny, nz;               //Dimension of the threads
  int gridx, gridy, gridz;      //Dimension of the 3D Grid
  /* setting the values into the object */
  pfts_Sizes *p_pfts_Sizes = new pfts_Sizes(npart,nrot,nk);
  mesh_3D *p_mesh_3D = new mesh_3D(s_polar);
  //TODO: fix the template argument when done
  //mesh_1D *p_mesh_1D = new mesh_1D(s_polar);
  /* device allocation */
  cuFloatComplex *d_A = NULL;   //device pointer for A (const in)
  cuFloatComplex *d_B = NULL;   //device pointer for B (const in)
  float          *d_C = NULL;   //device pointer for C (re[A*Bstar])
  /* checking arrays for debug true */
  float *h_C          = NULL;
  /* creating the working space on the GPU */
  //float      *d_chunk = NULL;   //device pointer for C (re[A*Bstar])
  float *h_partial_1D = NULL;
  float *d_partial_1D = NULL;
  /* summed values */
  double sum;                    // summed value of the partial sums
  double *sum_dble_vec = NULL;      // vector to check double precision sum
  float       *sum_vec = NULL;      // The CPU version
  float     *h_ptr_vec = NULL;      // The GPU version

  /* size of the element in consideration */
  int size_m = npts*sizeof(cuFloatComplex);
  int size_reC = npts*sizeof(float);
  int size_1D = npart * sizeof(float);
  int size_1D_dble = npart*sizeof(double);
  /* indexer */
  int i,j, ipart;
  /* the error handlers from the cuda library */
  cudaError_t err;
  /* printer index depth */
  float depth = 5;
  
  /*start of the execution commands */

  //nx=16; ny=16; nz=4; dim3 threads(16,16,4); */
  nx=s_polar->nx; ny=s_polar->ny; nz=s_polar->nz;
  dim3 threads(nx,ny,nz);
  gridx = npart/(float)nx+(npart%nx!=0);
  gridy =  nrot/(float)ny+( nrot%ny!=0);
  gridz =    nk/(float)nz+(   nk%nz!=0);
  dim3 grid(gridx,gridy,gridz);
  
  if (s_debug_gpu->debug_i == true ) {
    cuFloatComplex *C = (cuFloatComplex*)malloc(npts);

    rc = print_3D_mesh(0,p_mesh_3D,s_polar,p_pfts_Sizes,nx,ny,nz);
    
    /* printing input variables on screen */
    if ( s_debug_gpu->debug_high_i == true ) {
      rc = print_function_header_X_N_info(s_polar,
                                          r,
                                          C, A, B, 
                                          npart, nrot, nk,
                                          alpha);
    }
  }
  
  //allocating the memory on GPU
  err = cudaMalloc((void**)&d_A, size_m);
  err = cudaMalloc((void**)&d_B, size_m);
  //size_reC = npts*sizeof(float);
  err = cudaMalloc((void**)&d_C, size_reC);
  //uploading the matrices A, B to GPU
  err = cudaMemcpy(d_A, A, size_m, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_B, B, size_m, cudaMemcpyHostToDevice);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  //computing the r product
  xprod_3D_mat<<<grid,threads>>>( d_C, d_A, d_B, npart, nrot, nk, alpha);
  //synchronizing CPU with GPU
  cuCtxSynchronize();

  h_C = (float*)malloc(size_reC);
  sum_vec = (float*)malloc(size_1D);
  sum_dble_vec = (double*)malloc(size_1D_dble);
  if (s_debug_gpu->debug_i == true ) {
    err = cudaMemcpy(h_C, d_C, size_reC, cudaMemcpyDeviceToHost);

    //To be OpenMP if cannot get the double sum to work on GPU
#pragma omp parallel default(shared) private (i,j)
#pragma omp for schedule(auto)
    for (i=0; i<npart ; i++) {
      sum_dble_vec[i] = 0.0;
      sum_vec[i] = 0.0;
      for (j=0; j<nrot*nk ; j++) {
        sum_vec[i] += h_C[j*npart+i];
        sum_dble_vec[i] += (double)h_C[j*npart+i];
      }
    }
  }

  //now summing the over the chunk and mapping it to vector r  
  int chunk = nrot * nk;
  //  int size_chunk = chunk * sizeof(float);
  int threadsPerBlock = 256; //threads/per block=16*16 for more //lism
  int blocksPerGrid = (chunk/256)+(chunk%256!=0);//imin(32,(chunk+threadsPerBlock-1)/threadsPerBlock);
  int size_pB = blocksPerGrid * sizeof(float);
  if ( s_debug_gpu->debug_i == true ) {
    rc = print_summing_method(1);
    rc = print_1D_mesh(0, NULL, chunk, threadsPerBlock, blocksPerGrid); 
  }

  // first on GPU
  h_partial_1D = (float*)malloc(size_pB);
  err = cudaMalloc(&d_partial_1D,size_pB);

  //Start loop here
  for (ipart = 0 ; ipart < npart ; ipart++ ) {
    float *startptr = d_C + ipart;
    r[ipart] = 0.0;
    //computing on GPU
    mysum_23D<<<blocksPerGrid,threadsPerBlock>>>(d_partial_1D, startptr, chunk, npart);
    err = cudaMemcpy(h_partial_1D,d_partial_1D,size_pB,cudaMemcpyDeviceToHost);
    //looping over the partial blocks and summing.
    sum = 0.0;
    for ( int igrid=0 ; igrid < blocksPerGrid ; igrid++ ) {
      sum += (double)h_partial_1D[igrid];
    }
    r[ipart] = sum;
  }//finish the loop here.
  
  // second on cpu
  h_ptr_vec = (float*)malloc(npart*sizeof(float));
  if (s_debug_gpu->debug_i == true ) {
    for (i = 0 ; i < npart ; i++ ) {
      float *startptr = h_C + i;
      h_ptr_vec[i] = 0.0;
      for (int ichunk = 0 ; ichunk < chunk ; ichunk++ ) {
        float *endptr = startptr + ichunk*npart;
        h_ptr_vec[i] += *endptr;
      }
    }
  }

  //printing the resulting vec sums to check

  if(s_debug_gpu->debug_i == true){
    rc = print_sumVecs_X_N(npart, sum_vec, h_ptr_vec,r, sum_dble_vec, depth);}
  
  //freeing the temporary resources on CPU
  free(h_C);
  free(sum_vec);
  free(sum_dble_vec);
  free(h_ptr_vec);
  
  //freeing resources on device
  err = cudaFree(d_A);
  err = cudaFree(d_B);
  err = cudaFree(d_C);
  /* calling the destructor for the pfts_Sizes object */
  p_mesh_3D->~mesh_3D();
  //p_mesh_1D->~mesh_1D(); 

  return rc;
 
} /* End of polarft_corr_P_N */

#endif /* CUDA */
