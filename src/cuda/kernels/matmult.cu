/**
 * Tau (https://www.cs.uoregon.edu/research/tau/) Cuda example
 * Modified by Michael Eager (michael.eager@monash.edu) 2018
 */


#include <stdlib.h>
#include <unistd.h>
#include "cuda.h"
#include "cuda_runtime_api.h"

#define SIZE_OF_MATRIX 1000
#define SIZE_OF_BLOCK 16



#define idx(i,j,lda) ((j) + ((i)*(lda)))

__global__ void multiply_matrices(float *d_a, float *d_b, float *d_c, int lda, int M)
{
	unsigned int row = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int col = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int id  = idx(row,col,lda);

	float ctemp = 0.0;

	if (row < M && col < M)
	{
		for (unsigned int j=0; j<M; j++)
		{
			ctemp = ctemp + d_a[idx(row,j,lda)] * d_b[idx(j,col,lda)];
		}
		d_c[id] = ctemp;
	}
}

__global__ void multiply_matrices_shared_blocks(float *d_a, float *d_b, float *d_c,
                                                int lda, int M, int bs)
{

  //	int bs = SIZE_OF_BLOCK;
	unsigned int row = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int col = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int id  = idx(row,col,lda);

	//submatrices
	float *sub_a, *sub_b;

	//shared submatrices
	__shared__ float a[SIZE_OF_BLOCK][SIZE_OF_BLOCK], b[SIZE_OF_BLOCK][SIZE_OF_BLOCK];
	//temp element of d_c
	float c = 0;

	//top-level row,col of block
	int block_row = blockIdx.y * bs;
	int block_col = blockIdx.x * bs;

	//id inside each block
	int sub_row = threadIdx.y;
	int sub_col = threadIdx.x;

	//for each block
	for (int k = 0; k < (M / bs); k++)
	{

	  sub_a = &d_a[idx(block_row, bs*k, lda)];
		sub_b = &d_b[idx(bs*k, block_col, lda)];
		a[sub_row][sub_col] = sub_a[idx(sub_row, sub_col, lda)];
		b[sub_row][sub_col] = sub_b[idx(sub_row, sub_col, lda)];

		//wait for all threads to complete copy to shared memory.
		__syncthreads();

		//multiply each submatrix
		for (int j=0; j < bs; j++)
		{
			c = c + a[sub_row][j] * b[j][sub_col];
		}

		// move results to device memory.
		d_c[id] = c;

		// wait for multiplication to finish before moving onto the next submatrix.
		__syncthreads();

	}
}

   unsigned int cudaTileDim = 32;
   unsigned int blockRows = 8;
__global__ void transpose_matrix(float **idata, float **odata, int lda, int M)
{
	unsigned int row = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int col = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int id  = idx(row,col,lda);
  __shared__ float tile[cudaTileDim+1][cudaTileDim];
  unsigned int x,y,j;
	float ctemp = 0.0;

	if (row < M && col < M)
	{
    x = (blockIdx%x-1) * cudaTileDim + threadIdx%x;
    y = (blockIdx%y-1) * cudaTileDim + threadIdx%y;

   	for ( j = 0; j<cudaTileDim-1; j = j + blockRows)
      tile[threadIdx%x][threadIdx%y+j] = idata[x][y+j];

     __syncthreads();

     x = (blockIdx%y-1) * cudaTileDim + threadIdx%x;
     y = (blockIdx%x-1) * cudaTileDim + threadIdx%y;

    for ( j = 0; j<cudaTileDim-1; j = j + blockRows)
      odata[x][y+j] = tile[threadIdx%y+j][threadIdx%x];
  }
}



extern "C"
{
void multiply_by_element(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, cudaStream_t cStream)
{

	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	multiply_matrices<<< grid, threads, 0, cStream >>>(d_a, d_b, d_c, m);
  if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);

  if(cudaSuccess != err)
      printf("error in memcpy, #= %s\n", cudaGetErrorString(err));

}

void multiply_by_block(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, cudaStream_t cStream)
{
	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	multiply_matrices_shared_blocks<<< grid, threads, 0, cStream >>>(d_a, d_b, d_c, m);
  if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);
  if(cudaSuccess != err)
      printf("error in memcpy, #= %s\n", cudaGetErrorString(err));

}


void transpose_by_block(dim3 grid, dim3 threads, float *d_a, float *d_c, int m, cudaStream_t cStream)
{
	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	transpose_matrix<<< grid, threads, 0, cStream >>>(d_a, d_c, m);
  if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);
  if(cudaSuccess != err)
      printf("error in memcpy, #= %s\n", cudaGetErrorString(err));

}

}


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <cuComplex.h>
#include "cuda.h"


// Use cuCmulf from cuComplex.h

/* Define CUDA kernel that squares the input complex array */
__global__ void  mul2DComplex(cuFloatComplex *in1,cuFloatComplex *in2, cuFloatComplex *out, int Nsize, int Msize)
{

  unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

  if(n > Nsize || m > Msize) return;
  int grid_width = gridDim.x * blockDim.x;
  int index = m + (n * grid_width);

  out[index] = cuCmulf (in1[index] , in2[index]);

}
// math_functions in cuda.h
__global__ void  mul2DFloat(float *in1,float *in2, float *out, int Nsize, int Msize)
{

  unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m = blockIdx.y * blockDim.y + threadIdx.y;

  if(n > Nsize || m > Msize) return;
  int grid_width = gridDim.x * blockDim.x;
  int index = m + (n * grid_width);

  out[index] = (in1[index] * in2[index]);

}


/*
   Fortran subroutine arguments are passed by references.
   call fun( array_a, array_b, N) will be mapped to
   function (*a, *b, *N);
*/
extern "C"
{

void kernelmul2dfloat_(float *a, float *b, float *c, int *Np, int *Ns, int *Bsize)
{
  int block_size=*Bsize;
  cuFloat *a_d,*b_d,*c_d;
  int N=*Np;int M=*Ns;
  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloat)*N*M );
   cudaMalloc ((void **) &b_d , sizeof(cuFloat)*N*M );
   cudaMalloc ((void **) &c_d , sizeof(cuFloat)*N*M );
  // if(cudaGetLastError() != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloat)*N*M   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloat)*N*M   ,cudaMemcpyHostToDevice);
  // if( cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice) != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid ; //(N/dimBlock.x);
   dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


   if( N*N % block_size != 0 ) dimGrid.x+=1;

  /* Execute the kernel */
   mul2DFloat<<<dimGrid,dimBlock>>>(a_d,b_d,c_d,N,M);
  if( cudaGetLastError()!= cudaSuccess){
    printf("%s\n",cudaGetErrorString(cudaGetLastError()));
    exit(1);
  }

  /* Copy the result back */
  cudaMemcpy( c, c_d, sizeof(cuFloatComplex)*N*M,cudaMemcpyDeviceToHost);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Free memory on the device */
  cudaFree(a_d, b_d, c_d);

  return;
}

void kernelmul2dcomplex_(cuFloatComplex *a, cuFloatComplex *b, cuFloatComplex *c, int *Np, int *Ns, int *Bsize)
{
  int block_size=*Bsize;
  cuFloatComplex *a_d,*b_d;
  int N=*Np;int M=*Ns;
  cudaSetDevice(0);

  /* Allocate complex array on device */
   cudaMalloc ((void **) &a_d , sizeof(cuFloatComplex)*N*M*2 );
   cudaMalloc ((void **) &b_d , sizeof(cuFloatComplex)*N*M*2 );
  // if(cudaGetLastError() != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy array from host memory to device memory */
  cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
  cudaMemcpy( b_d, b,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice);
 // if( cudaMemcpy( a_d, a,  sizeof(cuFloatComplex)*N*M   ,cudaMemcpyHostToDevice) != cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Compute execution configuration */
   dim3 dimBlock(block_size, 8);
   dim3 dimGrid ; //(N/dimBlock.x);
   dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
   dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;


   if( N*N % block_size != 0 ) dimGrid.x+=1;

  /* Execute the kernel */
  mul2DComplex<<<dimGrid,dimBlock>>>(a_d,b_d,b_d,N,M);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Copy the result back */
  cudaMemcpy( c, b_d, sizeof(cuFloatComplex)*N*M,cudaMemcpyDeviceToHost);
  // if( cudaGetLastError()!= cudaSuccess){
  // printf("%s\n",cudaGetErrorString(cudaGetLastError()));
  // exit(1);
  // }

  /* Free memory on the device */
  cudaFree(a_d);

  return;
}

}
