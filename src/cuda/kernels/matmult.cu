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

__global__ void multiply_matrices(float *d_a, float *d_b, float *d_c, int M, int lda)
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

__global__ void multiply_matrices_shared_blocks(float *d_a, float *d_b, float *d_c, int M, int bs,
                                                int lda)
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

#define cudaTileDim  32
#define blockRows  8
__global__ void transpose_matrix_by_shared(float **idata, float **odata, int M, int lda)
{
	unsigned int row = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int col = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int id  = idx(row,col,lda);
  __shared__ float tile[cudaTileDim+1][cudaTileDim];
  unsigned int x,y,j;
	float ctemp = 0.0;

	if (row < M && col < M)
	{
    x = (blockIdx.x-1) * cudaTileDim + threadIdx.x;
    y = (blockIdx.y-1) * cudaTileDim + threadIdx.y;

   	for ( j = 0; j<cudaTileDim-1; j = j + blockRows)
      tile[threadIdx.x][threadIdx.y+j] = idata[x][y+j];

     __syncthreads();

     x = (blockIdx.y-1) * cudaTileDim + threadIdx.x;
     y = (blockIdx.x-1) * cudaTileDim + threadIdx.y;

    for ( j = 0; j<cudaTileDim-1; j = j + blockRows)
      odata[x][y+j] = tile[threadIdx.y+j][threadIdx.x];
  }
}

/*

extern "C"
{
  void multiply_by_element(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, cudaStream_t cStream)
{

	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	multiply_matrices<<< grid, threads, 0, cStream >>>(d_a, d_b, d_c, m, m);
  if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);

  if(cudaSuccess != err)
      printf("error in memcpy, #= %s\n", cudaGetErrorString(err));

}

  void multiply_by_block(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, int bs, cudaStream_t cStream)
{
	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	multiply_matrices_shared_blocks<<< grid, threads, 0, cStream >>>(d_a, d_b, d_c, m, bs, m);
  if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);
  if(cudaSuccess != err)
      printf("error in memcpy, #= %s\n", cudaGetErrorString(err));

}


void transpose_sharedblck(dim3 grid, dim3 threads, float *d_a, float *d_c, int m,  cudaStream_t cStream)
{
	cudaError err;
	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);
	float* c = (float*)malloc(matsize);

	transpose_matrix_by_shared<<< grid, threads, 0, cStream >>>(&d_a, &d_c, m,m);
  if(cudaSuccess != cudaGetLastError())
      printf("error in  transpose_sharedblck, #=%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaStreamSynchronize(cStream);
	err = cudaMemcpyAsync(c, d_c, matsize, cudaMemcpyDeviceToHost, cStream);
  if(cudaSuccess != err)
      printf("error in memcpy, #=%s\n", cudaGetErrorString(err));

}

void kernelmulbyelement(int nSize, int nBlock, float* a, float* b, float* c)
{
	unsigned int number_of_threads = min(nSize, nBlock);
	unsigned int number_of_blocks;
	if (nSize > nBlock)
		number_of_blocks = ceil(nSize / ((float) nBlock));
	else
		 number_of_blocks = 1;

	unsigned int matsize = nSize*nSize*sizeof(float);

	//cout << "blocks: " << number_of_blocks << " threads: " <<
	//number_of_threads << endl;

	//cout.flush();

	cudaError_t err;

	int count = 0;

	err = cudaGetDeviceCount(&count);

	printf("%d devices found.\n",count);

	int number_of_iterations = 1;

	int devices[count];
	int nDevices = 0;
	//default: use all the devices

		for (int d=0;d<count;d++)		{
			devices[d] = d;
		}
		nDevices = count;

	//cout << "finished mapping devices." << endl;
	float *d_a[nDevices], *d_b[nDevices], *d_c[nDevices];
	cudaStream_t streams[nDevices];
	for (int d=0;d<nDevices;d++)	{
		cudaSetDevice(devices[d]);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, devices[d]);
		//cout << "Using device " << devices[d] << ", name: " << deviceProp.name << endl;

		err = cudaSetDevice(devices[d]);
		if (err != cudaSuccess)
					printf("error setting device, #=%s\n", cudaGetErrorString(err) );

		err = cudaStreamCreate(&streams[d]);
		if (err != cudaSuccess)
					printf("error in stream creation, #=%s\n", cudaGetErrorString(err));




		err = cudaMalloc((void **) &d_a[d], matsize);
		if (err != cudaSuccess)
		      printf("error in malloc, #=%s\n", cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_b[d], matsize);
		if (err != cudaSuccess)
		      printf("error in malloc, #=%s\n", cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_c[d], matsize);
		if (err != cudaSuccess)
		      printf("error in malloc, #=%s\n", cudaGetErrorString(err)) ;


	}

	for (int i=0; i<number_of_iterations*nDevices; i++)
	{
		int cDevice = i%nDevices;
		cudaStream_t cStream = streams[cDevice];
		cudaSetDevice(devices[cDevice]);
		if (err != cudaSuccess)
		      printf("error setting device: %d  #=%s\n", devices[i%nDevices] , cudaGetErrorString(err) );


		err = cudaMemcpyAsync(d_a[cDevice], a, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
					 printf("error in memcpy, #=%s\n" << cudaGetErrorString(err) );

		err = cudaMemcpyAsync(d_b[cDevice], b, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
					 printf("error in memcpy, #=%s\n" , cudaGetErrorString(err) );


		//cout << "running on device " << cDevice << endl;

		dim3 grid(number_of_blocks, number_of_blocks);
		dim3 threads(number_of_threads, number_of_threads, 1);

		//multiply each element at a time.
		multiply_by_element(grid, threads, d_a[cDevice], d_b[cDevice], d_c[cDevice], m, cStream);

  }
	for (int d=0;d<nDevices;d++)	{
		cudaSetDevice(devices[d]);
		cudaStreamSynchronize(streams[d]);
	}
	for (int d=0;d<nDevices;d++)	{
		cudaStreamDestroy(streams[d]);
	}




	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	cudaThreadExit();
}

void kernelmulbyblock(int nSize, int nBlock, float* a, float* b, float* c)
{
	unsigned int number_of_threads = min(nSize, nBlock);
	unsigned int number_of_blocks;
	if (nSize > nBlock)
		number_of_blocks = ceil(nSize / ((float) nBlock));
	else
		 number_of_blocks = 1;

	unsigned int matsize = nSize*nSize*sizeof(float);


	cudaError_t err;

	int count = 0;

	err = cudaGetDeviceCount(&count);

	printf(" %d devices found." , count);

	int number_of_iterations = 1;
	int devices[count];
	int nDevices = 0;
		for (int d=0;d<count;d++)		{
			devices[d] = d;
		}
		nDevices = count;
	//printf("finished mapping devices." );
	float *d_a[nDevices], *d_b[nDevices], *d_c[nDevices];
	cudaStream_t streams[nDevices];
	for (int d=0;d<nDevices;d++)	{
		cudaSetDevice(devices[d]);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, devices[d]);
		printf( "Using device %d, name: %s\n" ,devices[d], deviceProp.name );

		err = cudaSetDevice(devices[d]);
		if (err != cudaSuccess)
					printf( "error setting device, #=%s\n" , cudaGetErrorString(err) );
				err = cudaStreamCreate(&streams[d]);
		if (err != cudaSuccess)
					printf( "error in stream creation, #=%s\n" , cudaGetErrorString(err) );




		err = cudaMalloc((void **) &d_a[d], matsize);
		if (err != cudaSuccess)
					printf( "error in malloc, #=%s\n" , cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_b[d], matsize);
		if (err != cudaSuccess)
					printf( "error in malloc, #=%s\n" , cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_c[d], matsize);
		if (err != cudaSuccess)
					printf( "error in malloc, #=%s\n" , cudaGetErrorString(err) );


	}

	for (int i=0; i<number_of_iterations*nDevices; i++)
	{
		int cDevice = i%nDevices;
		cudaStream_t cStream = streams[cDevice];
		cudaSetDevice(devices[cDevice]);
		if (err != cudaSuccess)
			printf( "error setting device: %d #=%s\n" ,devices[i%nDevices], cudaGetErrorString(err) );


		err = cudaMemcpyAsync(d_a[cDevice], a, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
      printf( "error in memcpy, #=%s\n" , cudaGetErrorString(err) );

		err = cudaMemcpyAsync(d_b[cDevice], b, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
					printf( "error in memcpy, #=%s\n" , cudaGetErrorString(err) );


		//cout << "running on device " << cDevice << endl;

		dim3 grid(number_of_blocks, number_of_blocks);
		dim3 threads(number_of_threads, number_of_threads, 1);

    //multiply by first load a 16x16 submatrix into shared memory.
		multiply_by_block(grid, threads, d_a[cDevice], d_b[cDevice], d_c[cDevice], m, cStream);
	}

	printf( "Finished %d iterations on %d devices. \n" , number_of_iterations, nDevices );

	for (int d=0;d<nDevices;d++)	{
		cudaSetDevice(devices[d]);
		cudaStreamSynchronize(streams[d]);
	}
	for (int d=0;d<nDevices;d++)	{
		cudaStreamDestroy(streams[d]);
	}
	//print c


	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	cudaThreadExit();
} // kernelmulbyblock

void kerneltransposebyblock(float *a, float *b, int nSize, int nBlock)
{
	unsigned int number_of_threads = min(nSize, nBlock);
	unsigned int number_of_blocks;
	if (nSize > nBlock)
		number_of_blocks = ceil(nSize / ((float) nBlock));
	else
		 number_of_blocks = 1;

	unsigned int matsize = nSize*nSize*sizeof(float);

	cudaError_t err;

	int count = 0;

	err = cudaGetDeviceCount(&count);

	printf("%d devices found.\n" , count);

	int number_of_iterations = 1;

	int devices[count];
	int nDevices = 0;
	//default: use all the devices

		for (int d=0;d<count;d++)
		{
			devices[d] = d;
		}
		nDevices = count;

	//cout << "finished mapping devices." << endl;
	float *d_a[nDevices], *d_b[nDevices];
	cudaStream_t streams[nDevices];
	for (int d=0;d<nDevices;d++) {
		cudaSetDevice(devices[d]);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, devices[d]);
		printf("Using device %d, name: %s\n",  devices[d], deviceProp.name );

		err = cudaSetDevice(devices[d]);
		if (err != cudaSuccess)
					printf( "error setting device, #=%s\n" , cudaGetErrorString(err) );
    err = cudaStreamCreate(&streams[d]);
		if (err != cudaSuccess)
      printf( "error in stream creation, #=%s\n" , cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_a[d], matsize);
		if (err != cudaSuccess)
      printf( "error in malloc, #=%s\n" , cudaGetErrorString(err) );

		err = cudaMalloc((void **) &d_b[d], matsize);
		if (err != cudaSuccess)
      printf( "error in malloc, #=%s\n" , cudaGetErrorString(err) );


	}

	for (int i=0; i<number_of_iterations*nDevices; i++)	{
		int cDevice = i%nDevices;
		cudaStream_t cStream = streams[cDevice];
		cudaSetDevice(devices[cDevice]);
		if (err != cudaSuccess)
		{
			printf( "error setting device: %d #=%s\n" ,  devices[i%nDevices], cudaGetErrorString(err) );
		}

		err = cudaMemcpyAsync(d_a[cDevice], a, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
		{
			printf( "error in memcpy, #=%s\n" , cudaGetErrorString(err) );
		}
		err = cudaMemcpyAsync(d_b[cDevice], b, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
		{
			printf( "error in memcpy, #=%s\n" , cudaGetErrorString(err) );
		}

		dim3 grid(number_of_blocks, number_of_blocks);
		dim3 threads(number_of_threads, number_of_threads, 1);


		//multiply by first load a 16x16 submatrix into shared memory.
		transpose_matrix_by_shared_block(grid, threads, d_a[cDevice], d_b[cDevice],  m, cStream);
	}

	//cout << "Finished " << number_of_iterations << " iterations on " << nDevices << " devices." << endl;

	for (int d=0;d<nDevices;d++)	{
		cudaSetDevice(devices[d]);
		cudaStreamSynchronize(streams[d]);
	}
	for (int d=0;d<nDevices;d++)	{
		cudaStreamDestroy(streams[d]);
	}

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	cudaThreadExit();
} // kerneltransposebyblock



} // extern C


*/
