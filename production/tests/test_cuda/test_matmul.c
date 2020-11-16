#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include "cuda_runtime_api.h"

#define SIZE_OF_MATRIX 1000
#define SIZE_OF_BLOCK 16



void multiply_by_element(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, cudaStream_t cStream);

void multiply_by_block(dim3 grid, dim3 threads, float *d_a, float *d_b, float *d_c, int m, cudaStream_t cStream);


int main(int argc, char** argv)
{
  int m = SIZE_OF_MATRIX;
	unsigned int number_of_threads = min(SIZE_OF_MATRIX, SIZE_OF_BLOCK);
	unsigned int number_of_blocks;
	if (SIZE_OF_MATRIX > SIZE_OF_BLOCK)
		number_of_blocks = ceil(SIZE_OF_MATRIX / ((float) SIZE_OF_BLOCK));
	else
		 number_of_blocks = 1;

	unsigned int matsize = SIZE_OF_MATRIX*SIZE_OF_MATRIX*sizeof(float);

	//cout << "blocks: " << number_of_blocks << " threads: " <<
	//number_of_threads );

	//cout.flush();

	float* a = (float*)malloc(matsize);
	float* b = (float*)malloc(matsize);
	float* c = (float*)malloc(matsize);

	//initalize matrices
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			//a[i*m+j] = i;
			//b[i*m+j] = i;
			a[i*m+j] = i-j*2 + i-j+1 + 1;
			b[i*m+j] = i-j*2 + i-j+1 + 1;
			c[i*m+j] = 0;
			//cout << a[i*m+j] << ", ";
		}
		//cout );
	}
	cudaError_t err;

	int count = 0;

	err = cudaGetDeviceCount(&count);
  if(cudaSuccess != err)
      printf("%s\n", cudaGetErrorString(err));
  printf("%d\t devices found.\n", count);
  if (count <= 0 ) exit(EXIT_FAILURE);

	int number_of_iterations = 1;

	int devices[count];
	int nDevices = 0;
	//default: use all the devices

		for (int d=0;d<count;d++)
		{
			devices[d] = d;
		}
		nDevices = count;

	//cout << "finnished mapping devices." );
	float *d_a[nDevices], *d_b[nDevices], *d_c[nDevices];
	cudaStream_t streams[nDevices];
	for (int d=0;d<nDevices;d++)
	{
		cudaSetDevice(devices[d]);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, devices[d]);
		printf( "Using device %d name: %s\n",  devices[d] ,deviceProp.name);

		err = cudaSetDevice(devices[d]);
		if (err != cudaSuccess)
		{
			printf( "error setting device, #=%s\n" , cudaGetErrorString(err) );
		}
		err = cudaStreamCreate(&streams[d]);
		if (err != cudaSuccess)
		{
			printf("error in stream creation, #=%s\n", cudaGetErrorString(err) );
		}



		err = cudaMalloc((void **) &d_a[d], matsize);
		if (err != cudaSuccess)
		{
			printf("error in malloc, %s\n", cudaGetErrorString(err) );
		}
		err = cudaMalloc((void **) &d_b[d], matsize);
		if (err != cudaSuccess)
		{
			printf("error in malloc, %s\n", cudaGetErrorString(err) );
		}
		err = cudaMalloc((void **) &d_c[d], matsize);
		if (err != cudaSuccess)
		{
			printf("error in malloc, %s\n", cudaGetErrorString(err) );
		}

	}

	for (int i=0; i<number_of_iterations*nDevices; i++)
	{
		int cDevice = i%nDevices;
		cudaStream_t cStream = streams[cDevice];
		cudaSetDevice(devices[cDevice]);
		if (err != cudaSuccess)
		{
			printf("error setting device: %d  %s\n", devices[i%nDevices], cudaGetErrorString(err) );
		}

		err = cudaMemcpyAsync(d_a[cDevice], a, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
		{
			printf("error in memcpy, %s\n", cudaGetErrorString(err) );
		}
		err = cudaMemcpyAsync(d_b[cDevice], b, matsize, cudaMemcpyHostToDevice, cStream);
		if (err != cudaSuccess)
		{
			printf("error in memcpy, %s\n", cudaGetErrorString(err) );
		}

		//printf("running on device " << cDevice );

		dim3 grid(number_of_blocks, number_of_blocks);
		dim3 threads(number_of_threads, number_of_threads, 1);

		//multiply each element at a time.
		multiply_by_element(grid, threads, d_a[cDevice], d_b[cDevice], d_c[cDevice], m, cStream);

		//multiply by first load a 16x16 submatrix into shared memory.
		multiply_by_block(grid, threads, d_a[cDevice], d_b[cDevice], d_c[cDevice], m, cStream);
	}

	printf("Finished %d iterations on %d devices.", number_of_iterations, nDevices );

	for (int d=0;d<nDevices;d++)
	{
		cudaSetDevice(devices[d]);
		cudaStreamSynchronize(streams[d]);
	}
	for (int d=0;d<nDevices;d++)
	{
		cudaStreamDestroy(streams[d]);
	}
	//print c

	printf(" results: " );
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			printf("%10.8f\t ", c[i*m+j] );
		}
        printf("\n");
	}



	//print c
	/*
	printf(" results: " );
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			printf(c[i*m+j] << ", ";
		}
		printf(endl;
	}
	*/
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	cudaThreadExit();
}
