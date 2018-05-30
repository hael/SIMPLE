/*  BAsed on Tau 2.27 example gpu/cuda/openmp */

#include <omp.h>
#include <stdio.h>      // stdio functions are used since C++ streams aren't necessarily thread safe
#include <unistd.h>

// a simple kernel that simply increments each array element by b
__global__ void kernelAddConstantInt(int *g_a, const int b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] += b;
}
// a simple kernel that simply increments each array element by b
__global__ void kernelMulConstantInt(int *g_a, const int b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] *= b;
}
// a simple kernel that simply increments each array element by b
__global__ void kernelAddConstantFloat(float *g_a, const float b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] += b;
}
// a simple kernel that multiplies each array element by b
__global__ void kernelMulConstantFloat(float *g_a, const float b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] *= b;
}
// a simple kernel that simply increments each array element by b
__global__ void kernelAddConstantComplex(float *g_a, const float b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] += b;
}
// a simple kernel that multiplies each array element by b
__global__ void kernelMulConstantComplex(complex *g_a, const complex b)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  g_a[idx] *= b;
}


extern "C"
{
  void ompVecAddConstantInt(int *A, int *B, int N)
                            //dim3 *dimGrid, dim3 *dimBlk,
                            //int N, cudaStream_t *stream)
  {
    int num_gpus = 0;       // number of CUDA GPUs
    unsigned int nbytes = sizeof(A)*4;
    if( N*4 != nbytes ){
      printf("ompVecAddConstantInt: incorrect N %d and sizeof(A) %d\n", N, sizeof(A) );
        return 1;
    }
    /////////////////////////////////////////////////////////////////
    // determine the number of CUDA capable GPUs
    //
    cudaGetDeviceCount(&num_gpus);
    if(num_gpus < 1)
      {
        printf("no CUDA capable devices were detected\n");
        return 1;
      }
    omp_set_num_threads(num_gpus);  // create as many CPU threads as there are CUDA devices
    //omp_set_num_threads(2*num_gpus);// create twice as many CPU threads as there are CUDA devices
#pragma omp parallel
    {
      unsigned int cpu_thread_id = omp_get_thread_num();
      unsigned int num_cpu_threads = omp_get_num_threads();
      // local_sleep();
      // set and check the CUDA device for this CPU thread
      int gpu_id = -1;
      cudaSetDevice(cpu_thread_id % num_gpus);        // "% num_gpus" allows more CPU threads than GPU devices
      cudaGetDevice(&gpu_id);

      printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);

      int *d_a = 0;   // pointer to memory on the device associated with this CPU thread
      int *sub_a = A + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data
      unsigned int nbytes_per_kernel = nbytes / num_cpu_threads;
      dim3 gpu_threads(128);  // 128 threads per block
      dim3 gpu_blocks(n / (gpu_threads.x * num_cpu_threads));

      cudaMalloc((void**)&d_a, nbytes_per_kernel);
      cudaMemset(d_a, 0, nbytes_per_kernel);
      cudaMemcpy(d_a, sub_a, nbytes_per_kernel, cudaMemcpyHostToDevice);
      kernelAddConstantInt<<<gpu_blocks, gpu_threads>>>(d_a, b);
      cudaMemcpy(sub_a, d_a, nbytes_per_kernel, cudaMemcpyDeviceToHost);
      cudaFree(d_a);


    }
    printf("---------------------------\n");

    if(cudaSuccess != cudaGetLastError())
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));

    cudaThreadExit();
    cudaDeviceSynchronize();
    cudaDeviceReset();
    return 0;
  }

}
