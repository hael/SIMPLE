/*  Based on Tau 2.27 example gpu/cuda/openmp
 *  Modified by Michael Eager (michael.eager@monash.edu) 2018
 */

#include <omp.h>
#include <stdio.h>      // stdio functions are used since C++ streams aren't necessarily thread safe
#include <unistd.h>
#include "cuda.h"

// a simple kernel that simply increments each array element by b
__global__ void kernelAddConstantInt(float *g_a, const int b)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    g_a[idx] += b;
}
// a simple kernel that simply increments each array element by b
__global__ void kernelMulConstantInt(float *g_a, const int b)
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
__global__ void kernelAddVecFloat(float *g_a, const float *b)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    g_a[idx] += b[idx];
}
// a simple kernel that multiplies each array element by b
__global__ void kernelMulVecFloat(float *g_a, const float* b)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    g_a[idx] *= b[idx];
}

extern "C"
{
    void multigpu_vecadd_(float *A, float *B, int N)
    //dim3 *dimGrid, dim3 *dimBlk,
    //int N, cudaStream_t *stream)
    {
        int num_gpus = 0;       // number of CUDA GPUs
        unsigned int nbytes = sizeof(A) * 4;
        if(N * 4 != nbytes) {
            printf(" multigpu_vecaddInt: incorrect N %d and sizeof(A) %lu\n", N, sizeof(A));
            return ;
        }
        /////////////////////////////////////////////////////////////////
        // determine the number of CUDA capable GPUs
        //
        cudaGetDeviceCount(&num_gpus);
        if(num_gpus < 1) {
            printf("no CUDA capable devices were detected\n");
            return ;
        }
        // run as many CPU threads as there are CUDA devices
        //   each CPU thread controls a different device, processing its
        //   portion of the data.  It's possible to use more CPU threads
        //   than there are CUDA devices, in which case several CPU
        //   threads will be allocating resources and launching kernels
        //   on the same device.  For example, try omp_set_num_threads(2*num_gpus);
        //   Recall that all variables declared inside an "omp parallel" scope are
        //   local to each CPU thread
        //
        // initialize data

        unsigned int n = num_gpus * 8192;
        omp_set_num_threads(num_gpus);  // create as many CPU threads as there are CUDA devices
        //omp_set_num_threads(2*num_gpus);// create twice as many CPU threads as there are CUDA devices
        printf("--------multigpu_vecaddInt Starting %d-GPUs\n", num_gpus);
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

            float *d_a = 0;   // pointer to memory on the device associated with this CPU thread
            float *sub_a = A + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data

            unsigned int nbytes_per_kernel = nbytes / num_cpu_threads;
            dim3 gpu_threads(128);  // 128 threads per block
            dim3 gpu_blocks(n / (gpu_threads.x * num_cpu_threads));

            cudaMalloc((void**)&d_a, nbytes_per_kernel);
            cudaMemset(d_a, 0, nbytes_per_kernel);
            cudaMemcpy(d_a, sub_a, nbytes_per_kernel, cudaMemcpyHostToDevice);

            kernelAddConstantFloat <<< gpu_blocks, gpu_threads>>>(d_a, *B);
            cudaMemcpy(sub_a, d_a, nbytes_per_kernel, cudaMemcpyDeviceToHost);
            cudaFree(d_a);


        }
        printf("--------multigpu_vecaddInt Done------------\n");

        if(cudaSuccess != cudaGetLastError())
            printf("%s\n", cudaGetErrorString(cudaGetLastError()));

        cudaThreadExit();
        cudaDeviceSynchronize();
        cudaDeviceReset();
        return ;
    }

    void multigpu_vecmulconstant_(float *A, float *B, int N)
    //dim3 *dimGrid, dim3 *dimBlk,
    //int N, cudaStream_t *stream)
    {
        int num_gpus = 0;       // number of CUDA GPUs


        unsigned int nbytes = sizeof(A) * 4;
        if(N * 4 != nbytes) {
            printf(" multigpu_vecaddInt: incorrect N %d and sizeof(A) %lu\n", N, sizeof(A));
            return ;
        }
        /////////////////////////////////////////////////////////////////
        // determine the number of CUDA capable GPUs
        //
        cudaGetDeviceCount(&num_gpus);
        if(num_gpus < 1) {
            printf("no CUDA capable devices were detected\n");
            return ;
        }
        // run as many CPU threads as there are CUDA devices
        //   each CPU thread controls a different device, processing its
        //   portion of the data.  It's possible to use more CPU threads
        //   than there are CUDA devices, in which case several CPU
        //   threads will be allocating resources and launching kernels
        //   on the same device.  For example, try omp_set_num_threads(2*num_gpus);
        //   Recall that all variables declared inside an "omp parallel" scope are
        //   local to each CPU thread
        //
        // initialize data

        unsigned int n = num_gpus * 8192;
        omp_set_num_threads(num_gpus);  // create as many CPU threads as there are CUDA devices
        //omp_set_num_threads(2*num_gpus);// create twice as many CPU threads as there are CUDA devices
        printf("--------multigpu_vecaddInt Starting %d-GPUs\n", num_gpus);
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

            float *d_a = 0;   // pointer to memory on the device associated with this CPU thread
            float *sub_a = A + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data

            unsigned int nbytes_per_kernel = nbytes / num_cpu_threads;
            dim3 gpu_threads(128);  // 128 threads per block
            dim3 gpu_blocks(n / (gpu_threads.x * num_cpu_threads));

            cudaMalloc((void**)&d_a, nbytes_per_kernel);
            cudaMemset(d_a, 0, nbytes_per_kernel);
            cudaMemcpy(d_a, sub_a, nbytes_per_kernel, cudaMemcpyHostToDevice);

            kernelMulConstantFloat <<< gpu_blocks, gpu_threads>>>(d_a, *B);
            cudaMemcpy(sub_a, d_a, nbytes_per_kernel, cudaMemcpyDeviceToHost);
            cudaFree(d_a);


        }
        printf("--------multigpu_vecaddInt Done------------\n");

        if(cudaSuccess != cudaGetLastError())
            printf("%s\n", cudaGetErrorString(cudaGetLastError()));

        cudaThreadExit();
        cudaDeviceSynchronize();
        cudaDeviceReset();
        return ;
    }


    void multigpu_vecaddvec_(float *A, float *B, int N)
    //dim3 *dimGrid, dim3 *dimBlk,
    //int N, cudaStream_t *stream)
    {
        int num_gpus = 0;       // number of CUDA GPUs
        unsigned int nbytes = sizeof(A) * 4;
        if(N * 4 != nbytes && sizeof(B) != N) {
            printf(" multigpu_vecaddInt: incorrect N %d and sizeof(A) %lu and sizeof(A) %lu\n", N, sizeof(A) , sizeof(B));
            return ;
        }
        /////////////////////////////////////////////////////////////////
        // determine the number of CUDA capable GPUs
        //
        cudaGetDeviceCount(&num_gpus);
        if(num_gpus < 1) {
            printf("no CUDA capable devices were detected\n");
            return ;
        }
        // run as many CPU threads as there are CUDA devices
        //   each CPU thread controls a different device, processing its
        //   portion of the data.  It's possible to use more CPU threads
        //   than there are CUDA devices, in which case several CPU
        //   threads will be allocating resources and launching kernels
        //   on the same device.  For example, try omp_set_num_threads(2*num_gpus);
        //   Recall that all variables declared inside an "omp parallel" scope are
        //   local to each CPU thread
        //
        // initialize data

        unsigned int n = num_gpus * 8192;
        omp_set_num_threads(num_gpus);  // create as many CPU threads as there are CUDA devices
        //omp_set_num_threads(2*num_gpus);// create twice as many CPU threads as there are CUDA devices
        printf("--------multigpu_vecaddInt Starting %d-GPUs\n", num_gpus);
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

            float *d_a = 0;   // pointer to memory on the device associated with this CPU thread
            float *sub_a = A + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data
            float *d_b = 0;   // pointer to memory on the device associated with this CPU thread
            float *sub_b = B + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data

            unsigned int nbytes_per_kernel = nbytes / num_cpu_threads;
            dim3 gpu_threads(128);  // 128 threads per block
            dim3 gpu_blocks(n / (gpu_threads.x * num_cpu_threads));
            // Set A on device
            cudaMalloc((void**)&d_a, nbytes_per_kernel);
            cudaMemset(d_a, 0, nbytes_per_kernel);
            cudaMemcpy(d_a, sub_a, nbytes_per_kernel, cudaMemcpyHostToDevice);
            // set B on device
            cudaMalloc((void**)&d_b, nbytes_per_kernel);
            cudaMemset(d_b, 0, nbytes_per_kernel);
            cudaMemcpy(d_b, sub_b, nbytes_per_kernel, cudaMemcpyHostToDevice);

            kernelAddVecFloat <<< gpu_blocks, gpu_threads>>>(d_a, d_b);
            cudaMemcpy(sub_a, d_a, nbytes_per_kernel, cudaMemcpyDeviceToHost);
            cudaFree(d_a);
            cudaFree(d_b);

        }
        printf("--------multigpu_vecaddvec Done------------\n");

        if(cudaSuccess != cudaGetLastError())
            printf("%s\n", cudaGetErrorString(cudaGetLastError()));

        cudaThreadExit();
        cudaDeviceSynchronize();
        cudaDeviceReset();
        return ;
    }
    void multigpu_vecmulvec_(float *A, float *B, int N)
    //dim3 *dimGrid, dim3 *dimBlk,
    //int N, cudaStream_t *stream)
    {
        int num_gpus = 0;       // number of CUDA GPUs
        unsigned int nbytes = sizeof(A) * 4;
        if(N * 4 != nbytes && sizeof(B) != N) {
            printf(" multigpu_vecaddInt: incorrect N %d and sizeof(A) %lu and sizeof(A) %lu\n", N, sizeof(A) , sizeof(B));
            return ;
        }
        /////////////////////////////////////////////////////////////////
        // determine the number of CUDA capable GPUs
        //
        cudaGetDeviceCount(&num_gpus);
        if(num_gpus < 1) {
            printf("no CUDA capable devices were detected\n");
            return ;
        }
        // run as many CPU threads as there are CUDA devices
        //   each CPU thread controls a different device, processing its
        //   portion of the data.  It's possible to use more CPU threads
        //   than there are CUDA devices, in which case several CPU
        //   threads will be allocating resources and launching kernels
        //   on the same device.  For example, try omp_set_num_threads(2*num_gpus);
        //   Recall that all variables declared inside an "omp parallel" scope are
        //   local to each CPU thread
        //
        // initialize data

        unsigned int n = num_gpus * 8192;
        omp_set_num_threads(num_gpus);  // create as many CPU threads as there are CUDA devices
        //omp_set_num_threads(2*num_gpus);// create twice as many CPU threads as there are CUDA devices
        printf("--------multigpu_vecaddInt Starting %d-GPUs\n", num_gpus);
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

            float *d_a = 0;   // pofloater to memory on the device associated with this CPU thread
            float *sub_a = A + cpu_thread_id * n / num_cpu_threads;   // pofloater to this CPU thread's portion of data
            float *d_b = 0;   // pofloater to memory on the device associated with this CPU thread
            float *sub_b = B + cpu_thread_id * n / num_cpu_threads;   // pointer to this CPU thread's portion of data

            unsigned int nbytes_per_kernel = nbytes / num_cpu_threads;
            dim3 gpu_threads(128);  // 128 threads per block
            dim3 gpu_blocks(n / (gpu_threads.x * num_cpu_threads));
            // Set A on device
            cudaMalloc((void**)&d_a, nbytes_per_kernel);
            cudaMemset(d_a, 0, nbytes_per_kernel);
            cudaMemcpy(d_a, sub_a, nbytes_per_kernel, cudaMemcpyHostToDevice);
            // set B on device
            cudaMalloc((void**)&d_b, nbytes_per_kernel);
            cudaMemset(d_b, 0, nbytes_per_kernel);
            cudaMemcpy(d_b, sub_b, nbytes_per_kernel, cudaMemcpyHostToDevice);

            kernelMulVecFloat <<< gpu_blocks, gpu_threads>>>(d_a, d_b);
            cudaMemcpy(sub_a, d_a, nbytes_per_kernel, cudaMemcpyDeviceToHost);
            cudaFree(d_a);
            cudaFree(d_b);

        }
        printf("--------multigpu_vecaddvec Done------------\n");

        if(cudaSuccess != cudaGetLastError())
            printf("%s\n", cudaGetErrorString(cudaGetLastError()));

        cudaThreadExit();
        cudaDeviceSynchronize();
        cudaDeviceReset();
        return ;
    }


}
