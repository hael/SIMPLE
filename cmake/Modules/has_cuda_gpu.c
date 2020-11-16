#include <stdio.h>
#include <cuda_runtime.h>
/* #include <../extras/Debugger/include/cudadebugger.h> */

/* int main() { */
/*   int deviceCount, device; */
/*   int gpuDeviceCount = 0; */
/*   struct cudaDeviceProp properties; */
/*   cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount); */
/*   if (cudaResultCode != cudaSuccess){ */
/*     printf(" CUDA error: %lu %s\n", cudaResultCode, cudaGetErrorString(cudaResultCode)); */
/*     deviceCount = 0; */
/*   } */
/*    printf("%d CUDA device(s) found\n", deviceCount); */
/*   /\* machines with no GPUs can still report one emulation device *\/ */
/*   for (device = 0; device < deviceCount; ++device) { */
/*     cudaGetDeviceProperties(&properties, device); */
/*     if (properties.major != 9999) /\* 9999 means emulation only *\/ */
/*       ++gpuDeviceCount; */
/*   } */
/*   printf("%d GPU CUDA device(s) found\n", gpuDeviceCount); */

/*   /\* don't just return the number of gpus, because other runtime cuda */
/*      errors can also yield non-zero return values *\/ */
/*   if (gpuDeviceCount > 0) */
/*     return 0; /\* success *\/ */
/*   else */
/*     return 1; /\* failure *\/ */
/* } */
#include <stdio.h>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;
#include <helper_functions.h>
#include <helper_cuda.h>

int main(int argc, char **argv)
{
    int nkernels = 8;               // number of concurrent kernels
    int nstreams = nkernels + 1;    // use one more stream than concurrent kernel
    int nbytes = nkernels * sizeof(clock_t);   // number of data bytes
    float kernel_time = 10; // time the kernel should run in ms
    float elapsed_time;   // timing variables
    int cuda_device = 0;

    printf("[%s] - Starting...\n", argv[0]);

    // get number of kernels if overridden on the command line
    if (checkCmdLineFlag(argc, (const char **)argv, "nkernels"))
    {
        nkernels = getCmdLineArgumentInt(argc, (const char **)argv, "nkernels");
        nstreams = nkernels + 1;
    }

    // use command-line specified CUDA device, otherwise use device with highest Gflops/s
    cuda_device = findCudaDevice(argc, (const char **)argv);

    cudaDeviceProp deviceProp;
    checkCudaErrors(cudaGetDevice(&cuda_device));

    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, cuda_device));

    if ((deviceProp.concurrentKernels == 0))
    {
        printf("> GPU does not support concurrent kernel execution\n");
        printf("  CUDA kernel runs will be serialized\n");
    }

    printf("> Detected Compute SM %d.%d hardware with %d multi-processors\n",
           deviceProp.major, deviceProp.minor, deviceProp.multiProcessorCount);

    return 0;
}
