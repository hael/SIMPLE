#include "cuda.h"
#include "cuda_runtime.h"
#include <stdio.h>      // stdio functions are used since C++ streams aren't necessarily thread safe
#include <unistd.h>
#include <errno.h>
// for cuda error checking
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            return ; \
        } \
    } while (0)



inline void GPUassert(cudaError_t code, char * file, int line, int Abort=1)
{
    if (code != 0) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),file,line);
        if (Abort !=0 ) exit(code);
    }
}

#define GPUerrchk(ans) { GPUassert((ans), __FILE__, __LINE__); }
