#include <omp.h>
#include <stdio.h>      // stdio functions are used since C++ streams aren't necessarily thread safe
#include <unistd.h>
#include "cuda.h"
#include "cuda_runtime_api.h"

// a simple kernel that simply increments each array element by b
void ompVecAddConstantInt(int *g_a, int *b, int N);

// a predicate that checks whether each array elemen is set to its index plus b
int correctResultAddInt(int *data, const int n, const int b)
{
  for(int i = 0; i < n; i++){
    //	printf("%d ne %d + %d\n", data[i],i,b);
    if(data[i] != i + b)
      {
        printf("%d ne %d + %d\n", data[i],i,b);
        return 0;
      }
  }
  return 1;
}

int local_sleep()
{
  sleep(5);
  return 0;
}


int main(int argc, char *argv[])
{
  int num_gpus = 0;       // number of CUDA GPUs

  /////////////////////////////////////////////////////////////////
  // determine the number of CUDA capable GPUs
  //
  cudaGetDeviceCount(&num_gpus);
  if(cudaSuccess != cudaGetLastError())
    printf("add.cu:  %s\n", cudaGetErrorString(cudaGetLastError()));
  if(num_gpus < 1)
    {
      printf("no CUDA capable devices were detected\n");
      return 1;
    }

  /////////////////////////////////////////////////////////////////
  // display CPU and GPU configuration
  //
  printf("number of host CPUs:\t%d\n", omp_get_num_procs());
  printf("number of CUDA devices:\t%d\n", num_gpus);
  for(int i = 0; i < num_gpus; i++)
    {
      cudaDeviceProp dprop;
      cudaGetDeviceProperties(&dprop, i);
      printf("   %d: %s\n", i, dprop.name);
    }
  printf("---------------------------\n");


  ////////////////////////////////////////////////////////////////r
  // initialize data
  //
  unsigned int n = num_gpus * 8192;
  unsigned int nbytes = n * sizeof(int);
  int *a = 0;             // pointer to data on the CPU
  int b = 3;              // value by which the array is incremented
  a = (int*)malloc(nbytes);
  if(0 == a)
    {
      printf("couldn't allocate CPU memory\n");
      return 1;
    }
  for(unsigned int i = 0; i < n; i++)
    a[i] = i;

  ompVecAddConstantInt(a, &b, n);

  printf("---------------------------\n");
  ////////////////////////////////////////////////////////////////
  // check the result
  //
  if(correctResultAddInt(a, n, b))
    printf("Test PASSED\n");
  else
    printf("Test FAILED\n");

  for(unsigned int i = 0; i < n; i++)
    a[i] = i;


  printf("---------------------------\n");
  free(a);    // free CPU memory


  return 0;
}
