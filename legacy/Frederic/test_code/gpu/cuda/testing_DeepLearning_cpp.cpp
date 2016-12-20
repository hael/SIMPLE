/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 02nd of May 2016
 *      Monash University
 *      May 2016
 *
 *      Routine which test the cuDNN deep learning library
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <iostream>

#if defined (CUDA) /* prepropressing for the cuda environment */
#include <cuda.h> // need CUDA_VERSION 7.0
#include <cudnn.h>
#endif

#define EXIT_WAIVED 0

#define TOSTR_(s)   #s
#define TOSTR(s)    TOSTR_(s)
#if defined(__GNUC__)
#define COMPILER_NAME "GCC"
#define COMPILER_VER  TOSTR(__GNUC__) "." TOSTR(__GNUC_MINOR__) "." TOSTR(__GNUC_PATCHLEVEL__)
#elif defined(__clang_major__)
#define COMPILER_NAME "CLANG"
#define COMPILER_VER  TOSTR(__clang_major__ ) "." TOSTR(__clang_minor__) "." TOSTR(__clang_patchlevel__)
#endif

#define STRNCASECMP strncasecmp
#define CUDNN_VERSION_STR  TOSTR(CUDNN_MAJOR) "." TOSTR (CUDNN_MINOR) "." TOSTR(CUDNN_PATCHLEVEL)

#define FatalError(s) {                                                \
    std::stringstream _where, _message;                                \
    _where << __FILE__ << ':' << __LINE__;                             \
    _message << std::string(s) + "\n" << __FILE__ << ':' << __LINE__;  \
    std::cerr << _message.str() << "\nAborting...\n";                  \
    cudaDeviceReset();                                                 \
    exit(EXIT_FAILURE);                                                \
}
#define checkCudaErrors(status) {                                      \
    std::stringstream _error;                                          \
    if (status != 0) {                                                 \
      _error << "Cuda failure\nError: " << cudaGetErrorString(status); \
      FatalError(_error.str());                                        \
    }                                                                  \
}
static void showDevices( void )
{
  int totalDevices;
  checkCudaErrors(cudaGetDeviceCount( &totalDevices ));
  printf("\nThere are %d CUDA capable devices on your machine :\n", totalDevices);
  for (int i=0; i< totalDevices; i++) {
    struct cudaDeviceProp prop;
    checkCudaErrors(cudaGetDeviceProperties( &prop, i ));
    printf( "device %d : sms %2d  Capabilities %d.%d, SmClock %.1f Mhz, MemSize (Mb) %d, MemClock %.1f Mhz, Ecc=%d, boardGroupID=%d\n",
            i, prop.multiProcessorCount, prop.major, prop.minor,
            (float)prop.clockRate*1e-3,
            (int)(prop.totalGlobalMem/(1024*1024)),
            (float)prop.memoryClockRate*1e-3,
            prop.ECCEnabled,
            prop.multiGpuBoardGroupID);
  }
}

inline int stringRemoveDelimiter(char delimiter, const char *string)
{
  int string_start = 0;
  while (string[string_start] == delimiter){string_start++;}
  if (string_start >= (int)strlen(string)-1){return 0;}
  return string_start;
}

inline bool checkCmdLineFlag(const int argc, const char **argv,
                             const char *string_ref)
{
  bool bFound = false;

  if (argc >= 1)
    {
      for (int i=1; i < argc; i++)
        {
          int string_start = stringRemoveDelimiter('-', argv[i]);
          const char *string_argv = &argv[i][string_start];

          const char *equal_pos = strchr(string_argv, '=');
          int argv_length = (int)(equal_pos == 0 ? strlen(string_argv) : equal_pos - string_argv);

          int length = (int)strlen(string_ref);

          if (length == argv_length && !STRNCASECMP(string_argv, string_ref, length))
            {
              bFound = true;
              continue;
            }
        }
    }

  return bFound;
}
void get_path(std::string& sFilename, const char *fname, const char *pname)
{sFilename = (std::string("data/") + std::string(fname));}
void displayUsage() {
  printf( "mnistCUDNN {<options>}\n");
  printf( "help                   : display this help\n");
  printf( "device=<int>           : set the device to run the sample\n");
  printf( "image=<name>           : classify specific image\n");
}

int main(int argc, char *argv[])
{   
  std::string image_path;
  //start of the execution commands
  printf("Ecxecutable= %s\n",argv[0]);

  if (checkCmdLineFlag(argc, (const char **)argv, "help"))
    {
      displayUsage();
      exit(EXIT_WAIVED);
    }

  int version = (int)cudnnGetVersion();
  printf("cudnnGetVersion() : %d , CUDNN_VERSION from cudnn.h : %d (%s)\n",
         version, CUDNN_VERSION, CUDNN_VERSION_STR);
  printf("Host compiler version : %s %s\r", COMPILER_NAME, COMPILER_VER);
  showDevices();

}
