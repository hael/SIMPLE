
/* 
 * File:   GetDevices.cpp
 * Author: Frederic D.R. Bonnet
 *
 * Created on April 12, 2011, 11:19 AM
 *
 * routine to get the devices cdDevices on the NVIDIA platform id
 * and set to the class
 * set_cpPlatformGPU on the GPU in OpenCL 
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

#include <simul.h>

//extern "C" void getGPU_interface( int* );

void getGPU_interface( int* nplat )
{
  // Program Vars
  int TheNumberOfPlatforms = 0;
  float alpha;
  float *vector1;
  //  float RET;
  unsigned int N;
  size_t szLocalWorkSize;               // 1D var for # of work items in the work group	

  //OpenCL varaibles
  /*  //  cl_platform_id cpPlatform = NULL;
  //  cl_uint ciDeviceCount = 0;
  cl_device_id *cdDevices = NULL;
  cl_context cxGPUContext;
  cl_program cpProgram;                // OpenCl program
  cl_int ciErrNum = CL_SUCCESS;
  cl_command_queue *commandQueue = NULL;
  */

  size_t szGlobalWorkSize_2;           // 1D var for Total # of work items
  size_t szLocalWorkSize_2;	       // 1D var for # of work items in the work group	
  //cl_mem cmDeVvector1;               // OpenCL device destination buffer 
  //cl_context cxGPUContext_2;

  //start of the execution commends

  //quick dummy initialization

  TheNumberOfPlatforms = *nplat;
  printf("hello getGPU_interface_\n");
  printf("the number of platforms is: %i\n",TheNumberOfPlatforms);
  printf("\n\nThe Big\n\n");
  shrLog ("\n");
  shrLog ("Starting...\n\n");

  //starting some log files
  StartLogs("matrixmul.log");

  //get the platform ID
  PlatformID();

  // Get the system info
  //  getplatformdetails ( cpPlatform );
  getPlatformDetails();

  //Get the devices 
  GetDevices();

  //create the context
  //cxGPUContext = CreateContext(ciDeviceCount, cpPlatform );
  CreateContext();

  //create the command queue
  //commandQueue = CreateCommandQueue(ciDeviceCount , cxGPUContext);
  CreateCommandQueue();

  //create the build of the source code for the GPU.
  //cpProgram = CreateBuild(cxGPUContext,cSourcePath);
  
  SourcePathGPU *path = new SourcePathGPU(cSourcePath);
  path -> set_SourcePathGPU(cSourcePath);              // setting the path
  SourceFileGPU *file = new SourceFileGPU(cSourceFile);
  file -> set_SourceFileGPU(cSourceFile);              // setting the file

  shrLog("the source path: %s and source files: %s used",path->get_SourcePathGPU(),file -> get_SourceFileGPU());

  //CreateBuild(cSourcePath);
  CreateBuild(path->get_SourcePathGPU(),file -> get_SourceFileGPU());

  //start of the routine to send to the GPU
  alpha = 6.0;
  N = 10;
  szLocalWorkSize = 32;

  worksizeGPU -> set_szVectorGPU(N);
  worksizeGPU -> set_szLocalWorkSizeGPU(szLocalWorkSize);
  worksizeGPU -> set_szGlobalWorkSizeGPU();

  vector1 = (float *)malloc(sizeof(cl_float) * worksizeGPU -> get_szVectorGPU());
  printf("vector1 = [");
  for (int i = 0 ; i < N ; i++)
    {
      vector1[i]=i;
      printf("%4.2f ",vector1[i]);
    }
  printf("]\n");

  //cxGPUContext_2 = context -> get_cxGPUContextGPU();
  //creating the buffer to send to the routine for the GPU
  //cmDeVvector1  = clCreateBuffer(cxGPUContext_2, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*N, vector1, &ciErrNum);

  CreateCmDeVvector(vector1);

  //cmDeVvector1 = DeVvectorGPU -> get_szVectorGPU();

  szLocalWorkSize_2 = worksizeGPU -> get_szLocalWorkSizeGPU();
  szGlobalWorkSize_2 = worksizeGPU -> get_szGlobalWorkSizeGPU();

  shrLog("Global Work Size \t\t= %u\nLocal Work Size \t\t= %u\n# of Work Groups \t\t= %u\n\n", 
  	 szGlobalWorkSize_2, szLocalWorkSize_2, (szGlobalWorkSize_2 % szLocalWorkSize_2 + szGlobalWorkSize_2/szLocalWorkSize_2));

  routine_saX(N, alpha, 1);

  shrLog("size of cl_mem: %d \n",sizeof(cl_mem));

  /*
  //releasing the ocl

  //ciErrNum = clReleaseMemObject(cmDeVvector1);
  //if( ciErrNum != CL_SUCCESS) shrLog("Error: Failed to release context: %d\n", ciErrNum);
  
  ciErrNum = clReleaseProgram(cpProgram);
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in getGPU_interface, to release cpProgram near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
  else
    {
      shrLog("Success %i: success in getGPU_interface, to release cpProgram near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
  ciErrNum = clReleaseContext(cxGPUContext);
  if( ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in getGPU_interface, to release cxGPUContext near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
  else
    {
      shrLog("Success %i: success in getGPU_interface, to release cxGPUContext near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }
  */

  //freeing the memory
  free(vector1);

  // finish

  shrLogEx(LOGBOTH | CLOSELOG, 0, "Exiting...\n"); 

}

extern "C" void getgpu_interface_gpu_(int*) __attribute__((weak,alias("getGPU_interface")));

#endif /* OPENCL */
