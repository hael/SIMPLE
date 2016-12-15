/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 12th of April 2011.15:01am
 *
 * Name:
 * CreateBuild.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to create the build and set to the class
 * set_cpProgramGPU on the GPU in OpenCL
 * classes to set and get the cqCommandQueueGPU on the GPU in OpenCL
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <global.h>
#include <openCL_GPU_funcDec.h>

//extern "C" void CreateBuild( char*);
//extern "C" void Cleanup(int);

void CreateBuild(char* filepath, char* filename)
{

  //local variables

  char* cPathAndName = NULL;      // var for full paths to data, src, etc.
  char* cSourceCL = NULL;         // Buffer to hold source for compilation 
  cl_context cxGPUContext_2;

  //OpenCl Variables
  size_t szKernelLength;	  // Byte size of kernel code
  cl_program cpProgram;           // OpenCl program

  cxGPUContext_2 = context -> get_cxGPUContextGPU();    

  //start of the execution commands

  shrLog("\n");
  shrLog("Now building the source code for the GPU...\n");
  shrLog("\n");

  //read the OpenCL kernel source code 

  //  cPathAndName = shrFindFilePath(cSourceFile,filepath);
  cPathAndName = shrFindFilePath(filename,filepath);
  cSourceCL = oclLoadProgSource(cPathAndName, "", &szKernelLength);
  //  shrLog("oclLoadProgSource (%s) (%s)...\n", filepath,cSourceFile);
  shrLog("oclLoadProgSource (%s) (%s)...\n", filepath,filename);

  //Now create the program 
  cpProgram = clCreateProgramWithSource(cxGPUContext_2, 1 , (const char **)&cSourceCL, &szKernelLength, &ciErrNum);
  shrLog("clCreateProgramWithSource...\n"); 
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in clCreateProgramWithSource, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      free(cPathAndName);                     //freeing thecPathAndName
      free(cSourceCL);                        //freeing the cSourceCL
      clReleaseProgram(cpProgram);            //freeing the cpProgram
      Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clCreateProgramWithSource, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      //      free(cSourceCL);
    }
  //------------------------
  //Now building the program
  //------------------------
  shrLog("clBuildProgram...\n");
  //ciErrNum = clBuildProgram(cpProgram, 0, NULL, "-cl-fast-relaxed-math", NULL, NULL);
  ciErrNum = clBuildProgram(cpProgram, 0, NULL, NULL, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
    {
      shrLog("Error %i: error in clBuildProgram, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
      shrLogEx(LOGBOTH | ERRORMSG, ciErrNum, STDERROR);
      oclLogBuildInfo(cpProgram, oclGetFirstDev(cxGPUContext_2));
      oclLogPtx(cpProgram, oclGetFirstDev(cxGPUContext_2), "error.ptx");
      Cleanup(EXIT_FAILURE);
    }
  else
    {
      shrLog("Success %i: success in clBuildProgram, near Line %i in file %s %s\n",ciErrNum,__LINE__,__FILE__,__FUNCTION__);
    }

  programGPU -> set_cpProgramGPU(cpProgram);//return cpProgram;

}

extern "C" void createbuild_gpu_(char*, char*) __attribute__((weak,alias("CreateBuild")));

#endif /* OPENCL */
