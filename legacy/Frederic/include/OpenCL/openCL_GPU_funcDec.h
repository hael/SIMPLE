/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 11:19am
 *
 * Name:
 * openCl_GPU_funcDec.h - basic definitions used in all modules.
 *
 * Description:
 * header file for all the class definitions and function.
 *******************************************************************************
 */
#ifndef _OPENCL_GPU_FUNCDEC_H
#define _OPENCL_GPU_FUNCDEC_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <complex.h>
#include <string.h>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

#include <complex.h>

typedef _Complex double doubleComplex;

//The extern cpp routine to interface with the f90

extern "C" void getGPU_interface( int* );

//the extern cpp routines

extern "C" void StartLogs(char* );
extern "C" void getPlatformDetails();
extern "C" void PlatformID();
extern "C" void GetDevices();
extern "C" void CreateContext();
extern "C" void CreateCommandQueue();
extern "C" void CreateBuild(char*, char*);
extern "C" void routine_saX(int, float, int);
extern "C" int routine_ZDgem(doubleComplex*, int, float, float, int);
extern "C" void CreateCmDeVvector(float *);
extern "C" void Cleanup(int);

extern "C" cl_kernel Createkernel(cl_program , char * );
extern "C" void LaunchKernel( cl_command_queue , cl_kernel , size_t, size_t );


//the classes declarations 

//the constructor
class cpPlatformGPU
{
 private:

 public:
  cl_platform_id cpPlatform;
  cpPlatformGPU(cl_platform_id Platform);
  void set_cpPlatformGPU(cl_platform_id Platform);
  cl_platform_id get_cpPlatformGPU();
  ~cpPlatformGPU();
};

//the constructor
class ciDeviceCountGPU
{
private:

public:
  cl_uint ciDeviceCount;
  ciDeviceCountGPU( cl_uint DeviceCount );
  void set_ciDeviceCountGPU(cl_uint DeviceCount);
  cl_uint get_ciDeviceCountGPU();
  ~ciDeviceCountGPU();
};

//the constructor
class cdDevicesGPU
{
private:

public:
  cl_device_id* cdDevices;
  cdDevicesGPU( cl_device_id* Devices );
  void set_cdDevicesGPU( cl_device_id* Devices);
  cl_device_id* get_cdDevicesGPU();
  ~cdDevicesGPU();
};

//the constructor
class cxGPUContextGPU
{
private:

public:
  cl_context cxGPUContext;
  cxGPUContextGPU(cl_context GPUContext );
  void set_cxGPUContextGPU(cl_context GPUContext);
  cl_context get_cxGPUContextGPU();
  ~cxGPUContextGPU();
};

//the constructor
class cqCommandQueueGPU
{
private:

public:
  cl_command_queue* cqCommandQueue;
  cqCommandQueueGPU(cl_command_queue* CommandQueue);
  void set_cqCommandQueueGPU(cl_command_queue* CommandQueue);
  cl_command_queue* get_cqCommandQueueGPU();
  ~cqCommandQueueGPU();
};

//the constructor
class cpProgramGPU
{
private:

public:
  cl_program cpProgram;
  cpProgramGPU(cl_program Prog );
  void set_cpProgramGPU(cl_program Prog);
  cl_program get_cpProgramGPU();
  ~cpProgramGPU();
};

//the constructor
class SourceFileGPU
{
private:

public:
  char* cSourceFile;
  SourceFileGPU(char* Filename);
  void set_SourceFileGPU(char* Filename);
  char* get_SourceFileGPU();
  ~SourceFileGPU();
};

//the constructor
class SourcePathGPU
{
private:

public:
  char* cSourcePath;
  SourcePathGPU(char* Path);
  void set_SourcePathGPU(char* Path);
  char* get_SourcePathGPU();
  ~SourcePathGPU();
};

//the constructor
class szWorkSizeGPU
{
 private:

 public:
  size_t szGlobalWorkSize;
  size_t szLocalWorkSize;
  unsigned int N;
  szWorkSizeGPU(size_t LocalWorkSize, unsigned int nn );
  void set_szVectorGPU(unsigned int nn);
  void set_szLocalWorkSizeGPU(size_t LocalWorkSize);
  void set_szGlobalWorkSizeGPU();
  unsigned int get_szVectorGPU();
  size_t get_szLocalWorkSizeGPU();
  size_t get_szGlobalWorkSizeGPU();
  ~szWorkSizeGPU();
};

//the constructor
class cmDeVvectorGPU
{
 private:

 public:
  cl_mem cmDeVvector;
  cmDeVvectorGPU(cl_mem Vvector);
  void set_szVectorGPU(cl_mem Vvector);
  cl_mem get_szVectorGPU();
  ~cmDeVvectorGPU();
};

#endif /* OPENCL */

#endif /* OPENCL_GPU_FUNCDEC_H */
