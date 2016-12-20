/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 21st May 2015
 *
 *      May 2015
 *
 *   Header files for cpp code to obtain information from the GPU cards
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

#ifndef _GET_DEVICEQUERY_GPU_H_
#define _GET_DEVICEQUERY_GPU_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- routines used to interface back to fortran
   *
   *  FORTRAN API -
   * 
   */
#if defined (CUDA) /*preprossing for the CUDA environment */

  /* sets the device for computation */
  int setDevice(int *idev);
  /* wrapper function to find the most suitable device on system */
  int get_findCudaDevice(int *idev);
  /* helper methods from CUDA */  
  inline int _ConvertSMVer2Cores(int major, int minor);
  /*error and warning handlers methods */
  int get_warning_message();
  int get_error_id (cudaError_t error_id);
  int get_cuda_cores_error(int ncores);
  int get_ecc_warning_message(deviceDetails_t *devD);
  /* quering handlers methods */
  int findCudaDevice(cudaDeviceProp deviceProp);
  inline int gpuGetMaxGflopsDeviceId();
  int get_dev_count(deviceDetails_t *devD);
  int get_dev_Name(int *idev, deviceDetails_t *devD);
  int deviceDetails_init(int *idev, deviceDetails_t *devD);
  int deviceDetails_free(deviceDetails_t *devD);
  int get_cuda_driverVersion(int *driverVersion, deviceDetails_t *devD);
  int get_cuda_runtimeVersion(int *runtimeVersion, deviceDetails_t *devD);
  int get_tot_global_mem_MB(float *totalGlobalMem_MB, deviceDetails_t *devD);
  int get_tot_global_mem_bytes(unsigned long long *totalGlobalMem_bytes,
			       deviceDetails_t *devD);
  int get_CUDA_cores(int *idev, int *gpu_major, int *gpu_minor, int *nmp,
		     int *cuda_cores_per_mp, int *ncuda_cores,
		     deviceDetails_t *devD);
  int get_thread_details(int *warpSize, int *max_threads_per_mp,
			 int *max_threads_per_blk,
			 deviceDetails_t *devD);
  int get_nregisters(int *nregisters, deviceDetails_t *devD);
  int get_eec_support(int *is_ecc, int *idev, deviceDetails_t *devD);
  inline bool checkCudaCapabilities(int major_version, int minor_version);
  int get_peer_to_peer_capabilities(deviceDetails_t *devD);

#endif /* CUDA */

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _GET_DEVICEQUERY_GPU_H_ */
