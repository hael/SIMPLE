/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 21st May 2015
 *
 *      May 2015
 *
 *      c code for the MACOSX code this wrapps around the the cpp for LINUX
 *      files methods for the cpp files to extract the information from the GPU
 *      devices.
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

#if defined (CUDA) /*preprossing for the CUDA environment */

/* sets the device for computation */
int setdevice_c_(int *idev) {
  int rc = RC_SUCCESS;

  rc = setDevice(idev);
  
  return rc;
}
void SETDEVICE_C_() __attribute__((weak,alias("setdevice_c_")));
void setdevice_c__() __attribute__((weak,alias("setdevice_c_")));
void SETDEVICE_C__() __attribute__((weak,alias("setdevice_c_")));
/* finds the best suitable CUDA device for computation */
int findcudadevice_c_(int *idev) {
  int rc = RC_SUCCESS;

  rc = get_findCudaDevice(idev);
  
  return rc;
  
}
void FINDCUDADEVICE_C_() __attribute__((weak,alias("findcudadevice_c_")));
void findcudadevice_c__() __attribute__((weak,alias("findcudadevice_c_")));
void FINDCUDADEVICE_C__() __attribute__((weak,alias("findcudadevice_c_")));
/* determine if peer to peer is allowed */
int get_peer_to_peer_capabilities_c_(deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_peer_to_peer_capabilities(devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_PEER_TO_PEER_CAPABILITIES_C_() __attribute__((weak,alias("get_peer_to_peer_capabilities_c_")));
void get_peer_to_peer_capabilities_c__() __attribute__((weak,alias("get_peer_to_peer_capabilities_c_")));
void GET_PEER_TO_PEER_CAPABILITIES_C__() __attribute__((weak,alias("get_peer_to_peer_capabilities_c_")));
/* getting the error correcting code memory device true or flase*/
int get_eec_support_c_(int *is_ecc, int *idev, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_eec_support(is_ecc, idev, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_EEC_SUPPORT_C_() __attribute__((weak,alias("get_eec_support_c_")));
void get_eec_support_c__() __attribute__((weak,alias("get_eec_support_c_")));
void GET_EEC_SUPPORT_C__() __attribute__((weak,alias("get_eec_support_c_")));
/* getting the threads details */
int get_thread_details_c_(int *warpSize, int *max_threads_per_mp,
			  int *max_threads_per_blk,
			  deviceDetails_t *devD) {
    int rc = RC_SUCCESS;
    rc = get_thread_details(warpSize, max_threads_per_mp,
			    max_threads_per_blk, devD);
    if (rc == RC_FAIL) { rc = get_error_c(); }
    return rc;
  }
void GET_THREAD_DETAILS_C_() __attribute__((weak,alias("get_thread_details_c_")));
void get_thread_details_c__() __attribute__((weak,alias("get_thread_details_c_")));
void GET_THREAD_DETAILS_C__() __attribute__((weak,alias("get_thread_details_c_")));
/* gettig the number of registers */
int get_nregisters_c_(int *nregisters_per_blk, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_nregisters(nregisters_per_blk, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_NREGISTERS_C_() __attribute__((weak,alias("get_nregisters_c_")));
void get_nregisters_c__() __attribute__((weak,alias("get_nregisters_c_")));
void GET_NREGISTERS_C__() __attribute__((weak,alias("get_nregisters_c_")));
/* getting # of Multiprocessors, CUDA cores/MP and total # of CUDA cores */ 
int get_CUDA_cores_c_(int *idev, int *gpu_major, int *gpu_minor, int *nmp,
		     int *cuda_cores_per_mp, int *ncuda_cores,
		     deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_CUDA_cores(idev, gpu_major, gpu_minor, nmp,
		      cuda_cores_per_mp, ncuda_cores, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_CUDA_CORES_C_() __attribute__((weak,alias("get_CUDA_cores_c_")));
void get_cuda_cores_c__() __attribute__((weak,alias("get_CUDA_cores_c_")));
void GET_CUDA_CORES_C__() __attribute__((weak,alias("get_CUDA_cores_c_")));
/*getting the total amount of memory available on the card in MB */
int get_tot_global_mem_MB_c_(float *totalGlobalMem_MB, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_tot_global_mem_MB(totalGlobalMem_MB, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
}
void GET_TOT_GLOBAL_MEM_MB_C_() __attribute__((weak,alias("get_tot_global_mem_MB_c_")));
void get_tot_global_mem_MB_c__() __attribute__((weak,alias("get_tot_global_mem_MB_c_")));
void GET_TOT_GLOBAL_MEM_MB_C__() __attribute__((weak,alias("get_tot_global_mem_MB_c_")));
/*getting the total amount of memory available on the card in bytes */
int get_tot_global_mem_bytes_c_(unsigned long long *totalGlobalMem_bytes,
				deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_tot_global_mem_bytes(totalGlobalMem_bytes, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
}
void GET_TOT_GLOBAL_MEM_BYTES_C_() __attribute__((weak,alias("get_tot_global_mem_bytes_c_")));
void get_tot_global_mem_bytes_c__() __attribute__((weak,alias("get_tot_global_mem_bytes_c_")));
void GET_TOT_GLOBAL_MEM_BYTES_C__() __attribute__((weak,alias("get_tot_global_mem_bytes_c_")));
/*getting the driver version  of the devices*/
int get_cuda_driverVersion_c_(int *driverVersion, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_cuda_driverVersion(driverVersion, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_CUDA_DRIVERVERSION_C_() __attribute__((weak,alias("get_cuda_driverVersion_c_")));
void get_cuda_driverversion_c__() __attribute__((weak,alias("get_cuda_driverVersion_c_")));
void GET_CUDA_DRIVERVERSION_C__() __attribute__((weak,alias("get_cuda_driverVersion_c_")));
/*getting the runtime version fo the devices*/
int get_cuda_runtimeVersion_c_(int *runtimeVersion, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_cuda_runtimeVersion(runtimeVersion, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_CUDA_RUNTIMEVERSION_C_() __attribute__((weak,alias("get_cuda_runtimeVersion_c_")));
void get_cuda_runtimeversion_c__() __attribute__((weak,alias("get_cuda_runtimeVersion_c_")));
void GET_CUDA_RUNTIMEVERSION_C__() __attribute__((weak,alias("get_cuda_runtimeVersion_c_")));
/*getting the name of the devices*/
int get_dev_Name_c_(int *idev, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = get_dev_Name(idev, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_DEV_NAME_C_() __attribute__((weak,alias("get_dev_Name_c_")));
void get_dev_name_c__() __attribute__((weak,alias("get_dev_Name_c_")));
void GET_DEV_NAME_C__() __attribute__((weak,alias("get_dev_Name_c_")));
/*device init*/
int devicedetails_init_c_(int *idev, deviceDetails_t *devD) {
  int rc = RC_SUCCESS;
  rc = deviceDetails_init(idev, devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void DEVICEDETAILS_INIT_C_() __attribute__((weak,alias("devicedetails_init_c_")));
void devicedetails_init_c__() __attribute__((weak,alias("devicedetails_init_c_")));
void DEVICEDETAILS_INIT_C__() __attribute__((weak,alias("devicedetails_init_c_")));
/*getting the device count*/
int get_dev_count_c_(deviceDetails_t *devD) {
  int rc = RC_SUCCESS; /* return code */
  int devCnt = 0;
  rc = get_dev_count(devD);
  if (rc == RC_FAIL) { rc = get_error_c(); }
  return rc;
}
void GET_DEV_COUNT_C_() __attribute__((weak,alias("get_dev_count_c_")));
void get_dev_count_c__() __attribute__((weak,alias("get_dev_count_c_")));
void GET_DEV_COUNT_C__() __attribute__((weak,alias("get_dev_count_c_")));
/*getting the error message and return code */
int get_error_c() {
  int rc = RC_SUCCESS;
  printf("rc  = %i, at Line %i in file %s %s\n",
	 rc,__LINE__,__FILE__,__FUNCTION__);
  printf("Result = FAIL\n");
  rc = RC_FAIL;
  exit(EXIT_FAILURE);
  return rc;
}

#endif /* CUDA */
