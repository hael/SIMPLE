! ============================================================================
! Name        : simple_deepLearning_defs
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 02nd of May 2016
! Description : Module to handle cuDNN librray within the cuda environemt
!               provides basic definitions for the types and declarations
!               used in gpu calculations in modules using cuda calls. Using
!               cuda-5.0+
! Path        : ./defs
! ============================================================================
!
module simple_deepLearning_defs
  use, intrinsic :: iso_c_binding
  
  use simple_defs
  use simple_cuda_defs

  implicit none
!******
!DESCRIPTION
!defining the types of the for the GPU environment
#define devptr_t integer*8
#define size_t integer*8
#if defined (OPENCL)
#define cl_mem integer*8
#endif
  ! Maximum supported number of tensor dimensions */
#define CUDNN_DIM_MAX 8

  ! CUDNN return codes
  integer, parameter :: CUDNN_STATUS_SUCCESS          = 0
  integer, parameter :: CUDNN_STATUS_NOT_INITIALIZED  = 1
  integer, parameter :: CUDNN_STATUS_ALLOC_FAILED     = 2
  integer, parameter :: CUDNN_STATUS_BAD_PARAM        = 3
  integer, parameter :: CUDNN_STATUS_INTERNAL_ERROR   = 4
  integer, parameter :: CUDNN_STATUS_INVALID_VALUE    = 5
  integer, parameter :: CUDNN_STATUS_ARCH_MISMATCH    = 6
  integer, parameter :: CUDNN_STATUS_MAPPING_ERROR    = 7
  integer, parameter :: CUDNN_STATUS_EXECUTION_FAILED = 8
  integer, parameter :: CUDNN_STATUS_NOT_SUPPORTED    = 9
  integer, parameter :: CUDNN_STATUS_LICENSE_ERROR    = 10

  !CUDNN data type
  integer, parameter :: CUDNN_DATA_FLOAT  = 0
  integer, parameter :: CUDNN_DATA_DOUBLE = 1
  integer, parameter :: CUDNN_DATA_HALF   = 2

  !CUDNN propagate Nan
  integer, parameter :: CUDNN_NOT_PROPAGATE_NAN  = 0
  integer, parameter :: CUDNN_PROPAGATE_NAN      = 1

  ! CUDNN data structure for return codes
  type, bind(c) :: cudnnStatus
     integer(c_int) :: CUDNN_STATUS_SUCCESS
     integer(c_int) :: CUDNN_STATUS_NOT_INITIALIZED
     integer(c_int) :: CUDNN_STATUS_ALLOC_FAILED
     integer(c_int) :: CUDNN_STATUS_BAD_PARAM
     integer(c_int) :: CUDNN_STATUS_INTERNAL_ERROR
     integer(c_int) :: CUDNN_STATUS_INVALID_VALUE
     integer(c_int) :: CUDNN_STATUS_ARCH_MISMATCH
     integer(c_int) :: CUDNN_STATUS_MAPPING_ERROR
     integer(c_int) :: CUDNN_STATUS_EXECUTION_FAILED
     integer(c_int) :: CUDNN_STATUS_NOT_SUPPORTED
     integer(c_int) :: CUDNN_STATUS_LICENSE_ERROR
  end type cudnnStatus

  ! CUDNN data type
  type, bind(c) :: cudnnDataType
     integer(c_int) :: CUDNN_DATA_FLOAT
     integer(c_int) :: CUDNN_DATA_DOUBLE
     integer(c_int) :: CUDNN_DATA_HALF
  end type cudnnDataType
  ! CUDNN propagate Nan
  type, bind(c) :: cudnnNanPropagation
     integer(c_int) :: CUDNN_NOT_PROPAGATE_NAN
     integer(c_int) :: CUDNN_PROPAGATE_NAN
  end type cudnnNanPropagation

#if defined (CUDA)
  !here set the external routines for the cuDNN environment
  integer, external :: simple_cudnn_Create   !(cudnnHandle_t *handle);
  integer, external :: simple_cudnn_Destroy  !(cudnnHandle_t handle);
  integer, external :: simple_cudnn_SetStream!(cudnnHandle_t handle,cudaStream_t streamId);
  integer, external :: simple_cudnn_GetStream!(cudnnHandle_t handle,cudaStream_t *streamId);

#endif

contains
! TODO: insert the error handlers
! Error message from the cuDNN calls methods
!

end module simple_deepLearning_defs
