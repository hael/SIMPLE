! ============================================================================
! Name        : testing_deviceQuery_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the functionality of the My_zgetri on cuda
!             :
! ============================================================================
!
program testing_deviceQuery_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use greeting_version
  use simple_timing
  use simple_deviceQuery_gpu

  implicit none
#define devptr_t integer*8

  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  integer                         :: err
  !local variables
  integer                         :: ndev

  !start of the execution commands
  !start of the greeting message
  call hello_deviceQuery_gpu(err)
  call timestamp()
  call start_Alltimers_cpu()

#if defined (CUDA)

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'

  call devQ%new_deviceQuery_gpu(devD)

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(devQ) and the data structure(devD) '
  write(*,*)'******************************************************************'
  ndev = devQ%get_devCnt()
  write(*,*)'cudaGetDeviceCount returned                    : ',ndev
  !writting the number of devices from data structure
  write(*,*)'cudaGetDeviceCount returned from data structure: ',devD%ndev
  if ( ndev /= devD%ndev) call devQ%get_warning_dataStruct_gpu()
  write(*,*)'******************************************************************'

  !shutting down the environment
  call simple_cuda_shutdown()

#else
  write(*,*)"**************************WARNING******************************"
  write(*,*)"You need to compile with -DCUDA                                "
  write(*,*)"to acces the CUDA environment computation using GPU            "
  write(*,*)"switching back to the CPU version of corr function             "
  write(*,*)"***************************************************************"
#endif

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_deviceQuery_gpu()

end program testing_deviceQuery_gpu
