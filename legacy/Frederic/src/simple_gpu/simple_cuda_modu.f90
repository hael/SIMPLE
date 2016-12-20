! ============================================================================
! Name        : simple_cuda
! Author      : Frederic Bonnet
! Version     : 2.0
! Date        : 05th of April 2016
! Description : Module to handle cuda environemt initialisation/shutdown
!             : and error
!             : 
! ============================================================================
!
module simple_cuda
  use simple_defs
  use simple_cuda_defs
  use simple_deviceQuery_gpu
  use simple_systemQuery_cpu

  implicit none
  type(deviceDetails),allocatable :: a_devD(:)
  
  private
  
  public :: simple_cuda_init
  public :: simple_cuda_shutdown
  public :: simple_cuda_get_a_devd
  public :: simple_cuda_deviceQuery_logic
  public :: Sanity_check_gpu
  
contains
!*******************************************************************************
! DESCRIPTION
! starts gpu initializing the cublas library
!
!*******************************************************************************
! SOURCE
subroutine simple_cuda_init(cuda_ierr,devQ,devD)
  implicit none
  !argument variable
  integer, intent(out), optional  :: cuda_ierr
  type(deviceQuery_gpu),optional  :: devQ
  type(deviceDetails),optional    :: devD
  !local variable
  !type(deviceDetails),allocatable :: a_devD(:)
  integer                         :: ndev
  integer                         :: rc !return code
  !indexers
  integer                         :: idev
  !counters
  integer                         :: counter
  !functions calls
  integer :: get_dev_count_c
  integer :: print_s_devD_struct_c
  !start of the execution commands

#if defined (CUDA)
  rc = cublas_init()
  if ( rc /= RC_SUCCESS ) then
     if ( present(cuda_ierr) ) then
        cuda_ierr = rc
        return
     end if
  end if
  if ( present(devQ) .and. present(devD) ) then
     rc = get_dev_count_c(devD)
     ndev = devD%ndev
     allocate(a_devD(0:devD%ndev-1))
     call simple_cuda_get_a_devd(ndev,devQ,devD,a_devD)
     call simple_cuda_deviceQuery_logic(counter,ndev,a_devD)
     if (counter /= 0 ) rc = print_s_devD_struct_c(a_devD)
     if (counter == 0 ) write(*,*) "No GPU devices found"
  else
     has_gpu = .false.
     write(*,*) "deviceQuery not executed assuming has_gpu = false"
  end if
#else
  if ( use_gpu .eqv. .true. ) then
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA to acces the CUDA environment   "
    write(*,*)"computation using GPU                                           "
    write(*,*)"****************************************************************"
  end if
#endif 

  if (present(cuda_ierr) ) cuda_ierr = 0
  
  return
end subroutine simple_cuda_init
!*******************************************************************************
! DESCRIPTION
! Flips the logic of the has_gpu
!
!*******************************************************************************
! SOURCE
subroutine simple_cuda_deviceQuery_logic(ncnt,ndev,a_devD)
  implicit none
  integer               :: ncnt  !usable devices
  integer               :: ndev  !devices on system
  type(deviceDetails)   :: a_devD(0:ndev-1)
  !local variables
  !indexers
  integer :: idev
  !counters
  integer :: counter
  !start of the execution commands
  counter = 0
  do idev = 0,ndev-1
     if (a_devD(idev)%is_SMsuitable == 1 ) then
        counter = counter + 1
        write(*,'(5x,i2,20x,i1,20x,i1)') idev,&
             a_devD(idev)%is_SMsuitable,&
             a_devD(idev)%is_ecc
        has_gpu = .true.
     end if
  end do
  ncnt = counter
  return
end subroutine simple_cuda_deviceQuery_logic
!*******************************************************************************
! DESCRIPTION
! Fills in the array of data (a_devD) from deviceQuery object (devQ)
!
!*******************************************************************************
! SOURCE
subroutine simple_cuda_get_a_devd(ndev,devQ,devD,a_devD)
  implicit none
  !argument variable
  type(deviceQuery_gpu) :: devQ
  type(deviceDetails)   :: devD
  type(deviceDetails)   :: a_devD(0:ndev-1)
  !local variables
  integer               :: ndev !devices on system
  integer               :: rc !return code
  !indexers
  integer               :: idev
 
  !functions calls
  integer :: get_dev_count_c
  
  write(*,*)'****************************************************************'
  write(*,*)' Device fills in the object(devQ) and the data structure(devD)  '
  write(*,*)'****************************************************************'
  do idev = 0, ndev-1
     call devQ%new_deviceQuery_gpu(devD,idev)
     !call Sanity_check_gpu(devQ, devD)
     !mapping the data structures into an array
     a_devD(idev) = devD
  end do
  write(*,*)'****************************************************************'
  write(*,'(a,3x,a,3x,a)') "GPU device","is suitable for SM>=3.5","is_ecc 0: no 1:yes"

  return
end subroutine simple_cuda_get_a_devd
!*******************************************************************************
! DESCRIPTION
! shut down gpu shuting down the cublas library
!
!*******************************************************************************
! SOURCE
subroutine simple_cuda_shutdown()
  implicit none
  !    type(deviceDetails),optional :: a_devD(*)
  !local variable
  integer :: rc !return code

#if defined (CUDA)
  rc = cublas_shutdown()
#else
  if ( use_gpu .eqv. .true. ) then
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA to acces the CUDA environment   "
    write(*,*)"computation using GPU                                           "
    write(*,*)"****************************************************************"
  end if
#endif 

  use_gpu = .false.
  has_gpu = .false.

  !freeing ressources on host for the data structure
  !Best to free the data structure at program termination
  !call simple_cuda_free_a_devD() !deallocate(a_devD)

  return
end subroutine simple_cuda_shutdown
!*******************************************************************************
! DESCRIPTION
! Deallocator of the a_devD array
!
!*******************************************************************************
! SOURCE
subroutine simple_cuda_free_a_devD()
  !argument variable
  if (allocated(a_devD)) deallocate(a_devD)
  return
end subroutine simple_cuda_free_a_devD
!*******************************************************************************
!    Subroutine to run sanity checks on the data structure passed GPU
!
!*******************************************************************************
!
subroutine Sanity_check_gpu(devQ, devD)
  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_deviceQuery_gpu
  implicit none
  
  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  !local variables
  integer                         :: ndev

  !start of the execution commands

  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(devQ) and the data structure(devD) '
  write(*,*)'******************************************************************'
  write(*,*)'CUDA driver Ver.              (devQ):',devQ%get_cuda_drv_version()
  write(*,*)'CUDA driver Ver.              (devD):',devD%d_ver
  if(devD%d_ver/=devQ%get_cuda_drv_version())call devQ%get_warning_dataStruct_gpu()
  ndev = devQ%get_devCnt()
  write(*,*)'CUDA driver run Ver.          (devD):',devD%d_runver
  write(*,*)'DeviceCount                   (devQ):',ndev
  write(*,*)'DeviceCount                   (devD):',devD%ndev
  if ( ndev /= devD%ndev) call devQ%get_warning_dataStruct_gpu()
  write(*,*)'Device name                   (devD): ',devQ%get_devname(0)
  write(*,*)'Total global Mem in MB        (devD):',devD%tot_global_mem_MB
  write(*,*)'Total global Mem in bytes     (devD):',devD%tot_global_mem_bytes
  write(*,*)'nMultiProc                    (devD):',devD%nMultiProc
  write(*,*)'ncc_per_mp                    (devD):',devD%ncc_per_mp
  write(*,*)'ncc                           (devD):',devD%ncc
  write(*,*)'is_SMsuitable                 (devD):',devD%is_SMsuitable
  write(*,*)'nregisters_per_blk            (devD):',devD%nregisters_per_blk
  write(*,*)'warpSze                       (devD):',devD%warpSze
  write(*,*)'maxthreads_per_mp             (devD):',devD%maxthreads_per_mp
  write(*,*)'maxthreads_per_blk            (devD):',devD%maxthreads_per_blk
  write(*,*)'is_ecc                        (devD):',devD%is_ecc
  write(*,*)'******************************************************************'
  
  return
end subroutine Sanity_check_gpu

end module simple_cuda
