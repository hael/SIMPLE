! ============================================================================
! Name        : simple_deepLearning
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 02nd of May 2016
! Description : Module to handle cuDNN librray within the cuda environemt
!               initialisation/shutdown and error
! Path        : ./src/simple_gpu/
! ============================================================================
!
module simple_deepLearning
  use simple_defs
  use simple_cuda_defs
  use simple_deepLearning_defs
  use simple_deviceQuery_gpu
  use simple_systemQuery_cpu
  
  implicit none
#define devptr_t integer*8
  type(deviceDetails),allocatable :: a_devD(:)
  type(cudnnStatus) :: s_cudnnStatus
  private

  public :: simple_dpl_init
  public :: simple_dpl_shutdown
  
contains
!*******************************************************************************
! DESCRIPTION
! starts cuDNN library on gpu, compulsory for any computation using the library
!
!*******************************************************************************
! SOURCE
  subroutine simple_dpl_init(dpl_ierr,devQ,devD)
    !use simple_cuda_defs
    implicit none
    !argument variable
    integer, intent(out), optional  :: dpl_ierr
    type(deviceQuery_gpu),optional  :: devQ
    type(deviceDetails),optional    :: devD
    !return code
    integer                         :: rc !return code
    !start of the execution commands

#if defined (CUDA)
    dpl_ierr = simple_cudnn_Create(s_cudnnStatus)
#else
    if ( use_gpu .eqv. .true. ) then
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA to acces the CUDA environment   "
    write(*,*)"computation using GPU in routine: simple_dpl_init               "
    write(*,*)"****************************************************************"
    end if
#endif 

  return
end subroutine simple_dpl_init
!*******************************************************************************
! DESCRIPTION
! shuts down cuDNN library on gpu, compulsory for any computation using the
! library
!
!*******************************************************************************
! SOURCE
subroutine simple_dpl_shutdown(dpl_ierr)
  implicit none
    integer, intent(out), optional  :: dpl_ierr

#if defined (CUDA)
    dpl_ierr = simple_cudnn_Destroy(s_cudnnStatus)
#else
    if ( use_gpu .eqv. .true. ) then
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA to acces the CUDA environment   "
    write(*,*)"computation using GPU in routine: simple_dpl_init               "
    write(*,*)"****************************************************************"
    end if
#endif
  return
end subroutine simple_dpl_shutdown
!*******************************************************************************
! DESCRIPTION
! shuts down cuDNN library on gpu, compulsory for any computation using the
! library
!
!*******************************************************************************
! SOURCE
  subroutine simple_cudnn_stat_return(dpl_ierr)
    implicit none

    !global varaiables
    integer :: dpl_ierr

    select case(dpl_ierr)
    case(CUDNN_STATUS_SUCCESS)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_SUCCESS'
    case(CUDNN_STATUS_NOT_INITIALIZED)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_NOT_INITIALIZED'
    case(CUDNN_STATUS_ALLOC_FAILED)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_ALLOC_FAILED'
    case(CUDNN_STATUS_BAD_PARAM)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_BAD_PARAM'
    case(CUDNN_STATUS_INTERNAL_ERROR)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_INTERNAL_ERROR'
    case(CUDNN_STATUS_INVALID_VALUE)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_INVALID_VALUE'
    case(CUDNN_STATUS_ARCH_MISMATCH)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_ARCH_MISMATCH'
    case(CUDNN_STATUS_MAPPING_ERROR)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_MAPPING_ERROR'
    case(CUDNN_STATUS_EXECUTION_FAILED)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_EXECUTION_FAILED'
    case(CUDNN_STATUS_NOT_SUPPORTED)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_NOT_SUPPORTED'
    case(CUDNN_STATUS_LICENSE_ERROR)
       write(*,*) 'Error=',dpl_ierr,': CUDNN_STATUS_LICENSE_ERROR'
    end select

    return
  end subroutine simple_cudnn_stat_return

end module simple_deepLearning
