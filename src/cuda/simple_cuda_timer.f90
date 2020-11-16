!> Simple CUDA timer module
!
!! TODO Working for pgi fortran but not gfortran/cuda8
!!
!! Michael Eager 2017-03-15
module simple_timer_cuda
include 'simple_lib.f08'
use, intrinsic :: ISO_C_BINDING
use CUDA
use simple_cuda
#include "simple_cuda_handle.inc"
implicit none
private
    !>
    type timer_cuda
        type(cudaEvent_t) :: start_point
        type(cudaEvent_t) :: end_point
        type (cudaStream_t) :: tStream
        logical :: exists = .false.
    contains
        procedure :: ticU
        procedure :: tocU
        procedure :: nowCU
        procedure :: setup
        procedure ::  kill_
        final :: destroy_cuda_timer
    end type timer_cuda

    interface timer_cuda
        module procedure constructor
    end interface timer_cuda

    public :: timer_cuda
#include "simple_local_flags.inc"
contains
    function constructor() result(this)
        type(timer_cuda) :: this
        call this%setup()
    end function constructor

    subroutine setup(this)
        class(timer_cuda), intent(inout) :: this
        type (cudaDeviceProp)            :: prop
        integer (KIND(cudaSuccess))      :: istat
        integer                          :: nDevices=0, ilen
        logical :: cuda_errors_found
        cuda_errors_found=.false.
        if(this%exists .eqv. .false.)then
            istat=cudaGetDeviceCount(nDevices)
            HANDLE_CUDAERROR( istat, cuda_errors_found, 'cudaGetDeviceCount failed')
            if (nDevices == 0 ) then
                write(*,"(/,'No CUDA devices found',/)")
                stop
            else
                write (*,"('Number of CUDA-capable devices: ', i0,/)") nDevices
            end if
            ! output device info and transfer size
            istat=cudaGetDeviceProperties(prop,0)
            HANDLE_CUDAERROR( istat, cuda_errors_found, 'cudaGetDeviceProperties failed')

            ! ilen=verify(prop%name,' ',.true.)
            ! call print_cuda_properties(prop)
            ! write(*,*) '  Device Name: ', prop%name
            ! write(*,"('  Compute Capability: ',i0,'.',i0)") &
            !     prop%major, prop%minor
            ! write(*,"('  Number of Multiprocessors: ',i0)") &
            !     prop%multiProcessorCount
            ! write(*,"('  Clock Rate: ',i0)") &
            !     prop%clockRate
            ! write(*,"('  Global Memory (GB): ',f9.3,/)") &
            !     prop%totalGlobalMem/1024.0**3


        !    call timer_kill
        !    if (this%start_point) this%istat = cudaEventDestroy(this%start_point)
        !    if (this%end_point ) this%istat = cudaEventDestroy(this%end_point)
            istat= cudaStreamCreate(this%tStream)
            HANDLE_CUDAERROR( istat, cuda_errors_found, 'cudaStreamCreate Failed')
            istat=cudaEventCreate(this%start_point)
            HANDLE_CUDAERROR( istat, cuda_errors_found, 'cudaEventCreate  Failed')
            istat=cudaEventCreate(this%end_point)
            HANDLE_CUDAERROR( istat, cuda_errors_found, 'cudaEventCreate Failed')
            this%exists = .true.
        end if
    end subroutine setup


    function ticU(this) result(res)
        class(timer_cuda), intent(inout) :: this
        type(cudaEvent_t)                :: res
        integer (KIND(cudaSuccess))      :: istat
        logical :: err
        if(this%exists)then

        istat=cudaEventRecord(this%start_point, this%tStream)
         HANDLE_CUDAERROR( istat, err, ' cudaEventRecord  Failed')
        res=this%start_point
        end if
    end function ticU

    function tocU(this,start_optional) result(elapsed)
        class(timer_cuda), intent(inout)        :: this
        type(cudaEvent_t), intent(in), optional ::  start_optional
        integer (KIND(cudaSuccess)) :: istat
        real                        :: elapsed
        logical :: err
        if (present(start_optional)) this%start_point=start_optional

        istat=cudaEventRecord(this%end_point,this%tStream)
       HANDLE_CUDAERROR( istat, err,' cudaEventRecord  Failed')
        istat=cudaEventSynchronize(this%end_point)
       HANDLE_CUDAERROR( istat, err,' cudaEventtSynchronize Failed')
        istat=cudaEventElapsedTime(elapsed,this%start_point,this%end_point)
       HANDLE_CUDAERROR( istat, err,' cudaEventElapsedTime  Failed')
        elapsed=elapsed/100.
    end function tocU

    !< Print Info and clock time
    subroutine nowCU(this)
        class(timer_cuda), intent(inout) :: this
        integer                    :: ilen
        character(len=8)           :: date
        character(len=10)          :: time
        character(len=33)          :: f_result
        integer(KIND(cudaSuccess)) :: istat
        type(cudaDeviceProp)       :: prop
        istat=cudaGetDeviceProperties(prop,0)
        if (istat .ne. cudaSuccess) THROW_HARD('GetDeviceProperties for device 0: Failed')
       ! ilen=verify(prop%name,' ',.true.)
        write (*,*) prop%name, &
            real(prop%clockRate)/1000.0, &
            real(prop%totalGlobalMem)/1024.0/1024.0
        call date_and_time(date,time)
        write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4)
        write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10)
    end subroutine nowCU

    subroutine kill_(this)
        class(timer_cuda) :: this
        call destroy_cuda_timer(this)
        this%exists = .false.
    end subroutine kill_


    !> \brief timer destructor
    subroutine destroy_cuda_timer(this)
        type(timer_cuda), intent(inout) :: this
        integer(KIND(cudaSuccess))      :: istat
        istat= cudaStreamDestroy(this%tStream)
        if (istat .ne. cudaSuccess) THROW_HARD('cudaStreamCreate for device 0: Failed')
        istat=cudaEventDestroy(this%start_point)
        if (istat .ne. cudaSuccess) THROW_HARD('cudaEventDestroyfor device 0: Failed')
        istat=cudaEventDestroy(this%end_point)
        if (istat .ne. cudaSuccess) THROW_HARD('cudaEventDestroyfor device 0: Failed')
        this%exists = .false.
    end subroutine destroy_cuda_timer

end module simple_timer_cuda
