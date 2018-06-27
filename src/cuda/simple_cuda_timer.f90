!> Simple CUDA timer module
!
!! TODO Working for pgi fortran but not gfortran/cuda8
!!
!! Michael Eager 2017-03-15
module simple_timer_cuda
include 'simple_lib.f08'
use, intrinsic :: ISO_C_BINDING
use CUDA

implicit none
private
    !>
    type timer_cuda
        type(cudaEvent_t) :: start_point
        type(cudaEvent_t) :: end_point
        type (cudaStream_t) :: tStream
    contains
        procedure :: ticU
        procedure :: tocU
        procedure :: nowCU
        procedure :: setup
        final :: destroy_cuda_timer
    end type timer_cuda

    interface timer_cuda
        module procedure constructor
    end interface timer_cuda

    public :: timer_cuda

contains
    function constructor() result(this)
        type(timer_cuda) :: this
        call this%setup()
    end function constructor

    subroutine setup(this)
        class(timer_cuda), intent(inout) :: this
        type (cudaDeviceProp) :: prop
        integer (KIND(cudaSuccess))  :: istat
        integer :: nDevices=0, ilen
        istat=cudaGetDeviceCount(nDevices)
        if (nDevices == 0) then
            write(*,"(/,'No CUDA devices found',/)")
            stop
        else
            write (*,"('Number of CUDA-capable devices: ', i0,/)") nDevices
        end if
        ! output device info and transfer size
        istat=cudaGetDeviceProperties(prop,0)
        if (istat .ne. cudaSuccess) call simple_stop('  GetDeviceProperties for device 0: Failed')
       ! ilen=verify(prop%name,' ',.true.)
        write(*,"('  Device Name: ',a)") prop%name
        write(*,"('  Compute Capability: ',i0,'.',i0)") &
            prop%major, prop%minor
        write(*,"('  Number of Multiprocessors: ',i0)") &
            prop%multiProcessorCount
        write(*,"('  Clock Rate: ',i0)") &
            prop%clockRate
        write(*,"('  Global Memory (GB): ',f9.3,/)") &
            prop%totalGlobalMem/1024.0**3

        ! Execution Configuration

        write(*,"('  Execution Configuration Limits')")
        write(*,"('    Max Grid Dims: ',2(i0,' x '),i0)") &
            prop%maxGridSize
        write(*,"('    Max Block Dims: ',2(i0,' x '),i0)") &
            prop%maxThreadsDim
        write(*,"('    Max Threads per Block: ',i0,/)") &
            prop%maxThreadsPerBlock
        !    call timer_kill
        !    if (this%start_point) this%istat = cudaEventDestroy(this%start_point)
        !    if (this%end_point ) this%istat = cudaEventDestroy(this%end_point)
        istat= cudaStreamCreate(this%tStream)
        if (istat .ne. cudaSuccess) call simple_stop(' cuStreamCreate for device 0: Failed')

        istat=cudaEventCreate(this%start_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventCreate for device 0: Failed')
        istat=cudaEventCreate(this%end_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventCreate for device 0: Failed')
    end subroutine setup


    function ticU(this) result(res)
        class(timer_cuda), intent(inout) :: this
        type(cudaEvent_t)  :: res
        integer (KIND(cudaSuccess)) :: istat
        istat=cudaEventRecord(this%start_point, this%tStream)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventRecord for device 0: Failed')
        res=this%start_point
    end function ticU

    function tocU(this,start_optional) result(elapsed)
        class(timer_cuda), intent(inout) :: this
        integer (KIND(cudaSuccess)) :: istat
        type(cudaEvent_t),intent(in),optional ::  start_optional
        real :: elapsed
        if (present(start_optional)) this%start_point=start_optional

        istat=cudaEventRecord(this%end_point,this%tStream)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventRecord for CUDA device : Failed')
        istat=cudaEventSynchronize(this%end_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventtSynchronize for device 0: Failed')
        istat=cudaEventElapsedTime(elapsed,this%start_point,this%end_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventElapsedTime for device 0: Failed')
        elapsed=elapsed/100.
    end function tocU

    !< Print Info and clock time
    subroutine nowCU(this)
        class(timer_cuda), intent(inout) :: this
        integer :: ilen
        character(len=8)  :: date
        character(len=10) :: time
        character(len=33) :: f_result
        integer(KIND(cudaSuccess))  :: istat
        type(cudaDeviceProp) :: prop
        istat=cudaGetDeviceProperties(prop,0)
        if (istat .ne. cudaSuccess) call simple_stop('  GetDeviceProperties for device 0: Failed')
       ! ilen=verify(prop%name,' ',.true.)
        write (*,*) prop%name, &
            real(prop%clockRate)/1000.0, &
            real(prop%totalGlobalMem)/1024.0/1024.0
        call date_and_time(date,time)
        write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4),'\n'
        write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10),'\n'
    end subroutine nowCU

    !> \brief timer destructor
    subroutine destroy_cuda_timer(this)
        type(timer_cuda), intent(inout) :: this
        integer(KIND(cudaSuccess))  :: istat
        istat= cudaStreamDestroy(this%tStream)
        if (istat .ne. cudaSuccess) call simple_stop(' cuStreamCreate for device 0: Failed')
        istat=cudaEventDestroy(this%start_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventDestroyfor device 0: Failed')
        istat=cudaEventDestroy(this%end_point)
        if (istat .ne. cudaSuccess) call simple_stop(' cudaEventDestroyfor device 0: Failed')
    end subroutine destroy_cuda_timer

end module simple_timer_cuda
