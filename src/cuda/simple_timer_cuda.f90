!> Simple CUDA timer module
!
!! TODO Working for pgi fortran but not gfortran/cuda8
!!
!! Michael Eager 2017-03-15
module simple_timer_cuda
include 'simple_lib.f08'
use cudafor
implicit none
    !> 
    type timer_cuda
#ifdef PGI
        type(cudaEvent) :: start_point
        type(cudaEvent) :: end_point
#else
        real(dp)  :: start_point=REAL(0.,dp)
        real(dp)  :: end_point=REAL(0.,dp)
#endif
        
    contains
        final :: CUtimer_kill
        procedure :: ticU
        procedure :: tocU
        procedure :: nowCU
        procedure :: CUtimer_setup
    end type timer_cuda

    interface timer_cuda
        module procedure constructor
    end interface timer_cuda

contains
    function constructor() result(this)
        class(timer_cuda),pointer :: this
        call this%CUtimer_setup()
    end function constructor

    subroutine CUtimer_setup(this)
#ifdef PGI
        use cudafor
        class(timer_cuda) :: this
        type (cudaDeviceProp) :: prop
        integer :: istat
#if defined(_DEBUG) || defined (DEBUG)
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
        if (istat .ne. 0) call simple_stop('  GetDeviceProperties for device 0: Failed')
        ilen=verify(prop%name,' ',.true.)
        write(*,"('  Device Name: ',a)") trim(prop%name)
        write(*,"('  Compute Capability: ',i0,'.',i0)") &
            prop%major, prop%minor
        write(*,"('  Number of Multiprocessors: ',i0)") &
            prop%multiProcessorCount
        write(*,"('  Max Threads per Multiprocessor: ',i0)") &
            prop%maxThreadsPerMultiprocessor
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
#endif
        !    call timer_kill
        !    if (this%start_point) this%istat = cudaEventDestroy(this%start_point)
        !    if (this%end_point ) this%istat = cudaEventDestroy(this%end_point)

        istat=cudaEventCreate(this%start_point)
        if (istat .ne. 0) call simple_stop(' cudaEventCreate for device 0: Failed')
        istat=cudaEventCreate(this%end_point)
        if (istat .ne. 0) call simple_stop(' cudaEventCreate for device 0: Failed')
#else
        this%start_point=0.0d0
        this%end_point=this%start_point
#endif
    end subroutine CUtimer_setup


    function ticU(this) result(res)
        class(timer_cuda) :: this
#ifdef PGI
        type(cudaEvent)  :: res
        integer :: istat
        istat=cudaEventRecord(this%start_point,0)
        if (istat .ne. 0) call simple_stop(' cudaEventRecord for device 0: Failed')
#else
        real(dp) :: res
        call cpu_time(this%start_point)
#endif
        ticU=this%start_point
    end function ticU

    real function tocU(this,start_optional)
        class(timer_cuda) :: this
#ifdef PGI
        integer :: istat
        type(cudaEvent),intent(in),optional ::  start_optional
#else
        real(dp),intent(in),optional :: start_optional
#endif
        real :: elapsed
        if (present(start_optional)) this%start_point=start_optional

#ifdef PGI
        istat=cudaEventRecord(this%end_point,0)
        if (istat .ne. 0) call simple_stop(' cudaEventRecord for device 0: Failed')
        istat=cudaEventSynchronize(this%end_point)
        if (istat .ne. 0) call simple_stop(' cudaEventtSynchronize for device 0: Failed')
        istat=cudaEventElapsedTime(elapsed,this%start_point,this%end_point)
        if (istat .ne. 0) call simple_stop(' cudaEventElapsedTime for device 0: Failed')
#else
        call cpu_time(this%end_point)
        elapsed=this%end_point-this%start_point
#endif
        tocU=elapsed
    end function tocU

    !< Print Info and clock time
    subroutine nowCU(this)
        class(timer_cuda) :: this
        integer :: ilen
        character(len=8)  :: date
        character(len=10) :: time
        character(len=33) :: f_result
        !***********************************************************************************
#ifdef PGI
        integer :: istat
        type(cudaDeviceProp) :: prop
        istat=cudaGetDeviceProperties(prop,0)
        if (istat .ne. 0) call simple_stop('  GetDeviceProperties for device 0: Failed')
        ilen=verify(prop%name,' ',.true.)
        write (*,*) prop%name(1:ilen), &
            real(prop%clockRate)/1000.0, &
            real(prop%totalGlobalMem)/1024.0/1024.0
#endif
        call date_and_time(date,time)
        write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4),'\n'
        write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10),'\n'
    end subroutine nowCU

    !> \brief timer destructor
    subroutine CUtimer_kill(this)
#ifdef PGI
        type(timer_cuda) :: this
        integer :: istat
        istat=cudaEventDestroy(this%start_point)
        if (istat .ne. 0) call simple_stop(' cudaEventDestroyfor device 0: Failed')
        istat=cudaEventDestroy(this%end_point)
        if (istat .ne. 0) call simple_stop(' cudaEventDestroyfor device 0: Failed')
#endif
    end subroutine CUtimer_kill

end module simple_timer_cuda

