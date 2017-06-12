!= Module simple_timer_cuda
!
! TODO Working for pgi fortran but not gfortran/cuda8
! Michael Eager 2017-03-15

module simple_timer_cuda
  use simple_defs
  use cudafor

  implicit none

  type timer_cuda
#ifdef PGI
    type(cudaEvent) :: start_point=>null()
    type(cudaEvent) :: end_point=>null()
#else
    real(dp)  :: start_point=REAL(0.,dp)
    real(dp)  :: end_point=REAL(0.,dp)
#endif
    integer :: istat
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
    class(timer_cuda) :: this
#ifdef _DEBUG
    integer :: nDevices
    this%istat=cudaGetDeviceCount(nDevices)
    write (*,"('Number of CUDA-capable devices: ', i0,/)") nDevices
    ! output device info and transfer size
    this%istat=cudaGetDeviceProperties(prop,0)
    ilen=verify(prop%name,' ',.true.)
    write (*,900) prop%name(1:ilen), &
      real(prop%clockRate)/1000.0, &
      real(prop%totalGlobalMem)/1024.0/1024.0
#endif
!    call timer_kill
!    if (this%start_point) this%istat = cudaEventDestroy(this%start_point)
    !    if (this%end_point ) this%istat = cudaEventDestroy(this%end_point)

    this%istat=cudaEventCreate(this%start_point)
    this%istat=cudaEventCreate(this%end_point)
#else
    this%start_point=0.0d0
    this%end_point=this%start_point
#endif
  end subroutine CUtimer_setup

#ifdef PGI
  type(cudaEvent) function ticU(this)
    class(timer_cuda) :: this
    this%istat=cudaEventRecord(this%start_point,0)
#else
    real(dp) function ticU(this)
      call cpu_time(this%start_point)
#endif
      ticU=this%start_point
    end function ticU

    real(dp) function tocU(this,start_optional)
#ifdef PGI
      class(timer_cuda) :: this
      type(cudaEvent),intent(in),optional ::  start_optional
#else
      real(dp),intent(in),optional :: start_optional
#endif
      real(fp_kind) :: elapsed
      if (present(start_optional)) this%start_point=start_optional

#ifdef PGI
      this%istat=cudaEventRecord(this%end_point,0)
      this%istat=cudaEventSynchronize(this%end_point)
      this%istat=cudaEventElapsedTime(elapsed,this%start_point,this%end_point)
#else
      call cpu_time(this%end_point)
      elapsed=this%end_point-this%start_point
#endif
      tocU=elapsed
    end function tocU

    !< Print Info and clock time
    subroutine nowCU(this)
      !integer :: ilen
      character(len=8)  :: date
      character(len=10) :: time
      character(len=33) :: f_result
      !***********************************************************************************
#ifdef PGI
      class(timer_cuda) :: this
      type(cudaDeviceProp) :: prop
      this%istat=cudaGetDeviceProperties(prop,0)
      ilen=verify(prop%name,' ',.true.)
      write (*,900) prop%name(1:ilen), &
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
      this5istat=cudaEventDestroy(this%start_point)
      this%istat=cudaEventDestroy(this%end_point)
#endif
    end subroutine CUtimer_kill

  end module simple_timer_cuda

