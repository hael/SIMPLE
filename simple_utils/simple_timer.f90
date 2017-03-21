!= Module simple_timer
!
!
! Michael Eager 2017-03-15

module simple_timer
  use precision_m
#ifdef _CUDA
  use cudafor
#endif
  implicit none


  public :: get_clock_rate,tic, toc, now, timer_setup, timer_kill
  private
  !type timer ! system or CUDA timer


  ! contains
  ! declare a static constructor        
  ! initial  :: timer_setup  ! this will be evoked before any
  ! calls to timer variables or
  ! functions
  !     final :: timer_destructor
  ! methods
  !    procedure :: kill
  !    procedure :: tic
  !    procedure :: toc
  !   procedure :: now
  !    procedure :: setup
  !     procedure :: assign
  !     procedure :: add
  !     procedure :: sub
  !    procedure :: writef
  !     generic :: assignment(=) => assign
  !     generic :: operator(+) => add
  !     generic :: operator(-) => sub
  !    generic :: write(formatted) => writef

  !end type timer

  ! interface timer
  !       module procedure  constructor
  !       procedure :: constructFromArray
  !  end interface timer

#ifdef CUDA
  integer :: istat
  type (cudaEvent) :: &
#else
       integer(fp_kind) :: &
#endif
       clock_start, clock_end, clock_rate
  real(fp_kind) :: elapsed

contains

  !  function constructor() result(self)
  !    call setup
  !  end function constructor


  !< Timer constructor
  subroutine timer_setup()
    call timer_kill
#ifdef CUDA
    integer :: nDevices
    istat = cudaGetDeviceCount(nDevices)
    write(*,"('Number of CUDA-capable devices: ', i0,/)") &
         nDevices 
    istat = cudaEventCreate(clock_start)
    istat = cudaEventCreate(clock_end)
    ! output device info and transfer size
    !  istat = cudaGetDeviceProperties(prop, 0)

#else
    intrinsic system_clock
#ifdef GF_CLOCK_MONOTONIC
    ! 2010 fortran intrinsic module ISO FORTRAN ENV  
    ! build with -lrt to ensure high precision clock is used
    call system_clock_8(count, count_rate, count_max)
#else
    call system_clock(count_rate=clock_rate) ! Find the rate
#endif

#endif
  end subroutine timer_setup

  !< Get the clock_rate variable, make sure the timer is setup 
  function get_clock_rate() result(rate)
    real(fp_kind) :: rate
    call timer_setup
    rate = REAL(clock_rate,fp_kind)/1e9
  end function get_clock_rate

    function tic() result(itic)
#ifdef CUDA
      integer(fp_kind) :: itic
      istat = cudaEventRecord(clock_start,0) ! start event 
#else
      integer(fp_kind) :: itic
      intrinsic system_clock
#ifdef GF_CLOCK_MONOTONIC
      call system_clock_8(count=clock_start) ! Start timing
#else
      call system_clock(count=clock_start) ! Start timing
#endif
      itic=clock_start !REAL(clock_start,fp_kind)
#endif

    end function tic

    function toc (arg) result (elapsed_time)
      real(fp_kind) :: elapsed_time
#ifdef CUDA
      type (cudaEvent), optional, intent(in) :: arg
      if( present(arg)) then
         clock_start = arg
      end if
      istat = cudaEventRecord(clock_end,0)
      istat = cudaEventSynchronize(clock_end)
      istat = cudaEventElapsedTime(elapsed, &
           clock_start, clock_end)
#else
      real(fp_kind) :: itoc
      integer(8), optional, intent(in) :: arg
      intrinsic system_clock
      if( present(arg)) then
         clock_start = arg
      end if
#ifdef GF_CLOCK_MONOTONIC
      call system_clock_8(count=clock_end) ! End timing
#else
      call system_clock(count=clock_end) ! End timing
#endif
      itoc = REAL(clock_end,fp_kind)
      elapsed=REAL(clock_end-clock_start,fp_kind)/REAL(clock_rate,fp_kind)
#endif
      elapsed_time=elapsed
    end function toc

    !< Timer destructor
    subroutine timer_kill () !timer_destructor(this)
#ifdef _CUDA
      ! nothing to do here
      istat = cudaEventDestroy(clock_start)
      istat = cudaEventDestroy(clock_end)
#else
      ! nothing much to do here

      clock_end=0.
      clock_start=0.
#endif

    end subroutine timer_kill !timer_destructor

    subroutine now()

    end subroutine now

  end module simple_timer



