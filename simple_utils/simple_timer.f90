!==Module simple_timer
!
!>\brief High resolution timer in fortran
!
!  64 bit INT implementation of system_clock
!  present in gfortran and pgfortran
!<------------------------------------
module simple_timer
!     use simple_jiffys ! singleton
!     use simple_defs   ! singleton
  use precision_m
  implicit none

!  private :: raise_sys_error
  private
  integer(dp),save   :: clock_ticks_per_second=INT(0,dp) !< Number of counts per second
  integer(dp),save   :: last_time_point=INT(0,dp) !< State timesamp
  integer(dp),save   :: end_point=INT(0,dp) !< End timestamp
  integer,save       :: idx_elapsed,num_elapsed=3
  logical,save       :: inloop=.false.
  real(dp),allocatable :: elapsed_times(:)
  public :: tic,tickrate
  public :: toc,tdiff,tocprint
  public :: now,reset_timer,timer_loop_start,in_timer_loop,timer_loop_end
  public :: last_time_point,clock_ticks_per_second
contains

  !< Force timestamps and clock rate to zero
  subroutine reset_timer()
    last_time_point=INT(0,dp)
    end_point=INT(0,dp)
    if (allocated(elapsed_times)) deallocate (elapsed_times)
  end subroutine reset_timer

  !< Get system_clock timestamp
  integer(dp) function tic()
    call system_clock(count=tic)
    last_time_point=tic
  end function tic

  !< Get the clock tick count per second
  integer(dp) function tickrate()
    tickrate=INT(0,dp)
    if (clock_ticks_per_second.eq.0) call system_clock(count_rate=tickrate)
    clock_ticks_per_second=tickrate
#ifdef _DEBUG
    write (*,'(A,1d20.10)') " CLOCK_RATE(ticks/sec) ",REAL(clock_ticks_per_second,dp)
#endif
  end function tickrate

  !< Calculate the time from two timestamps
  real(dp) function tdiff(tfinal,tstart)
    integer(dp),intent(in) :: tfinal
    integer(dp),intent(in) :: tstart
    ! integer(dp)                       ::  end_point
    ! if(present(tstart)) last_time_point = tstart
    ! if(.not. present(tfinal)) then
    !    call system_clock(count=end_point)
    !    tfinal = end_point
    ! end if
    if (clock_ticks_per_second.eq.0) call system_clock(count_rate=clock_ticks_per_second)
    ! Calulate the time difference
    tdiff=REAL(tfinal-tstart,dp)/REAL(clock_ticks_per_second,dp)
    !      last_time_point = tfinal
  end function tdiff

  !< Complete the timing regime using a reference timestamp or the one
  !  in last_time_point
  real(dp) function toc(start_optional)
    integer(dp),intent(in),optional ::  start_optional
    integer(dp)                       ::  end_point
    if (present(start_optional)) last_time_point=start_optional
    call system_clock(count=end_point)
    toc=tdiff(end_point,last_time_point)
    last_time_point=end_point
  end function toc

  !< Complete the timing regime using a reference timestamp or the one
  !  in last_time_point
  real(dp) function tocprint(start_optional)
    integer(dp),intent(in),optional ::  start_optional
    integer(dp)                       ::  end_point
    if (present(start_optional)) last_time_point=start_optional
    call system_clock(count=end_point)
#ifdef _DEBUG
    write (*,'(A,1d20.10)') " TOC Time stamp ",REAL(end_point,dp)
#endif
    tocprint=tdiff(end_point,last_time_point)
    last_time_point=end_point
    write (*,'(A,1d20.10)') " Elapsed time ",tocprint
  end function tocprint

  !> print current time and date
  subroutine now()
    character(len=8)  :: date
    character(len=10) :: time

    print*,"System_clock: ",tic()
    !***********************************************************************************
    call date_and_time(date,time)
    write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4),'\n'
    write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10),'\n'
  end subroutine now

  !< in_timer_loop checks the time within a timer loop
  ! It does not start the timer or set the inloop variable
  ! It returns false on the final loop- so that timer_loop_end can finish
  logical function in_timer_loop()
    in_timer_loop=.false.
    if (.not.inloop) then
      print*,"Failed timer_loop: Timer loop did not start"
    else
      if (idx_elapsed.lt.num_elapsed) then
        elapsed_times(idx_elapsed)=toc()
        idx_elapsed=idx_elapsed+1
        in_timer_loop=.true.
      end if
    end if
  end function in_timer_loop

  subroutine timer_loop_start(num)
    integer,intent(in),optional :: num
    integer(dp)  :: dummytimestamp
    num_elapsed=3
    if (present(num)) num_elapsed=num
    call reset_timer()
    if (num_elapsed.gt.1) then
      allocate (elapsed_times(num_elapsed))
      idx_elapsed=1
      dummytimestamp=tic()
      inloop=.true.
    end if
#ifdef _DEBUG
    print*,'Size of elapsed array ',size(elapsed_times)
#endif
  end subroutine timer_loop_start

  subroutine timer_loop_end()
    if (.not.inloop) then
      print*,"Failed timer_loop_end: Timer loop did not start"
    else
      write (*,'(A)') "******* TIMER LOOP **************"
      write (*,'(A,A,A,1i4)') '*** FILE:LINE: ',__FILE__,":",__LINE__
      write (*,'(A,1i8,A)') '*** Iterations: ',num_elapsed,' timed loops'

      if (idx_elapsed.eq.num_elapsed) then
        elapsed_times(idx_elapsed)=toc()
        write (*,'(A,1d20.10)') "*** Average (sec):",SUM(elapsed_times,DIM=1)/REAL(num_elapsed,dp)
        write (*,'(A,1d20.10,A,1i3)') "*** Longest run(sec) ",MAXVAL(elapsed_times,DIM=1), &
          '    at ',MAXLOC(elapsed_times,DIM=1)
        write (*,'(A,1d20.10,A,1i3)') "*** Shortest run(sec) ",MINVAL(elapsed_times,DIM=1), &
          '   at ',MINLOC(elapsed_times,DIM=1)
      else
        write (*,'(A,1i8)') '*** Failed at iteration ',idx_elapsed
      end if
      write (*,'(A)') "******* TIMER LOOP **************"
      inloop=.false.
      if (allocated(elapsed_times)) then
        deallocate (elapsed_times)
      end if
    end if
  end subroutine timer_loop_end
end module simple_timer
