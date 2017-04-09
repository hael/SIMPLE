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
   logical :: debug = .false.
   logical :: warn = .false.
   integer(dp), save :: clock_rate = INT(0, dp) !< Number of counts per second
   integer(dp), save :: last_time_point = INT(0, dp) !< State timesamp
   integer(dp), save :: end_point = INT(0, dp) !< End timestamp

   public ::  tic, toc, tdiff, now, get_clock_rate, reset_timer, tocprint

   public :: last_time_point, end_point, clock_rate
contains

   !< Force timestamps and clock rate to zero
   subroutine reset_timer()
      last_time_point = INT(0, dp)
      end_point = INT(0, dp)
      clock_rate = INT(0, dp)
   end subroutine reset_timer

   !< Get system_clock timestamp
   integer(dp) function tic() ! result (timestamp)
      call system_clock(count=tic)
#ifdef DEBUG
      write (*, '(A,1d20.10)') " TIC Time stamp", REAL(tic, dp)
#endif
      last_time_point = tic
      end_point = INT(0, dp)
   end function tic

   !< Get the clock tick count per second
   real(dp) function get_clock_rate()
      if (clock_rate .eq. 0) call system_clock(count_rate=clock_rate)
      get_clock_rate = REAL(clock_rate, dp)
#ifdef DEBUG
      write (*, '(A,1d20.10)') " CLOCK_RATE (ticks/sec) ", REAL(clock_rate, dp)
#endif
   end function get_clock_rate

   !< Calculate the time from two timestamps
   real(dp) function tdiff(tend, tstart)
      integer(dp), intent(in), optional :: tend
      integer(dp), intent(in), optional :: tstart
      if (present(tend)) end_point = tend
      if (present(tstart)) last_time_point = tstart
      if (clock_rate .eq. 0) call system_clock(count_rate=clock_rate)
      ! Calulate the time difference
      tdiff = REAL(end_point - last_time_point, dp)/REAL(clock_rate, dp)
   end function tdiff

   !< Complete the timing regime using a reference timestamp or the one
   !  in last_time_point
   real(dp) function toc(start_optional)
      integer(dp), intent(in), optional ::  start_optional
      if (present(start_optional)) last_time_point = start_optional
      call system_clock(count=end_point)
      toc = tdiff(end_point, last_time_point)
      last_time_point = end_point
   end function toc

   !< Complete the timing regime using a reference timestamp or the one
   !  in last_time_point
   real(dp) function tocprint(start_optional)
      integer(dp), intent(in), optional ::  start_optional
      if (present(start_optional)) last_time_point = start_optional
      call system_clock(count=end_point)
#ifdef DEBUG
      write (*, '(A,1d20.10)') " TOC Time stamp ", REAL(end_point, dp)
#endif
      tocprint = tdiff(end_point, last_time_point)
      last_time_point = end_point
      write (*, '(A,1d20.10)') " Elapsed time ", tocprint
   end function tocprint

   !> print current time and date
   subroutine now()
      character(len=8)  :: date
      character(len=10) :: time
      character(len=33) :: f_result
      !***********************************************************************************
      call date_and_time(date, time)
      write (*, '(A,A,A,A,A,A,A)') 'Date: ', date(7:8), '-', date(5:6), '-', date(1:4), '\n'
      write (*, '(A,A,A,A,A,A,A)') 'Time: ', time(1:2), ':', time(3:4), ':', time(5:10), '\n'
   end subroutine now

end module simple_timer
