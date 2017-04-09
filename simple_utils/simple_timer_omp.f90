!=Module simple_timer_omp
!
!
! Michael Eager 2017-03-15

module simple_timer_omp
   use precision_m
   use omp_lib
   implicit none

   real(dp), save :: start_point = REAL(0., dp) !< Starting timesamp
   public ::  ticOMP, tocOMP, tdiffOMP, nowOMP, reset_timerOMP
   public :: start_point
contains
   !< Force timestamps and clock rate to zero
   subroutine reset_timerOMP()
      start_point = REAL(0., dp)
   end subroutine reset_timerOMP

   !< Get system_clock timestamp
   real(dp) function ticOMP()
#ifdef OPENMP
     ticOMP = OMP_get_wtime()
#else
     call cpu_time(ticOMP)
#endif
#ifdef DEBUG
      write (*, '(A,1d20.10)') " TIC Time stamp", REAL(ticOMP, dp)
#endif
      start_point = ticOMP
   end function ticOMP

   !< Calculate the time from two timestamps
   real(dp) function tdiffOMP(tend, tstart)
      real(dp), intent(in) :: tend
      real(dp), intent(in), optional :: tstart
      if (present(tstart)) start_point = tstart
      ! Calulate the time difference
      tdiffOMP = REAL(tend - start_point, dp)
   end function tdiffOMP

   !< Complete the timing regime using a reference timestamp or the one
   !  in start_point
   real(dp) function tocOMP(start_optional)
      real(dp), intent(in), optional ::  start_optional
      real(dp) :: end_point = REAL(0.0,dp)
      if (present(start_optional)) start_point = start_optional
#ifdef OPENMP
      end_point = OMP_get_wtime()
#else
      call cpu_time(end_point)
#endif
      tocOMP = tdiffOMP(end_point, start_point)
      start_point = end_point
   end function tocOMP

   !> print current time and date
   subroutine nowOMP()
      character(len=8)  :: date
      character(len=10) :: time
      character(len=33) :: f_result
      !***********************************************************************************
      call date_and_time(date, time)
      write (*, '(A,A,A,A,A,A,A)') 'Date: ', date(7:8), '-', date(5:6), '-', date(1:4), '\n'
      write (*, '(A,A,A,A,A,A,A)') 'Time: ', time(1:2), ':', time(3:4), ':', time(5:10), '\n'
   end subroutine nowOMP
end module simple_timer_omp



