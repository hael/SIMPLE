!=Module simple_timer_omp
!
!
! Michael Eager 2017-03-15

module simple_timer_omp
  use precision_m
  use omp_lib
  implicit none

  real(dp),save :: last_time_point=REAL(0.,dp) !< Starting timesamp
  public :: ticOMP,tickrateOMP
  public :: tocOMP,tdiffOMP,tocprintOMP
  public :: nowOMP,reset_timerOMP
  public :: last_time_point
contains
!< Force timestamps and clock rate to zero
  subroutine reset_timerOMP()
    last_time_point=REAL(0.,dp)
  end subroutine reset_timerOMP
!< Get the clock tick count per second
  real(dp) function tickrateOMP()
    integer(dp)  :: sysclock
#ifdef _DEBUG
    write (*,'(A)') " OMP timer doesn't have a CLOCK_RATE "
#endif
    call system_clock(count_rate=sysclock)
    tickrateOMP=REAL(sysclock,dp)
  end function tickrateOMP
!< Get system_clock timestamp
  real(dp) function ticOMP()
#ifdef OPENMP
    ticOMP=OMP_get_wtime()
#else
    call cpu_time(ticOMP)
#endif
    last_time_point=ticOMP
  end function ticOMP

!< Calculate the time from two timestamps
  real(dp) function tdiffOMP(tend,tstart)
    real(dp),intent(in) :: tend
    real(dp),intent(in),optional :: tstart
    if (present(tstart)) last_time_point=tstart
! Calulate the time difference
    tdiffOMP=REAL(tend-last_time_point,dp)
  end function tdiffOMP

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point
  real(dp) function tocOMP(start_optional)
    real(dp),intent(in),optional ::  start_optional
    real(dp) :: end_point=REAL(0.0,dp)
    if (present(start_optional)) last_time_point=start_optional
#ifdef OPENMP
    end_point=OMP_get_wtime()
#else
    call cpu_time(end_point)
#endif
    tocOMP=tdiffOMP(end_point,last_time_point)
    last_time_point=end_point
  end function tocOMP

!< Complete the timing regime using a reference timestamp or the one
!  in last_time_point then print the results
  real(dp) function tocprintOMP(start_optional)
    real(dp),intent(in),optional ::  start_optional
    real(dp) :: end_point=REAL(0.0,dp)
    if (present(start_optional)) last_time_point=start_optional
#ifdef OPENMP
    end_point=OMP_get_wtime()
#else
    call cpu_time(end_point)
#endif
    tocprintOMP=tdiffOMP(end_point,last_time_point)
    last_time_point=end_point
    write (*,'(A,1d20.10)') " Elapsed time ",tocprintOMP
  end function tocprintOMP

!> print current time and date
  subroutine nowOMP()
    character(len=8)  :: date
    character(len=10) :: time
    print*,"OpenMP time: ",ticOMP()
!***********************************************************************************
    call date_and_time(date,time)
    write (*,'(A,A,A,A,A,A,A)') 'Date: ',date(7:8),'-',date(5:6),'-',date(1:4),'\n'
    write (*,'(A,A,A,A,A,A,A)') 'Time: ',time(1:2),':',time(3:4),':',time(5:10),'\n'
  endsubroutine nowOMP
end module simple_timer_omp

