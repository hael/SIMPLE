
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

#define NREP_MAX INT(10000000,8)

#include "simple_timer.h"

module simple_timer_basic_test
  use simple_defs
  use simple_timer
  implicit none
  public:: exec_timertest
contains
  subroutine exec_timertest(be_verbose)
    logical,optional,intent(in)    :: be_verbose
    integer(dp),parameter :: nrep=10000000_dp
    real(dp)    :: xx,c,cfac,b
    real(dp)    :: etime,sysclockrate
    integer(dp) ::  t1,t2
    integer(dp) :: i

    c=.1
    cfac=.25
    b=1.
    xx=12.0_dp
    if (be_verbose) print*,'High resolution Fortran Timer'

    sysclockrate=REAL(tickrate(),dp)
    if (be_verbose) print*,'Clock rate (ticks per second): ',sysclockrate
    if (be_verbose) print*,"   "

    if (be_verbose) write (*,'(A)') "1. Simple timestamp and diff "
    t1=tic()
    if (be_verbose) write (*,'(A,1d20.10)') "    t1 = ",real(t1,dp)
    do i=1,nrep
      c=cfac*c+b
    end do
    t2=tic()
    if (be_verbose) write (*,'(A,1d20.10)') "    t2 = ",real(t2,dp)
    etime=tdiff(t2,t1)
    if (be_verbose) write (*,'(A,1d20.10)') '    Time for simple evaluation (s) = ',etime

    call reset_timer()
    if (be_verbose) write (*,"(A)") "  "
    if (be_verbose) write (*,'(A)') "2. Simple tic toc usage "
    if (be_verbose) write (*,'(A)') "    t1=TIC etime=TOC(T1)"

    t1=tic()
    if (be_verbose) write (*,'(A,1d20.10)') "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    etime=toc(t1)
    if (be_verbose) write (*,'(A,1d20.10)') '    toc(t1) = ',etime
    if (be_verbose) write (*,"(A)") "  "

    call reset_timer()
    if (be_verbose) write (*,'(A)') "3. Simple tic toc "
    t1=tic()
    if (be_verbose) write (*,'(A,1d20.10)') "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    etime=toc()
    if (be_verbose) write (*,'(A,1d20.10)') '    toc() = ',etime
    if (be_verbose) write (*,"(A)") ' '

    call reset_timer()
    if (be_verbose) write (*,'(A)') "4.  Testing return of toc in write cmd "
    t1=tic()
    if (be_verbose) write (*,'(A,1d20.10)') "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    if (be_verbose) write (*,'(A,1d20.10)') '    toc in write ',toc(t1)
    if (be_verbose) write (*,"(A)") "  "

    call reset_timer()
    if (be_verbose) write (*,'(A)') "5.  Testing empty tic() lhs "
    open (10,file=DEV_NULL)
    if (be_verbose) write (10,'(A,1i20)') "    Inline tic ",tic()
    ! if(be_verbose) write (*, '(A,1d20.10)') "2. t1 = ", real(t1, dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    if (be_verbose) write (*,'(A,1d20.10)') "    tic/toc in write ",toc()

    call reset_timer()
    if (be_verbose) write (*,"(A)") ' '
    if (be_verbose) write (*,'(A)') '6.  Testing tic toc in preprocessor macro  '
    TBLOCK()
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    TSTOP()
    if (be_verbose) write (*,"(A)") ' '

    call reset_timer()
    if (be_verbose) write (*,"(A)") ' '
    if (be_verbose) write (*,'(A)') '7.  Testing tic toc in preprocessor macro with subroutine  '
    TBLOCK()
    c=saxy(c)
    TSTOP()
    if (be_verbose) write (*,"(A)") ' '

    call reset_timer()
    if (be_verbose) write (*,"(A)") ' '
    if (be_verbose) write (*,'(A)') '8.  Testing timed_loop in subroutine '
    c=test_loop()

    call reset_timer()
    if (be_verbose) write (*,"(A)") ' '
    if (be_verbose) write (*,'(A)') '9.  Testing timed loop macro '
    START_TIMER_LOOP(10)
    c=.1
    c=saxy(c)
    STOP_TIMER_LOOP_( 'my test loop comment')

    call reset_timer()
    if (be_verbose) write (*,"(A)") ' '
    if (be_verbose) write (*,'(A)') '10.  Testing timed block macro '
    TIMER_BLOCK(
    c=.1;
    c=saxy(c)
    , ' my block comment ')

    end subroutine exec_timertest

  function saxy(c_in) result(c)
    real(dp),intent(in) :: c_in
    integer(dp),parameter :: nrep=NREP_MAX
    real(dp)              :: c,cfac,b
    integer(dp)           :: i
    b=1.
    c=c_in
    c=.1
    cfac=.25
    do i=1,nrep
      c=cfac*c+b
    end do
  end function saxy

  ! Loop macro has to declare an array and an int
  ! -- this means it has to be at the top or in
  !    its own function
  function test_loop() result(cloop)
    real(dp)                :: cloop
    integer(dp)             :: timestamp
    integer, parameter     :: N=100
    real(dp),dimension(N)  :: elapsed
    integer :: ii
    do ii=1,N
      timestamp=tic()
      cloop=REAL(.1,dp)
      cloop=saxy(cloop)
      elapsed(ii)=toc(timestamp)
    end do
    write (*,'(A)') '    Timed loop in subroutine *** '
    write (*,'(A,1d20.10)') "    Average over 3 (sec):", SUM(elapsed,DIM=1)/REAL(N,dp)
    write (*,'(A,1d20.10)') "    Min time (sec) ", MINVAL(elapsed,DIM=1)
    write (*,"(A)") ' '
  end function test_loop
end module simple_timer_basic_test
