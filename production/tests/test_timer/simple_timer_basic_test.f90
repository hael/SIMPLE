!------------------------------------------------------------------------------!
! SIMPLE , Elmlund & Elmlund Lab,     simplecryoem.com                         !
!------------------------------------------------------------------------------!
!> test basic routines for simple_timer
!!
!! Test the basic timing functions in the SIMPLE library.
!!
!! @author
!! Michael Eager 2017
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
! -----------------------------------------------------------------------------!
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

#define NREP_MAX INT(100000,dp)

#include "simple_timer.h"

module simple_timer_basic_test
  use simple_defs
  use simple_timer
  implicit none
  public:: exec_timertest
  private
#include "simple_local_flags.inc"
contains
  subroutine exec_timertest(be_verbose)
    logical,optional,intent(in)    :: be_verbose
    integer(dp),parameter :: nrep=NREP_MAX
    real(dp)    :: xx,c,cfac,b
    real(dp)    :: etime,sysclockrate
    integer(dp) ::  t1,t2
    integer(dp) :: i
    integer :: io_stat

    if(present(be_verbose)) verbose=be_verbose
    c=.1
    cfac=.25
    b=1.
    xx=12.0_dp
    VerbosePrint 'High resolution Fortran Timer'

    sysclockrate=REAL(tickrate(),dp)
    VerbosePrint 'Clock rate (ticks per second): ',sysclockrate
    VerbosePrint "   "

    VerbosePrint "1. Simple timestamp and diff "
    t1=tic()
    VerbosePrint "    t1 = ",real(t1,dp)
    do i=1,nrep
      c=cfac*c+b
    end do
    t2=tic()
    VerbosePrint "    t2 = ",real(t2,dp)
    etime=tdiff(t2,t1)
    VerbosePrint '    Time for simple evaluation (s) = ',etime

    call reset_timer()
    VerbosePrint "  "
    VerbosePrint "2. Simple tic toc usage "
    VerbosePrint "    t1=TIC etime=TOC(T1)"

    t1=tic()
    VerbosePrint "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    etime=toc(t1)
    VerbosePrint '    toc(t1) = ',etime
    VerbosePrint "  "

    call reset_timer()
    VerbosePrint "3. Simple tic toc "
    t1=tic()
    VerbosePrint "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    etime=toc()
    VerbosePrint '    toc() = ',etime
    VerbosePrint ' '

    call reset_timer()
    VerbosePrint "4.  Testing return of toc in write cmd "
    t1=tic()
    VerbosePrint "    t1 = ",real(t1,dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    VerbosePrint '    toc in write ',toc(t1)
    VerbosePrint "  "

    call reset_timer()
    VerbosePrint "5.  Testing empty tic() lhs "
    open (10,file=DEV_NULL, IOSTAT=io_stat)
    if (be_verbose) write (10,'(A,1i20)') "    Inline tic ",tic()
    ! if(be_verbose) write (*, '(A,1d20.10)') "2. t1 = ", real(t1, dp)
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    VerbosePrint "    tic/toc in write ",toc()

    call reset_timer()
    VerbosePrint ' '
    VerbosePrint '6.  Testing tic toc in preprocessor macro  '
    TBLOCK()
    c=.1
    do i=1,nrep
      c=cfac*c+b
    end do
    TSTOP()
    VerbosePrint ' '

    call reset_timer()
    VerbosePrint ' '
    VerbosePrint '7.  Testing tic toc in preprocessor macro with subroutine  '
    TBLOCK()
    c=saxy(c)
    TSTOP()
    VerbosePrint ' '

    call reset_timer()
    VerbosePrint ' '
    VerbosePrint '8.  Testing timed_loop in subroutine '
    c=test_loop()

#if ! defined(PGI)
    call reset_timer()
    VerbosePrint ' '
    VerbosePrint '9.  Testing timed loop macro '

    START_TIMER_LOOP(10)
    c=.1
    c=saxy(c)
    STOP_TIMER_LOOP_( 'my test loop comment')

    call reset_timer()
    VerbosePrint ' '
    VerbosePrint '10.  Testing timed block macro '
    DebugPrint 'Testing line num'
#if 0
    ! defined(GNU)
    TIMER_BLOCK(
    c=.1;
    c=saxy(c)
    , ' my multiline timer block')
#else
    TIMER_BLOCK(   c=.1;  c=saxy(c) , ' my block comment (Intel does not like multi line in macro)')
#endif

#endif
    end subroutine exec_timertest

  function saxy(c_in) result(c)
    real(dp),intent(in) :: c_in
    integer(dp),parameter :: nrep=NREP_MAX
    real(dp)              :: c,cfac,b
    integer(dp)           :: i
    b=1.0_dp
    c=c_in
    c=0.1_dp
    cfac=0.25_dp
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
    real(dp)  :: elapsed(N)
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

  ! from GCC test suite -- F2008
   subroutine test_precision_timer ( )
       integer(1)    :: count1, irate1, mymax1
       integer(2)    :: count2, irate2, mymax2
       integer(4)    :: count4, irate4, mymax4
       integer(8)    :: count8, irate8, mymax8
       real(4)       :: rrate4
#if defined(GNU) && __GNUC__ < 5
#define NOREAL4CLOCK
#endif
#ifdef NOREAL4CLOCK
       if (FC_COMPILER_VERSION(1) >= 5)then
       call system_clock(count=count4, count_rate=irate4, count_max=mymax4)
       if (count4.ne.-127.or.irate4.ne.0.or.mymax4.ne.0) call abort
       !call system_clock(count=count4, count_rate=rrate4, count_max=mymax1)
       !if (count4.ne.-127.or.rrate4.ne.0.0.or.mymax4.ne.0) call abort
       call system_clock(count=count2, count_rate=irate2, count_max=mymax2)
       if (count2.ne.-32767.or.irate2.ne.0.or.mymax2.ne.0) call abort
       !call system_clock(count=count4, count_rate=rrate4, count_max=mymax2)
       !if (count4.ne.-32767.or.rrate4.ne.0.0.or.mymax2.ne.0) call abort
       call system_clock(count=count4, count_rate=irate4, count_max=mymax4)
       if (irate4.ne.1000.or.mymax4.ne.huge(0_4)) call abort
       !call system_clock(count=count4, count_rate=rrate4, count_max=mymax4)
       !if (rrate4.ne.1000.0.or.mymax4.ne.huge(0_4)) call abort
       call system_clock(count=count8, count_rate=irate8, count_max=mymax8)
       if (irate8.ne.1000000000.or.mymax4.ne.huge(0_8)) call abort
    end if
#endif
  end subroutine test_precision_timer

end module simple_timer_basic_test
