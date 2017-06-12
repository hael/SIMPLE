
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

#define NREP_MAX INT(10000000,8)

#include "simple_timer.h"

module simple_timer_omp_test
   use simple_defs
   use simple_timer_omp
   implicit none

   public:: exec_OMPtimertest
contains

   subroutine exec_OMPtimertest(be_verbose)
      logical, optional, intent(in)    :: be_verbose
      integer(dp), parameter :: nrep = 10000000_dp
      real(dp)    :: xx, c, cfac, b
      real(dp)    :: etime, sysclockrate
      real(dp) ::  t1, t2
      integer(dp) :: i
      c = .1
      cfac = .25
      b = 1.
      xx = 12.0_dp
      if (be_verbose) print *, 'OpenMP Fortran Timer'
      if (be_verbose) print *, 'Note: in debug, OpenMP may not be present, timer defaults to cpu_time'
      if (be_verbose) write (*, '(A)') "1. Simple timestamp and diff "
      t1 = ticOMP()
      if (be_verbose) write (*, '(A,1d20.10)') "    t1 = ", real(t1, dp)
      do i = 1, nrep
         c = cfac*c + b
      end do
      t2 = ticOMP()
      if (be_verbose) write (*, '(A,1d20.10)') "    t2 = ", real(t2, dp)
      etime = tdiffOMP(t2, t1)
      if (be_verbose) write (*, '(A,1d20.10)') '    Time for simple evaluation (s) = ', etime

      call reset_timerOMP()
      if (be_verbose) write (*, "(A)") "  "
      if (be_verbose) write (*, '(A)') "2. Simple tic toc usage "
      if (be_verbose) write (*, '(A)') "    t1=TIC; etime=TOC(T1)"

      t1 = ticOMP()
      if (be_verbose) write (*, '(A,1d20.10)') "    t1 = ", real(t1, dp)
      c = .1
      do i = 1, nrep
         c = cfac*c + b
      end do
      etime = tocOMP(t1)
      if (be_verbose) write (*, '(A,1d20.10)') '    toc(t1) = ', etime
      if (be_verbose) write (*, "(A)") "  "

      call reset_timerOMP()
      if (be_verbose) write (*, '(A)') "3. Simple tic toc "
      t1 = ticOMP()
      if (be_verbose) write (*, '(A,1d20.10)') "    t1 = ", real(t1, dp)
      c = .1
      do i = 1, nrep
         c = cfac*c + b
      end do
      etime = tocOMP()
      if (be_verbose) write (*, '(A,1d20.10)') '    toc() = ', etime
      if (be_verbose) write (*, "(A)") ' '

      call reset_timerOMP()
      if (be_verbose) write (*, '(A)') "4.  Testing return of toc in write cmd "
      t1 = ticOMP()
      if (be_verbose) write (*, '(A,1d20.10)') "     t1 = ", real(t1, dp)
      c = .1
      do i = 1, nrep
         c = cfac*c + b
      end do
      if (be_verbose) write (*, '(A,1d20.10)') '    toc in write ', tocOMP(t1)
      if (be_verbose) write (*, "(A)") "  "

      call reset_timerOMP()
      if (be_verbose) write (*, '(A)') "5.  Testing empty tic() lhs "
      open (10, file=DEV_NULL)
      if (be_verbose) write (10, '(A,1d20.10)') "     inline tic ", ticOMP()
      ! if(be_verbose) write (*, '(A,1d20.10)') "2. t1 = ", real(t1, dp)
      c = .1
      do i = 1, nrep
         c = cfac*c + b
      end do
      if (be_verbose) write (*, '(A,1d20.10)') "    tic/toc in write ", tocOMP()

      call reset_timerOMP()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '6.  Testing tic toc in preprocessor macro  '
      TBLOCKOMP()
      c = .1
      do i = 1, nrep
         c = cfac*c + b
      end do
      TSTOPOMP()
      if (be_verbose) write (*, "(A)") ' '

      call reset_timerOMP()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '7.  Testing tic toc in preprocessor macro with subroutine  '
      TBLOCKOMP()
      c = saxy(c)
      TSTOPOMP()
      if (be_verbose) write (*, "(A)") ' '

      call reset_timerOMP()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '8.  Testing timed_loop macro - must be in own subroutine '
      c = test_loop()

   end subroutine exec_OMPtimertest

   function saxy(c_in) result(c)
      real(dp), intent(in) :: c_in
      integer(dp), parameter :: nrep = NREP_MAX
      real(dp)              :: c, cfac, b
      integer(dp)           :: i
      b = 1.
      c = c_in
      c = .1
      cfac = .25
      do i = 1, nrep
         c = cfac*c + b
      end do
   end function saxy

   ! Loop macro has to declare an array and an int
   ! -- this means it has to be at the top or in
   !    its own function
   function test_loop() result(cloop)
      real(dp) :: cloop
      real(dp) :: timestamp
      real(dp), dimension(3):: elapsed
      integer :: ii
      do ii = 1, 3
         timestamp = ticOMP(); 
         cloop = REAL(.1, dp)
         cloop = saxy(cloop)
         elapsed(ii) = tocOMP(timestamp)
      end do
      write (*, '(A,A,1i4,A)') __FILE__, ":", __LINE__, ' *** Timed loop *** '
      write (*, '(A,1d20.10)') "    Average over 3 (sec):", SUM(elapsed, DIM=1)/REAL(3., dp)
      write (*, '(A,1d20.10)') "    Min time (sec) ", MINVAL(elapsed, DIM=1)
      write (*, "(A)") ' '

   end function test_loop

end module simple_timer_omp_test
