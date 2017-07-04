
#if defined  _WIN32
#define DEV_NULL "nul"
#else
#define DEV_NULL "/dev/null"
#endif

#define NREP_MAX INT(10000000,8)


#include "simple_timer.h"

module simple_timer_profile_test
   use simple_defs
   use simple_timer
   implicit none
   public:: exec_profiletest
contains
   subroutine exec_profiletest(be_verbose)
      logical, optional, intent(in)    :: be_verbose
      integer(dp), parameter :: nrep = 100000_dp
      real(dp)    :: c, cfac, b
      common b, c, cfac
      real(dp)    :: xx, etime, sysclockrate
      integer(dp) ::  t1, t2
      integer(dp) :: i, j

      !    The following is a Statement function
      ! slaxpy(x,y,z)= x*y+z

      c = .1
      cfac = .25
      b = 1.5
      xx = 12.0_dp
      if (be_verbose) print *, 'Fortran Timer and Profiler'
#ifdef PGI
      print *, 'PGI cannot process block statements'
      stop
#else

      call reset_timer()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '1.  Testing profiler using macros inside loop'
      c = .1
#ifndef Intel
      TPROFILER(nrep, i, foo, bar)
#else
      TPROFILER(nrep, i, 2, foo bar)
#endif
      do i = 1, nrep
         TBEG(foo)
         c = cfac*c + b
         TEND(foo)
         TBEG(bar)
         c = cfac*c + b
         TEND(bar)
         if (mod(nrep, nrep/1000) .eq. 1) print *, 'Repetition ', i

      end do
      TREPORT( Testing profiler using macros )

      call reset_timer()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '2.  Testing profiler using macros and seperate loops'
      c = .1
#ifndef Intel
      TPROFILER(1, i, standard, subrout, common, empty)
#else
      TPROFILER(nrep, i, 2, foo bar)
#endif
      do i = 1, nrep
         TBEG(foo)
         c = cfac*c + b
         TEND(foo)
         TBEG(bar)
         c = cfac*c + b
         TEND(bar)
         if (mod(nrep, nrep/1000) .eq. 1) print *, 'Repetition ', i

      end do
      TREPORT( Testing profiler using macros )

      call reset_timer()
      if (be_verbose) write (*, "(A)") ' '
      if (be_verbose) write (*, '(A)') '2.  Testing profiler using macros and seperate loops'
      c = .1
#ifndef Intel
      TPROFILER(1, i, standard, subrout, common, empty)
#else
      TPROFILER(1, i, 4, standard subrout common empty)
#endif
      TBEG(standard)
      do i = 1, nrep
         c = cfac*c + b
      end do
      TEND(standard)
      TBEG(subrout)
      c = .1
      c = saxy(c)
      TEND(subrout)
      ! TBEG(statement)
      ! do i=1,nrep
      !     c = slaxpy(cfac,c,b)
      ! end do
      ! TEND(statement)
      TBEG(common)
      do i = 1, nrep
         call scsaxpy
      end do
      TEND(common)
      TBEG(empty)
      do i = 1, nrep

      end do
      TEND(empty)

      TREPORT( Testing different loop unrolling methods)
#endif
   end subroutine exec_profiletest

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

   subroutine scsaxpy
      !   Multply all contents of  "c" by the scalar "cfac"
      !   then add the result to  "b"
      !   John Mahaffy    4/3/96

      implicit none
      real(dp) b, c, cfac
      common b, c, cfac
      c = cfac*c + b
      return
   end subroutine scsaxpy

end module simple_timer_profile_test
