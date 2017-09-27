!------------------------------------------------------------------------------!
! SIMPLE , Elmlund & Elmlund Lab,     simplecryoem.com                         !
!------------------------------------------------------------------------------!
!> test module for profiling routines in simple_timer
!!
!! Test the profile timing functions in the SIMPLE library.
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



module simple_timer_profile_test
#include "simple_lib.f08"
   use simple_timer
   implicit none
   public:: exec_profiletest
private
#include "simple_local_flags.inc"
#include "simple_timer.h"   
contains

    subroutine exec_profiletest(be_verbose)
      logical, optional, intent(in)    :: be_verbose
      integer(dp), parameter :: nrep = INT(100000,dp)
      real(dp)    :: c, cfac, b
      common b, c, cfac
      real(dp)    :: xx

      if(present(be_verbose)) verbose=be_verbose
      !    The following is a Statement function
      ! slaxpy(x,y,z)= x*y+z

      c = .1
      cfac = .25
      b = 1.5
      xx = 12.0_dp
      if (verbose) print *, 'Fortran Timer and Profiler'
#ifdef PGI
      print *, 'PGI cannot process block statements'
      stop
#else

      call reset_timer()
      if (verbose) write (*, "(A)") ' '
      if (verbose) write (*, '(A)') '1.  Testing profiler using macros inside loop'
      c = .1

      TPROFILER(nrep, i, foo bar)
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
      if (verbose) write (*, "(A)") ' '
      if (verbose) write (*, '(A)') '2.  Testing profiler using macros and seperate loops'
      c = .1

      TPROFILER(nrep, i,  foo bar)
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
      if (verbose) write (*, "(A)") ' '
      if (verbose) write (*, '(A)') '2.  Testing profiler using macros and seperate loops'
      c = .1

      TPROFILER(1, i, standard subrout common empty)
      TBEG(standard)
      do i = 1, nrep
         c = cfac*c + b
      end do
      TEND(standard)
      TBEG(subrout)
      c = .1
      c = saxy(c)
      TEND(subrout)
      TBEG(common)
      do i = 1, nrep
         call scsaxpy
      end do
      TEND(common)
      TBEG(empty)
      do i = 1, nrep

      end do
      TEND(empty)

      TREPORT( Testing different loop unrolling methods )
#endif

   end subroutine exec_profiletest

   function saxy(c_in) result(c)
      real(dp), intent(in) :: c_in
      integer(dp), parameter :: nrep = INT(100000,dp)
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

