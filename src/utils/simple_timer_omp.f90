!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple timer module using OpenMP functions
!!
!! \author
!!     Michael Eager
!! Original version: 2017-03-15
!
! The SIMPLE code is distributed with the hope that it will be
! useful, but WITHOUT ANY WARRANTY. Redistribution and modification is regulated
! by the GNU General Public License.
! -----------------------------------------------------------------------------!
module simple_timer_omp
   use simple_defs
!$   use omp_lib
   implicit none
   real(dp), save :: last_time_point_mp = REAL(0., dp) !< Starting timesamp
   public :: tic_omp, tickrate_omp
   public :: toc_omp, tdiff_omp, tocprint_omp
   public :: now_omp, reset_timer_omp
   public :: last_time_point_mp
private
#include "simple_local_flags.inc"

contains

    !> Force timestamps and clock rate to zero
   subroutine reset_timer_omp()
      last_time_point_mp = REAL(0., dp)
   end subroutine reset_timer_omp
   !< Get the clock tick count per second
   real(dp) function tickrate_omp()
       tickrate_omp = REAL(1.0, dp)
   end function tickrate_omp
   !< Get system_clock timestamp
   real(dp) function tic_omp()
#ifdef OPENMP
      tic_omp = OMP_get_wtime()
#else
      call cpu_time(tic_omp)
#endif
      last_time_point_mp = tic_omp
   end function tic_omp

   !> Calculate the time from two timestamps
   real(dp) function tdiff_omp(tend, tstart)
      real(dp), intent(in) :: tend
      real(dp), intent(in), optional :: tstart
      if (present(tstart)) last_time_point_mp = tstart
      ! Calulate the time difference
      tdiff_omp = REAL(tend - last_time_point_mp, dp)
   end function tdiff_omp

   !> Complete the timing regime using a reference timestamp or the one
   !  in last_time_point_mp
   real(dp) function toc_omp(start_optional)
      real(dp), intent(in), optional ::  start_optional
      real(dp) :: end_point_mp = REAL(0.0, dp)
      if (present(start_optional)) last_time_point_mp = start_optional
#ifdef OPENMP
      end_point_mp = OMP_get_wtime()
#else
      call cpu_time(end_point_mp)
#endif
      toc_omp = tdiff_omp(end_point_mp, last_time_point_mp)
      last_time_point_mp = end_point_mp
   end function toc_omp

   !> Complete the timing regime using a reference timestamp or the one
   !  in last_time_point_mp then print the results
   real(dp) function tocprint_omp(start_optional)
      real(dp), intent(in), optional ::  start_optional
      real(dp) :: end_point = REAL(0.0, dp)
      if (present(start_optional)) last_time_point_mp = start_optional
#ifdef OPENMP
      end_point = OMP_get_wtime()
#else
      call cpu_time(end_point)
#endif
      tocprint_omp = tdiff_omp(end_point, last_time_point_mp)
      last_time_point_mp = end_point
      write (*, '(A,1d20.10)') " Elapsed time ", tocprint_omp
   end function tocprint_omp

   !> print current time and date
   subroutine now_omp()
      character(len=8)  :: date
      character(len=10) :: time
      write(logfhandle,*) "OpenMP time: ", tic_omp()
      call date_and_time(date, time)
      write (*, '(A,A,A,A,A,A,A)') 'Date: ', date(7:8), '-', date(5:6), '-', date(1:4), '\n'
      write (*, '(A,A,A,A,A,A,A)') 'Time: ', time(1:2), ':', time(3:4), ':', time(5:10), '\n'
   endsubroutine now_omp
end module simple_timer_omp
