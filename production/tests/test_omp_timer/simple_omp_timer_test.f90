!------------------------------------------------------------------------------!
! SIMPLE , Elmlund & Elmlund Lab,     simplecryoem.com                         !
!------------------------------------------------------------------------------!
!> test program for simple_timer_omp
!!
!! Test the OpenMP timing functions in the SIMPLE library.
!!
!! @author
!! Michael Eager 2017
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
! -----------------------------------------------------------------------------!
program simple_omp_timer_test
include 'simple_lib.f08'
use simple_timer_omp_test
use simple_test_omp
implicit none

real              :: starttime, stoptime
logical           :: be_verbose=.false.
character(STDLEN) :: timestr
call date_and_time(TIME=timestr)
starttime = str2real(timestr)

be_verbose = .true.

call exec_OpenMP_timer_test(be_verbose)
call date_and_time(TIME=timestr)
stoptime = str2real(timestr)
write(*,'(a,1x,f9.2)') '<<< intrinsic date_and_time elapsed (s): ', stoptime - starttime

call test_omp_basics(10000)
call test_internal_omp
call test_parallel_omp
call test_shared_race_condition


end program simple_omp_timer_test
