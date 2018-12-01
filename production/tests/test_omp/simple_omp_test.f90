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
program simple_omp_test
!$ use omp_lib
include 'simple_lib.f08'

use simple_timer_omp_test
use simple_test_omp_basics

#ifdef OPENMP_VERSION
#if OPENMP_VERSION >= 201511
use simple_test_omp45
#endif
#endif

implicit none
! include 'omp_lib.h'

real              :: starttime, stoptime
logical           :: be_verbose=.false.
character(STDLEN) :: timestr
call date_and_time(TIME=timestr)
starttime = str2real(timestr)
be_verbose = .true.

#ifdef _OPENMP
    print *, " Preprocessor macro _OPENMP ", _OPENMP
#else
    print *, " Preprocessor macro _OPENMP not defined"
#endif
write(logfhandle,'(a,i0)') 'OpenMP version: ', openmp_version
print *, ' Test omp flush '
!$omp flush
call test_omp_basics(10000)
print *, ' Test internal openMP '
call test_internal_omp
print *, ' Test parallel sections openMP '
call test_parallel_sections_omp
print *, ' Test parallel loops openMP '
call test_parallel_loops_omp
print *, ' Test shared race-conditions openMP '
call test_shared_race_condition
print *, ' Test OpenMP stand-alone '
call test_standalone_ok
print *, ' Test OpenMP collapse order '
call test_collapse_order
print *, ' Test OpenMP ordered example '
call test_ordered_example
print *, ' Test OpenMP firstprivate example '
call test_first_private
print *, ' Test OpenMP firstprivate example 2'
call test_first_private2
print *, ' Test OpenMP SIMD example 2'
call test_simd_example2
print *, ' Test OpenMP SIMD example 3'
call test_simd_example3
print *, ' Test OpenMP SIMD example 4'
call test_simd_example4
print *, ' Test OpenMP timer '

call exec_OpenMP_timer_test(be_verbose)
call date_and_time(TIME=timestr)
stoptime = str2real(timestr)
write(logfhandle,'(a,1x,f9.2)') '<<< intrinsic date_and_time elapsed (s): ', stoptime - starttime


#ifdef OPENMP_VERSION
#if OPENMP_VERSION >= 201511
print *, ' Test OpenMP SIMD example 5'
call test_simd_example5
print *, ' Test OpenMP affinity example'
call test_omp_affinity
print *, ' Test OpenMP cancellation example '
call test_omp_cancellation
#endif
#endif

end program simple_omp_test
