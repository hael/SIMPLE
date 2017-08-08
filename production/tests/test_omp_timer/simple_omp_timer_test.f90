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
use simple_defs
use simple_timer
use simple_timer_omp_test
use simple_cmdline, only: cmdline
use simple_strings, only: str2real
use simple_syscalls
implicit none
type(cmdline)     :: cline
real              :: starttime, stoptime
logical           :: be_verbose=.false.
character(STDLEN) :: time
call date_and_time(TIME=time)
starttime = str2real(time)
!if( command_argument_count() < 0 )then
!    write(*,'(a)') 'simple_test_timer [verbose=<yes|no{no}>]'
!    stop
!endif
!call cline%parse
! call cline%checkvar('nthr', 1)
!call cline%check
be_verbose = .true.
!if( cline%defined('verbose') )then
!    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
!        be_verbose = .true.
!    endif
!endif
call exec_OpenMP_timer_test(be_verbose)
call date_and_time(TIME=time)
stoptime = str2real(time)
write(*,'(a,1x,f9.2)') '<<< intrinsic date_and_time elapsed (s): ', stoptime - starttime
end program simple_omp_timer_test
