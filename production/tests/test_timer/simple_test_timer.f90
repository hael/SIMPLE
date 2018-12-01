!------------------------------------------------------------------------------!
! SIMPLE , Elmlund & Elmlund Lab,     simplecryoem.com                         !
!------------------------------------------------------------------------------!
!> test program for simple_timer
!!
!! Test the timing functions in the SIMPLE library.
!!
!! @author
!! Michael Eager 2017
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
! -----------------------------------------------------------------------------!
program simple_test_timer
include 'simple_lib.f08'
use simple_timer_basic_test
use simple_timer_profile_test
use simple_cmdline, only: cmdline
implicit none
type(cmdline)     :: cline
real              :: starttime, stoptime
logical           :: be_verbose=.false.
character(STDLEN) :: timestr
call date_and_time(TIME=timestr)
starttime = str2real(timestr)
!if( command_argument_count() < 0 )then
!    write(logfhandle,'(a)') 'simple_test_timer [verbose=<yes|no{no}>]'
!    stop
!endif
!call cline%parse
! call cline%checkvar('nthr', 1)
! call cline%check
be_verbose = .true.
!if( cline%defined('verbose') )then
!    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
!        be_verbose = .true.
!    endif
!endif
write(logfhandle,'(a)') '------------- Basic timer test ----------------- '
call exec_timertest(be_verbose)
call date_and_time(TIME=timestr)
stoptime = str2real(timestr)
write(logfhandle,'(a,1x,f9.2)') 'simple_test_timer_basic total runtime (s): ', stoptime - starttime
starttime=stoptime
write(logfhandle,'(a)') '------------- Profiling timer test ----------------- '
call exec_profiletest(be_verbose)
call date_and_time(TIME=timestr)
stoptime = str2real(timestr)
write(logfhandle,'(a,1x,f9.2)') 'simple_test_timer_basic total runtime (s): ', stoptime - starttime
end program simple_test_timer
