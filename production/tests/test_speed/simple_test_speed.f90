program simple_test_speed
use simple_defs
use simple_speedtester
use simple_cmdline, only: cmdline
use simple_strings, only: str2real
use simple_syslib
implicit none
type(cmdline)     :: cline
real              :: starttime, stoptime
logical           :: be_verbose=.false.
character(STDLEN) :: time
call date_and_time(TIME=time)
starttime = str2real(time)
if( command_argument_count() < 1 )then
    write(*,'(a)') 'simple_test_speed nthr=<number of threads> [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('nthr', 1)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
call exec_speedtest(cline, be_verbose)
call date_and_time(TIME=time)
stoptime = str2real(time)
write(*,'(a,1x,f9.2)') 'time elapsed (s): ', stoptime - starttime
end program simple_test_speed
