


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

character(STDLEN) :: time
call date_and_time(TIME=time)
starttime = str2real(time)
if( command_argument_count() < 0 )then
    write(*,'(a)') 'simple_test_timer [verbose=<yes|no{no}>]'
    stop
endif

call exec_OMPtimertest()
call date_and_time(TIME=time)
stoptime = str2real(time)
write(*,'(a,1x,f9.2)') 'simple_omp_timer_test:  date_and_time elapsed (s): ', stoptime - starttime
end program simple_omp_timer_test
