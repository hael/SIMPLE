program simple_test_classrefine
use simple_classrefine_tester
use simple_cmdline, only: cmdline
implicit none
type(cmdline) :: cline
logical       :: be_verbose=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'simple_test_classrefine stk=<particle.mrc> msk=<mask radius(in pixels)>'
    write(*,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('msk',  2)
call cline%checkvar('smpd', 3)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
call exec_classrefine_test( cline, be_verbose )
end program simple_test_classrefine