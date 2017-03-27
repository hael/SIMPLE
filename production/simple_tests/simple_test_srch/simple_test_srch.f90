program simple_test_srch
use simple_prime2D_srch_tester
use simple_wiener2D_tester
use simple_optimiser_tester
use simple_prime3D_srch_tester
use simple_volpft_srch_tester
use simple_cmdline, only: cmdline
implicit none
type(cmdline) :: cline
logical       :: be_verbose=.false.
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'simple_test_srch vol1=<volume.mrc> msk=<mask radius(in pixels)>'
    write(*,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('vol1', 1)
call cline%checkvar('msk',  2)
call cline%checkvar('smpd', 3)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
call exec_prime2D_srch_test( cline, be_verbose )
! call exec_prime3D_srch_test( cline, be_verbose )
! call exec_wiener2D_test    ( cline, be_verbose )
! call exec_optimiser_test   (        be_verbose )
! call exec_volpft_srch_test ( cline, be_verbose )
end program simple_test_srch