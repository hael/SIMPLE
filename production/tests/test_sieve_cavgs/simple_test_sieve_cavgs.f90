program simple_test_sieve_cavgs
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_commanders_stream2D, only:stream_test_sieve_cavgs
implicit none
type(cmdline)                 :: cline
character(len=:), allocatable :: last_prev_dir, fname
integer                       :: idir
if( command_argument_count() < 8 )then
    write(logfhandle,'(a)') 'ERROR! Required parameters:'
    write(logfhandle,'(a)') 'smpd,cs,kv,fraca,filetab,pickrefs,nthr,ncls'
    write(logfhandle,'(a)') 'Optional parameters:'
    write(logfhandle,'(a)') 'nptcls_per_cls,box,mskdiam,maxpop'
    call exit(-1)
endif
call cline%parse_oldschool
call cline%checkvar('smpd',            1)
call cline%checkvar('cs',              2)
call cline%checkvar('kv',              3)
call cline%checkvar('fraca',           4)
call cline%checkvar('pickrefs  ',      5)
call cline%checkvar('filetab',         6)
call cline%checkvar('nthr',            7)
call cline%checkvar('ncls',            8)
call cline%check
! For test program convenience
idir = find_next_int_dir_prefix(PATH_HERE, last_prev_dir)
call cline%set('exec_dir', int2str(idir)//'_test_sieve_cavgs')
call cline%set('filetab',  simple_abspath(cline%get_carg('filetab')))
call cline%set('pickrefs', simple_abspath(cline%get_carg('pickrefs')))
call simple_mkdir( filepath(PATH_HERE, trim(cline%get_carg('exec_dir'))))
call simple_chdir( filepath(PATH_HERE, trim(cline%get_carg('exec_dir'))))
! execution
call stream_test_sieve_cavgs(cline)
end program simple_test_sieve_cavgs