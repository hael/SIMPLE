program simple_test_cmdline
use simple_cmdline
use simple_chash
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"
type(cmdline)         :: cline
type(chash)           :: job_descr
character(len=STDLEN) :: xarg, line
type(string)          :: fname, str_prg
logical               :: test_passed
integer               :: cmdstat, cmdlen, pos
test_passed  = .true.
fname        = 'file_command.txt'
xarg         = "prg=symmetrize_map"
cmdlen       = len(trim(xarg))
line         = "projname=system_name smpd=1.3 cs=2.7 kv=300 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies &
              & gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200"
cmdstat      = 1
pos  = index(xarg, '=')
call cmdline_err(cmdstat, cmdlen, xarg, pos)
call cline%read(line)
print *, '>>> CLINE'
call cline%printline()
call cline%writeline(fname)
print *, '>>> CLINE'
call cline%printline()
call cline%checkvar('projname',        1)
call cline%checkvar('smpd',            2)
call cline%checkvar('cs',              3)
call cline%checkvar('kv',              4)
call cline%checkvar('fraca',           5)
call cline%checkvar('total_dose',      6)
call cline%checkvar('dir_movies',      7)
call cline%checkvar('gainref',         8)
call cline%checkvar('nparts',          9)
call cline%checkvar('nthr',           10)
call cline%checkvar('moldiam_max',    11)
if( .not. cline%defined('projname')   ) test_passed=.false.
if( .not. cline%defined('smpd')       ) test_passed=.false.
if( .not. cline%defined('dir_movies') ) test_passed=.false.
if( .not. cline%defined('gainref')    ) test_passed=.false.
if( .not. cline%defined('moldiam_max')) test_passed=.false.
call cline%check()
call cline%gen_job_descr(job_descr)
call job_descr%set('prg',      'scale')
call job_descr%set('autoscale','no')
print *,'>>> JOB_DESCR'
call job_descr%print_key_val_pairs(logfhandle) 
call cline%set('prg', 'list')
str_prg = cline%get_carg('prg')
print *,'>>> PROGRAM ', str_prg%to_char()
if( .not. (cline%get_carg('prg').eq.'list')) test_passed=.false.
call cline%delete('kv')
print *,'>>> CS ',cline%get_rarg('cs')
print *,'>>> NPARTS ',cline%get_iarg('nparts')
if( test_passed )then
    print *, '>>> TEST PASSED'
else
    THROW_HARD('>>> TEST FAILED')
endif
call cline%kill()
call job_descr%kill()
end program simple_test_cmdline
