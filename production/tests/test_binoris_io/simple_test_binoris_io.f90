program simple_test_binoris_io
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_commanders_project, only: commander_new_project
use simple_sp_project
use simple_binoris_io
implicit none
type(commander_new_project) :: xnew_project
type(cmdline)               :: cline
type(oris)                  :: o
type(sp_project)            :: spproj 
logical                     :: test_passed
type(string)                :: fname, projname
integer                     :: fromto(2)
#include "simple_local_flags.inc"
test_passed = .true.
fname       = 'spproject.simple'
projname    = 'spproject'
fromto      = [ 1, 2 ]
call cline%set('projname', projname)
call cline%set('mkdir',        'no')
call cline%check()
call xnew_project%execute_safe(cline)
print *, '>>> BINREAD ORITAB'
call binread_oritab(fname, spproj, o, fromto)
print *, '>>> CTFPARAMS STATE EO'
call binread_ctfparams_state_eo(fname, spproj, o, fromto)
print *, '>>> READ NLINES ',binread_nlines(fname)
print *, '>>> WRITE ORITAB'
call binwrite_oritab(fname, spproj, o, fromto)
if( test_passed )then
    print *, '>>> TEST PASSED'
else
    THROW_HARD('>>> TEST FAILED')
endif
end program simple_test_binoris_io

