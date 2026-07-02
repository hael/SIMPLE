program simple_test_cavg_compatibility
use simple_core_module_api
use simple_cavg_compatibility_analysis, only: cavg_compatibility_analysis
use simple_parameters,               only: parameters
use simple_cmdline, only: cmdline
use simple_sp_project,              only: sp_project
use simple_string,                  only: string

implicit none
#include "simple_local_flags.inc"

type(cavg_compatibility_analysis) :: analysis
type(parameters) :: params
type(sp_project) :: spproj
type(cmdline) :: cline
character(len=STDLEN)             :: projfile
character(len=STDLEN)             :: refstk

if( command_argument_count() < 1 )then
  write(logfhandle,'(A)') 'Usage: simple_test_cavg_compatibility <project.simple> [reference_stack.mrcs]'
    stop
end if

call get_command_argument(1, projfile)
if( len_trim(projfile) == 0 ) THROW_HARD('simple_test_cavg_compatibility: empty project path')
refstk = ''
if( command_argument_count() >= 2 )then
  call get_command_argument(2, refstk)
  if( len_trim(refstk) == 0 ) THROW_HARD('simple_test_cavg_compatibility: empty reference stack path')
end if

call cline%set('projfile', trim(projfile))
call params%new(cline)
call spproj%read(params%projfile)
call analysis%new()
if( len_trim(refstk) > 0 )then
  call analysis%train(string(trim(refstk)))
end if
call analysis%train(spproj)
call analysis%infer(spproj)
call analysis%kill()

write(logfhandle,'(A,A)') 'simple_test_cavg_compatibility complete for project: ', trim(projfile)

end program simple_test_cavg_compatibility
