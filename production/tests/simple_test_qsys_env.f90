program simple_test_qsys_env
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none
#include "simple_local_flags.inc"
character(len=*), parameter :: LEGACY_SIMPLE_PATH = '/remote/install/that/is/not/local'
type(cmdline)    :: cline
type(parameters) :: params
type(qsys_env)   :: qenv
type(sp_project) :: project
type(string)     :: exec_bin, expected_exec_bin, local_simple_path, projfile
integer          :: iostat

projfile = 'test_qsys_env_path_policy.simple'
call del_file(projfile)

! Updating project metadata must strip installation paths left by older projects.
call project%projinfo%new(1, is_ptcl=.false.)
call project%projinfo%set(1, 'projname', 'test_qsys_env_path_policy')
call project%compenv%new(1, is_ptcl=.false.)
call project%compenv%set(1, 'simple_path', LEGACY_SIMPLE_PATH)
call cline%set('qsys_name', 'local')
call project%update_compenv(cline)
if( project%compenv%isthere('simple_path') )then
    THROW_HARD('update_compenv retained the runtime-local simple_path')
endif

! Emulate a legacy project and verify that distributed execution ignores its path.
call project%compenv%set(1, 'simple_path', LEGACY_SIMPLE_PATH)
call project%write(projfile)
local_simple_path = simple_getenv('SIMPLE_PATH', iostat)
if( iostat /= 0 ) THROW_HARD('SIMPLE_PATH is required for qsys path-policy test')

params%prg             = 'qsys_env_path_policy'
params%projfile        = projfile
params%qsys_name       = 'local'
params%nparts          = 1
params%ncunits         = 1
params%nptcls          = 1
params%nthr            = 1
params%worker_priority = 'normal'
params%worker_server   = ''
call qenv%new(params, 1, exec_bin=string('simple_exec'))
exec_bin          = qenv%get_exec_bin()
expected_exec_bin = filepath(local_simple_path, string('bin'), string('simple_exec'))
if( exec_bin /= expected_exec_bin )then
    write(logfhandle,'(A)') 'Expected local executable: '//expected_exec_bin%to_char()
    write(logfhandle,'(A)') 'Resolved executable:       '//exec_bin%to_char()
    THROW_HARD('qsys_env used a project-stored installation path')
endif
if( qenv%qdescr%isthere('simple_path') )then
    THROW_HARD('qsys_env retained the legacy simple_path in its queue description')
endif

call qenv%kill
call project%kill
call cline%kill
call del_file(projfile)
write(logfhandle,'(A)') 'QSYS LOCAL EXECUTABLE PATH POLICY TEST PASSED'
end program simple_test_qsys_env
