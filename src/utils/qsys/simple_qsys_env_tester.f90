!@descr: unit tests for the queue-system environment's installation-path policy (simple_qsys_env, sp_project%update_compenv)
! An installation path stored in a project by an older SIMPLE is runtime-local: update_compenv
! drops it, and a distributed run takes the local SIMPLE_PATH (set by the CTest environment) for
! its queue description and executable, never the project's (the policy of 9a87c12a8: executables
! come from the environment), and derives the job time of an empty project (0-0:1:40). The project
! file is written in the suite's directory and deleted.
module simple_qsys_env_tester
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
use simple_test_utils
implicit none
private
public :: run_all_qsys_env_tests

character(len=*), parameter :: LEGACY_SIMPLE_PATH = '/remote/install/that/is/not/local'

contains

    subroutine run_all_qsys_env_tests()
        write(*,'(A)') '**** running all qsys environment tests ****'
        call test_installation_path_policy()
    end subroutine run_all_qsys_env_tests

    subroutine test_installation_path_policy()
        type(cmdline)    :: cline
        type(parameters) :: params
        type(qsys_env)   :: qenv
        type(sp_project) :: project
        type(string)     :: exec_bin, expected_exec_bin, job_time, local_simple_path, projfile
        integer          :: iostat
        write(*,'(A)') 'test_installation_path_policy'
        projfile = 'test_qsys_env_path_policy.simple'
        call del_file(projfile)
        ! updating project metadata strips the installation path left by an older project
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'test_qsys_env_path_policy')
        call project%compenv%new(1, is_ptcl=.false.)
        call project%compenv%set(1, 'simple_path', LEGACY_SIMPLE_PATH)
        call cline%set('qsys_name', 'local')
        call project%update_compenv(cline)
        call assert_false(project%compenv%isthere('simple_path'), 'update_compenv drops the runtime-local simple_path')
        ! a legacy project: distributed execution must ignore its stored path
        call project%compenv%set(1, 'simple_path', LEGACY_SIMPLE_PATH)
        call project%write(projfile)
        local_simple_path = simple_getenv('SIMPLE_PATH', iostat)
        call assert_int(0, iostat, 'SIMPLE_PATH is set in the test environment')
        if( iostat == 0 )then
            params%prg             = 'qsys_env_path_policy'
            params%projfile        = projfile
            params%qsys_name       = 'local'
            params%nparts          = 1
            params%ncunits         = 1
            params%nptcls          = 0
            params%nthr            = 1
            params%worker_priority = 'normal'
            params%worker_server   = ''
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
            exec_bin          = qenv%get_exec_bin()
            expected_exec_bin = filepath(local_simple_path, string('bin'), string('simple_exec'))
            call assert_string_eq(expected_exec_bin%to_char(), exec_bin, &
                &'the executable resolves from the local SIMPLE_PATH, not the project')
            call assert_true(qenv%qdescr%isthere('simple_path'), 'the queue description carries a simple_path')
            if( qenv%qdescr%isthere('simple_path') )then
                call assert_string_eq(local_simple_path%to_char(), qenv%qdescr%get('simple_path'), &
                    &'the queue description carries the local SIMPLE_PATH, not the project one')
            endif
            call assert_true(qenv%qdescr%isthere('job_time'), 'the queue description carries a derived job time')
            if( qenv%qdescr%isthere('job_time') )then
                job_time = qenv%qdescr%get('job_time')
                call assert_string_eq('0-0:1:40', job_time, 'an empty project gets the minimum job time')
            endif
            call qenv%kill
        endif
        call project%kill
        call cline%kill
        call del_file(projfile)
    end subroutine test_installation_path_policy

end module simple_qsys_env_tester
