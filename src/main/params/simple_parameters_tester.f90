!@descr: parameter-constructor bindings used by focused high-level test selection
module simple_parameters_tester
use simple_string,         only: string
use simple_test_utils,     only: assert_char, assert_int, assert_real, assert_true, enter_fixture, leave_fixture, &
                                 tests_failed
use simple_forked_process, only: forked_process, FORK_STATUS_FAILED, FORK_STATUS_STOPPED
use simple_cmdline,        only: cmdline
implicit none
private
public :: run_all_parameters_tests
real, parameter :: SMPD = 1.5

!> a process that builds parameters from its command line (a yes/no key's value decides whether it stops)
type, extends(forked_process) :: parameters_child
contains
    procedure :: execute => build_parameters
end type parameters_child

contains

    subroutine run_all_parameters_tests()
        use simple_defs, only: CWD_GLOB, CWD_GLOB_ORIG
        type(string)                  :: saved, root
        character(len=:), allocatable :: cwd_before, origin_before
        integer                       :: nfail0
        write(*,'(A)') '**** running all parameters tests ****'
        nfail0 = tests_failed
        if( allocated(CWD_GLOB) ) cwd_before = CWD_GLOB
        if( allocated(CWD_GLOB_ORIG) ) origin_before = CWD_GLOB_ORIG
        call enter_fixture('parameters', saved, root)
        call test_cont_prior_binding()
        call test_cont_prior_rejections()
        call test_suite_binding()
        call test_accepted_keys()
        call test_ui_keys_bound()
        call test_yes_no_refused()
        call test_odd_box_refused()
        call test_worker_project_recorded()
        call test_mirr_proj_binding()
        call leave_fixture(saved, root, nfail0)
        ! The constructor must not leave directory globals pointing into the removed fixture.
        if( allocated(CWD_GLOB) ) deallocate(CWD_GLOB)
        if( allocated(CWD_GLOB_ORIG) ) deallocate(CWD_GLOB_ORIG)
        if( allocated(cwd_before) ) CWD_GLOB = cwd_before
        if( allocated(origin_before) ) CWD_GLOB_ORIG = origin_before
    end subroutine run_all_parameters_tests

    !> The explicit prior survives the same constructor used by serialized worker commands.
    subroutine test_cont_prior_binding()
        use simple_parameters, only: parameters
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline
        allocate(params, cline)
        call base_cline(cline)
        call params%new(cline, silent=.true.)
        call assert_char('none', trim(params%cont_prior), 'Cartesian prior defaults off')
        call cline%set('objfun', 'euclid')
        call cline%set('refine', 'cont')
        call cline%set('cont_prior', 'shift')
        call params%new(cline, silent=.true.)
        call assert_char('shift', trim(params%cont_prior), 'Cartesian prior survives parameter construction')
        call cline%set('cont_prior', 'pose')
        call cline%set('cont_prior_rot_sigma', 5.)
        call cline%set('cont_prior_shift_sigma', 2.)
        call params%new(cline, silent=.true.)
        call assert_char('pose', trim(params%cont_prior), 'pose prior survives parameter construction')
        call assert_real(5., params%cont_prior_rot_sigma, 0., 'rotation prior width survives construction')
        call assert_real(2., params%cont_prior_shift_sigma, 0., 'fixed shift prior width stays in Angstroms through construction')
        call cline%kill()
    end subroutine test_cont_prior_binding

    !> Invalid controls and the retiring split route fail before matching.
    subroutine test_cont_prior_rejections()
        type(parameters_child) :: child
        class(cmdline), allocatable :: cline
        allocate(cline)
        call base_cline(cline)
        call cline%set('objfun', 'euclid')
        call cline%set('refine', 'cont')
        call cline%set('cont_prior', 'pose')
        call child%start(name=string('PRIOR_RESERVED'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'pose prior without a width is rejected')
        call child%destroy
        call cline%set('cont_prior', 'shift')
        call cline%set('cont_route', 'shift_then_joint')
        call child%start(name=string('PRIOR_ROUTE'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'shift prior rejects split route')
        call child%destroy
        call cline%set('cont_route', 'joint')
        call cline%set('objfun', 'cc')
        call child%start(name=string('PRIOR_OBJECTIVE'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'shift prior rejects correlation scaling')
        call child%destroy
        call cline%kill()
    end subroutine test_cont_prior_rejections

    !> the command line every test here starts from
    subroutine base_cline( cline )
        class(cmdline), intent(inout) :: cline
        call cline%set('prg', 'cont_refine3D_1jxy')
        call cline%set('mkdir', 'no')
        call cline%set('box', 32)
        call cline%set('smpd', SMPD)
        call cline%set('mskdiam', 12.)
        call cline%set('nthr', 1)
    end subroutine base_cline

    subroutine test_mirr_proj_binding()
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline
        allocate(params, cline)
        call cline%set('box',32)
        call cline%set('smpd',SMPD)
        call cline%set('mskdiam',12.)
        call cline%set('nthr',1)
        call params%new(cline,silent=.true.)
        call assert_char('no',params%mirr_proj,'projection mirroring defaults to no')
        call assert_char('no',params%sigma2_frozen,'frozen bootstrap defaults to no')
        call cline%set('mirr_proj','yes')
        call cline%set('sigma2_frozen','yes')
        call params%new(cline,silent=.true.)
        call assert_char('yes',params%mirr_proj,'projection mirroring survives parameter construction')
        call assert_char('yes',params%sigma2_frozen,'frozen bootstrap survives parameter construction')
        call cline%kill()
    end subroutine test_mirr_proj_binding

    ! The parameter constructor must retain the selector before the expensive simulation.
    subroutine test_suite_binding()
        use simple_parameters, only: parameters
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: cline
        write(*,'(A)') 'test_suite_binding'
        allocate(params, cline)
        call base_cline(cline)
        call params%new(cline, silent=.true.)
        call assert_char('', trim(params%suite), 'suite parser: default selects the full gate')
        call cline%set('suite', 'nu_handoff')
        call params%new(cline, silent=.true.)
        call assert_char('nu_handoff', trim(params%suite), 'suite parser: focused NU handoff survives construction')
        call cline%set('suite', 'all')
        call params%new(cline, silent=.true.)
        call assert_char('all', trim(params%suite), 'suite parser: explicit full gate survives construction')
        call cline%kill()
    end subroutine test_suite_binding

    !> the accepted keys are the components of type(parameters) the generator reads, up to its end:
    !! no dummy argument of an interface after it, no fragment of an array constructor, no derived field
    subroutine test_accepted_keys()
        use simple_args, only: args
        type(args) :: accepted
        write(*,'(A)') 'test_accepted_keys'
        accepted = args()
        call assert_true(accepted%is_present('mul'),         'mul is an accepted key')
        call assert_true(accepted%is_present('import_type'), 'import_type is an accepted key')
        call assert_true(accepted%is_present('gold_std_stage'), 'gold_std_stage is an accepted key')
        call assert_true(.not. accepted%is_present('cntfile'),       'cntfile (a dummy argument after the type) is no key')
        call assert_true(.not. accepted%is_present('silent'),        'silent (a dummy argument after the type) is no key')
        call assert_true(.not. accepted%is_present('final_planned'), 'final_planned (a function result) is no key')
        call assert_true(.not. accepted%is_present('1.0]'),          'an array constructor''s fragment is no key')
        call assert_true(.not. accepted%is_present('l_update_frac'), 'a derived flag is no key')
    end subroutine test_accepted_keys

    !> UI inputs that were accepted but never bound (the field kept its default)
    subroutine test_ui_keys_bound()
        use simple_parameters, only: parameters
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: cline
        write(*,'(A)') 'test_ui_keys_bound'
        allocate(params, cline)
        call base_cline(cline)
        call cline%set('mul', 2.5)
        call cline%set('import_type', 'mic')
        call cline%set('tiltgroupmax', 3.)
        call assert_int(0, params%gold_std_stage, 'gold-standard start stage defaults to workflow policy')
        call cline%set('gold_std_stage', 4)
        call params%new(cline, silent=.true.)
        call assert_real(2.5, params%mul, 1.e-6, 'mul binds (oriops, volops)')
        call assert_char('mic', trim(params%import_type), 'import_type binds (import_starproject)')
        call assert_real(3., params%tiltgroupmax, 1.e-6, 'tiltgroupmax binds')
        call assert_int(4, params%gold_std_stage, 'gold_std_stage binds')
        call cline%kill()
    end subroutine test_ui_keys_bound

    !> a yes/no key refuses any other value; a len=3 field would otherwise keep 'yes' of 'yess'
    subroutine test_yes_no_refused()
        type(parameters_child)      :: child
        class(cmdline), allocatable :: cline
        write(*,'(A)') 'test_yes_no_refused'
        allocate(cline)
        call base_cline(cline)
        call cline%set('async', 'yes')
        call child%start(name=string('PARAMETERS_YES'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_STOPPED, child%status(), 'a yes/no key takes yes')
        call child%destroy
        call cline%set('async', 'yess')
        call child%start(name=string('PARAMETERS_YESS'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a yes/no key refuses yess')
        call child%destroy
        call cline%set('async', 'true')
        call child%start(name=string('PARAMETERS_TRUE'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a yes/no key refuses true')
        call child%destroy
        call cline%delete('async')
        call cline%set('sigma2_frozen', 'yess')
        call child%start(name=string('PARAMETERS_FROZEN_INVALID'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'sigma2_frozen refuses misspelled yes')
        call child%destroy
        call cline%kill()
    end subroutine test_yes_no_refused

    !> working grids are even: an odd box or box_crop stops at new (doc/policies/even_working_grid_policy.md)
    subroutine test_odd_box_refused()
        type(parameters_child)      :: child
        class(cmdline), allocatable :: cline
        write(*,'(A)') 'test_odd_box_refused'
        allocate(cline)
        call base_cline(cline)
        call cline%set('box', 33)
        call child%start(name=string('PARAMETERS_ODD_BOX'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'an odd box stops')
        call child%destroy
        call cline%set('box', 32)
        call cline%set('box_crop', 31)
        call child%start(name=string('PARAMETERS_ODD_BOX_CROP'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'an odd box_crop stops')
        call child%destroy
        call cline%set('box_crop', 30)
        call child%start(name=string('PARAMETERS_EVEN_BOX_CROP'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_STOPPED, child%status(), 'an even box_crop passes')
        call child%destroy
        call cline%kill()
    end subroutine test_odd_box_refused

    !> a partition worker's project is recorded at new (the project writers refuse it); a stream job's is not
    subroutine test_worker_project_recorded()
        use simple_parameters, only: parameters
        use simple_defs,       only: l_distr_worker_glob, worker_projfile_glob
        use simple_syslib,     only: simple_abspath
        use simple_sp_project, only: sp_project
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: cline
        class(sp_project), allocatable :: proj
        type(string) :: absname
        write(*,'(A)') 'test_worker_project_recorded'
        allocate(params, cline, proj)
        call proj%update_projinfo(string('worker_recorded.simple'))
        call proj%write(string('worker_recorded.simple'))
        call proj%kill
        absname = simple_abspath(string('worker_recorded.simple'))
        call base_cline(cline)
        call cline%set('projfile', 'worker_recorded.simple')
        call cline%set('part',     1)
        call cline%set('nparts',   2)
        call params%new(cline, silent=.true.)
        call assert_true(allocated(worker_projfile_glob), 'a partition worker records its project')
        if( allocated(worker_projfile_glob) ) call assert_char(absname%to_char(), worker_projfile_glob, &
            &'the recorded project is the worker''s input')
        call cline%set('stream', 'yes')
        call params%new(cline, silent=.true.)
        call assert_true(.not. allocated(worker_projfile_glob), 'a stream job owns its job project')
        call cline%delete('part')
        call cline%delete('stream')
        call params%new(cline, silent=.true.)
        call assert_true(.not. l_distr_worker_glob, 'the next new is no worker')
        call cline%kill()
    end subroutine test_worker_project_recorded

    subroutine build_parameters( self, cline )
        use simple_parameters, only: parameters
        class(parameters_child), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        class(parameters), allocatable :: params
        allocate(params)
        call params%new(cline, silent=.true.)
    end subroutine build_parameters

end module simple_parameters_tester
