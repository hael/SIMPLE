!@descr: parameter-constructor bindings and execution-directory project ownership
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
        call test_execution_project_copy()
        call test_mirr_proj_binding()
        call test_sigma2_anneal_defaults()
        call test_cont_view()
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

    !> Invalid controls fail before matching.
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
        call cline%set('objfun', 'cc')
        call child%start(name=string('PRIOR_OBJECTIVE'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'shift prior rejects correlation scaling')
        call child%destroy
        call cline%kill()
    end subroutine test_cont_prior_rejections

    !> cont_view (polar_continuous_view.md): accepted on the refine3D search of one state with refine=prob_neigh,
    !! objfun=euclid and nspace_sub = nspace; any other value, band of directions, refinement, state count or objective
    !! stops; another program (the master's children inherit its command line) ignores it
    subroutine test_cont_view()
        use simple_parameters, only: parameters
        type(parameters_child)         :: child
        class(cmdline),    allocatable :: cline
        class(parameters), allocatable :: params
        write(*,'(A)') 'test_cont_view'
        allocate(cline, params)
        call base_cline(cline)
        call cline%set('cont_view', 'lm')
        call params%new(cline, silent=.true.)
        call assert_true(.not. params%l_cont_view, 'cont_view is ignored by a program other than the refine3D search')
        call cline%set('prg',        'refine3D')
        call cline%set('refine',     'prob_neigh')
        call cline%set('objfun',     'euclid')
        call cline%set('nspace',     1000)
        call cline%set('nspace_sub', 1000)
        call view_case('PARAMETERS_CONT_VIEW', FORK_STATUS_STOPPED, 'cont_view=lm on refine3D refine=prob_neigh is accepted')
        call cline%set('prob_neigh_mode', 'geom')
        call cline%set('topk', 1)
        call view_case('PARAMETERS_CONT_VIEW_GEOM', FORK_STATUS_STOPPED, 'geom view seeding accepts the existing topk=1')
        call params%new(cline, silent=.true.)
        call assert_int(1, params%topk, 'geometry seed budget is the parsed topk')
        call assert_char('geom', trim(params%prob_neigh_mode), 'geometry view seed mode survives parsing')
        call cline%set('topk', 0)
        call view_case('PARAMETERS_CONT_VIEW_TOPK', FORK_STATUS_FAILED, 'geometry view seeding refuses topk=0')
        call cline%set('topk', 3)
        call cline%set('prob_neigh_mode', 'shc')
        call view_case('PARAMETERS_CONT_VIEW_SHC', FORK_STATUS_FAILED, 'cont_view refuses stochastic seed modes')
        call cline%set('prob_neigh_mode', 'state')
        call cline%set('cont_view', 'lbfgsb')
        call view_case('PARAMETERS_CONT_VIEW_LBFGSB', FORK_STATUS_STOPPED, 'cont_view=lbfgsb is accepted')
        call cline%set('cont_view', 'yes')
        call view_case('PARAMETERS_CONT_VIEW_VALUE', FORK_STATUS_FAILED, 'cont_view refuses a value other than no|lbfgsb|lm')
        call cline%set('cont_view', 'lm')
        call cline%set('nspace_sub', 500)
        call view_case('PARAMETERS_CONT_VIEW_SUB', FORK_STATUS_FAILED, 'cont_view refuses nspace_sub different from nspace')
        call cline%set('nspace_sub', 1000)
        call cline%set('refine', 'cont')
        call view_case('PARAMETERS_CONT_VIEW_CONT', FORK_STATUS_FAILED, 'cont_view refuses refine=cont')
        call cline%set('refine', 'prob_neigh')
        call cline%set('nstates', 2)
        call view_case('PARAMETERS_CONT_VIEW_STATES', FORK_STATUS_FAILED, 'cont_view refuses two states')
        call cline%delete('nstates')
        call cline%set('objfun', 'cc')
        call view_case('PARAMETERS_CONT_VIEW_CC', FORK_STATUS_FAILED, 'cont_view refuses objfun=cc')
        call cline%kill
    contains
        subroutine view_case( name, expected, label )
            character(len=*), intent(in) :: name, label
            integer,          intent(in) :: expected
            call child%start(name=string(name), cline=cline)
            call child%await_final_status()
            call assert_int(expected, child%status(), label)
            call child%destroy
        end subroutine view_case
    end subroutine test_cont_view

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
        call assert_char('no',params%sigma2_anneal,'sigma annealing defaults to no')
        call cline%set('mirr_proj','yes')
        call cline%set('sigma2_anneal','yes')
        call cline%set('sigma2_weight_start',0.5)
        call cline%set('sigma2_weight_end',1.)
        call cline%set('sigma2_weight',0.75)
        call params%new(cline,silent=.true.)
        call assert_char('yes',params%mirr_proj,'projection mirroring survives parameter construction')
        call assert_char('yes',params%sigma2_anneal,'sigma annealing survives parameter construction')
        call assert_real(0.5,params%sigma2_weight_start,1.e-6,'sigma ramp start survives construction')
        call assert_real(1.,params%sigma2_weight_end,1.e-6,'sigma ramp end survives construction')
        call assert_real(0.75,params%sigma2_weight,1.e-6,'worker sigma weight survives construction')
        call cline%set('sigma2_weight_start',-1.)
        call cline%set('sigma2_weight',-0.5)
        call params%new(cline,silent=.true.)
        call assert_char('yes',params%sigma2_anneal,'enabled sigma annealing accepts the white start')
        call assert_real(-1.,params%sigma2_weight_start,1.e-6,'sigma ramp starts at the white ceiling')
        call assert_real(-0.5,params%sigma2_weight,1.e-6,'negative worker sigma weight survives construction')
        call cline%kill()
    end subroutine test_mirr_proj_binding

    subroutine test_sigma2_anneal_defaults()
        use simple_default_clines, only: set_refine2D_defaults
        class(cmdline), allocatable :: cline
        write(*,'(A)') 'test_sigma2_anneal_defaults'
        allocate(cline)
        call cline%set('mskdiam',12.)
        call set_refine2D_defaults(cline)
        call assert_true(.not. cline%defined('sigma2_weight_start'),'default refine2D adds no annealing ramp')
        call cline%set('sigma2_anneal','no')
        call set_refine2D_defaults(cline)
        call assert_true(.not. cline%defined('sigma2_weight_start'),'disabled annealing adds no ramp')
        call cline%set('sigma2_anneal','yes')
        call set_refine2D_defaults(cline)
        call assert_real(-1.,cline%get_rarg('sigma2_weight_start'),0.,'standalone annealing starts white')
        call cline%set('sigma2_weight_start',0.5)
        call set_refine2D_defaults(cline)
        call assert_real(0.5,cline%get_rarg('sigma2_weight_start'),0.,'explicit stage start overrides the white default')
        call cline%kill()
    end subroutine test_sigma2_anneal_defaults

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
        call cline%set('sigma2_anneal', 'yess')
        call child%start(name=string('PARAMETERS_SIGMA2_ANNEAL_INVALID'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'sigma2_anneal refuses misspelled yes')
        call child%destroy
        call cline%set('sigma2_anneal', 'flat')
        call child%start(name=string('PARAMETERS_SIGMA2_FLAT'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'sigma2_anneal refuses the retired flat value')
        call child%destroy
        call cline%set('sigma2_anneal', 'no')
        call cline%set('sigma2_weight_start', -1.)
        call child%start(name=string('PARAMETERS_SIGMA2_NEGATIVE_NO'), cline=cline)
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a weight below 0 needs sigma2_anneal=yes')
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

    !> A public workflow owns a copy of its project in its numbered run directory. The executable is identified
    !! by its file name alone: a build path containing 'private' made a public run behave as a worker (no copy).
    subroutine test_execution_project_copy()
        use simple_parameters, only: parameters
        use simple_sp_project, only: sp_project
        use simple_fileio,     only: filepath
        use simple_syslib,     only: file_exists, simple_getcwd
        use simple_ui,         only: make_ui
        class(parameters), allocatable :: params
        class(cmdline), allocatable    :: cline
        class(sp_project), allocatable :: proj
        type(string) :: saved, root, input_project, run_project, cwd, cline_project
        integer :: nfail0
        write(*,'(A)') 'test_execution_project_copy'
        allocate(params, cline, proj)
        call make_ui
        nfail0 = tests_failed
        call enter_fixture('execution_context', saved, root)
        input_project = filepath(root, 'input.simple')
        run_project   = filepath(filepath(root, '1_context'), 'input.simple')
        call proj%update_projinfo(input_project)
        call proj%os_ptcl2D%new(1, is_ptcl=.true.)
        call proj%os_ptcl2D%set(1, 'state', 1)
        call proj%write(input_project)
        call proj%kill
        call base_cline(cline)
        call cline%set('prg',      'make_cavgs')
        call cline%set('projfile', input_project)
        call cline%set('oritype',  'ptcl2D')
        call cline%set('mkdir',    'yes')
        call cline%set('dir_exec', 'context')
        call params%new(cline, silent=.true.)
        call assert_char('simple_test_exec', params%executable%to_char(), 'executable identity excludes its directory')
        call assert_true(params%sp_required, 'a public project workflow requires its project')
        call assert_true(file_exists(run_project), 'mkdir copies the input project to the run directory')
        call simple_getcwd(cwd)
        call assert_char(root%to_char()//'/1_context', cwd%to_char(), 'a public run enters its numbered directory')
        call assert_char(run_project%to_char(), params%projfile%to_char(), 'parameters point to the run project')
        cline_project = cline%get_carg('projfile')
        call assert_char(run_project%to_char(), cline_project%to_char(), 'child commands inherit the run project')
        call proj%read_segment('ptcl2D', params%projfile)
        call proj%os_ptcl2D%set(1, 'state', 0)
        call proj%write_segment_inside('ptcl2D', params%projfile)
        call proj%kill
        call proj%read_segment('ptcl2D', input_project)
        call assert_int(1, proj%os_ptcl2D%get_state(1), 'writing the run project leaves the input state unchanged')
        call proj%kill
        if( file_exists(run_project) )then
            call proj%read_segment('ptcl2D', run_project)
            call assert_int(0, proj%os_ptcl2D%get_state(1), 'the changed state belongs to the run project')
            call proj%kill
        endif
        call cline%kill
        call leave_fixture(saved, root, nfail0)
    end subroutine test_execution_project_copy

    subroutine build_parameters( self, cline )
        use simple_parameters, only: parameters
        class(parameters_child), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        class(parameters), allocatable :: params
        allocate(params)
        call params%new(cline, silent=.true.)
    end subroutine build_parameters

end module simple_parameters_tester
