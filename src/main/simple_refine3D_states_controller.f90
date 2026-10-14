!@descr: the refine3D_states frequency march and its flex_pca state initialization
! The commander keeps its defaults and its ending; the settings that belong to its lineage arrive in a
! refine3D_states_policy.
module simple_refine3D_states_controller
use simple_class_sample_io,  only: deallocate_class_samples
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_defs,             only: STDLEN, logfhandle
use simple_error,            only: simple_exception
use simple_gui_communicator, only: gui_communicator
use simple_parameters,       only: parameters
use simple_sp_project,       only: sp_project
use simple_string,           only: string
use simple_string_utils,     only: int2str
use simple_syslib,           only: file_exists
use simple_type_defs,        only: class_sample, lp_crop_inf
implicit none
private
#include "simple_local_flags.inc"

public :: refine3D_states_controller, refine3D_states_policy

integer, parameter :: NSAMPLE_PER_STATE_REFINE3D_STATES           = 10000
integer, parameter :: NSAMPLE_REFINE3D_STATES_CAP                 = 100000
integer, parameter :: STAGE2_NSPACE                               = 5000
integer, parameter :: STAGE2_NSPACE_SUB                           = 500
integer, parameter :: FREQUENCY_BLOCK_NITS                        = 3
integer, parameter :: STAGE2_MINITS                               = 5
integer, parameter :: MINITS_REFINE3D_STATES                      = 10
integer, parameter :: MAXITS_REFINE3D_STATES_CAP                  = 50
real,    parameter :: TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES = 4.0
real,    parameter :: MIN_STATE_FRAC_FLEX                         = 0.1

!> the settings that differ between the lineages the controller serves; the defaults are the particle lineage's
type :: refine3D_states_policy
    character(len=STDLEN) :: label                   = ''         !< workflow label of the log lines and errors
    character(len=STDLEN) :: flex_rec_backend        = 'gridding' !< rec_backend of the flex_pca initialization
    character(len=STDLEN) :: flex_rec_states_backend = ''         !< its state-map backend; blank keeps flex's default
    integer               :: flex_min_neff           = 0          !< its minimum effective state size; 0 keeps flex's default
    logical               :: l_flex_native_box       = .false.    !< flex works at the native box
    logical               :: l_flex_lp_lpstart       = .false.    !< flex's covariance band is the march's lpstart
    logical               :: l_flex_shmem            = .false.    !< flex runs in shared memory
    logical               :: l_full_batch            = .false.    !< every active row in every iteration, no sample
    logical               :: l_gui_stages            = .true.     !< every stage map is published to the GUI
end type refine3D_states_policy

type :: refine3D_states_controller
    private
    type(refine3D_states_policy)    :: policy
    type(lp_crop_inf)               :: lpinfo_multi(2)
    type(string),       allocatable :: init_vols(:)
    type(class_sample), allocatable :: clssmp_units(:)
    integer :: nstates_project    = 0
    integer :: nptcls_eff         = 0
    integer :: nptcls_per_iter    = 0
    integer :: local_nspace_sub   = STAGE2_NSPACE_SUB
    integer :: maxits_user        = 0
    integer :: stage_cap          = 0
    integer :: total_iter         = 0
    integer :: gui_stage          = 0
    real    :: local_ang_bound    = -1.
    real    :: local_inpl_bound   = -1.
    real    :: local_shift_bound  = -1.
    logical :: l_maxits_defined   = .false.
    logical :: l_nstates_on_cline = .false.
    logical :: l_flex_init        = .false.
    logical :: l_nsample_auto     = .false.
  contains
    procedure          :: new
    procedure          :: normalize_inputs
    procedure          :: resolve_states
    procedure          :: resolve_sampling
    procedure          :: read_maxits
    procedure          :: validate
    procedure          :: initialize_states
    procedure          :: plan
    procedure          :: startup_maps
    procedure          :: run_march
    procedure          :: ensure_all_updated
    procedure          :: kill
    procedure, private :: label
    procedure, private :: default_nsample
    procedure, private :: reject_input_volumes
    procedure, private :: project_has_multistate_labels
    procedure, private :: set_nstates
    procedure, private :: set_sampling
    procedure, private :: prepare_sampling_units
    procedure, private :: report_coverage
    procedure, private :: set_downscaling
    procedure, private :: initialize_state_volumes
    procedure, private :: project_state_volumes_compatible
    procedure, private :: run_flex_pca
    procedure, private :: run_frequency_march
    procedure, private :: run_stage
    procedure, private :: read_update_coverage
    procedure, private :: run_missing_update
    procedure, private :: next_iteration
    procedure, private :: read_state_overlap
end type refine3D_states_controller

contains

    subroutine new( self, policy )
        class(refine3D_states_controller), intent(inout) :: self
        type(refine3D_states_policy),      intent(in)    :: policy
        call self%kill
        self%policy = policy
    end subroutine new

    !> Before params%new: refuses volume inputs and derived keys, and translates pose_policy and the
    !! local pose bounds into the matcher's keys
    subroutine normalize_inputs( self, cline )
        use simple_linalg, only: deg2rad
        use simple_math,   only: round2even
        class(refine3D_states_controller), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(string) :: pose_policy_arg
        call self%reject_input_volumes(cline)
        self%l_nstates_on_cline = cline%defined('nstates')
        if( cline%defined('multivol_mode') ) THROW_HARD(self%label()//' uses pose_policy, not multivol_mode')
        if( cline%defined('prob_neigh_mode') ) THROW_HARD(self%label()//' derives prob_neigh_mode from pose_policy')
        if( .not. cline%defined('pose_policy') ) call cline%set('pose_policy', 'global')
        pose_policy_arg = cline%get_carg('pose_policy')
        select case(trim(pose_policy_arg%to_char()))
            case('local')
                call cline%set('multivol_mode',   'input_oris_refine')
                call cline%set('prob_neigh_mode', 'geom')
            case('global')
                call cline%set('multivol_mode',   'input_oris_refine')
                call cline%set('prob_neigh_mode', 'state')
            case default
                THROW_HARD(self%label()//' supports pose_policy=local|global')
        end select
        if( cline%defined('local_ang_bound') ) self%local_ang_bound = cline%get_rarg('local_ang_bound')
        if( cline%defined('local_inpl_bound') ) self%local_inpl_bound = cline%get_rarg('local_inpl_bound')
        if( cline%defined('local_shift_bound') ) self%local_shift_bound = cline%get_rarg('local_shift_bound')
        if( trim(pose_policy_arg%to_char()).ne.'local' .and. &
            &any([self%local_ang_bound, self%local_inpl_bound, self%local_shift_bound] >= 0.) )then
            THROW_HARD('local pose-bound overrides require pose_policy=local')
        endif
        if( self%local_ang_bound >= 0. )then
            if( self%local_ang_bound > 180. ) THROW_HARD('local_ang_bound must be between 0 and 180 degrees')
            if( self%local_ang_bound < 0.01 )then
                self%local_nspace_sub = STAGE2_NSPACE
            else
                self%local_nspace_sub = min(STAGE2_NSPACE, max(2, round2even(2. / &
                    &max(1.e-6, 1. - cos(deg2rad(self%local_ang_bound))))))
            endif
        endif
        if( self%local_inpl_bound >= 0. ) call cline%set('prob_athres', self%local_inpl_bound)
        if( self%local_shift_bound >= 0. ) call cline%set('trs', self%local_shift_bound)
        call pose_policy_arg%kill
    end subroutine normalize_inputs

    !> Before params%new: state=0/1 input is initialized by flex_pca; a project that already carries
    !! multi-state labels refines those states instead
    subroutine resolve_states( self, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        self%l_flex_init = .not. self%project_has_multistate_labels(cline)
        if( self%l_flex_init )then
            write(logfhandle,'(A)') '>>> '//self%label()//' STATE INITIALIZATION BY FLEX PCA'
            self%nstates_project = 0
            if( self%l_nstates_on_cline ) self%nstates_project = cline%get_iarg('nstates')
            if( self%nstates_project < 3 )then
                THROW_HARD(self%label()//' state=0/1 input requires nstates >= 3 for flex_pca initialization')
            endif
            ! population floor for the flex states: no under-populated cluster
            ! enters the volume refinement (flex_pca min_state_frac)
            if( .not. cline%defined('min_state_frac') ) call cline%set('min_state_frac', MIN_STATE_FRAC_FLEX)
        else
            call self%set_nstates(cline)
        endif
    end subroutine resolve_states

    !> Before params%new: the default sample size and the planner-driven downscaling
    subroutine resolve_sampling( self, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        if( .not. self%policy%l_full_batch )then
            self%l_nsample_auto = .not. cline%defined('nsample')
            if( .not. self%l_nsample_auto ) self%l_nsample_auto = cline%get_iarg('nsample') <= 0
            if( self%l_nsample_auto )then
                call cline%set('nsample', self%default_nsample())
            endif
        endif
        ! staged frequency marching always uses planner-driven downscaling
        if( cline%defined('autoscale') ) write(logfhandle,'(A)') &
            &'>>> '//self%label()//' IGNORES autoscale: frequency marching owns downscaling'
        call cline%set('autoscale', 'yes')
    end subroutine resolve_sampling

    !> Before params%new: an explicit maxits caps the march
    subroutine read_maxits( self, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        self%l_maxits_defined = cline%defined('maxits')
        if( self%l_maxits_defined )then
            self%maxits_user = cline%get_iarg('maxits')
            if( self%maxits_user < STAGE2_MINITS )then
                THROW_HARD('maxits must be >= '//int2str(STAGE2_MINITS)//' for '//self%label())
            endif
        endif
    end subroutine read_maxits

    subroutine validate( self, params )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        write(logfhandle,'(A,A)') '>>> '//self%label()//' POSE_POLICY: ', trim(params%pose_policy)
        select case(trim(params%multivol_mode))
            case('input_oris_refine')
                ! supported
            case default
                THROW_HARD('Unsupported multivol_mode for '//self%label()//': '//trim(params%multivol_mode))
        end select
        write(logfhandle,'(A,A)') '>>> '//self%label()//' MULTIVOL_MODE: ', trim(params%multivol_mode)
        if( trim(params%combine_eo).ne.'no' )then
            THROW_HARD(self%label()//' does not support combine_eo; it is an LP-set multi-state workflow')
        endif
        select case(trim(params%filt_mode))
            case('fsc', 'nonuniform_lpset', 'none')
                ! supported
            case default
                THROW_HARD(self%label()//' supports filt_mode=fsc|nonuniform_lpset|none')
        end select
        select case(trim(params%prob_neigh_mode))
            case('geom', 'state')
                ! supported
            case default
                THROW_HARD(self%label()//' supports prob_neigh_mode=geom|state')
        end select
    end subroutine validate

    !> state=0/1 input: flex_pca labels the particles and delivers the state maps; flex_mean replaces the
    !! project consensus map flex would otherwise pick up
    subroutine initialize_states( self, params, cline, flex_mean )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        class(string),           optional, intent(in)    :: flex_mean
        if( .not. self%l_flex_init ) return
        call self%run_flex_pca(params, cline, flex_mean)
        call self%set_nstates(cline)
        params%nstates = self%nstates_project
        if( self%l_nsample_auto )then
            params%nsample = self%default_nsample()
            call cline%set('nsample', params%nsample)
        endif
    end subroutine initialize_states

    !> the sample per iteration, the stage cap, the downscaling and the coverage report
    subroutine plan( self, params, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        call self%set_sampling(params, cline)
        call self%set_downscaling(params, cline)
        call self%report_coverage(params)
    end subroutine plan

    !> labelled input: state maps from the project or reconstructed from the labels
    subroutine startup_maps( self, params, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        if( .not. self%l_flex_init ) call self%initialize_state_volumes(params, cline)
    end subroutine startup_maps

    !> iters_done: refinement iterations the caller already ran, so the march numbers its own after them
    subroutine run_march( self, params, cline, xrefine3D, gui_comm, iters_done )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(commander_base),             intent(inout) :: xrefine3D
        class(gui_communicator),           intent(inout) :: gui_comm
        integer,                 optional, intent(in)    :: iters_done
        integer :: nspace_sub, stage_cap, stage2_niters
        call cline%set('prg', 'refine3D')
        self%total_iter = 0
        if( present(iters_done) ) self%total_iter = iters_done
        call cline%set('maxits_glob', max(1, self%stage_cap))
        ! plain copies: the march updates self, so none of its components is passed as an argument
        nspace_sub = self%local_nspace_sub
        stage_cap  = self%stage_cap
        call self%run_frequency_march(params, cline, xrefine3D, gui_comm, 'prob_neigh', STAGE2_NSPACE, &
            &nspace_sub, STAGE2_MINITS, stage2_niters, stage_cap, params%overlap)
        write(logfhandle,'(A,I0,A,I0)') '>>> '//self%label()//' STAGE ITERATIONS PROB_NEIGH/TOTAL: ', &
            &stage2_niters, '/', self%total_iter
    end subroutine run_march

    subroutine ensure_all_updated( self, params, cline, xrefine3D )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(commander_base),             intent(inout) :: xrefine3D
        integer :: nactive, nupdated, nmissing
        call self%read_update_coverage(params, nactive, nupdated, nmissing)
        if( nactive < 1 )then
            THROW_HARD(self%label()//' has no active particles after staged refinement')
        endif
        if( nmissing > 0 )then
            call self%run_missing_update(params, cline, xrefine3D, nmissing, nactive)
            call self%read_update_coverage(params, nactive, nupdated, nmissing)
            if( nmissing > 0 )then
                THROW_HARD(self%label()//' final missing-update pass failed to update every active particle')
            endif
        endif
    end subroutine ensure_all_updated

    subroutine kill( self )
        class(refine3D_states_controller), intent(inout) :: self
        integer :: state
        if( allocated(self%init_vols) )then
            do state = 1,size(self%init_vols)
                call self%init_vols(state)%kill
            enddo
            deallocate(self%init_vols)
        endif
        if( allocated(self%clssmp_units) ) call deallocate_class_samples(self%clssmp_units)
        self%policy             = refine3D_states_policy()
        self%lpinfo_multi       = lp_crop_inf()
        self%nstates_project    = 0
        self%nptcls_eff         = 0
        self%nptcls_per_iter    = 0
        self%local_nspace_sub   = STAGE2_NSPACE_SUB
        self%maxits_user        = 0
        self%stage_cap          = 0
        self%total_iter         = 0
        self%gui_stage          = 0
        self%local_ang_bound    = -1.
        self%local_inpl_bound   = -1.
        self%local_shift_bound  = -1.
        self%l_maxits_defined   = .false.
        self%l_nstates_on_cline = .false.
        self%l_flex_init        = .false.
        self%l_nsample_auto     = .false.
    end subroutine kill

    ! ---------- private ----------

    function label( self ) result( str )
        class(refine3D_states_controller), intent(in) :: self
        character(len=:), allocatable :: str
        str = trim(self%policy%label)
    end function label

    integer function default_nsample( self ) result( nsample_default )
        class(refine3D_states_controller), intent(in) :: self
        nsample_default = min(NSAMPLE_REFINE3D_STATES_CAP, &
            &NSAMPLE_PER_STATE_REFINE3D_STATES * self%nstates_project)
    end function default_nsample

    subroutine reject_input_volumes( self, cline )
        use simple_defs, only: MAXS
        class(refine3D_states_controller), intent(in)    :: self
        class(cmdline),                    intent(inout) :: cline
        integer :: state
        do state = 1,MAXS
            if( cline%defined('vol'//int2str(state)) )then
                THROW_HARD(self%label()//' does not accept vol1..volN inputs; use classify3D_refs for supplied references')
            endif
        enddo
    end subroutine reject_input_volumes

    !> whether the input project already carries populated multi-state labels
    logical function project_has_multistate_labels( self, cline ) result( l_multi )
        class(refine3D_states_controller), intent(in)    :: self
        class(cmdline),                    intent(inout) :: cline
        type(sp_project) :: proj
        type(string)     :: projfile
        integer :: nactive
        l_multi = .false.
        if( .not. cline%defined('projfile') ) THROW_HARD('projfile is required for '//self%label())
        projfile = cline%get_carg('projfile')
        call proj%read_segment('ptcl3D', projfile)
        if( proj%os_ptcl3D%isthere('state') )then
            nactive = proj%os_ptcl3D%count_state_gt_zero()
            if( nactive > 0 ) l_multi = proj%os_ptcl3D%get_n('state') > 1
        endif
        call proj%kill
        call projfile%kill
    end function project_has_multistate_labels

    subroutine set_nstates( self, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(sp_project) :: state_proj
        type(string)     :: projfile
        integer, allocatable :: pops(:)
        integer :: state, nactive_labels, nstates_cline, nstates_labels
        logical :: l_has_state_labels, l_has_project_multistates
        if( .not. cline%defined('projfile') )then
            THROW_HARD('projfile is required for '//self%label())
        endif
        projfile = cline%get_carg('projfile')
        call state_proj%read_segment('ptcl3D', projfile)
        nactive_labels    = 0
        nstates_labels    = 1
        l_has_state_labels = state_proj%os_ptcl3D%isthere('state')
        if( l_has_state_labels )then
            nactive_labels = state_proj%os_ptcl3D%count_state_gt_zero()
            if( nactive_labels > 0 ) nstates_labels = state_proj%os_ptcl3D%get_n('state')
        endif
        l_has_project_multistates = nactive_labels > 0 .and. nstates_labels > 1
        if( l_has_project_multistates )then
            self%nstates_project = nstates_labels
            if( self%l_nstates_on_cline )then
                nstates_cline = cline%get_iarg('nstates')
                if( nstates_cline /= self%nstates_project )then
                    THROW_HARD('command-line nstates does not match project state labels for '//self%label())
                endif
            endif
            call state_proj%os_ptcl3D%get_pops(pops, 'state', maxn=self%nstates_project)
            do state = 1,self%nstates_project
                if( pops(state) < 1 )then
                    write(logfhandle,*) 'state, population: ', state, pops(state)
                    THROW_HARD(self%label()//' requires every state label to have at least one active particle')
                endif
            enddo
            write(logfhandle,'(A,I0)') '>>> '//self%label()//' NSTATES FROM PROJECT: ', self%nstates_project
        else
            ! state=0/1 input reaches this point only labelled by flex_pca
            THROW_HARD(self%label()//' found no multi-state labels in the project')
        endif
        call cline%set('nstates', self%nstates_project)
        if( allocated(pops) ) deallocate(pops)
        call state_proj%kill
        call projfile%kill
    end subroutine set_nstates

    subroutine set_sampling( self, params, cline )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        type(sp_project) :: sampling_proj
        integer :: maxits_auto, nsample_target, sweep_iters
        real    :: update_frac_auto
        nsample_target = params%nsample
        if( .not. self%policy%l_full_batch .and. nsample_target < 1 ) THROW_HARD('nsample must be >= 1 for '//self%label())
        call sampling_proj%read(params%projfile)
        ! the ptcl3D states are this workflow's own (flex drops particles to state 0 there) and the
        ! sampling units leave those rows out; the project-level count reads the ptcl2D states
        self%nptcls_eff = sampling_proj%os_ptcl3D%count_state_gt_zero()
        call sampling_proj%kill
        if( self%nptcls_eff < 1 ) THROW_HARD('no active particles available for '//self%label())
        if( self%policy%l_full_batch ) nsample_target = self%nptcls_eff
        self%nptcls_per_iter = min(self%nptcls_eff, nsample_target)
        if( self%nptcls_eff <= nsample_target )then
            params%update_frac   = 1.0
            params%l_update_frac = .false.
            params%l_trail_rec   = .false.
            call cline%delete('update_frac')
            write(logfhandle,'(A,I0,A,I0,A)') &
                &'>>> '//self%label()//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                &self%nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
        else
            update_frac_auto = real(nsample_target) / real(self%nptcls_eff)
            if( update_frac_auto <= 0.99 )then
                params%update_frac   = update_frac_auto
                params%l_update_frac = .true.
                params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                call cline%set('update_frac', update_frac_auto)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//self%label()//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                    &self%nptcls_eff, '/', nsample_target, '/', update_frac_auto
            else
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)') &
                    &'>>> '//self%label()//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &self%nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
            endif
        endif
        ! draws for one full sweep of the active particles, from the sampling-unit table: the
        ! sampler draws lowest-updatecnt particles first inside every unit's quota
        if( params%l_update_frac )then
            sweep_iters = self%prepare_sampling_units(params, cline)
        else
            sweep_iters = 1
        endif
        if( self%l_maxits_defined )then
            self%stage_cap = self%maxits_user
            params%maxits = self%stage_cap
            write(logfhandle,'(A,I0)') &
                &'>>> '//self%label()//' STAGE MAXITS COMMAND-LINE OVERRIDE: ', self%stage_cap
        else
            maxits_auto = ceiling(TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES * real(sweep_iters))
            self%stage_cap = max(MINITS_REFINE3D_STATES, min(MAXITS_REFINE3D_STATES_CAP, max(STAGE2_MINITS, maxits_auto)))
            params%maxits = self%stage_cap
            call cline%set('maxits', params%maxits)
            write(logfhandle,'(A,I0,A,F5.1,A,I0,A)') '>>> '//self%label()//' STAGE MAXITS: ', &
                &self%stage_cap, ' FOR ~', TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES, ' UPDATES/PARTICLE (SWEEP ', &
                &sweep_iters, ' DRAWS)'
        endif
    end subroutine set_sampling

    !> The sampling units of balance=cavg (none without selected class averages) on CLASS_SAMPLING_FILE;
    !! returns the draws one full sweep of the active particles takes
    integer function prepare_sampling_units( self, params, cline ) result( sweep )
        use simple_class_sample_io,         only: write_class_samples
        use simple_defs_fname,              only: CLASS_SAMPLING_FILE
        use simple_oris,                    only: class_sample_sweep
        use simple_view_partition_sampling, only: make_class_samples, has_selected_cavgs
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        type(sp_project) :: sampling_proj
        sweep = max(1, ceiling(real(self%nptcls_eff) / real(max(1, self%nptcls_per_iter))))
        call sampling_proj%read(params%projfile)
        if( trim(params%balance) /= 'none' .and. .not. has_selected_cavgs(sampling_proj) )then
            write(logfhandle,'(A)') '>>> '//self%label()//' NO SELECTED CLASS AVERAGES: balance=none'
            params%balance = 'none'
        endif
        call cline%set('balance', params%balance)
        if( trim(params%balance) /= 'none' )then
            call make_class_samples(params, sampling_proj, params%balance, params%nthr, .true., self%clssmp_units)
            call write_class_samples(self%clssmp_units, string(CLASS_SAMPLING_FILE))
            sweep = class_sample_sweep(self%clssmp_units, self%nptcls_per_iter)
        endif
        call sampling_proj%kill
    end function prepare_sampling_units

    !> the unit table against the draws planned: the frequency blocks of the prob_neigh march, one
    !! cohort each
    subroutine report_coverage( self, params )
        use simple_view_partition_sampling, only: report_class_sample_coverage
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        if( .not. allocated(self%clssmp_units) ) return
        call report_class_sample_coverage(self%clssmp_units, params%balance, &
            &self%nptcls_per_iter, ceiling(real(self%stage_cap) / real(FREQUENCY_BLOCK_NITS)), 'frequency block')
        call deallocate_class_samples(self%clssmp_units)
    end subroutine report_coverage

    subroutine set_downscaling( self, params, cline )
        use simple_estimate_ssnr, only: lpstages_setlims
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        call lpstages_setlims(params%box, 2, params%smpd, params%lpstart, params%lpstop, self%lpinfo_multi)
        if( self%local_shift_bound >= 0. .and. trim(params%pose_policy).eq.'local' )then
            params%trs = self%local_shift_bound
            call cline%set('trs', params%trs)
        else
            params%trs = self%lpinfo_multi(2)%trslim
            call cline%set('trs', params%trs)
        endif
        if( self%lpinfo_multi(2)%l_autoscale )then
            params%box_crop  = self%lpinfo_multi(2)%box_crop
            params%smpd_crop = self%lpinfo_multi(2)%smpd_crop
            call cline%set('box_crop', params%box_crop)
            call cline%set('smpd_crop', params%smpd_crop)
            write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                &'>>> '//self%label()//' DOWNSCALING BOX/SMPD_CROP: ', &
                &params%box, '/', params%box_crop, '/', params%smpd_crop
        else
            call cline%delete('box_crop')
            call cline%delete('smpd_crop')
            write(logfhandle,'(A)') '>>> '//self%label()//' DOWNSCALING: native sampling retained'
        endif
    end subroutine set_downscaling

    subroutine initialize_state_volumes( self, params, cline )
        use simple_commanders_rec,  only: commander_rec3D
        use simple_refine3D_fnames, only: refine3D_state_vol_fname
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(cmdline)         :: cline_rec3D
        integer :: state
        allocate(self%init_vols(self%nstates_project))
        if( self%project_state_volumes_compatible(params) )then
            do state = 1,self%nstates_project
                call cline%set('vol'//int2str(state), self%init_vols(state))
            enddo
            write(logfhandle,'(A)') '>>> '//self%label()//' USING PROJECT STATE VOLUMES'
            return
        endif
        cline_rec3D = cline
        call cline_rec3D%set('prg', 'reconstruct3D')
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('refine')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%delete('ufrac_trec')
        call cline_rec3D%delete('endit')
        call cline_rec3D%delete('box_crop')
        call cline_rec3D%delete('smpd_crop')
        call cline_rec3D%set('objfun', 'cc')
        call cline_rec3D%set('postprocess', 'no')
        call xrec3D%execute(cline_rec3D)
        do state = 1,self%nstates_project
            call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
        enddo
        write(logfhandle,'(A)') '>>> '//self%label()//' INITIALIZED STATE VOLUMES BY RECONSTRUCTION'
        call cline_rec3D%kill
    end subroutine initialize_state_volumes

    logical function project_state_volumes_compatible( self, params ) result( l_compatible )
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        type(sp_project) :: spproj
        real    :: init_smpd
        integer :: state, init_box
        l_compatible = .false.
        call spproj%read_segment('out', params%projfile)
        do state = 1,self%nstates_project
            if( .not. spproj%isthere_in_osout('vol', state) )then
                call spproj%kill
                return
            endif
            call spproj%get_vol('vol', state, self%init_vols(state), init_smpd, init_box)
            if( .not. file_exists(self%init_vols(state)) )then
                call spproj%kill
                return
            endif
            if( init_smpd <= 0. )then
                call spproj%kill
                return
            endif
            ! any downscaled sampling of the native grid is acceptable: base
            ! refine3D rescales references to the stage crop
            if( init_box > params%box .or. &
                &abs(real(init_box)*init_smpd - real(params%box)*params%smpd) > &
                &0.01 * real(params%box)*params%smpd )then
                call spproj%kill
                return
            endif
        enddo
        l_compatible = .true.
        call spproj%kill
    end function project_state_volumes_compatible

    subroutine run_flex_pca( self, params, cline, flex_mean )
        use simple_commanders_flex_pca, only: commander_flex_pca
        use simple_state_weight_set,    only: discard_state_weight_set
        use simple_syslib,              only: simple_abspath
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(string),           optional, intent(in)    :: flex_mean
        type(commander_flex_pca) :: xflex
        type(cmdline)            :: cline_flex
        type(sp_project)         :: flex_proj
        integer :: state, flex_box, nstates_requested, nstates_flex
        integer :: nactive_labels, nstates_labels
        real    :: flex_smpd
        ! validate inputs
        nstates_requested = params%nstates
        if( nstates_requested < 3 ) THROW_HARD(self%label()//' flex_pca initialization requires nstates >= 3')
        ! validate project
        call flex_proj%read_segment('ptcl3D', params%projfile)
        nactive_labels = 0
        nstates_labels = 1
        if( flex_proj%os_ptcl3D%isthere('state') )then
            nactive_labels = flex_proj%os_ptcl3D%count_state_gt_zero()
            if( nactive_labels > 0 ) nstates_labels = flex_proj%os_ptcl3D%get_n('state')
        endif
        call flex_proj%kill
        if( nactive_labels > 0 .and. nstates_labels > 1 )then
            THROW_HARD(self%label()//' flex_pca initialization requires an input project with a single state')
        endif
        ! execution prg=flex_pca; without flex_mean it picks up the project consensus map
        ! (out segment, state 1) itself and validates it at native sampling
        cline_flex = cline
        call cline_flex%set('prg',        'flex_pca')
        call cline_flex%set('mkdir',      'no')
        call cline_flex%set('npreimages', nstates_requested)
        call cline_flex%set('rec_backend', trim(self%policy%flex_rec_backend))
        call cline_flex%delete('nstates')
        if( len_trim(self%policy%flex_rec_states_backend) > 0 ) &
            &call cline_flex%set('rec_states_backend', trim(self%policy%flex_rec_states_backend))
        if( self%policy%flex_min_neff > 0 ) call cline_flex%set('min_neff', self%policy%flex_min_neff)
        if( self%policy%l_flex_native_box ) call cline_flex%set('box_crop', params%box)
        if( self%policy%l_flex_lp_lpstart ) call cline_flex%set('lp', params%lpstart)
        if( self%policy%l_flex_shmem )then
            call cline_flex%delete('nparts')
            call cline_flex%delete('numlen')
        endif
        if( present(flex_mean) ) call cline_flex%set('vol1', simple_abspath(flex_mean))
        call xflex%execute(cline_flex)
        ! output parsing
        call flex_proj%read_segment('ptcl3D', params%projfile)
        if( .not. flex_proj%os_ptcl3D%isthere('state') )then
            THROW_HARD(self%label()//' flex_pca did not write particle state labels')
        endif
        nstates_flex = flex_proj%os_ptcl3D%get_n('state')
        call flex_proj%kill
        if( nstates_flex <= 1 ) THROW_HARD(self%label()//' flex_pca produced fewer than two states')
        call cline%set('nstates', nstates_flex)
        allocate(self%init_vols(nstates_flex))
        call flex_proj%read_segment('out', params%projfile)
        do state = 1,nstates_flex
            if( .not. flex_proj%isthere_in_osout('vol', state) )then
                THROW_HARD(self%label()//' flex_pca did not write every state volume')
            endif
            call flex_proj%get_vol('vol', state, self%init_vols(state), flex_smpd, flex_box)
            if( .not. file_exists(self%init_vols(state)) )then
                THROW_HARD(self%label()//' flex_pca state volume does not exist')
            endif
            if( flex_box /= params%box .or. flex_smpd <= 0. .or. abs(flex_smpd - params%smpd) > 1.e-6 )then
                THROW_HARD(self%label()//' flex_pca state volumes must match the project particle sampling')
            endif
            ! update command line with flex volumes
            self%init_vols(state) = simple_abspath(self%init_vols(state))
            call cline%set('vol'//int2str(state), self%init_vols(state))
        enddo
        call cline%delete('neigs')
        ! the flex solution enters as hard labels only: its state weight set is not kept
        call discard_state_weight_set(flex_proj, params%projfile)
        call flex_proj%kill
        write(logfhandle,'(A,I0)') '>>> '//self%label()//' FLEX_PCA INITIALIZED NSTATES: ', nstates_flex
        call cline_flex%kill
    end subroutine run_flex_pca

    subroutine run_frequency_march( self, params, cline, xrefine3D, gui_comm, refine_mode, nspace_stage, &
        &nspace_sub_stage, min_total_iters, niters, max_total_iters, overlap_target )
        use simple_refine3D_stage_plan, only: refine3D_stage_plan_entry, plan_refine3D_frequency_stages
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(commander_base),             intent(inout) :: xrefine3D
        class(gui_communicator),           intent(inout) :: gui_comm
        character(len=*),                  intent(in)    :: refine_mode
        integer,                           intent(in)    :: nspace_stage, nspace_sub_stage, min_total_iters, max_total_iters
        integer,                           intent(out)   :: niters
        real,                              intent(in)    :: overlap_target
        type(refine3D_stage_plan_entry), allocatable :: stages(:)
        integer :: block, block_niters, block_minits
        niters = 0
        if( max_total_iters < min_total_iters )then
            THROW_HARD(self%label()//' frequency-stage cap is below its required minimum')
        endif
        call plan_refine3D_frequency_stages(params%box, params%smpd, params%lpstart, params%lpstop, &
            &max_total_iters, FREQUENCY_BLOCK_NITS, self%total_iter + 1, self%lpinfo_multi(2), stages)
        write(logfhandle,'(A,I0,A,F7.2,A,F7.2)') &
            &'>>> '//self%label()//' FREQUENCY MARCHING NSTAGES/LPSTART/LPSTOP: ', &
            &size(stages), '/', params%lpstart, '/', params%lpstop
        do block = 1,size(stages)
            call cline%set('trs',    stages(block)%lpinfo%trslim)
            if( self%local_shift_bound >= 0. .and. trim(params%pose_policy).eq.'local' )then
                call cline%set('trs', self%local_shift_bound)
            endif
            call cline%set('lp',     stages(block)%lpinfo%lp)
            call cline%set('lpstop', stages(block)%lpinfo%lp)
            if( stages(block)%lpinfo%l_autoscale )then
                call cline%set('box_crop',  stages(block)%lpinfo%box_crop)
                call cline%set('smpd_crop', stages(block)%lpinfo%smpd_crop)
            else
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
            endif
            ! a block may stop early on the state-overlap target once the march-wide minimum is met;
            ! every block is still visited, so the low-pass schedule completes
            block_minits = max(1, min(stages(block)%nits, min_total_iters - niters))
            call self%run_stage(params, cline, xrefine3D, gui_comm, block, refine_mode, nspace_stage, &
                &nspace_sub_stage, block_minits, block_niters, stages(block)%nits, overlap_target)
            niters = niters + block_niters
        enddo
        call cline%delete('cohort_sampling')
        deallocate(stages)
    end subroutine run_frequency_march

    subroutine run_stage( self, params, cline, xrefine3D, gui_comm, stage, refine_mode, nspace_stage, &
        &nspace_sub_stage, min_stage_iters, niters, max_stage_iters, stage_overlap_target )
        use simple_defs_fname,      only: LP_SUFFIX, MRC_EXT
        use simple_fileio,          only: add2fbody
        use simple_refine3D_fnames, only: refine3D_state_vol_fname
        use simple_solve3D_utils,   only: gen_ortho_reprojs4viz, write_solve3D_lowpass_snapshot
        use simple_string_utils,    only: int2str_pad
        class(refine3D_states_controller), intent(inout) :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(commander_base),             intent(inout) :: xrefine3D
        class(gui_communicator),           intent(inout) :: gui_comm
        integer,                           intent(in)    :: stage, nspace_stage, nspace_sub_stage, min_stage_iters
        character(len=*),                  intent(in)    :: refine_mode
        integer,                           intent(out)   :: niters
        integer,                           intent(in)    :: max_stage_iters
        real,                              intent(in)    :: stage_overlap_target
        type(sp_project) :: gui_proj
        type(string)     :: gui_stage_suffix, gui_vol, gui_stage_vol, gui_stage_lp
        integer :: stage_start, stage_limit, gui_box, state
        real    :: gui_smpd, gui_lp, state_overlap
        logical :: l_gui_comm_active
        l_gui_comm_active = self%policy%l_gui_stages .and. params%niceserver .ne. "" .and. params%niceprocid > 0
        niters            = 0
        state_overlap     = 0.
        stage_limit       = max_stage_iters
        stage_start       = self%total_iter + 1
        write(logfhandle,'(A,I0,A,A,A,I0,A,I0,A,F7.4)') '>>> '//self%label()//' ENTERING STAGE ', stage, &
            &' REFINE=', trim(refine_mode), ' NSPACE/NSPACE_SUB=', nspace_stage, '/', nspace_sub_stage, &
            &' OVERLAP_TARGET=', stage_overlap_target
        call cline%set('refine', refine_mode)
        call cline%set('nspace', nspace_stage)
        ! a prob_neigh frequency block refines one particle cohort, drawn at its first iteration
        if( trim(refine_mode).eq.'prob_neigh' .and. params%l_update_frac )then
            call cline%set('cohort_sampling', 'yes')
        else
            call cline%delete('cohort_sampling')
        endif
        if( nspace_sub_stage > 0 )then
            call cline%set('nspace_sub', nspace_sub_stage)
        else
            call cline%delete('nspace_sub')
        endif
        call cline%set('maxits', stage_limit)
        call cline%set('minits', min_stage_iters)
        call cline%set('startit', stage_start)
        call cline%set('which_iter', stage_start)
        call cline%set('extr_iter', stage_start)
        call cline%set('overlap', stage_overlap_target)
        call cline%delete('endit')
        call xrefine3D%execute(cline)
        if( cline%defined('endit') )then
            self%total_iter = nint(cline%get_rarg('endit'))
        else
            self%total_iter = stage_start + stage_limit - 1
        endif
        niters = self%total_iter - stage_start + 1
        state_overlap = self%read_state_overlap(params)
        write(logfhandle,'(A,I0,A,I0,A,F7.4,A,F7.4)') '>>> '//self%label()//' STAGE ', stage, &
            &' NITERS ', niters, ' FINAL STATE_OVERLAP: ', state_overlap, ' TARGET: ', stage_overlap_target
        if( niters >= stage_limit .and. state_overlap < stage_overlap_target )then
            write(logfhandle,'(A,I0,A,F7.4)') '>>> '//self%label()//' STAGE ', stage, &
                &' REACHED STAGE CAP BEFORE STATE_OVERLAP TARGET: ', state_overlap
        endif
        if( l_gui_comm_active )then
            ! Freeze the current stage map before the next block overwrites
            ! recvol_stateNN.mrc, then publish matching project metadata.
            self%gui_stage = self%gui_stage + 1
            gui_box   = params%box
            gui_smpd  = params%smpd
            gui_lp    = cline%get_rarg('lp')
            if( cline%defined('box_crop')  ) gui_box  = cline%get_iarg('box_crop')
            if( cline%defined('smpd_crop') ) gui_smpd = cline%get_rarg('smpd_crop')
            gui_stage_suffix = '_stage'//int2str_pad(self%gui_stage,2)
            do state = 1,params%nstates
                gui_vol       = refine3D_state_vol_fname(state)
                if( .not. file_exists(gui_vol) ) cycle
                gui_stage_vol = add2fbody(gui_vol, string(MRC_EXT), gui_stage_suffix)
                gui_stage_lp  = add2fbody(gui_stage_vol, MRC_EXT, LP_SUFFIX)
                call write_solve3D_lowpass_snapshot(gui_vol, gui_lp, gui_stage_lp, gui_smpd, box=gui_box)
            enddo
            call gui_proj%read(params%projfile)
            call gen_ortho_reprojs4viz(params, gui_proj)
            call gui_comm%add_metadata(gui_proj, oritype='cls3D', stage=self%gui_stage)
            call gui_proj%kill
            call gui_stage_suffix%kill
            call gui_vol%kill
            call gui_stage_vol%kill
            call gui_stage_lp%kill
        endif
    end subroutine run_stage

    subroutine read_update_coverage( self, params, nactive, nupdated, nmissing )
        class(refine3D_states_controller), intent(in)  :: self
        class(parameters),                 intent(in)  :: params
        integer,                           intent(out) :: nactive, nupdated, nmissing
        type(sp_project) :: update_proj
        integer, allocatable :: states(:), updatecnts(:)
        call update_proj%read_segment('ptcl3D', params%projfile)
        if( .not. update_proj%os_ptcl3D%isthere('updatecnt') )then
            call update_proj%kill
            THROW_HARD(self%label()//' cannot finish before active particles are updated')
        endif
        states     = update_proj%os_ptcl3D%get_all_asint('state')
        updatecnts = update_proj%os_ptcl3D%get_all_asint('updatecnt')
        nactive    = count(states > 0)
        nupdated   = count(states > 0 .and. updatecnts > 0)
        nmissing   = nactive - nupdated
        write(logfhandle,'(A,I0,A,I0,A,I0)') &
            &'>>> '//self%label()//' ACTIVE PARTICLE UPDATE COVERAGE UPDATED/ACTIVE/MISSING: ', &
            &nupdated, '/', nactive, '/', nmissing
        if( allocated(states)     ) deallocate(states)
        if( allocated(updatecnts) ) deallocate(updatecnts)
        call update_proj%kill
    end subroutine read_update_coverage

    subroutine run_missing_update( self, params, cline, xrefine3D, nmissing, nactive )
        use simple_defs_fname, only: ASSIGNMENT_FBODY, DIST_FBODY
        use simple_fileio,     only: del_files
        use simple_syslib,     only: del_file
        class(refine3D_states_controller), intent(in)    :: self
        class(parameters),                 intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        class(commander_base),             intent(inout) :: xrefine3D
        integer,                           intent(in)    :: nmissing, nactive
        type(cmdline) :: cline_missing
        integer       :: iter_missing
        iter_missing = self%next_iteration(cline)
        write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
            &'>>> '//self%label()//' FINAL MISSING-UPDATE ASSIGNMENT', &
            &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_missing
        call flush(logfhandle)
        cline_missing = cline
        call cline_missing%set('prg',           'refine3D')
        call cline_missing%set('mkdir',              'no')
        call cline_missing%set('balance',          'none')
        call cline_missing%set('frac_best',           1.0)
        call cline_missing%set('fillin',             'no')
        call cline_missing%set('update_frac',         1.0)
        call cline_missing%set('trail_rec',          'no')
        call cline_missing%set('volrec',             'no')
        call cline_missing%set('maxits',                1)
        call cline_missing%set('startit',    iter_missing)
        call cline_missing%set('which_iter', iter_missing)
        call cline_missing%set('extr_iter',  iter_missing)
        call cline_missing%delete('endit')
        call cline_missing%delete('cohort_sampling')
        call cline_missing%set('refine',          'greedy')
        call cline_missing%set('greedy_sampling',   'yes')
        call cline_missing%set('update_missing',    'yes')
        call xrefine3D%execute(cline_missing)
        call del_files(DIST_FBODY,       params%nparts, ext='.dat')
        call del_files(ASSIGNMENT_FBODY, params%nparts, ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        call cline_missing%kill
    end subroutine run_missing_update

    integer function next_iteration( self, cline ) result( iter )
        class(refine3D_states_controller), intent(in)    :: self
        class(cmdline),                    intent(inout) :: cline
        iter = 1
        if( cline%defined('endit') )then
            iter = cline%get_iarg('endit') + 1
        else if( cline%defined('which_iter') )then
            iter = cline%get_iarg('which_iter') + 1
        endif
        iter = max(1, iter)
    end function next_iteration

    real function read_state_overlap( self, params ) result( overlap )
        class(refine3D_states_controller), intent(in) :: self
        class(parameters),                 intent(in) :: params
        type(sp_project) :: state_proj
        real,    allocatable :: states(:), mi_state(:), sampled(:)
        logical, allocatable :: mask(:)
        real :: sampled_lb
        overlap = 0.
        call state_proj%read_segment('ptcl3D', params%projfile)
        if( .not. state_proj%os_ptcl3D%isthere('mi_state') )then
            call state_proj%kill
            return
        endif
        states   = state_proj%os_ptcl3D%get_all('state')
        mi_state = state_proj%os_ptcl3D%get_all('mi_state')
        if( state_proj%os_ptcl3D%isthere('sampled') )then
            sampled = state_proj%os_ptcl3D%get_all('sampled')
            sampled_lb = maxval(sampled) - 0.5
            allocate(mask(size(states)), source=sampled > sampled_lb .and. states > 0.5)
        else
            allocate(mask(size(states)), source=states > 0.5)
        endif
        if( count(mask) > 0 ) overlap = sum(mi_state, mask=mask) / real(count(mask))
        if( allocated(states)   ) deallocate(states)
        if( allocated(mi_state) ) deallocate(mi_state)
        if( allocated(sampled)  ) deallocate(sampled)
        if( allocated(mask)     ) deallocate(mask)
        call state_proj%kill
    end function read_state_overlap

end module simple_refine3D_states_controller
