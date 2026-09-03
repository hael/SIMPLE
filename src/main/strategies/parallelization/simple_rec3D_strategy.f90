module simple_rec3D_strategy
use simple_core_module_api
use simple_builder,              only: builder
use simple_parameters,           only: parameters
use simple_cmdline,              only: cmdline
use simple_qsys_env,             only: qsys_env
use simple_matcher_3Drec,        only: calc_3Drec, calc_projdir3Drec
use simple_commanders_rec_distr, only: commander_volassemble, filter_pcg_nonuniform_maps
use simple_refine3D_fnames,      only: refine3D_fsc_fname, refine3D_state_vol_fname, &
    &refine3D_pcg_raw_accum_fname
use simple_rec3D_pcg_strategy,   only: execute_rec3D_pcg_shared, execute_rec3D_pcg_distributed_master
use simple_sigma2_files,         only: load_sigma2_groups
implicit none

public :: rec3D_strategy, rec3D_inmem_strategy, rec3D_distr_strategy, create_rec3D_strategy
public :: rec3D_pcg_inmem_strategy
public :: rec3D_backend_id, rec3D_backend_is_wired
public :: REC3D_BACKEND_INVALID, REC3D_BACKEND_GRIDDING, REC3D_BACKEND_PCG
private
#include "simple_local_flags.inc"

integer, parameter :: REC3D_BACKEND_INVALID  = 0
integer, parameter :: REC3D_BACKEND_GRIDDING = 1
integer, parameter :: REC3D_BACKEND_PCG      = 2

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: rec3D_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type rec3D_strategy

! Shared-memory
type, extends(rec3D_strategy) :: rec3D_inmem_strategy
contains
    procedure :: initialize   => inmem_initialize
    procedure :: execute      => inmem_execute
    procedure :: finalize_run => inmem_finalize_run
    procedure :: cleanup      => inmem_cleanup
end type rec3D_inmem_strategy

! Shared-memory kernel PCG. Initialization, final postprocessing and cleanup
! reuse the established in-memory reconstruct3D lifecycle; only execution is
! replaced because PCG produces dense halfmaps rather than gridding partials.
type, extends(rec3D_inmem_strategy) :: rec3D_pcg_inmem_strategy
contains
    procedure :: execute => pcg_inmem_execute
end type rec3D_pcg_inmem_strategy

! Distributed-memory
type, extends(rec3D_strategy) :: rec3D_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    integer        :: nthr_master = 1
contains
    procedure :: initialize   => distr_initialize
    procedure :: execute      => distr_execute
    procedure :: finalize_run => distr_finalize_run
    procedure :: cleanup      => distr_cleanup
end type rec3D_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(in)    :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, build, cline)
        import :: rec3D_strategy, parameters, builder, cmdline
        class(rec3D_strategy), intent(inout) :: self
        type(parameters),      intent(in)    :: params
        type(builder),         intent(inout) :: build
        class(cmdline),        intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    pure integer function rec3D_backend_id(name) result(backend_id)
        character(len=*), intent(in) :: name
        select case(trim(name))
            case('gridding')
                backend_id = REC3D_BACKEND_GRIDDING
            case('pcg')
                backend_id = REC3D_BACKEND_PCG
            case DEFAULT
                backend_id = REC3D_BACKEND_INVALID
        end select
    end function rec3D_backend_id

    pure logical function rec3D_backend_is_wired(backend_id) result(l_wired)
        integer, intent(in) :: backend_id
        l_wired = backend_id == REC3D_BACKEND_GRIDDING .or. backend_id == REC3D_BACKEND_PCG
    end function rec3D_backend_is_wired

    ! --------------------------------------------------------------------
    ! Strategy selection
    ! --------------------------------------------------------------------

    function create_rec3D_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(rec3D_strategy), allocatable :: strategy
        type(string) :: backend
        backend = string('gridding')
        if( cline%defined('rec_backend') ) backend = cline%get_carg('rec_backend')
        select case(rec3D_backend_id(backend%to_char()))
            case(REC3D_BACKEND_GRIDDING)
                ! Distributed master iff: nparts defined AND part not defined.
                ! Keep this branch identical to the pre-selector strategy choice.
                if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
                    allocate(rec3D_distr_strategy :: strategy)
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                        &'>>> DISTRIBUTED-MEMORY REC3D EXECUTION (gridding)'
                else
                    allocate(rec3D_inmem_strategy :: strategy)
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                        &'>>> SHARED-MEMORY REC3D EXECUTION (gridding)'
                endif
            case(REC3D_BACKEND_PCG)
                if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
                    if( cline%get_iarg('nparts') > 1 )then
                        allocate(rec3D_distr_strategy :: strategy)
                        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') &
                            &'>>> DISTRIBUTED-MEMORY REC3D EXECUTION (kernel PCG)'
                    else
                        allocate(rec3D_pcg_inmem_strategy :: strategy)
                    endif
                else
                    allocate(rec3D_pcg_inmem_strategy :: strategy)
                endif
                if( L_VERBOSE_GLOB .and. .not. cline%defined('nparts') )then
                    write(logfhandle,'(A)') '>>> SHARED-MEMORY REC3D EXECUTION (kernel PCG)'
                endif
            case DEFAULT
                THROW_HARD('rec_backend must be gridding or pcg')
        end select
        call backend%kill
    end function create_rec3D_strategy

    ! =====================================================================
    ! SHARED-MEMORY IMPLEMENTATION
    ! =====================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),                   intent(inout) :: params
        type(builder),                      intent(inout) :: build
        class(cmdline),                     intent(inout) :: cline
        call build%init_params_and_build_general_tbox(cline, params)
        call sync_resolved_rec_params(params, cline)
        call build%build_strategy3D_tbox(params)
        ! Even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! Update eo flags in project
        call build%spproj%write_segment_inside(params%oritype)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, build, cline)
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        type(commander_volassemble) :: xvolassemble
        type(cmdline)               :: cline_volassemble
        type(string)                :: volname, vol_in
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update, state
        logical              :: l_sigma_loaded
        ! Sampling
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            ! Sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        ! Sigma weighting belongs to the Euclidean data objective. ML
        ! regularization is a separate FSC/SSNR prior applied by volassemble.
        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, &
                &cline, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('gridding objfun=euclid requires sigma2 files')
        endif
        ! Legacy handshake for rec-writing helpers that still inspect this key.
        ! The strategy owns the actual assembly dispatch decision.
        call cline%set('force_volassemble', 'yes')
        if( trim(params%projrec) == 'yes' )then
            call calc_projdir3Drec(params, build, cline, nptcls2update, pinds)
        else
            call calc_3Drec(params, build, cline, nptcls2update, pinds)
        endif
        cline_volassemble = cline
        call cline_volassemble%set('prg',  'volassemble')
        call cline_volassemble%set('nthr', params%nthr)
        call cline_volassemble%delete('objfun_den')
        call cline_volassemble%delete('objfun_den_w')
        do state = 1, params%nstates
            volname = refine3D_state_vol_fname(state)
            if( cline_volassemble%defined('vol'//int2str(state)) )then
                vol_in = cline_volassemble%get_carg('vol'//int2str(state))
                if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                    if( .not. file_exists(volname) ) call cline_volassemble%delete('vol'//int2str(state))
                endif
            endif
        end do
        call xvolassemble%execute(cline_volassemble)
        do state = 1, params%nstates
            volname = refine3D_state_vol_fname(state)
            params%vols(state) = volname
            call cline%set('vol'//int2str(state), volname)
        end do
        call cline%delete('force_volassemble')
        call cline_volassemble%kill
        if( allocated(pinds) ) deallocate(pinds)
    end subroutine inmem_execute

    subroutine pcg_inmem_execute(self, params, build, cline)
        class(rec3D_pcg_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call execute_rec3D_pcg_shared(params, build, cline)
    end subroutine pcg_inmem_execute

    subroutine inmem_finalize_run(self, params, build, cline)
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        call maybe_postprocess_reconstruct3D(params, cline)
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, build, cline)
        use simple_qsys_funs, only: qsys_job_finished
        class(rec3D_inmem_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_rec3D_strategy :: exec_rec3D'))
    end subroutine inmem_cleanup

    ! =====================================================================
    ! DISTRIBUTED-MEMORY IMPLEMENTATION
    ! =====================================================================

    subroutine distr_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        logical :: fall_over
        ! master thread count
        call set_master_num_threads(self%nthr_master, string('rec3D'))
        ! parse parameters and project
        call build%init_params_and_build_spproj(cline, params)
        call sync_resolved_rec_params(params, cline)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
            case DEFAULT
                THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        ! avoid nested directory structure for jobs
        call cline%set('mkdir', 'no')
        ! Even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
        endif
        ! Update eo flags in project
        call build%spproj%write_segment_inside(params%oritype)
        ! setup distributed execution
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, build, cline)
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        type(commander_volassemble) :: xvolassemble
        type(cmdline)               :: cline_volassemble
        type(string)                :: volname, vol_in, raw_fname
        logical, allocatable        :: l_trail_bootstrap(:)
        real,    allocatable        :: nu_replay_lps(:)
        integer                     :: state, part, eo
        if( trim(params%rec_backend) == 'pcg' )then
            ! A stale final filename must never masquerade as a completed worker
            ! from this launch. Workers publish through .tmp + atomic rename.
            do state = 1, params%nstates
                do eo = 0, 1
                    do part = 1, params%nparts
                        raw_fname = refine3D_pcg_raw_accum_fname(state, part, params%numlen, &
                            &merge('odd ', 'even', eo == 1))
                        call del_file(raw_fname)
                        call del_file(raw_fname//'.tmp')
                    enddo
                enddo
            enddo
            call raw_fname%kill
        endif
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        if( trim(params%rec_backend) == 'pcg' )then
            if( params%l_nonuniform )then
                ! reconstruct3D must leave behind the same Q_NU-regularized
                ! primary references and evidence-derived matching low-pass
                ! handoff as a refinement iteration.
                allocate(l_trail_bootstrap(params%nstates), source=.false.)
                allocate(nu_replay_lps(params%nstates),     source=0.0)
                call execute_rec3D_pcg_distributed_master(params, build, cline, &
                    &trail_bootstrap_states=l_trail_bootstrap, nu_replay_finest_lps=nu_replay_lps)
                call filter_pcg_nonuniform_maps(params, build, l_trail_bootstrap, nu_replay_lps)
                deallocate(l_trail_bootstrap, nu_replay_lps)
            else
                call execute_rec3D_pcg_distributed_master(params, build, cline)
            endif
            return
        endif
        ! Assemble volumes on master
        cline_volassemble = cline
        call cline_volassemble%set('prg',  'volassemble')
        call cline_volassemble%set('nthr', self%nthr_master)
        call cline_volassemble%delete('objfun_den')
        call cline_volassemble%delete('objfun_den_w')
        do state = 1, params%nstates
            volname = refine3D_state_vol_fname(state)
            if( cline_volassemble%defined('vol'//int2str(state)) )then
                vol_in = cline_volassemble%get_carg('vol'//int2str(state))
                if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                    if( .not. file_exists(volname) ) call cline_volassemble%delete('vol'//int2str(state))
                endif
            endif
        end do
        call xvolassemble%execute(cline_volassemble)
        call cline_volassemble%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, build, cline)
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        type(string) :: fsc_file
        integer      :: state
        ! updates project file only if mkdir is set to yes
        if( params%mkdir.eq.'yes' )then
            do state = 1, params%nstates
                fsc_file = refine3D_fsc_fname(state)
                call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%add_vol2os_out(refine3D_state_vol_fname(state), &
                        &params%smpd_crop, state, 'vol_cavg')
                else
                    call build%spproj%add_vol2os_out(refine3D_state_vol_fname(state), &
                        &params%smpd_crop, state, 'vol')
                endif
                call fsc_file%kill
            enddo
            call build%spproj%write_segment_inside('out', params%projfile)
        endif
        call maybe_postprocess_reconstruct3D(params, cline)
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, build, cline)
        use simple_qsys_funs, only: qsys_cleanup
        class(rec3D_distr_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        type(builder),               intent(inout) :: build
        class(cmdline),              intent(inout) :: cline
        call qsys_cleanup(params)
        call build%spproj_field%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    subroutine sync_resolved_rec_params(params, cline)
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        if( params%box > 0 ) call cline%set('box', params%box)
        if( params%smpd > TINY ) call cline%set('smpd', params%smpd)
        if( params%box_crop > 0 ) call cline%set('box_crop', params%box_crop)
        if( params%smpd_crop > TINY ) call cline%set('smpd_crop', params%smpd_crop)
        if( params%mskdiam > 0. )then
            call cline%set('mskdiam', params%mskdiam)
        else
            call cline%delete('mskdiam')
        endif
    end subroutine sync_resolved_rec_params

    subroutine maybe_postprocess_reconstruct3D(params, cline)
        use simple_commanders_volops, only: postprocess_volume_from_files
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        type(parameters) :: params_pp
        type(string)     :: fname_vol, fname_fsc
        real             :: smpd
        integer          :: state, nptcls, ldim(3)
        if( trim(params%postprocess) /= 'yes' )return
        if( cline%defined('part') )return
        if( params%l_nonuniform )then
            write(logfhandle,'(A)') &
                &'>>> reconstruct3D postprocess: using classical postprocessing'
        endif
        do state = 1, params%nstates
            params_pp = params
            fname_vol = refine3D_state_vol_fname(state)
            fname_fsc = refine3D_fsc_fname(state)
            if( .not. file_exists(fname_vol) )then
                call fname_vol%kill
                call fname_fsc%kill
                cycle
            endif
            call find_ldim_nptcls(fname_vol, ldim, nptcls)
            smpd = params%smpd_crop
            call postprocess_volume_from_files(fname_vol, fname_fsc, ldim(1), smpd, params_pp, cline, state, &
                &density_window_bfac=trim(params%rec_backend) == 'pcg')
            call fname_vol%kill
            call fname_fsc%kill
        enddo
    end subroutine maybe_postprocess_reconstruct3D

end module simple_rec3D_strategy
