module simple_refine3D_strategy
use simple_core_module_api
use simple_refine3D_fnames
use simple_matcher_refvol_utils
use simple_builder,                 only: builder
use simple_parameters,              only: parameters
use simple_cmdline,                 only: cmdline
use simple_qsys_env,                only: qsys_env
use simple_convergence,             only: convergence
use simple_decay_funs,              only: inv_cos_decay, cos_decay
use simple_cluster_seed,            only: gen_labelling
use simple_euclid_sigma2,           only: sigma2_star_from_iter
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
implicit none

public :: refine3D_strategy, refine3D_inmem_strategy, refine3D_distr_strategy
public :: create_refine3D_strategy

private
#include "simple_local_flags.inc"

!> Minimal strategy interface - only divergent operations
type, abstract :: refine3D_strategy
contains
    procedure(init_interface),          deferred :: initialize
    procedure(exec_iter_interface),     deferred :: execute_iteration
    procedure(finalize_iter_interface), deferred :: finalize_iteration
    procedure(finalize_run_interface),  deferred :: finalize_run
    procedure(cleanup_interface),       deferred :: cleanup
end type refine3D_strategy

type :: refine3D_bench_state
    integer(timer_int_kind) :: t_init     = 0
    integer(timer_int_kind) :: t_prob     = 0
    integer(timer_int_kind) :: t_sched    = 0
    integer(timer_int_kind) :: t_assemble = 0
    integer(timer_int_kind) :: t_tot      = 0
    real(timer_int_kind)    :: rt_init     = 0.
    real(timer_int_kind)    :: rt_prob     = 0.
    real(timer_int_kind)    :: rt_sched    = 0.
    real(timer_int_kind)    :: rt_assemble = 0.
    real(timer_int_kind)    :: rt_tot      = 0.
end type refine3D_bench_state

! ======================================================================
! SHARED-MEMORY IMPLEMENTATION
! ======================================================================

type, extends(refine3D_strategy) :: refine3D_inmem_strategy
    logical :: l_sigma
    type(cmdline) :: cline_calc_group_sigmas
    type(convergence) :: conv
contains
    procedure :: initialize         => inmem_initialize
    procedure :: execute_iteration  => inmem_execute_iteration
    procedure :: finalize_iteration => inmem_finalize_iteration
    procedure :: finalize_run       => inmem_finalize_run
    procedure :: cleanup            => inmem_cleanup
end type refine3D_inmem_strategy

! ======================================================================
! DISTRIBUTED MASTER IMPLEMENTATION
! ======================================================================

type, extends(refine3D_strategy) :: refine3D_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    type(refine3D_bench_state), private :: bench
    integer        :: nthr_master
    logical        :: have_oris
    logical        :: l_multistates
    logical        :: l_combine_eo
    ! Prototypes / persistent command lines
    type(cmdline) :: cline_rec3D
    type(cmdline) :: cline_calc_pspec_distr
    type(cmdline) :: cline_prob_align_distr
    type(cmdline) :: cline_calc_group_sigmas
    type(cmdline) :: cline_volassemble
    type(cmdline) :: cline_postprocess
    type(convergence) :: conv
contains
    procedure :: initialize         => distr_initialize
    procedure :: execute_iteration  => distr_execute_iteration
    procedure :: finalize_iteration => distr_finalize_iteration
    procedure :: finalize_run       => distr_finalize_run
    procedure :: cleanup            => distr_cleanup
end type refine3D_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: refine3D_strategy, parameters, builder, cmdline
        class(refine3D_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_iter_interface(self, params, build, cline, converged)
        import :: refine3D_strategy, parameters, builder, cmdline
        class(refine3D_strategy), intent(inout) :: self
        type(parameters),         intent(inout) :: params
        type(builder),            intent(inout) :: build
        type(cmdline),            intent(inout) :: cline
        logical,                  intent(out)   :: converged
    end subroutine exec_iter_interface

    subroutine finalize_iter_interface(self, params, build)
        import :: refine3D_strategy, parameters, builder
        class(refine3D_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        type(builder),            intent(inout) :: build
    end subroutine finalize_iter_interface

    subroutine finalize_run_interface(self, params, build, cline)
        import :: refine3D_strategy, parameters, builder, cmdline
        class(refine3D_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
        type(builder),            intent(inout) :: build
        type(cmdline),            intent(inout) :: cline
    end subroutine finalize_run_interface

    subroutine cleanup_interface(self, params)
        import :: refine3D_strategy, parameters
        class(refine3D_strategy), intent(inout) :: self
        type(parameters),         intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    !> Strategy selection based on command-line shape.
    function create_refine3D_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(refine3D_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            allocate(refine3D_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(refine3D_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_refine3D_strategy

    subroutine prepare_assembly_cline( cline, params, nthr, cline_assembly )
        type(cmdline),    intent(in)    :: cline
        type(parameters), intent(in)    :: params
        integer,          intent(in)    :: nthr
        type(cmdline),    intent(inout) :: cline_assembly
        type(string) :: volname, vol_in
        integer      :: state
        cline_assembly = cline
        call cline_assembly%set('which_iter', params%which_iter)
        call cline_assembly%set('nthr',       nthr)
        call cline_assembly%set('combine_eo', params%combine_eo)
        if( params%l_update_frac ) call cline_assembly%set('update_frac', params%update_frac)
        do state = 1, params%nstates
            volname = refine3D_state_vol_fname(state)
            if( cline_assembly%defined('vol'//int2str(state)) )then
                vol_in = cline_assembly%get_carg('vol'//int2str(state))
                if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                    if( .not. file_exists(volname) ) call cline_assembly%delete('vol'//int2str(state))
                endif
            endif
        end do
        call volname%kill
        call vol_in%kill
    end subroutine prepare_assembly_cline

    subroutine promote_assembly_nspace_if_needed( params, cline_assembly, force )
        type(parameters), intent(in)    :: params
        type(cmdline),    intent(inout) :: cline_assembly
        logical, optional, intent(in)   :: force
        logical :: l_force
        l_force = .false.
        if( present(force) ) l_force = force
        ! Cartesian and polar=obsfield can emit the next iteration's reference grid.
        if( l_force .and. params%can_promote_assembly_ref_nspace() .and. params%nspace_next > params%nspace )then
            call cline_assembly%set('nspace', params%nspace_next)
            call cline_assembly%delete('nspace_next')
            if( params%pftsz_next > 0 ) call cline_assembly%set('pftsz', params%pftsz_next)
            call cline_assembly%delete('pftsz_next')
        else if( params%uses_next_assembly_ref_nspace() )then
            call cline_assembly%set('nspace', params%assembly_ref_nspace())
            call cline_assembly%delete('nspace_next')
            if( params%pftsz_next > 0 ) call cline_assembly%set('pftsz', params%pftsz_next)
            call cline_assembly%delete('pftsz_next')
        endif
    end subroutine promote_assembly_nspace_if_needed

    subroutine refresh_resolution_fields_from_fsc( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(string) :: fsc_file
        real,    allocatable :: fsc(:), res(:), res0143s(:)
        logical, allocatable :: has_fsc(:)
        real    :: fsc05, fsc0143
        integer :: state, iptcl, istate
        allocate(res0143s(params%nstates), source=0.)
        allocate(has_fsc(params%nstates), source=.false.)
        res = get_resarr(params%box_crop, params%smpd_crop)
        do state = 1, params%nstates
            fsc_file = refine3D_fsc_fname(state)
            if( .not. file_exists(fsc_file) ) cycle
            fsc = file2rarr(fsc_file)
            call get_resolution(fsc, res, fsc05, fsc0143)
            res0143s(state) = fsc0143
            has_fsc(state)  = .true.
            deallocate(fsc)
        end do
        if( any(has_fsc) )then
            if( params%nstates == 1 )then
                if( has_fsc(1) ) call build%spproj_field%set_all2single('res', res0143s(1))
            else
                do iptcl = 1, build%spproj_field%get_noris()
                    istate = build%spproj_field%get_state(iptcl)
                    if( istate > 0 .and. istate <= params%nstates )then
                        if( has_fsc(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                    endif
                end do
            endif
        endif
        call fsc_file%kill
        if( allocated(fsc)      ) deallocate(fsc)
        if( allocated(res)      ) deallocate(res)
        if( allocated(res0143s) ) deallocate(res0143s)
        if( allocated(has_fsc)  ) deallocate(has_fsc)
    end subroutine refresh_resolution_fields_from_fsc

    subroutine remove_polar_partial_sum_files( params )
        type(parameters), intent(in) :: params
        type(string) :: fname
        integer :: ipart, numlen_part, state
        numlen_part = max(1, params%numlen)
        do ipart = 1,max(1, params%nparts)
            fname = refine3D_polar_cavgs_part_fname('even', ipart, numlen_part)
            if( file_exists(fname) ) call del_file(fname)
            fname = refine3D_polar_cavgs_part_fname('odd', ipart, numlen_part)
            if( file_exists(fname) ) call del_file(fname)
            fname = refine3D_polar_ctfsqsums_part_fname('even', ipart, numlen_part)
            if( file_exists(fname) ) call del_file(fname)
            fname = refine3D_polar_ctfsqsums_part_fname('odd', ipart, numlen_part)
            if( file_exists(fname) ) call del_file(fname)
            do state = 1,max(1, params%nstates)
                fname = refine3D_obsfield_part_fname(state, ipart, numlen_part)
                if( file_exists(fname) ) call del_file(fname)
            enddo
        end do
        call fname%kill
    end subroutine remove_polar_partial_sum_files

    subroutine remove_cartesian_partial_rec_files( params )
        type(parameters), intent(in) :: params
        type(string) :: fname
        integer :: state, ipart, numlen_part
        numlen_part = max(1, params%numlen)
        do state = 1,max(1, params%nstates)
            do ipart = 1,max(1, params%nparts)
                fname = refine3D_partial_rec_fname(state, ipart, numlen_part, 'even')
                if( file_exists(fname) ) call del_file(fname)
                fname = refine3D_partial_rec_fname(state, ipart, numlen_part, 'odd')
                if( file_exists(fname) ) call del_file(fname)
                fname = refine3D_partial_rho_fname(state, ipart, numlen_part, 'even')
                if( file_exists(fname) ) call del_file(fname)
                fname = refine3D_partial_rho_fname(state, ipart, numlen_part, 'odd')
                if( file_exists(fname) ) call del_file(fname)
            enddo
        enddo
        call fname%kill
    end subroutine remove_cartesian_partial_rec_files

    subroutine remove_partial_assembly_input_files( params )
        type(parameters), intent(in) :: params
        if( params%l_polar )then
            call remove_polar_partial_sum_files(params)
        else
            call remove_cartesian_partial_rec_files(params)
        endif
    end subroutine remove_partial_assembly_input_files

    subroutine delete_volume_source_keys( cline, nstates )
        type(cmdline), intent(inout) :: cline
        integer,       intent(in)    :: nstates
        integer :: state
        do state = 1,nstates
            call cline%delete('vol'//int2str(state))
        enddo
    end subroutine delete_volume_source_keys

    subroutine delete_volume_source_job_keys( job_descr, nstates )
        type(chash), intent(inout) :: job_descr
        integer,     intent(in)    :: nstates
        integer :: state
        do state = 1,nstates
            call job_descr%delete('vol'//int2str(state))
        enddo
    end subroutine delete_volume_source_job_keys

    ! ======================================================================
    ! SHARED-MEMORY STRATEGY METHODS
    ! ======================================================================

    subroutine inmem_initialize(self, params, build, cline)
        use simple_commanders_euclid,       only: commander_calc_group_sigmas, commander_calc_pspec
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(inout) :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        type(commander_calc_pspec)            :: xcalc_pspec
        type(cmdline)                         :: cline_calc_pspec
        integer                               :: startit, state
        logical                               :: l_proj_dirty
        ! Full in-memory toolbox build (required for refine3D_exec)
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        ! startit
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        l_proj_dirty = .false.
        ! Some refine modes manage sampling/updatecnt internally
        if( params%l_prob_align_mode )then
            ! random sampling and updatecnt dealt with in prob_align
        else
            if( startit == 1 )then
                call build%spproj_field%clean_entry('updatecnt', 'sampled')
                l_proj_dirty = .true.
            endif
        endif
        if( trim(params%continue) == 'yes' )then
            THROW_HARD('shared-memory implementation of refine3D does not support continue=yes')
        endif
        ! Input reference validation
        if( params%l_polar )then
            if( any_volume_source_defined(cline, params%nstates) &
                &.and. (.not. complete_volume_source_defined(cline, params%nstates)) )then
                THROW_HARD('incomplete multi-state volume source; provide vol1..volN or use POLAR_REFS')
            endif
            if( complete_volume_source_defined(cline, params%nstates) )then
                do state = 1, params%nstates
                    if( .not. file_exists(params%vols(state)) ) then
                        THROW_HARD('shared-memory implementation of refine3D needs starting volume input')
                    endif
                end do
            else
                if( .not. polar_ref_sections_available(params) ) then
                    call set_bp_range3D(params, build, cline)
                    call ensure_polar_refs_on_disk(params, build, cline, 1, 'refine3D shared-memory initialization')
                endif
            endif
        else
            if( .not. file_exists(params%vols(1)) ) then
                THROW_HARD('shared-memory implementation of refine3D needs starting volume input')
            endif
        endif
        ! Initial orientation parameters
        if( build%spproj%is_virgin_field(params%oritype) )then
            call build%spproj_field%rnd_oris
            l_proj_dirty = .true.
        endif
        if( l_proj_dirty ) call build%spproj%write_segment_inside(params%oritype)
        ! objfun=euclid initialisation
        self%l_sigma = (params%cc_objfun == OBJFUN_EUCLID)
        self%cline_calc_group_sigmas = cline
        call self%cline_calc_group_sigmas%set('prg', 'calc_group_sigmas')
        if( self%l_sigma )then
            if( file_exists(SIGMA2_GROUP_FBODY//int2str(params%which_iter)//STAR_EXT) )then
                ! assume sigmas2 already available and all corresponding flags have been set
            else
                ! Ensure e/o partitioning prior to calc_pspec
                if( build%spproj_field%get_nevenodd() == 0 )then
                    call build%spproj_field%partition_eo
                    call build%spproj%write_segment_inside(params%oritype)
                endif
                if( startit == 1 )then
                    ! make sure we have weights for the initial noise power estimation
                    call build%spproj_field%set_all2single('w', 1.0)
                    call build%spproj%write_segment_inside(params%oritype)
                endif
                cline_calc_pspec   = cline
                call cline_calc_pspec%set('prg', 'calc_pspec')
                call xcalc_pspec%execute( cline_calc_pspec )
            endif
        endif
        ! Keep run-time counters consistent with the new high-level loop
        params%startit    = startit
        params%which_iter = params%startit
        if( .not.cline%defined('extr_iter') ) params%extr_iter = params%startit
        params%outfile    = 'algndoc'//METADATA_EXT
    end subroutine inmem_initialize

    subroutine inmem_execute_iteration(self, params, build, cline, converged)
        use simple_strategy3D_matcher, only: refine3D_exec
        use simple_commanders_euclid,  only: commander_calc_group_sigmas
        use simple_commanders_prob,    only: commander_prob_align, commander_prob_align_neigh
        use simple_commanders_rec_distr, only: commander_cartesian_volassemble, commander_polar_volassemble
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(inout) :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        logical,                        intent(out)   :: converged
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_prob_align)        :: xprob_align
        type(commander_prob_align_neigh)  :: xprob_align_neigh
        type(commander_cartesian_volassemble) :: xvolassemble
        type(commander_polar_volassemble) :: xpolar_volassemble
        type(cmdline)                     :: cline_prob_align
        type(cmdline)                     :: cline_volassemble
        type(cmdline)                     :: cline_build
        integer                           :: state, iter, extr_iter
        logical                           :: l_prob_state_mode, l_prob_neigh_mode
        logical                           :: l_write_partial_recs
        type(string)                      :: volname
        601 format(A,1X,F12.3)
        iter      = params%which_iter
        extr_iter = params%extr_iter
        call self%conv%print_iteration(params%which_iter)
        ! In shared-memory mode, rebuild the toolbox from the settled
        ! per-iteration settings. Keep the stage startit in params so
        ! final-iteration planning still sees the whole refine3D stage.
        cline_build = cline
        call cline_build%set('startit',    params%startit)
        call cline_build%set('which_iter', iter)
        call cline_build%set('extr_iter',  extr_iter)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call build%init_params_and_build_strategy3D_tbox(cline_build, params)
        params%which_iter = iter
        params%extr_iter  = extr_iter
        params%outfile    = 'algndoc'//METADATA_EXT
        ! communicate iteration counters
        call cline%set('startit',    params%startit)
        call cline%set('which_iter', params%which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        ! noise regularization / annealing
        if( params%l_noise_reg )then
            params%eps = inv_cos_decay(params%which_iter, params%maxits_glob, params%eps_bounds)
            write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
        endif
        if( params%l_lam_anneal )then
            params%lambda = cos_decay(params%which_iter, params%maxits_glob, params%lam_bounds)
            write(logfhandle,601) '>>> LAMBDA, MAP CONNECTIVITY ANNEALING        ', params%lambda
        endif
        ! Per-iteration sigma update (euclid)
        if( self%l_sigma )then
            call self%cline_calc_group_sigmas%set('which_iter', params%which_iter)
            call xcalc_group_sigmas%execute(self%cline_calc_group_sigmas)
        endif
        l_prob_state_mode = trim(params%refine) == 'prob_state'
        l_prob_neigh_mode = trim(params%refine) == 'prob_neigh'
        if( params%l_polar .and. (.not. params%l_prob_align_mode) )then
            if( .not. complete_volume_source_defined(cline, params%nstates) )then
                call set_bp_range3D(params, build, cline)
                call ensure_polar_refs_on_disk(params, build, cline, 1, 'refine3D shared-memory iteration')
            endif
        endif
        ! refine=prob* pre-step
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call cline_prob_align%set('prg', 'prob_align_neigh')
            else
                call cline_prob_align%set('prg', 'prob_align')
            endif
            call cline_prob_align%set('which_iter', params%which_iter)
            if( .not.params%l_polar )then
                do state = 1, params%nstates
                    call cline_prob_align%set('vol'//int2str(state), params%vols(state))
                enddo
            endif
            ! communicate changes to probabilistic alignment
            call build%spproj%write_segment_inside(params%oritype)
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call xprob_align_neigh%execute( cline_prob_align )
            else
                call xprob_align%execute( cline_prob_align )
            endif
            ! communicate back changes made by probabilistic alignment including sampling
            call build%spproj%read_segment(params%oritype, params%projfile)
            ! make sure prob_align and refine see the same information
            if( cline%defined('lp') ) params%lp = cline%get_rarg('lp')
        endif
        ! main refinement step
        l_write_partial_recs = trim(params%volrec) .eq. 'yes' .or. params%l_polar
        if( l_write_partial_recs )then
            ! Legacy handshake for rec-writing helpers that still inspect this key.
            ! The strategy owns the actual assembly dispatch decision.
            call cline%set('force_volassemble', 'yes')
            call remove_partial_assembly_input_files(params)
        endif
        call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs)
        if( l_write_partial_recs )then
            call prepare_assembly_cline(cline, params, params%nthr, cline_volassemble)
            call promote_assembly_nspace_if_needed(params, cline_volassemble)
            if( params%l_polar )then
                call xpolar_volassemble%execute(cline_volassemble)
            else
                call xvolassemble%execute(cline_volassemble)
            endif
            if( trim(params%volrec) .eq. 'yes' )then
                do state = 1, params%nstates
                    volname = refine3D_state_vol_fname(state)
                    params%vols(state) = volname
                    call cline%set('vol'//int2str(state), volname)
                end do
            endif
            if( params%l_polar ) params%refs = refine3D_iter_refs_fname(params%which_iter)
            call cline%delete('force_volassemble')
        endif
        if( l_write_partial_recs ) call refresh_resolution_fields_from_fsc(params, build)
        converged = .false.
        select case(trim(params%refine))
            case('sigma')
                ! no convergence report for sigma-only updates
            case default
                converged = self%conv%check_conv3D(params, cline, build%spproj_field, params%msk)
        end select
        if( converged .and. params%l_polar .and. trim(params%polar) == 'obsfield' .and. &
            &(.not. params%uses_next_assembly_ref_nspace()) .and. &
            &params%can_promote_assembly_ref_nspace() .and. params%nspace_next > params%nspace )then
            call prepare_assembly_cline(cline, params, params%nthr, cline_volassemble)
            call promote_assembly_nspace_if_needed(params, cline_volassemble, force=.true.)
            call xpolar_volassemble%execute(cline_volassemble)
        endif
        ! input volume should only be used once in polar mode
        if( params%l_polar ) call delete_volume_source_keys(cline, params%nstates)
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration(self, params, build)
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        ! Shared-memory rebuilds the builder every iteration, so persist the
        ! updated orientations before the next rebuild reads the project file.
        if( trim(params%refine) /= 'sigma' ) call build%spproj%write_segment_inside(params%oritype)
    end subroutine inmem_finalize_iteration

    subroutine inmem_finalize_run(self, params, build, cline)
        use simple_commanders_euclid, only: commander_calc_group_sigmas
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(string) :: fsc_file, vol, vol_iter
        integer      :: state
        ! report last iteration
        call cline%delete( 'startit' )
        call cline%set('endit', real(params%which_iter))
        ! update project with new orientations
        call build%spproj%write_segment_inside(params%oritype)
        call del_file(params%outfile)
        if( trim(params%volrec) .eq. 'yes' )then
            do state = 1, params%nstates
                if( build%spproj_field%get_pop(state, 'state') == 0 )then
                    ! cleanup empty state
                    if( trim(params%oritype).eq.'cls3D' )then
                        call build%spproj%remove_entry_from_osout('vol_cavg', state)
                    else
                        call build%spproj%remove_entry_from_osout('vol', state)
                    endif
                    call build%spproj%remove_entry_from_osout('fsc', state)
                else
                    ! add state volume, fsc to os_out
                    fsc_file  = refine3D_fsc_fname(state)
                    call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                    vol       = refine3D_state_vol_fname(state)
                    vol_iter  = vol
                    if( trim(params%oritype).eq.'cls3D' )then
                        call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                    else
                        call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                    endif
                endif
            end do
            call build%spproj%write_segment_inside('out')
        endif
        if( self%l_sigma )then
            ! so final sigma2 can be used for a subsequent refine3D run
            call self%cline_calc_group_sigmas%set('which_iter',params%which_iter+1)
            call xcalc_group_sigmas%execute(self%cline_calc_group_sigmas)
        endif
        call simple_touch(JOB_FINISHED_FBODY)
        call fsc_file%kill
        call vol%kill
        call vol_iter%kill
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params)
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        ! no-op
    end subroutine inmem_cleanup

    ! ======================================================================
    ! DISTRIBUTED STRATEGY METHODS
    ! ======================================================================

    subroutine distr_initialize(self, params, build, cline)
        use simple_exec_helpers,      only: set_master_num_threads
        use simple_commanders_rec,    only: commander_rec3D
        use simple_commanders_euclid, only: commander_calc_pspec
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(inout) :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        type(commander_rec3D)      :: xrec3D
        type(commander_calc_pspec) :: xcalc_pspec_distr
        type(cmdline) :: cline_tmp
        type(string)  :: prev_refine_path, target_name, fname_vol, vol, fsc_file
        type(string), allocatable :: list(:)
        real    :: smpd
        integer :: state, box, nfiles, i
        logical :: err, fall_over, vol_defined, l_prob_state_mode, l_prob_neigh_mode
        ! deal with #threads for the master process
        call set_master_num_threads(self%nthr_master, string('REFINE3D'))
        ! Local options / flags
        self%l_multistates = cline%defined('nstates')
        ! init project
        call build%init_params_and_build_spproj(cline, params)
        ! Keep worker command lines explicit; assembly commanders own polar reference assembly.
        if( params%l_polar )then
            call cline%set('nspace', params%nspace)
            ! build%fsc is normally allocated in build_general_tbox (not called for the distr master).
            ! Allocate it here so that set_bp_range3D can read on-disk FSCs in distr_execute_iteration.
            if( .not. allocated(build%fsc) .and. params%box_crop > 0 )then
                allocate( build%fsc(params%nstates, fdim(params%box_crop)-1), source=0.0 )
            endif
        endif
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
            case DEFAULT
                write(logfhandle,*)'Unsupported ORITYPE; simple_refine3D_strategy :: distr_initialize'
        end select
        if( fall_over ) THROW_HARD('no particles found! refine3D distributed master')

        if( .not. self%l_multistates .and. params%nstates > 1 )then
            THROW_HARD('nstates > 1 but refine mode is single')
        endif
        ! final iteration with combined e/o
        self%l_combine_eo = .false.
        if( trim(params%combine_eo).eq.'yes' )then
            self%l_combine_eo = .true.
            call cline%set('combine_eo','no')
            params%combine_eo = 'no'
        endif
        ! set mkdir to no (to avoid nested directory structure in scheduled parts)
        call cline%set('mkdir', 'no')
        ! distributed environment
        call self%qenv%new(params, params%nparts)
        ! prepare prototype command lines
        self%cline_rec3D = cline
        self%cline_calc_pspec_distr    = cline
        self%cline_prob_align_distr    = cline
        self%cline_postprocess         = cline
        self%cline_calc_group_sigmas   = cline
        call self%cline_rec3D%set( 'prg', 'reconstruct3D' )
        call self%cline_rec3D%delete( 'nspace_next' )
        call self%cline_calc_pspec_distr%set(    'prg', 'calc_pspec' )
        l_prob_state_mode = trim(params%refine) == 'prob_state'
        l_prob_neigh_mode = trim(params%refine) == 'prob_neigh'
        if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
            call self%cline_prob_align_distr%set( 'prg', 'prob_align_neigh' )
        else
            call self%cline_prob_align_distr%set( 'prg', 'prob_align' )
        endif
        call self%cline_postprocess%set(         'prg', 'postprocess' )
        call self%cline_calc_group_sigmas%set(   'prg', 'calc_group_sigmas' )
        call self%cline_postprocess%set('mirr',    'no')
        call self%cline_postprocess%set('mkdir',   'no')
        call self%cline_postprocess%set('imgkind', 'vol')
        if( trim(params%oritype).eq.'cls3D' ) call self%cline_postprocess%set('imgkind', 'vol_cavg')
        ! remove unnecessary volume keys from postprocess
        do state = 1, params%nstates
            vol = 'vol'//int2str(state)
            call self%cline_postprocess%delete(vol%to_char())
        enddo
        ! E/O PARTITIONING (required for consistent e/o reconstructions)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! STATE LABEL INIT
        if( self%l_multistates )then
            if( build%spproj_field%get_n('state') /= params%nstates )then
                call gen_labelling(build%spproj_field, params%nstates, 'squared_uniform')
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! CONTINUE / STARTUP
        self%have_oris = .true.
        if( params%continue .eq. 'yes' )then
            ! Continuing from previous refinement round
            do state=1,params%nstates
                vol = 'vol' // int2str(state)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%get_vol('vol_cavg', state, fname_vol, smpd, box)
                else
                    call build%spproj%get_vol('vol', state, fname_vol, smpd, box)
                endif
                call cline%set(vol%to_char(), fname_vol)
                params%vols(state) = fname_vol
            end do
            prev_refine_path = get_fpath(fname_vol)
            if( simple_abspath(prev_refine_path,check_exists=.false.) .eq. CWD_GLOB )then
                do state=1,params%nstates
                    fsc_file  = refine3D_fsc_fname(state)
                    if( .not.file_exists(fsc_file)) THROW_HARD('Missing file: '//fsc_file%to_char())
                end do
                if( params%l_update_frac )then
                    call simple_list_files(refine3D_partial_rec_glob(prev_refine_path%to_char()), list)
                    nfiles = size(list)
                    err = params%nparts * 4 /= nfiles
                    if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
                    deallocate(list)
                endif
                if( params%cc_objfun==OBJFUN_EUCLID )then
                    call simple_list_files(prev_refine_path%to_char()//SIGMA2_FBODY//'*', list)
                    nfiles = size(list)
                    if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                    deallocate(list)
                endif
            else
                ! carry over FSCs
                do state=1,params%nstates
                    fsc_file  = refine3D_fsc_fname(state)
                    call simple_copy_file(prev_refine_path//fsc_file, fsc_file)
                end do
                if( params%cc_objfun==OBJFUN_EUCLID )then
                    call simple_list_files(prev_refine_path%to_char()//SIGMA2_FBODY//'*', list)
                    nfiles = size(list)
                    if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                    do i=1,nfiles
                        target_name = string(PATH_HERE)//basename(list(i))
                        call simple_copy_file(list(i), target_name)
                    end do
                    deallocate(list)
                endif
            endif
        else
            ! generate initial noise power estimates
            if( .not.file_exists(sigma2_star_from_iter(params%startit)) )then
                call build%spproj_field%set_all2single('w', 1.0)
                call build%spproj%write_segment_inside(params%oritype)
                call xcalc_pspec_distr%execute(self%cline_calc_pspec_distr)
            endif
            ! check if we have input volume(s) and/or 3D orientations
            vol_defined = .true.
            do state = 1,params%nstates
                if( .not. cline%defined('vol'//int2str(state)) ) vol_defined = .false.
            enddo
            self%have_oris = .not. build%spproj%is_virgin_field(params%oritype)
            if( .not. self%have_oris )then
                call build%spproj_field%rnd_oris
                self%have_oris = .true.
                call build%spproj%write_segment_inside(params%oritype)
            endif
            if( params%l_polar )then
                if( .not.vol_defined )then
                    if( polar_ref_sections_available(params) ) vol_defined = .true.
                endif
            endif
            if( .not. vol_defined )then
                if( params%l_polar .and. trim(params%polar) == 'obsfield' )then
                    THROW_HARD('obsfield refine3D initialization requires obsfield-derived POLAR_REFS; Cartesian bootstrap is not allowed')
                endif
                ! reconstructions needed
                cline_tmp = self%cline_rec3D
                call cline_tmp%delete('trail_rec')
                call cline_tmp%delete('objfun')
                call cline_tmp%delete('sigma_est')
                call cline_tmp%set('objfun', 'cc')
                call xrec3D%execute( cline_tmp )
                do state = 1,params%nstates
                    ! rename volumes and update cline/params
                    call simple_rename(refine3D_state_vol_fname(state), refine3D_startvol_fname(state))
                    params%vols(state) = refine3D_startvol_fname(state)
                    vol = 'vol'//int2str(state)
                    call cline%set(vol%to_char(), params%vols(state))
                    ! keep unfiltered copies
                    call simple_copy_file(refine3D_state_halfvol_fname(state, 'even'), &
                        &refine3D_startvol_half_fname(state, 'even', unfil=.true.))
                    call simple_rename(refine3D_state_halfvol_fname(state, 'even'), &
                        &refine3D_startvol_half_fname(state, 'even'))
                    call simple_copy_file(refine3D_state_halfvol_fname(state, 'odd'), &
                        &refine3D_startvol_half_fname(state, 'odd', unfil=.true.))
                    call simple_rename(refine3D_state_halfvol_fname(state, 'odd'), &
                        &refine3D_startvol_half_fname(state, 'odd'))
                enddo
                vol_defined = .true.
            endif
            ! euclid first-sigmas
            if( params%cc_objfun==OBJFUN_EUCLID )then
                call self%cline_calc_group_sigmas%set('nthr', self%nthr_master)
            endif
        endif
        ! prepare job description
        call cline%gen_job_descr(self%job_descr)
        call self%job_descr%set('prg', 'refine3D')
        ! Keep consistent iteration counters
        if( .not.cline%defined('extr_iter') ) params%extr_iter = params%startit
        call prev_refine_path%kill
        call target_name%kill
        call fname_vol%kill
        call vol%kill
        call fsc_file%kill
    end subroutine distr_initialize

    subroutine distr_execute_iteration(self, params, build, cline, converged)
        use simple_commanders_rec_distr, only: commander_cartesian_volassemble, commander_polar_volassemble
        use simple_commanders_volops, only: commander_postprocess
        use simple_commanders_euclid, only: commander_calc_group_sigmas
        use simple_commanders_prob,   only: commander_prob_align, commander_prob_align_neigh
        use simple_fsc,               only: plot_fsc
        use simple_image,             only: image
        use simple_image_msk,         only: image_msk
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(inout) :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        logical,                        intent(out)   :: converged
        type(commander_postprocess)       :: xpostprocess
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_prob_align)        :: xprob_align_distr
        type(commander_prob_align_neigh)  :: xprob_align_neigh_distr
        type(commander_cartesian_volassemble) :: xvolassemble
        type(commander_polar_volassemble) :: xpolar_volassemble
        type(cmdline) :: cline_prob_align, cline_volassemble
        type(string)  :: str
        type(string)  :: vol, vol_iter, fsc_templ, fsc_file
        type(string)  :: fname_vol, volpproc, vollp
        real, allocatable :: res(:), fsc(:)
        integer, allocatable :: state_pops(:)
        integer :: state, iter
        logical :: l_prob_state_mode, l_prob_neigh_mode
        if( L_BENCH_GLOB )then
            self%bench%t_init = tic()
            self%bench%t_tot  = tic()
        endif
        601 format(A,1X,F12.3)
        iter     = params%which_iter
        call self%conv%print_iteration(iter)
        ! annealing
        if( params%l_noise_reg )then
            params%eps = inv_cos_decay(iter, params%maxits_glob, params%eps_bounds)
            write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
        endif
        if( params%l_lam_anneal )then
            params%lambda = cos_decay(iter, params%maxits_glob, params%lam_bounds)
            write(logfhandle,601) '>>> LAMBDA, MAP CONNECTIVITY ANNEALING        ', params%lambda
        endif
        ! per-iteration group sigmas (euclid)
        if( trim(params%objfun).eq.'euclid' )then
            call self%cline_calc_group_sigmas%set('which_iter', iter)
            call xcalc_group_sigmas%execute(self%cline_calc_group_sigmas)
        endif
        ! ensure spproj is current
        if( self%have_oris .or. iter > params%startit )then
            call build%spproj%read(params%projfile)
        endif
        ! prob refinement
        if( L_BENCH_GLOB )then
            self%bench%rt_init = toc(self%bench%t_init)
            self%bench%t_prob  = tic()
        endif
        l_prob_state_mode = trim(params%refine) == 'prob_state'
        l_prob_neigh_mode = trim(params%refine) == 'prob_neigh'
        if( params%l_polar .and. (.not. params%l_prob_align_mode) )then
            if( .not. complete_volume_source_defined(cline, params%nstates) )then
                call set_bp_range3D(params, build, cline)
                call ensure_polar_refs_on_disk(params, build, cline, 1, 'refine3D distributed iteration')
            endif
        endif
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call cline_prob_align%set('prg', 'prob_align_neigh')
            else
                call cline_prob_align%set('prg', 'prob_align')
            endif
            call cline_prob_align%set('which_iter', iter)
            call build%spproj%write_segment_inside(params%oritype)
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call xprob_align_neigh_distr%execute( cline_prob_align )
            else
                call xprob_align_distr%execute( cline_prob_align )
            endif
        endif
        if( L_BENCH_GLOB )then
            self%bench%rt_prob = toc(self%bench%t_prob)
            self%bench%t_sched = tic()
        endif
        ! update job description
        call self%job_descr%set( 'which_iter', int2str(iter))
        call cline%set(          'which_iter', iter)
        call self%job_descr%set( 'extr_iter',  int2str(params%extr_iter))
        call cline%set(          'extr_iter',  params%extr_iter)
        call self%job_descr%set( 'startit',    int2str(params%startit))
        call cline%set(          'startit',    params%startit)
        if( (trim(params%volrec).eq.'yes') .or. params%l_polar )then
            call remove_partial_assembly_input_files(params)
        endif
        ! schedule distributed jobs
        call self%qenv%gen_scripts_and_schedule_jobs( self%job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
        ! merge alignment docs
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY)
        ! assemble volumes, postprocess, automask
        if( L_BENCH_GLOB )then
            self%bench%rt_sched   = toc(self%bench%t_sched)
            self%bench%t_assemble = tic()
        endif
        if( (trim(params%volrec).eq.'yes') .or. params%l_polar )then
            select case(trim(params%refine))
                case('eval')
                    ! nothing
                case DEFAULT
                    call prepare_assembly_cline(cline, params, self%nthr_master, cline_volassemble)
                    call promote_assembly_nspace_if_needed(params, cline_volassemble)
                    if( params%l_polar )then
                        call xpolar_volassemble%execute(cline_volassemble)
                        params%refs = refine3D_iter_refs_fname(iter)
                    else
                        call xvolassemble%execute(cline_volassemble)
                    endif
                    if( trim(params%volrec).eq.'yes' )then
                        ! rename & add volumes to project & update job_descr
                        call build%spproj_field%get_pops(state_pops, 'state')
                        do state = 1,params%nstates
                            if( state_pops(state) == 0 )then
                                vol = 'vol'//int2str(state)
                                call cline%delete(vol%to_char())
                                call self%job_descr%delete(vol%to_char() )
                                if( trim(params%oritype).eq.'cls3D' )then
                                    call build%spproj%remove_entry_from_osout('vol_cavg', state)
                                else
                                    call build%spproj%remove_entry_from_osout('vol', state)
                                endif
                                call build%spproj%remove_entry_from_osout('fsc', state)
                            else
                                vol_iter  = refine3D_state_vol_fname(state)
                                fsc_file  = refine3D_fsc_fname(state)
                                call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                                ! FSC plot
                                res       = get_resarr(params%box_crop, params%smpd_crop)
                                fsc       = file2rarr(fsc_file)
                                fsc_templ = refine3D_fsc_plot_fbody(state, iter)
                                call plot_fsc(size(fsc), fsc, res, params%smpd_crop, fsc_templ%to_char())
                                if( trim(params%oritype).eq.'cls3D' )then
                                    call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                                else
                                    call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                                endif
                                vol = 'vol'//int2str(state)
                                call self%job_descr%set( vol, vol_iter )
                                call cline%set(vol, vol_iter)
                            endif
                        enddo
                        call build%spproj%write_segment_inside('out')
                        ! per-state postprocess (and optional automask)
                        do state = 1,params%nstates
                            if( state_pops(state) == 0 ) cycle
                            call self%cline_postprocess%set('state', state)
                            call self%cline_postprocess%set('nthr',  self%nthr_master)
                            if( cline%defined('lp') ) call self%cline_postprocess%set('lp', params%lp)
                            call xpostprocess%execute(self%cline_postprocess)
                            volpproc = refine3D_state_vol_suffix_fname(state, PPROC_SUFFIX)
                            vollp    = refine3D_state_vol_suffix_fname(state, LP_SUFFIX)
                            ! keep per-iteration postprocessed copies
                            vol_iter = refine3D_iter_vol_fname(state, iter, PPROC_SUFFIX)
                            call simple_copy_file(volpproc, vol_iter)
                            vol_iter = refine3D_iter_vol_fname(state, iter, LP_SUFFIX)
                            call simple_copy_file(vollp, vol_iter)
                            if( iter > 1 .and. params%keepvol.eq.'no' )then
                                call del_file(refine3D_iter_vol_fname(state, iter-1, PPROC_SUFFIX))
                                call del_file(refine3D_iter_vol_fname(state, iter-1, LP_SUFFIX))
                            endif
                        enddo
                    endif
            end select
            if( L_BENCH_GLOB ) self%bench%rt_assemble = toc(self%bench%t_assemble)
        endif
        ! convergence
        ! For strict same-iteration metrics, evaluate convergence after assembly
        ! and refresh the FSC-derived resolution fields in the active table.
        converged = .false.
        select case(trim(params%refine))
            case('eval')
                ! nothing
            case DEFAULT
                if( (trim(params%volrec).eq.'yes') .or. params%l_polar ) call refresh_resolution_fields_from_fsc(params, build)
                converged = self%conv%check_conv3D(params, cline, build%spproj_field, params%msk)
                converged = converged .and. (iter >= params%startit + 2)
        end select
        ! Force termination at requested number of iterations (maxits is run-length)
        if( (iter - params%startit + 1) >= params%maxits ) converged = .true.
        if( converged .and. params%l_polar .and. trim(params%polar) == 'obsfield' .and. &
            &(.not. params%uses_next_assembly_ref_nspace()) .and. &
            &params%can_promote_assembly_ref_nspace() .and. params%nspace_next > params%nspace )then
            call prepare_assembly_cline(cline, params, self%nthr_master, cline_volassemble)
            call promote_assembly_nspace_if_needed(params, cline_volassemble, force=.true.)
            call xpolar_volassemble%execute(cline_volassemble)
            params%refs = refine3D_iter_refs_fname(iter)
        endif
        ! polar reference bookkeeping; the polar assembly commander owns the actual reduction.
        if( params%l_polar )then
            call delete_volume_source_job_keys(self%job_descr, params%nstates)
            call delete_volume_source_keys(cline, params%nstates)
            if( iter > 1 .and. params%keepvol.eq.'no' )then
                call del_file(refine3D_iter_refs_fname(iter-1))
            endif
        endif
        ! combine even/odd final iteration
        if ( self%l_combine_eo .and. converged )then
            converged            = .false.
            self%l_combine_eo    = .false.
            params%maxits = (iter - params%startit + 1) + 1
            call cline%set("maxits", params%maxits)
            call self%job_descr%set("maxits", int2str(params%maxits))
            params%combine_eo    = 'yes'
            params%l_update_frac = .false.
            params%update_frac   = 1.0
            params%lplim_crit    = min(0.143,params%lplim_crit)
            call cline%set('combine_eo',  'yes')
            call cline%set('lplim_crit',  params%lplim_crit)
            call cline%set('update_frac', 1.0)
            call self%job_descr%set('combine_eo',  'yes')
            call self%job_descr%set('lplim_crit',  real2str(params%lplim_crit))
            call self%job_descr%set('update_frac', real2str(1.0))
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A)')'>>> PERFORMING FINAL ITERATION WITH COMBINED EVEN/ODD VOLUMES'
        endif
        ! iteration-dependent updates
        if( params%l_doshift .and. .not.self%job_descr%isthere('trs') )then
            str = real2str(params%trs)
            call self%job_descr%set( 'trs', str )
            call cline%set( 'trs', params%trs )
        endif
        call str%kill
        call vol%kill
        call vol_iter%kill
        call fsc_templ%kill
        call fsc_file%kill
        call fname_vol%kill
        call volpproc%kill
        call vollp%kill
        if( allocated(res) ) deallocate(res)
        if( allocated(fsc) ) deallocate(fsc)
        if( allocated(state_pops) ) deallocate(state_pops)
        if( L_BENCH_GLOB ) self%bench%rt_tot = toc(self%bench%t_tot)
    end subroutine distr_execute_iteration

    subroutine distr_finalize_iteration(self, params, build)
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        type(string) :: benchfname
        integer :: fnr
        if( L_BENCH_GLOB )then
            benchfname = refine3D_distr_bench_fname(params%which_iter)
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'refine3D init, distributed              : ', self%bench%rt_init
            write(fnr,'(a,t52,f9.2)') 'refine3D prob tab, distributed       : ', self%bench%rt_prob
            write(fnr,'(a,t52,f9.2)') 'refine3D 3D align & rec, distributed : ', self%bench%rt_sched
            write(fnr,'(a,t52,f9.2)') 'refine3D assemble parts              : ', self%bench%rt_assemble
            write(fnr,'(a,t52,f9.2)') 'refine3D total time                  : ', self%bench%rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,t52,f9.2)') 'refine3D init, distributed              : ', &
                &(self%bench%rt_init/self%bench%rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'refine3D prob tab, distributed       : ', &
                &(self%bench%rt_prob/self%bench%rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'refine3D 3D align & rec, distributed : ', &
                &(self%bench%rt_sched/self%bench%rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'refine3D assemble parts              : ', &
                &(self%bench%rt_assemble/self%bench%rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'refine3D % accounted for             : ',&
                &((self%bench%rt_init + self%bench%rt_prob + self%bench%rt_sched + &
                &self%bench%rt_assemble)/self%bench%rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine distr_finalize_iteration

    subroutine distr_finalize_run(self, params, build, cline)
        use simple_commanders_euclid, only: commander_calc_group_sigmas
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        ! assemble sigma2 for next run
        if( trim(params%objfun).eq.'euclid' )then
            call self%cline_calc_group_sigmas%set('which_iter', params%which_iter + 1)
            call self%cline_calc_group_sigmas%set('nthr',       self%nthr_master)
            call xcalc_group_sigmas%execute(self%cline_calc_group_sigmas)
        endif
        if(trim(params%oritype).eq.'cls3D') call build%spproj%map2ptcls
        ! safest to write the whole thing here as multiple fields updated
        call build%spproj%write
        ! report last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(params%which_iter))
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params)
        use simple_qsys_funs, only: qsys_cleanup
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        call qsys_cleanup(params)
        call self%job_descr%kill
    end subroutine distr_cleanup

end module simple_refine3D_strategy
