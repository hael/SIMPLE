module simple_refine3D_strategy
use simple_core_module_api
use simple_builder,       only: builder
use simple_parameters,    only: parameters
use simple_cmdline,       only: cmdline
use simple_qsys_env,      only: qsys_env
use simple_convergence,   only: convergence
use simple_decay_funs,    only: inv_cos_decay, cos_decay
use simple_cluster_seed,  only: gen_labelling
use simple_euclid_sigma2, only: sigma2_star_from_iter
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

! BENCHMARKING VARIABLES
! initialization
integer(timer_int_kind) ::  t_init
real(timer_int_kind)    :: rt_init
! prob, scheduled jobs assemble
integer(timer_int_kind) ::  t_prob,   t_sched,  t_assemble
real(timer_int_kind)    :: rt_prob,  rt_sched, rt_assemble
! total
integer(timer_int_kind) ::  t_tot
real(timer_int_kind)    :: rt_tot

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
        integer                               :: startit
        ! Full in-memory toolbox build (required for refine3D_exec)
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        ! startit
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        ! Some refine modes manage sampling/updatecnt internally
        if( params%l_prob_align_mode )then
            ! random sampling and updatecnt dealt with in prob_align
        else
            if( startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        if( trim(params%continue) == 'yes' )then
            THROW_HARD('shared-memory implementation of refine3D does not support continue=yes')
        endif
        ! Input reference validation
        if( params%l_polar )then
            if( cline%defined('vol1') )then
                if( .not. file_exists(params%vols(1)) ) then
                    THROW_HARD('shared-memory implementation of refine3D needs starting volume input')
                endif
            else
                if( .not. file_exists(POLAR_REFS_FBODY//BIN_EXT) ) then
                    THROW_HARD('polar references are required when VOL1 not provided')
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
        endif
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
        use simple_commanders_rec_distr, only: commander_volassemble
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(inout) :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        logical,                        intent(out)   :: converged
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_prob_align)        :: xprob_align
        type(commander_prob_align_neigh)  :: xprob_align_neigh
        type(commander_volassemble)       :: xvolassemble
        type(cmdline)                     :: cline_prob_align
        type(cmdline)                     :: cline_volassemble
        integer                           :: state
        logical                           :: l_prob_state_mode, l_prob_neigh_mode
        type(string)                      :: volname, vol_in
        601 format(A,1X,F12.3)
        call self%conv%print_iteration(params%which_iter)
        ! communicate iteration counters
        call cline%set('startit',    params%which_iter)
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
        ! refine=prob* pre-step (except ptree, which runs direct tree-guided search)
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
        if( trim(params%volrec) .eq. 'yes' ) call cline%set('force_volassemble', 'yes')
        call refine3D_exec(params, build, cline, params%which_iter, converged)
        if( trim(params%volrec) .eq. 'yes' )then
            cline_volassemble = cline
            call cline_volassemble%set('which_iter', params%which_iter)
            call cline_volassemble%set('nthr',       params%nthr)
            call cline_volassemble%set('combine_eo', params%combine_eo)
            if( params%l_update_frac ) call cline_volassemble%set('update_frac', params%update_frac)
            do state = 1, params%nstates
                volname = string(VOL_FBODY)//int2str_pad(state,2)//params%ext
                if( cline_volassemble%defined('vol'//int2str(state)) )then
                    vol_in = cline_volassemble%get_carg('vol'//int2str(state))
                    if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                        if( .not. file_exists(volname) ) call cline_volassemble%delete('vol'//int2str(state))
                    endif
                endif
            end do
            call xvolassemble%execute(cline_volassemble)
            do state = 1, params%nstates
                volname = string(VOL_FBODY)//int2str_pad(state,2)//params%ext
                params%vols(state) = volname
                call cline%set('vol'//int2str(state), volname)
            end do
            call cline%delete('force_volassemble')
        endif
        ! input volume should only be used once in polar mode
        if( params%l_polar ) call cline%delete('vol1')
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration(self, params, build)
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        ! no-op (kept for symmetry)
    end subroutine inmem_finalize_iteration

    subroutine inmem_finalize_run(self, params, build, cline)
        use simple_commanders_euclid, only: commander_calc_group_sigmas
        class(refine3D_inmem_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        type(cmdline),                  intent(inout) :: cline
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(string) :: str_state, fsc_file, vol, vol_iter
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
                    str_state = int2str_pad(state,2)
                    fsc_file  = string(FSC_FBODY)//str_state//BIN_EXT
                    call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                    vol       = string(VOL_FBODY)//str_state//params%ext
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
        call str_state%kill
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
        type(string)  :: prev_refine_path, target_name, fname_vol, vol, str_state, fsc_file
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
        ! polar toolboxes needed on master for reference assembly
        if( params%l_polar )then
            ! Ensure scheduled worker command lines carry explicit angular-space size.
            ! This avoids any fallback to parser defaults in downstream private workers.
            call cline%set('nspace', params%nspace)
            call build%build_general_tbox(params, cline, do3d=.true.)
            call build%pftc%new(params, 1, [1,1], params%kfromto)
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
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' ) call build%spproj%split_stk(params%nparts, dir=string(PATH_PARENT))
        ! prepare prototype command lines
        self%cline_rec3D = cline
        self%cline_calc_pspec_distr    = cline
        self%cline_prob_align_distr    = cline
        self%cline_postprocess         = cline
        self%cline_calc_group_sigmas   = cline
        call self%cline_rec3D%set( 'prg', 'reconstruct3D' )
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
                    str_state = int2str_pad(state,2)
                    fsc_file  = FSC_FBODY//str_state%to_char()//BIN_EXT
                    if( .not.file_exists(fsc_file)) THROW_HARD('Missing file: '//fsc_file%to_char())
                end do
                if( params%l_update_frac )then
                    call simple_list_files(prev_refine_path%to_char()//'*recvol_state*part*', list)
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
                    str_state = int2str_pad(state,2)
                    fsc_file  = FSC_FBODY//str_state%to_char()//BIN_EXT
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
            vol_defined = .false.
            do state = 1,params%nstates
                vol_defined = cline%defined('vol'//int2str(state))
            enddo
            self%have_oris = .not. build%spproj%is_virgin_field(params%oritype)
            if( .not. self%have_oris )then
                call build%spproj_field%rnd_oris
                self%have_oris = .true.
                call build%spproj%write_segment_inside(params%oritype)
            endif
            if( params%l_polar )then
                if( .not.vol_defined )then
                    if( file_exists(POLAR_REFS_FBODY//BIN_EXT) ) vol_defined = .true.
                endif
            endif
            if( .not. vol_defined )then
                ! reconstructions needed
                cline_tmp = self%cline_rec3D
                call cline_tmp%delete('trail_rec')
                call cline_tmp%delete('objfun')
                call cline_tmp%delete('sigma_est')
                call cline_tmp%set('objfun', 'cc')
                call xrec3D%execute( cline_tmp )
                do state = 1,params%nstates
                    str_state = int2str_pad(state,2)
                    ! rename volumes and update cline/params
                    call simple_rename(string(VOL_FBODY)//str_state//params%ext, string(STARTVOL_FBODY)//str_state//params%ext)
                    params%vols(state) = string(STARTVOL_FBODY)//str_state//params%ext
                    vol = 'vol'//int2str(state)
                    call cline%set(vol%to_char(), params%vols(state))
                    ! keep unfiltered copies
                    call simple_copy_file(string(VOL_FBODY)//str_state//'_even'//params%ext, string(STARTVOL_FBODY)//str_state//'_even_unfil'//params%ext)
                    call simple_rename(  string(VOL_FBODY)//str_state//'_even'//params%ext, string(STARTVOL_FBODY)//str_state//'_even'//params%ext)
                    call simple_copy_file(string(VOL_FBODY)//str_state//'_odd' //params%ext, string(STARTVOL_FBODY)//str_state//'_odd_unfil' //params%ext)
                    call simple_rename(  string(VOL_FBODY)//str_state//'_odd' //params%ext, string(STARTVOL_FBODY)//str_state//'_odd' //params%ext)
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
        call str_state%kill
        call fsc_file%kill
    end subroutine distr_initialize

    subroutine distr_execute_iteration(self, params, build, cline, converged)
        use simple_commanders_rec_distr, only: commander_volassemble
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
        type(commander_volassemble)       :: xvolassemble
        type(cmdline) :: cline_prob_align, cline_volassemble
        type(string)  :: str, str_iter, str_state
        type(string)  :: vol, vol_iter, fsc_templ, fsc_file
        type(string)  :: fname_vol, volpproc, vollp, volname, vol_in
        real, allocatable :: res(:), fsc(:)
        integer, allocatable :: state_pops(:)
        integer :: state, iter
        logical :: l_prob_state_mode, l_prob_neigh_mode
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = tic()
        endif
        601 format(A,1X,F12.3)
        iter     = params%which_iter
        call self%conv%print_iteration(iter)
        str_iter = int2str_pad(iter,3)
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
            rt_init = toc(t_init)
            t_prob = tic()
        endif
        l_prob_state_mode = trim(params%refine) == 'prob_state'
        l_prob_neigh_mode = trim(params%refine) == 'prob_neigh'
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call cline_prob_align%set('prg', 'prob_align_neigh')
            else
                call cline_prob_align%set('prg', 'prob_align')
            endif
            call cline_prob_align%set('which_iter', iter)
            call cline_prob_align%set('startit',    iter)
            call build%spproj%write_segment_inside(params%oritype)
            if( l_prob_neigh_mode .and. (.not. l_prob_state_mode) )then
                call xprob_align_neigh_distr%execute( cline_prob_align )
            else
                call xprob_align_distr%execute( cline_prob_align )
            endif
        endif
        if( L_BENCH_GLOB )then
            rt_prob = toc(t_prob)
            t_sched = tic()
        endif
        ! update job description
        call self%job_descr%set( 'which_iter', int2str(iter))
        call cline%set(          'which_iter', iter)
        call self%job_descr%set( 'extr_iter',  int2str(params%extr_iter))
        call cline%set(          'extr_iter',  params%extr_iter)
        call self%job_descr%set( 'startit',    int2str(iter))
        call cline%set(          'startit',    iter)
        ! schedule distributed jobs
        call self%qenv%gen_scripts_and_schedule_jobs( self%job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
        ! merge alignment docs
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY)
        ! convergence
        converged = .false.
        select case(trim(params%refine))
            case('eval')
                ! nothing
            case DEFAULT
                converged = self%conv%check_conv3D(params, cline, build%spproj_field, params%msk)
                converged = converged .and. (iter >= params%startit + 2)
        end select
        ! Force termination at requested number of iterations (maxits is run-length)
        if( (iter - params%startit + 1) >= params%maxits ) converged = .true.
        ! assemble volumes, postprocess, automask
        if( L_BENCH_GLOB )then
            rt_sched   = toc(t_sched)
            t_assemble = tic()
        endif
        if( (trim(params%volrec).eq.'yes') )then
            select case(trim(params%refine))
                case('eval')
                    ! nothing
                case DEFAULT
                    cline_volassemble = cline
                    call cline_volassemble%set('which_iter', iter)
                    call cline_volassemble%set('nthr',       self%nthr_master)
                    call cline_volassemble%set('combine_eo', params%combine_eo)
                    if( params%l_update_frac ) call cline_volassemble%set('update_frac', params%update_frac)
                    do state = 1, params%nstates
                        volname = string(VOL_FBODY)//int2str_pad(state,2)//params%ext
                        if( cline_volassemble%defined('vol'//int2str(state)) )then
                            vol_in = cline_volassemble%get_carg('vol'//int2str(state))
                            if( trim(vol_in%to_char()) == trim(volname%to_char()) )then
                                if( .not. file_exists(volname) ) call cline_volassemble%delete('vol'//int2str(state))
                            endif
                        endif
                    end do
                    call xvolassemble%execute(cline_volassemble)
                    ! rename & add volumes to project & update job_descr
                    call build%spproj_field%get_pops(state_pops, 'state')
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
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
                            vol_iter  = string(VOL_FBODY)//str_state//params%ext
                            fsc_file  = string(FSC_FBODY)//str_state//BIN_EXT
                            call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                            ! FSC plot
                            res       = get_resarr(params%box_crop, params%smpd_crop)
                            fsc       = file2rarr(fsc_file)
                            fsc_templ = 'fsc_state'//str_state%to_char()//'_iter'//str_iter%to_char()
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
                        str_state = int2str_pad(state,2)
                        if( state_pops(state) == 0 ) cycle
                        call self%cline_postprocess%set('state', state)
                        call self%cline_postprocess%set('nthr',  self%nthr_master)
                        if( cline%defined('lp') ) call self%cline_postprocess%set('lp', params%lp)
                        call xpostprocess%execute(self%cline_postprocess)
                        volpproc = string(VOL_FBODY)//str_state//PPROC_SUFFIX//params%ext%to_char()
                        vollp    = string(VOL_FBODY)//str_state//LP_SUFFIX//params%ext%to_char()
                        ! keep per-iteration postprocessed copies
                        vol_iter = string(VOL_FBODY)//str_state//'_iter'//int2str_pad(iter,3)//PPROC_SUFFIX//params%ext%to_char()
                        call simple_copy_file(volpproc, vol_iter)
                        vol_iter = string(VOL_FBODY)//str_state//'_iter'//int2str_pad(iter,3)//LP_SUFFIX//params%ext%to_char()
                        call simple_copy_file(vollp, vol_iter)
                        if( iter > 1 .and. params%keepvol.eq.'no' )then
                            call del_file(string(VOL_FBODY)//str_state//'_iter'//int2str_pad(iter-1,3)//PPROC_SUFFIX//params%ext%to_char())
                            call del_file(string(VOL_FBODY)//str_state//'_iter'//int2str_pad(iter-1,3)//LP_SUFFIX//params%ext%to_char())
                        endif
                    enddo
            end select
            if( L_BENCH_GLOB ) rt_assemble = toc(t_assemble)
        endif
        ! polar references assembly
        if( params%l_polar )then
            params%refs = string(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//params%ext%to_char()
            call build%pftc%polar_cavger_new(.true., nrefs=params%nspace)
            call build%pftc%polar_cavger_calc_pops(build%spproj)
            call build%pftc%polar_cavger_assemble_sums_from_parts
            call build%pftc%polar_cavger_merge_eos_and_norm(build%eulspace, build%pgrpsyms, cline, build%spproj_field%get_update_frac())
            call build%pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
            call build%pftc%polar_cavger_kill
            call self%job_descr%delete('vol1')
            call cline%delete('vol1')
            if( iter > 1 .and. params%keepvol.eq.'no' )then
                call del_file(string(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext%to_char())
            endif
            if( L_BENCH_GLOB ) rt_assemble = toc(t_assemble)
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
        call str_iter%kill
        call str_state%kill
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
        if( L_BENCH_GLOB ) rt_tot = toc(t_tot)
    end subroutine distr_execute_iteration

    subroutine distr_finalize_iteration(self, params, build)
        class(refine3D_distr_strategy), intent(inout) :: self
        type(parameters),               intent(in)    :: params
        type(builder),                  intent(inout) :: build
        type(string) :: benchfname
        integer :: fnr
        if( L_BENCH_GLOB )then
            benchfname = 'DISTR_REFINE3D_BENCH_ITER'//int2str_pad(params%which_iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation              : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'prob tab, distributed       : ', rt_prob
            write(fnr,'(a,1x,f9.2)') '3D align & rec, distributed : ', rt_sched
            write(fnr,'(a,1x,f9.2)') 'assemble parts              : ', rt_assemble
            write(fnr,'(a,1x,f9.2)') 'total time                  : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation              : ', (rt_init/rt_tot)     * 100.
            write(fnr,'(a,1x,f9.2)') 'prob tab, distributed       : ', (rt_prob/rt_tot)     * 100.
            write(fnr,'(a,1x,f9.2)') '3D align & rec, distributed : ', (rt_sched/rt_tot)    * 100.
            write(fnr,'(a,1x,f9.2)') 'assemble parts           : ', (rt_assemble/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for             : ',&
                &((rt_init+rt_prob+rt_sched+rt_assemble)/rt_tot) * 100.
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
