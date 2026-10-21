!@descr: refine2D execution strategies: shared memory and distributed master
module simple_refine2D_strategy
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_builder,     only: builder
use simple_parameters,  only: parameters
use simple_cmdline,     only: cmdline
use simple_eo_partition, only: ensure_field_halves
use simple_qsys_env,    only: qsys_env
use simple_convergence, only: convergence
use simple_matcher_smpl_and_lplims, only: refine2D_requires_full_assignment, all_active_ptcls_2D_assigned, &
    &refine2D_blends_carryover
use simple_cavg_sums,   only: cavg_sums, cavg_contrib_fname, CAVG_SUMS_STATE, CAVG_SUMS_OK
use simple_classaverager, only: cavger_commit_carryover, cavger_restore_cavgs, cavger_gen2Dclassdoc, &
    &cavger_write_merged, cavger_write_eo, cavger_kill
use simple_ptcl_cache,  only: ptcl_cache_ensure
use simple_gui_utils,   only: mrc2jpeg_tiled
use simple_progress,    only: progressfile_update
use simple_euclid_sigma2, only: sigma2_commit_iter
use simple_strategy_sigma2, only: prepare_canonical_sigma_update, read_sigma2_bench_counts, commit_sigma2_exec
use simple_refine2D_fnames, only: refine2D_cavgs_fname, refine2D_cavgs_half_fname
implicit none

public :: refine2D_strategy, refine2D_inmem_strategy, refine2D_distr_strategy, create_refine2D_strategy
private
#include "simple_local_flags.inc"

!> Minimal strategy interface - only divergent operations
type, abstract :: refine2D_strategy
contains
    procedure(init_interface),          deferred :: initialize
    procedure(exec_iter_interface),     deferred :: execute_iteration
    procedure(finalize_iter_interface), deferred :: finalize_iteration
    procedure(finalize_run_interface),  deferred :: finalize_run
    procedure(cleanup_interface),       deferred :: cleanup
end type refine2D_strategy

!> Timings of an iteration, from its start to the end of finalize_iteration (the REFINE2D_STRATEGY bench record)
type :: refine2D_bench_state
    integer(timer_int_kind) :: t_tot        = 0
    integer(timer_int_kind) :: t_sigma_prep = 0
    integer(timer_int_kind) :: t_sigma      = 0
    integer(timer_int_kind) :: t            = 0   !< start of the phase being timed
    real(timer_int_kind)    :: rt_tot        = 0.
    real(timer_int_kind)    :: rt_sigma_prep = 0.
    real(timer_int_kind)    :: rt_sigma      = 0.
    real(timer_int_kind)    :: rt_setup      = 0.  !< job description, stale artifacts, the carried-sums check
    real(timer_int_kind)    :: rt_prestep    = 0.  !< prob_align2D
    real(timer_int_kind)    :: rt_search     = 0.  !< shared memory: the matcher; distributed: scheduling and search
    real(timer_int_kind)    :: rt_merge      = 0.  !< merging the partitions' alignment documents
    real(timer_int_kind)    :: rt_cavg       = 0.  !< distributed class-average assembly
    real(timer_int_kind)    :: rt_proj_io    = 0.  !< the master's project segment writes and rereads
    real(timer_int_kind)    :: rt_conv       = 0.  !< convergence
    real(timer_int_kind)    :: rt_star       = 0.  !< the cls2D star file
    real(timer_int_kind)    :: rt_jpeg       = 0.  !< the class-average JPEG
    ! the canonical sigma2 store this iteration updated, and the rows its update replaced
    integer                 :: sigma2_nptcls  = 0
    integer                 :: sigma2_nshells = 0
    integer(int64)          :: sigma2_nbytes  = 0_int64
    integer                 :: sigma2_nrows_updated = 0
end type refine2D_bench_state

type, extends(refine2D_strategy) :: refine2D_inmem_strategy
    type(convergence) :: conv
    type(refine2D_bench_state), private :: bench
contains
    procedure :: initialize         => inmem_initialize
    procedure :: execute_iteration  => inmem_execute_iteration
    procedure :: finalize_iteration => inmem_finalize_iteration
    procedure :: finalize_run       => inmem_finalize_run
    procedure :: cleanup            => inmem_cleanup
end type refine2D_inmem_strategy

type, extends(refine2D_strategy) :: refine2D_distr_strategy
    type(qsys_env)     :: qenv
    type(chash)        :: job_descr
    type(convergence)  :: conv
    integer            :: nthr_master
    type(refine2D_bench_state), private :: bench
contains
    procedure :: initialize         => distr_initialize
    procedure :: execute_iteration  => distr_execute_iteration
    procedure :: finalize_iteration => distr_finalize_iteration
    procedure :: finalize_run       => distr_finalize_run
    procedure :: cleanup            => distr_cleanup
end type refine2D_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: refine2D_strategy, parameters, builder, cmdline
        class(refine2D_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_iter_interface(self, params, build, cline, converged)
        import :: refine2D_strategy, parameters, builder, cmdline
        class(refine2D_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
        logical,                   intent(out)   :: converged
    end subroutine exec_iter_interface

    subroutine finalize_iter_interface(self, params, build)
        import :: refine2D_strategy, parameters, builder
        class(refine2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
    end subroutine finalize_iter_interface

    subroutine finalize_run_interface(self, params, build, cline)
        import :: refine2D_strategy, parameters, builder, cmdline
        class(refine2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
    end subroutine finalize_run_interface

    subroutine cleanup_interface(self, params)
        import :: refine2D_strategy, parameters
        class(refine2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    !> Strategy selection based on command-line shape.
    function create_refine2D_strategy(cline) result(strategy)
        class(cmdline),   intent(in) :: cline
        class(refine2D_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
            allocate(refine2D_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(refine2D_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_refine2D_strategy

    ! ======================================================================
    ! SHARED-MEMORY STRATEGY METHODS
    ! ======================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(refine2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer :: startit
        call cline%set('outfile', ALGN_FBODY//int2str_pad(params%part,params%numlen)//METADATA_EXT)
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call cline%set('mkdir', 'no')
        params%which_iter = max(1, params%startit)
        ! the halves are on disk before the references, whose class averages may seed the sigma2 state
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call init_refine2D_refs(cline, params, build)
        if( .not. cline%defined('part') ) call ensure_canonical_sigma2(params, build, cline)
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! once per refine2D invocation, before any iteration reads particles; may
        ! fall back to cache=no on cline, which child commands inherit via their copy
        call ptcl_cache_ensure(params, build, cline)
    end subroutine inmem_initialize

    subroutine inmem_execute_iteration( self, params, build, cline, converged)
        use simple_strategy2D_matcher, only: refine2D_exec
        use simple_starproject,        only: starproject
        use simple_commanders_euclid,  only: commander_commit_sigma2
        class(refine2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        logical,                         intent(out)   :: converged
        type(commander_commit_sigma2) :: xcommit_sigma2
        type(starproject) :: starproj
        type(cmdline)     :: cline_prob_align
        type(convergence) :: conv     !< a fresh convergence object every iteration
        real              :: update_frac_req
        logical           :: l_full_update, l_update_frac_req, l_frac
        if( L_BENCH_GLOB )then
            self%bench    = refine2D_bench_state()
            self%bench%t_tot = tic()
        endif
        call self%conv%print_iteration(params%which_iter)
        call cline%set('startit',    params%startit)
        call cline%set('which_iter', params%which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        ! without usable carried class sums this iteration updates every particle
        l_full_update = carryover_needs_full_update(params)
        if( l_full_update )then
            update_frac_req     = params%update_frac
            l_update_frac_req   = params%l_update_frac
            params%update_frac   = 1.0
            params%l_update_frac = .false.
        endif
        if( L_BENCH_GLOB ) self%bench%rt_setup = toc(self%bench%t_tot)
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            call cline_prob_align%set('prg', 'prob_align2D')
            call cline_prob_align%set('which_iter', params%which_iter)
            call cline_prob_align%set('startit',    params%startit)
            if( l_full_update ) call cline_prob_align%set('update_frac', 1.0)
            call run_prob_align2D(params, build, cline_prob_align, self%bench)
        endif
        ! main clustering/alignment step
        if( params%cc_objfun==OBJFUN_EUCLID )then
            if( L_BENCH_GLOB ) self%bench%t_sigma_prep = tic()
            call prepare_canonical_sigma_update(params, build)
            if( L_BENCH_GLOB ) self%bench%rt_sigma_prep = toc(self%bench%t_sigma_prep)
        endif
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call refine2D_exec(params, build, cline, params%which_iter, l_frac)
        if( L_BENCH_GLOB )then
            self%bench%rt_search = toc(self%bench%t)
            self%bench%t         = tic()
        endif
        call restore_cavgs_inmem(params, build, cline, l_frac)
        if( L_BENCH_GLOB )then
            self%bench%rt_cavg = toc(self%bench%t)
            self%bench%t       = tic()
        endif
        ! before update_frac is restored: check_conv2D reads params%l_update_frac
        call check_refine2D_convergence(params, build, cline, conv, converged)
        if( L_BENCH_GLOB ) self%bench%rt_conv = toc(self%bench%t)
        if( l_full_update )then
            params%update_frac   = update_frac_req
            params%l_update_frac = l_update_frac_req
        endif
        ! Euclid sigma2 consolidation for next iteration
        if( params%cc_objfun==OBJFUN_EUCLID )then
            call cline%set('which_iter', sigma2_commit_iter(params%which_iter, matcher_completed=.true.))
            if( L_BENCH_GLOB ) self%bench%t_sigma = tic()
            call xcommit_sigma2%execute(cline)
            if( L_BENCH_GLOB )then
                self%bench%rt_sigma = toc(self%bench%t_sigma)
                call read_sigma2_bench_counts(build, self%bench%sigma2_nptcls, self%bench%sigma2_nshells,&
                &self%bench%sigma2_nbytes, self%bench%sigma2_nrows_updated)
            endif
            call cline%set('which_iter', params%which_iter)
        endif
        ! Write starfile
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call starproj%export_cls2D(build%spproj, params%which_iter)
        if( L_BENCH_GLOB ) self%bench%rt_star = toc(self%bench%t)
        ! cleanup
        call starproj%kill
        call cline_prob_align%kill
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration( self, params, build)
        class(refine2D_inmem_strategy), intent(inout)   :: self
        type(parameters),                intent(in)     :: params
        type(builder),                   intent(inout)  :: build
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call gen_jpeg(params%which_iter)
        if( L_BENCH_GLOB )then
            self%bench%rt_jpeg = toc(self%bench%t)
            self%bench%rt_tot  = toc(self%bench%t_tot)
        endif
        call write_strategy_bench_report(params, self%bench, 'shared-memory')
    end subroutine inmem_finalize_iteration

    subroutine inmem_cleanup( self, params)
        class(refine2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        ! No cleanup needed for in-memory strategy, but could add here if needed in the future
    end subroutine inmem_cleanup

    subroutine inmem_finalize_run( self, params, build, cline )
        class(refine2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        type(string) :: finalcavgs
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            finalcavgs = refine2D_cavgs_fname(params%which_iter, MRC_EXT)
            call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
            call finalcavgs%kill
            call build%spproj%write_segment_inside('out', params%projfile)
        endif
        call cline%set('endit', params%which_iter)
    end subroutine inmem_finalize_run

    ! ======================================================================
    ! DISTRIBUTED STRATEGY METHODS
    ! ======================================================================

    subroutine distr_initialize( self, params, build, cline )
        use simple_exec_helpers, only: set_master_num_threads
        class(refine2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call cline%set('mkdir', 'no')
        params%which_iter = max(1, params%startit)
        ! the halves are on disk before the references, whose class averages may seed the sigma2 state
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call init_refine2D_refs(cline, params, build)
        call ensure_canonical_sigma2(params, build, cline)
        if( params%startit == 1 )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! master-side, before any job is scheduled, so the workers find it ready;
        ! must precede gen_job_descr below so a fallback to cache=no reaches the
        ! worker command lines
        call ptcl_cache_ensure(params, build, cline)
        call set_master_num_threads(self%nthr_master, string('REFINE2D'))
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute_iteration( self, params, build, cline, converged )
        use simple_stream_utils,         only: terminate_stream
        !$ use omp_lib,                  only: omp_get_max_threads, omp_set_num_threads
        class(refine2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        logical,                         intent(out)   :: converged
        type(cmdline)                     :: cline_prob_align
        integer                           :: nthr_here, nthr_glob_here
        type(string)                      :: update_frac_req
        real                              :: frac_srch_space
        logical                           :: l_full_update, l_update_frac_key
        if( L_BENCH_GLOB )then
            self%bench    = refine2D_bench_state()
            self%bench%t_tot = tic()
        endif
        call self%conv%print_iteration(params%which_iter)
        ! Update job description
        call cline%set('nparts',     params%nparts)
        call cline%set('startit',    params%startit)
        call cline%set('which_iter', params%which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        call self%job_descr%set('refs',       params%refs)
        call self%job_descr%set('nparts',     int2str(params%nparts))
        call self%job_descr%set('startit',    int2str(params%startit))
        call self%job_descr%set('which_iter', int2str(params%which_iter))
        call self%job_descr%set('extr_iter',  int2str(params%extr_iter))
        call self%job_descr%set('frcs',       FRCS_FILE)
        call cleanup_distributed_iteration_artifacts(params)
        ! without usable carried class sums this iteration updates every particle
        l_full_update     = carryover_needs_full_update(params)
        l_update_frac_key = self%job_descr%isthere('update_frac')
        if( l_full_update )then
            if( l_update_frac_key ) update_frac_req = self%job_descr%get('update_frac')
            call self%job_descr%set('update_frac', '1.0')
        endif
        if( L_BENCH_GLOB ) self%bench%rt_setup = toc(self%bench%t_tot)
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            call cline_prob_align%set('prg', 'prob_align2D')
            call cline_prob_align%set('which_iter', params%which_iter)
            call cline_prob_align%set('startit',    params%startit)
            if( l_full_update ) call cline_prob_align%set('update_frac', 1.0)
            call prob_align2D_master(params, build, cline_prob_align, l_full_update, self%bench)
        endif
        ! The master creates one generation-scoped candidate after every
        ! project-state update and before workers read the committed state.
        ! Each worker emits only its exclusive global particle range; the
        ! commit_sigma2 barrier below validates exact coverage and is the
        ! sole owner of the global reduction and atomic publication.
        if( params%cc_objfun==OBJFUN_EUCLID )then
            if( L_BENCH_GLOB ) self%bench%t_sigma_prep = tic()
            call prepare_canonical_sigma_update(params, build)
            if( L_BENCH_GLOB ) self%bench%rt_sigma_prep = toc(self%bench%t_sigma_prep)
        endif
        ! Schedule distributed jobs
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, &
                                                     algnfbody=string(ALGN_FBODY), &
                                                     array=L_USE_SLURM_ARR, &
                                                     extra_params=params)
        if( L_BENCH_GLOB ) self%bench%rt_search = toc(self%bench%t)
        call terminate_stream(params, 'SIMPLE_DISTR_REFINE2D HARD STOP 1')
        if( l_full_update )then
            if( l_update_frac_key )then
                call self%job_descr%set('update_frac', update_frac_req%to_char())
            else
                call self%job_descr%delete('update_frac')
            endif
        endif
        ! Merge alignment docs
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
        if( L_BENCH_GLOB ) self%bench%rt_merge = toc(self%bench%t)
        ! Assemble class averages
        if( trim(params%restore_cavgs) .eq. 'yes' )then
            if( L_BENCH_GLOB ) self%bench%t = tic()
            call run_distributed_cavg_assembly(params, build, cline, self%nthr_master)
            if( L_BENCH_GLOB ) self%bench%rt_cavg = toc(self%bench%t)
        endif
        ! Sigma2 consolidation
        if( params%cc_objfun==OBJFUN_EUCLID )then
            if( L_BENCH_GLOB ) self%bench%t_sigma = tic()
            ! on the master's parameters and builder, at the master's thread count; which_iter is the iteration the
            ! committed store feeds
            nthr_here      = self%nthr_master
            !$ nthr_here   = omp_get_max_threads()
            nthr_glob_here = nthr_glob
            call commit_sigma2_exec(params, build, sigma2_commit_iter(params%which_iter, matcher_completed=.true.), &
                &self%nthr_master)
            !$ call omp_set_num_threads(nthr_here)
            nthr_glob = nthr_glob_here
            if( L_BENCH_GLOB )then
                self%bench%rt_sigma = toc(self%bench%t_sigma)
                call read_sigma2_bench_counts(build, self%bench%sigma2_nptcls, self%bench%sigma2_nshells,&
                &self%bench%sigma2_nbytes, self%bench%sigma2_nrows_updated)
            endif
        endif
        ! Check convergence
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call check_refine2D_convergence(params, build, cline, self%conv, converged)
        frac_srch_space = 0.
        if( params%which_iter > 1 ) frac_srch_space = self%conv%get('frac_srch')
        ! Activate shift search if needed
        if( params%which_iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. params%l_doshift) )then
            if( .not. self%job_descr%isthere('trs') )then
                call self%job_descr%set('trs', real2str(params%trs))
            endif
        endif
        if( L_BENCH_GLOB ) self%bench%rt_conv = toc(self%bench%t)
    end subroutine distr_execute_iteration

    subroutine distr_finalize_iteration( self, params, build )
        class(refine2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        if( L_BENCH_GLOB ) self%bench%t = tic()
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        if( L_BENCH_GLOB )then
            self%bench%rt_proj_io = self%bench%rt_proj_io + toc(self%bench%t)
            self%bench%t          = tic()
        endif
        call gen_jpeg(params%which_iter)
        if( L_BENCH_GLOB )then
            self%bench%rt_jpeg = toc(self%bench%t)
            self%bench%rt_tot  = toc(self%bench%t_tot)
        endif
        call write_strategy_bench_report(params, self%bench, 'distributed')
    end subroutine distr_finalize_iteration

    subroutine distr_cleanup( self, params )
        use simple_qsys_funs, only: qsys_cleanup
        class(refine2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        call self%qenv%kill()
        call qsys_cleanup(params)
        call self%job_descr%kill
    end subroutine distr_cleanup

    subroutine distr_finalize_run( self, params, build, cline )
        class(refine2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        type(string) :: finalcavgs
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            finalcavgs = refine2D_cavgs_fname(params%which_iter, MRC_EXT)
            call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
            call finalcavgs%kill
            call build%spproj%write_segment_inside('out', params%projfile)
        endif
        call cline%set('endit', params%which_iter)
    end subroutine distr_finalize_run

    ! private helpers

    !> The canonical sigma2 state at refinement entry, through the one rebuild entry; a rebuild
    !! registers the state in the project file, so projinfo is re-read
    subroutine ensure_canonical_sigma2( params, build, cline )
        use simple_sigma2_bootstrap, only: ensure_sigma2_for_iteration
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(in)    :: cline
        logical :: l_rebuilt
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, params%startit, params%box, &
            &params%smpd, 'REFINE2D SIGMA2', l_rebuilt)
        if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
    end subroutine ensure_canonical_sigma2

    !> Shared-memory restoration of the class averages from the sums the matcher accumulated: the carried sums
    !! committed (l_frac: the iteration blends them), the averages restored and written, the class document, and the
    !! cls2D and cls3D segments written to the project
    subroutine restore_cavgs_inmem( params, build, cline, l_frac )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline
        logical,          intent(in)    :: l_frac
        real, allocatable :: states(:)
        if( trim(params%restore_cavgs) /= 'yes' ) return
        if( cline%defined('which_iter') )then
            params%refs      = refine2D_cavgs_fname(params%which_iter, MRC_EXT)
            params%refs_even = refine2D_cavgs_half_fname(params%which_iter, 'even', MRC_EXT)
            params%refs_odd  = refine2D_cavgs_half_fname(params%which_iter, 'odd', MRC_EXT)
        else
            THROW_HARD('which_iter expected to be part of command line in shared-memory execution')
        endif
        call cavger_commit_carryover(l_frac)
        call cavger_restore_cavgs( params%frcs )
        call cavger_gen2Dclassdoc
        call cavger_write_merged( params%refs )
        call cavger_write_eo( params%refs_even, params%refs_odd )
        call cavger_kill
        call cline%set('refs', params%refs)
        call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
        states = build%spproj%os_cls2D%get_all('state')
        call build%spproj%os_cls3D%set_all('state',states)
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        call build%spproj%write_segment_inside('cls3D', params%projfile)
        deallocate(states)
    end subroutine restore_cavgs_inmem

    !> The convergence decision of a refine2D iteration in both modes, after the class averages are restored
    !! (cavger_gen2Dclassdoc writes the particles' res that check_conv2D averages): check_conv2D's criteria, minits,
    !! maxits and the full-assignment coverage; prints the verdict and updates the progress file
    subroutine check_refine2D_convergence( params, build, cline, conv, converged )
        type(parameters),  intent(inout) :: params
        type(builder),     intent(inout) :: build
        type(cmdline),     intent(inout) :: cline
        type(convergence), intent(inout) :: conv
        logical,           intent(out)   :: converged
        integer :: n_unassigned
        logical :: l_full_assignment
        converged = conv%check_conv2D(params, cline, build%spproj_field, build%spproj_field%get_n('class'), params%msk)
        converged = converged .and. (params%which_iter >= params%minits)
        converged = converged .or.  (params%which_iter >= params%maxits)
        if( refine2D_requires_full_assignment(params) )then
            l_full_assignment = all_active_ptcls_2D_assigned(build%spproj_field, [params%fromp,params%top], n_unassigned)
            if( .not. l_full_assignment )then
                write(logfhandle,'(A,I8)') &
                    '>>> REFINE2D FULL-ASSIGNMENT COVERAGE: UNASSIGNED ACTIVE PARTICLES =', n_unassigned
            endif
            converged = converged .and. l_full_assignment
        endif
        call conv%print_verdict(converged)
        if( trim(params%stream2d) .eq. 'no' ) call progressfile_update(conv%get('progress'))
    end subroutine check_refine2D_convergence

    subroutine write_strategy_bench_report( params, bench, execution_mode )
        use simple_matcher_ptcl_io, only: write_mem_bench_context
        use simple_bench_io,        only: open_bench_report
        type(parameters),           intent(in) :: params
        type(refine2D_bench_state), intent(in) :: bench
        character(len=*),           intent(in) :: execution_mode
        integer :: fnr
        real(timer_int_kind) :: rt_accounted
        if( .not. L_BENCH_GLOB ) return
        rt_accounted = 0.
        if( bench%rt_tot > TINY ) rt_accounted = 100. * (bench%rt_setup + bench%rt_prestep + bench%rt_sigma_prep +&
            &bench%rt_search + bench%rt_merge + bench%rt_cavg + bench%rt_proj_io + bench%rt_sigma + bench%rt_conv +&
            &bench%rt_star + bench%rt_jpeg) / bench%rt_tot
        call open_bench_report(fnr, 'REFINE2D_STRATEGY', iter=params%which_iter)
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,a)')  'refine2D execution mode             : ', trim(execution_mode)
        write(fnr,'(a,i0)') 'refine2D ncls                       : ', params%ncls
        write(fnr,'(a,i0)') 'refine2D sigma2 store particles     : ', bench%sigma2_nptcls
        write(fnr,'(a,i0)') 'refine2D sigma2 store shells        : ', bench%sigma2_nshells
        write(fnr,'(a,i0)') 'refine2D sigma2 store bytes         : ', bench%sigma2_nbytes
        write(fnr,'(a,i0)') 'refine2D sigma2 rows updated        : ', bench%sigma2_nrows_updated
        call write_mem_bench_context(fnr, 'refine2D')
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f0.3)') 'refine2D iteration setup            :', bench%rt_setup
        write(fnr,'(a,1x,f0.3)') 'refine2D probabilistic pre-step     :', bench%rt_prestep
        write(fnr,'(a,1x,f0.3)') 'refine2D sigma2 candidate prep      :', bench%rt_sigma_prep
        write(fnr,'(a,1x,f0.3)') 'refine2D matcher/scheduler          :', bench%rt_search
        write(fnr,'(a,1x,f0.3)') 'refine2D merge alignment documents  :', bench%rt_merge
        write(fnr,'(a,1x,f0.3)') 'refine2D class-average assembly     :', bench%rt_cavg
        write(fnr,'(a,1x,f0.3)') 'refine2D project writes and rereads :', bench%rt_proj_io
        write(fnr,'(a,1x,f0.3)') 'refine2D sigma2 commit              :', bench%rt_sigma
        write(fnr,'(a,1x,f0.3)') 'refine2D convergence                :', bench%rt_conv
        write(fnr,'(a,1x,f0.3)') 'refine2D cls2D star file            :', bench%rt_star
        write(fnr,'(a,1x,f0.3)') 'refine2D class-average JPEG         :', bench%rt_jpeg
        write(fnr,'(a,1x,f0.2)') 'refine2D total time                 :', bench%rt_tot
        write(fnr,'(a,1x,f0.2)') 'refine2D % accounted for            :', rt_accounted
        call fclose(fnr)
    end subroutine write_strategy_bench_report

    !> prob_align2D with the project segment written before and reread after it (both master project I/O)
    subroutine run_prob_align2D( params, build, cline_prob_align, bench )
        use simple_commanders_prob, only: commander_prob_align2D
        type(parameters),           intent(in)    :: params
        type(builder),              intent(inout) :: build
        type(cmdline),              intent(inout) :: cline_prob_align
        type(refine2D_bench_state), intent(inout) :: bench
        type(commander_prob_align2D) :: xprob_align2D
        if( L_BENCH_GLOB ) bench%t = tic()
        call build%spproj%write_segment_inside(params%oritype)
        if( L_BENCH_GLOB )then
            bench%rt_proj_io = bench%rt_proj_io + toc(bench%t)
            bench%t          = tic()
        endif
        call xprob_align2D%execute(cline_prob_align)
        if( L_BENCH_GLOB )then
            bench%rt_prestep = toc(bench%t)
            bench%t          = tic()
        endif
        call build%spproj%read_segment(params%oritype, params%projfile)
        if( L_BENCH_GLOB ) bench%rt_proj_io = bench%rt_proj_io + toc(bench%t)
    end subroutine run_prob_align2D

    !> prob_align2D on the distributed master's parameters and builder: the sample is drawn in its particle field (a
    !! full update draws every particle; the master's own update_frac is restored after), at the workers' thread count
    !! (the master's own is restored after)
    subroutine prob_align2D_master( params, build, cline_prob_align, l_full_update, bench )
        use simple_prob_align_strategy, only: prob_align2D_exec
        !$ use omp_lib,                 only: omp_get_max_threads, omp_set_num_threads
        type(parameters),           intent(inout) :: params
        type(builder),              intent(inout) :: build
        type(cmdline),              intent(inout) :: cline_prob_align
        logical,                    intent(in)    :: l_full_update
        type(refine2D_bench_state), intent(inout) :: bench
        real    :: update_frac_req
        integer :: nthr_here, nthr_glob_here
        logical :: l_update_frac_req
        if( L_BENCH_GLOB ) bench%t = tic()
        call cline_prob_align%set('mkdir',  'no')
        call cline_prob_align%set('stream', 'no')
        update_frac_req   = params%update_frac
        l_update_frac_req = params%l_update_frac
        if( l_full_update )then
            params%update_frac   = 1.0
            params%l_update_frac = .false.
        endif
        nthr_here      = params%nthr
        !$ nthr_here   = omp_get_max_threads()
        nthr_glob_here = nthr_glob
        call prob_align2D_exec(params, build, cline_prob_align, params%nthr)
        !$ call omp_set_num_threads(nthr_here)
        nthr_glob            = nthr_glob_here
        params%update_frac   = update_frac_req
        params%l_update_frac = l_update_frac_req
        if( L_BENCH_GLOB ) bench%rt_prestep = toc(bench%t)
    end subroutine prob_align2D_master

    !> Initialize references for refine2D (used by inmem and distr modes).
    subroutine init_refine2D_refs( cline, params, build )
        use simple_procimgstk,         only: copy_imgfile
        use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr
        use simple_commanders_imgops,  only: commander_scale
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(commander_scale) :: xscale
        type(cmdline)         :: cline_make_cavgs, cline_scalerefs
        type(string)          :: refs_sc
        logical               :: l_scale_inirefs
        if( cline%defined('refs') ) return
        cline_make_cavgs = cline
        params%refs      = 'start2Drefs'     //MRC_EXT
        params%refs_even = 'start2Drefs_even'//MRC_EXT
        params%refs_odd  = 'start2Drefs_odd' //MRC_EXT
        l_scale_inirefs  = .false.
        if( build%spproj%is_virgin_field('ptcl2D') .or. params%which_iter <= 1 )then
            if( params%tseries .eq. 'yes' )then
                call init_tseries_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            else
                call init_standard_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            endif
        else
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        endif
        if( l_scale_inirefs )then
            refs_sc = 'refs'//SCALE_SUFFIX//MRC_EXT
            call cline_scalerefs%set('stk',    params%refs)
            call cline_scalerefs%set('outstk', refs_sc)
            call cline_scalerefs%set('smpd',   params%smpd)
            call cline_scalerefs%set('newbox', params%box_crop)
            call xscale%execute(cline_scalerefs)
            call simple_rename(refs_sc, params%refs)
        endif
        call copy_imgfile(params%refs, params%refs_even, params%smpd_crop, [1,params%ncls])
        call copy_imgfile(params%refs, params%refs_odd,  params%smpd_crop, [1,params%ncls])
        call cline%set('refs', params%refs)
        call cline_make_cavgs%kill
        call cline_scalerefs%kill
        call refs_sc%kill
    end subroutine init_refine2D_refs

    subroutine cleanup_distributed_iteration_artifacts( params )
        type(parameters), intent(in) :: params
        integer :: ipart
        call del_files(DIST_FBODY,      params%nparts, ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params%nparts, ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        ! stale current-iteration class-sum contributions; the carried set stays
        do ipart = 1, params%nparts
            call del_file(cavg_contrib_fname(ipart))
        enddo
    end subroutine cleanup_distributed_iteration_artifacts

    !> A refine2D iteration that blends carried class sums needs a usable carried set: one
    !! that exists, reads back whole and matches the run's classes, box and sampling. Without
    !! one the iteration runs as a full update (every particle, no carry-over), as a fresh start
    !! does. Owner-side check: workers never read carried state.
    logical function carryover_needs_full_update( params ) result( l_full )
        type(parameters), intent(in) :: params
        type(cavg_sums) :: carried
        integer         :: status
        l_full = .false.
        if( trim(params%restore_cavgs) /= 'yes' ) return
        if( .not. refine2D_blends_carryover(params, params%which_iter) ) return
        call carried%read(string(CAVG_STATE_FILE), CAVG_SUMS_STATE, status)
        if( status == CAVG_SUMS_OK ) l_full = .not. carried%matches(params%ncls, params%box_crop, params%smpd_crop)
        if( status /= CAVG_SUMS_OK ) l_full = .true.
        call carried%kill
        if( l_full ) write(logfhandle,'(A)') '>>> REFINE2D: NO USABLE CARRIED CLASS SUMS; THIS ITERATION IS A FULL UPDATE'
    end function carryover_needs_full_update

    !> class-average assembly on the distributed master's parameters and builder: the master's cls2D, cls3D, out and
    !! particle field are updated in memory and those segments written; the master's thread count
    subroutine run_distributed_cavg_assembly( params, build, cline, nthr_master )
        use simple_cavgassemble_strategy, only: cavgassemble_exec
        use simple_stream_utils,          only: terminate_stream
        !$ use omp_lib,                   only: omp_get_max_threads, omp_set_num_threads
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline
        integer,          intent(in)    :: nthr_master
        type(cmdline)                :: cline_cavgassemble
        integer                      :: nthr_here, nthr_glob_here
        params%refs      = refine2D_cavgs_fname(params%which_iter, MRC_EXT)
        params%refs_even = refine2D_cavgs_half_fname(params%which_iter, 'even', MRC_EXT)
        params%refs_odd  = refine2D_cavgs_half_fname(params%which_iter, 'odd', MRC_EXT)
        call cline%set('refs', params%refs)
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg',  'cavgassemble')
        call cline_cavgassemble%delete('which_iter')
        call cline_cavgassemble%set('refs', params%refs)
        call cline_cavgassemble%set('nthr', nthr_master)
        call terminate_stream(params, 'SIMPLE_DISTR_REFINE2D HARD STOP 2')
        call cline_cavgassemble%set('oritype', 'ptcl2D')
        nthr_here      = nthr_master
        !$ nthr_here   = omp_get_max_threads()
        nthr_glob_here = nthr_glob
        ! the previous iteration's class FRCs for the even/odd low-resolution averaging, as in shared memory
        call cavgassemble_exec(params, build, cline_cavgassemble, 0._timer_int_kind, nthr_master, prev_frcs=.true.)
        !$ call omp_set_num_threads(nthr_here)
        nthr_glob = nthr_glob_here
        ! the segments assembly changed: the class documents, the outputs and the particles' res
        call build%spproj%write_segment_inside('cls2D',        params%projfile)
        call build%spproj%write_segment_inside('cls3D',        params%projfile)
        call build%spproj%write_segment_inside('out',          params%projfile)
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call cline_cavgassemble%kill
    end subroutine run_distributed_cavg_assembly

    subroutine init_tseries_refs( cline, params, build, cline_make_cavgs, l_scale_inirefs )
        use simple_procimgstk, only: selection_from_tseries_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: cnt, iptcl, ptclind
        if( cline%defined('nptcls_per_cls') )then
            if( build%spproj%os_ptcl2D%any_state_zero() )then
                THROW_HARD('refine2D_nano does not allow state=0 particles, prune project before execution')
            endif
            cnt = 0
            do iptcl=1,params%nptcls,params%nptcls_per_cls
                cnt = cnt + 1
                params%ncls = cnt
                do ptclind=iptcl,min(params%nptcls, iptcl + params%nptcls_per_cls - 1)
                    call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                end do
            end do
            call cline%set('ncls', params%ncls)
            call cline_make_cavgs%set('ncls', params%ncls)
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        else
            if( trim(params%refine).eq.'inpl' )then
                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                call cline%set('ncls', params%ncls)
                call cline_make_cavgs%set('ncls', params%ncls)
                call cline_make_cavgs%delete('tseries')
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            else
                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            endif
        endif
    end subroutine init_tseries_refs

    subroutine init_standard_refs( cline, params, build, cline_make_cavgs, l_scale_inirefs )
        use simple_procimgstk, only: random_selection_from_imgfile, noise_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: iptcl
        select case(trim(params%cls_init))
            case('ptcl')
                call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            case('rand')
                call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                l_scale_inirefs = .false.
            case('randcls')
                if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS')
                do iptcl=1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                    call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                    call build%spproj_field%e3set(iptcl,ran3()*360.0)
                end do
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            case('prev')
                THROW_HARD('CLS_INIT=PREV is a solve2D seeded restart; refine2D does not support it')
            case DEFAULT
                THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
        end select
    end subroutine init_standard_refs

    subroutine execute_make_cavgs( cline_make_cavgs, cline, params)
        use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr
        type(cmdline),    intent(inout) :: cline_make_cavgs
        class(cmdline),   intent(in)    :: cline
        type(parameters), intent(in)    :: params
        type(commander_make_cavgs)       :: xmake_cavgs
        type(commander_make_cavgs_distr) :: xmake_cavgs_distr
        call cline_make_cavgs%set('prg', 'make_cavgs')
        if( (params%nparts > 1) .and. (.not.cline%defined('part')) )then
            call xmake_cavgs_distr%execute(cline_make_cavgs)
        else
            call xmake_cavgs%execute(cline_make_cavgs)
        endif
    end subroutine execute_make_cavgs

    subroutine gen_jpeg( which_iter ) 
        integer, intent(in) :: which_iter
        type(string) :: mrc, jpg
        ! Generate JPEG
        mrc = string(CWD_GLOB)//'/'//refine2D_cavgs_fname(which_iter, MRC_EXT)
        jpg = string(CWD_GLOB)//'/'//refine2D_cavgs_fname(which_iter, JPG_EXT)
        if( file_exists(mrc) )then
            call mrc2jpeg_tiled(mrc, jpg)
            write(logfhandle,'(A,A)') '>>> JPEG ', jpg%to_char()
        endif
        call mrc%kill
        call jpg%kill
    end subroutine gen_jpeg

end module simple_refine2D_strategy
