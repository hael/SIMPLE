!@descr: high-level search routines for the refine2D and solve2D applications
module simple_strategy2D_matcher
use simple_pftc_srch_api
use simple_classaverager
use simple_binoris_io,               only: binwrite_oritab
use simple_eo_partition,             only: ensure_field_halves
use simple_strategy2D_alloc,         only: clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob, &
                                           s2D, set_strategy2D_stoch_bound, is_fresh_2D_start
use simple_builder,                  only: builder
use simple_qsys_funs,                only: qsys_declare_part_finished
use simple_syslib,                   only: get_peak_rss_bytes, get_node_mem_in_use_bytes
use simple_matcher_ptcl_io,          only: ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, &
    &write_ptcl_read_bench_context, ptcl_read_session, release_pass_buffer
use simple_strategy2D,               only: strategy2D, strategy2D_per_ptcl
use simple_matcher_pftc_prep,        only: prep_pftc4align2D
use simple_matcher_smpl_and_lplims,  only: set_bp_range2d, sample_ptcls4update2D, refine2D_requires_full_assignment, &
                                           refine2D_blends_carryover
use simple_matcher_ptcl_batch,       only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D, &
    &write_sigma2_bench_context, write_sigma2_bench_timings, crop_particles2D
use simple_ptcl_cache,               only: ptcl_cache_in_use, ptcl_cache_assert_ready
use simple_imgarr_utils,             only: alloc_imgarr
use simple_strategy2D_greedy,        only: strategy2D_greedy
use simple_strategy2D_greedy_smpl,   only: strategy2D_greedy_smpl
use simple_strategy2D_inpl,          only: strategy2D_inpl
use simple_strategy2D_inpl_smpl,     only: strategy2D_inpl_smpl
use simple_strategy2D_snhc,          only: strategy2D_snhc
use simple_strategy2D_snhc_smpl,     only: strategy2D_snhc_smpl
use simple_strategy2D_snhc_smpl_many,only: strategy2D_snhc_smpl_many
use simple_strategy2D_prob,          only: strategy2D_prob
use simple_strategy2D_srch,          only: strategy2D_spec
use simple_strategy2D_tseries,       only: strategy2D_tseries
use simple_eul_prob_tab2D,           only: eul_prob_tab2D
implicit none

public :: refine2D_exec
public :: set_b_p_ptrs2D
public :: ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
private
#include "simple_local_flags.inc"

type(image),   allocatable :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
class(builder),    pointer :: b_ptr => null()
class(parameters), pointer :: p_ptr => null()
real(timer_int_kind)       :: rt_startup, rt_alloc_ptcl_imgs2D, rt_prep_pftc_refs2D
real(timer_int_kind)       :: rt_build_batch_particles2D, rt_align, rt_cavg, rt_tot
real(timer_int_kind)       :: rt_cavg_interp_splat, rt_sigma2_write
real(dp)                   :: ts_search, wall_search, ts_sums, wall_sums  !< thread-seconds and wall of the regions
integer(timer_int_kind)    :: t_startup, t_alloc_ptcl_imgs2D, t_prep_pftc_refs2D
integer(timer_int_kind)    :: t_build_batch_particles2D, t_align, t_cavg, t_tot, t_sigma2_write

type :: refine2D_ctrl
    character(len=:), allocatable :: refine_flag
    logical :: l_sample_updates
    logical :: l_frac_restore
    logical :: l_ctf
    logical :: l_snhc
    logical :: l_stream
    logical :: l_greedy
    logical :: l_np_cls_defined
    logical :: l_prob_align
    logical :: l_restore_cavgs
    logical :: l_require_full_assignment
    logical :: l_cached
    logical :: l_crop_sums
    logical :: do_bench
end type refine2D_ctrl

contains

    subroutine set_b_p_ptrs2D( params, build )
        class(parameters), target, intent(in) :: params
        class(builder),    target, intent(in) :: build
        p_ptr => params
        b_ptr => build
    end subroutine set_b_p_ptrs2D

    !>  \brief  is the prime2D algorithm: samples, searches, accumulates the class sums and writes the orientations
    !!  and sigma2 rows; a distributed worker writes its class-sum contribution, shared memory returns with the sums in
    !!  the class averager for the strategy to restore. l_frac: whether this iteration blends carried class sums
    subroutine refine2D_exec( params, build, cline, which_iter, l_frac )
        use simple_decay_funs,  only: extremal_decay2D
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        class(cmdline),            intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(out)   :: l_frac
        type(strategy2D_per_ptcl), allocatable   :: strategy2Dsrch(:)
        real,                      allocatable   :: incr_shifts(:,:)
        integer,                   allocatable   :: pinds(:), batches(:,:)
        type(eul_prob_tab2D),      target        :: eulprob_obj_part
        type(refine2D_ctrl)                      :: ctrl
        type(ori)             :: orientation
        type(strategy2D_spec) :: strategy2Dspec
        logical :: attempted, improved, no_improvement, invalid
        real    :: frac_srch_space, neigh_frac
        integer :: iptcl, fnr, iptcl_map, iptcl_batch, ibatch, nptcls2update
        integer :: batchsz_max, batchsz, nbatches, batch_start, batch_end, ithr, nthr_team
        real(dp), allocatable  :: ts_thr(:)
        real(dp)               :: t_region, t_ptcl
        type(ptcl_read_counts) :: reads0, reads_srch
        type(ptcl_read_session) :: reads
        p_ptr => params
        b_ptr => build
        call init_ctrl()
        if( trim(p_ptr%inpl_cont) == 'yes' )then
            call b_ptr%spproj_field%set_all2single('cont_inpl_attempted', 0.)
            call b_ptr%spproj_field%set_all2single('cont_inpl_improved',  0.)
        endif
        if( ctrl%do_bench )then
            t_startup = tic()
            t_tot     = t_startup
        endif
        frac_srch_space = b_ptr%spproj_field%get_avg('frac')
        call sample_particles_for_update()
        call compute_neigh_frac( neigh_frac )
        if( file_exists(p_ptr%frcs) ) call b_ptr%clsfrcs%read(p_ptr%frcs)
        call prepare_batches()
        call ensure_even_odd_partition()
        call prepare_class_averages_and_restoration()
        call set_bp_range2D(p_ptr, b_ptr, cline, which_iter, frac_srch_space)
        if( ctrl%do_bench )then
            rt_startup                 = toc(t_startup)
            rt_build_batch_particles2D = 0.0
            t_alloc_ptcl_imgs2D        = tic()
        endif
        ! ptcl_imgs holds the images handed to class-average restoration; cropped (from the cache, or made in
        ! memory) they are at box_crop, and from the cache the read buffer shrinks to box_crop too
        if( ctrl%l_cached )then
            call alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max, &
                &imgbatch_box=p_ptr%box_crop, imgbatch_smpd=p_ptr%smpd_crop)
            call alloc_imgarr(batchsz_max, [p_ptr%box_crop, p_ptr%box_crop, 1], p_ptr%smpd_crop, ptcl_imgs)
        else if( ctrl%l_crop_sums )then
            call alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
            call alloc_imgarr(batchsz_max, [p_ptr%box_crop, p_ptr%box_crop, 1], p_ptr%smpd_crop, ptcl_imgs)
        else
            call alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
            call alloc_imgarr(batchsz_max, [p_ptr%box, p_ptr%box, 1], p_ptr%smpd, ptcl_imgs)
        endif
        if( ctrl%do_bench )then
            rt_alloc_ptcl_imgs2D = toc(t_alloc_ptcl_imgs2D)
            t_prep_pftc_refs2D   = tic()
        endif
        call set_strategy2D_stoch_bound(params%ncls, neigh_frac)
        call prepare_alignment_references(batchsz_max)
        if( ctrl%do_bench ) rt_prep_pftc_refs2D = toc(t_prep_pftc_refs2D)
        call prep_strategy2D_glob(p_ptr, b_ptr%spproj, b_ptr%pftc%get_nrots(), neigh_frac)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'
        if( ctrl%l_prob_align )then
            call eulprob_obj_part%new_assignment(p_ptr,b_ptr,pinds)
            call eulprob_obj_part%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        endif
        allocate(strategy2Dsrch(batchsz_max))
        if( ctrl%do_bench ) reads0 = ptcl_read_counts_now()
        call reads%new(p_ptr)
        call reads%expect(pinds(1:nptcls2update))
        rt_align = 0.0
        rt_cavg  = 0.0
        rt_cavg_interp_splat = 0.0
        ! busy time per thread of the team the search runs on
        nthr_team = omp_get_max_threads()
        allocate(ts_thr(nthr_team), source=0.d0)
        ts_search   = 0.d0
        wall_search = 0.d0
        ts_sums     = 0.d0
        wall_sums   = 0.d0
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles_local()
            call prep_strategy2D_batch( p_ptr, b_ptr%spproj, which_iter, batchsz, pinds(batch_start:batch_end) )
            if( ctrl%do_bench ) t_align = tic()
            t_region = omp_get_wtime()
            ! Preserve initialized lifecycle flags and pointers in each thread's working objects. A particle's cost
            ! varies, so the schedule is dynamic: a particle is searched whole by one thread, so the order changes
            ! only which thread, and so which random draws, it meets
            !$omp parallel do private(iptcl,iptcl_batch,iptcl_map,attempted,improved,no_improvement,invalid,ithr,t_ptcl)&
            !$omp firstprivate(orientation,strategy2Dspec)&
            !$omp default(shared) schedule(dynamic) proc_bind(close)
            do iptcl_batch = 1, batchsz
                t_ptcl     = omp_get_wtime()
                iptcl_map  = batch_start + iptcl_batch - 1
                iptcl      = pinds(iptcl_map)
                call allocate_strategy_for_particle(iptcl, iptcl_batch)
                strategy2Dspec%iptcl       = iptcl
                strategy2Dspec%iptcl_batch = iptcl_batch
                strategy2Dspec%iptcl_map   = iptcl_map
                strategy2Dspec%stoch_bound = neigh_frac
                if( ctrl%l_prob_align ) strategy2Dspec%eulprob_obj_part2D => eulprob_obj_part
                call strategy2Dsrch(iptcl_batch)%ptr%new(p_ptr, strategy2Dspec, b_ptr)
                call strategy2Dsrch(iptcl_batch)%ptr%srch(b_ptr%spproj_field)
                incr_shifts(:,iptcl_batch) = strategy2Dsrch(iptcl_batch)%ptr%s%best_shvec
                if ( p_ptr%cc_objfun == OBJFUN_EUCLID ) then
                    call b_ptr%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch))
                    call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'class')
                end if
                if( trim(p_ptr%inpl_cont) == 'yes' )then
                    call strategy2Dsrch(iptcl_batch)%ptr%s%get_continuous_route_status( &
                        &attempted, improved, no_improvement, invalid)
                    call b_ptr%spproj_field%set(iptcl, 'cont_inpl_attempted', merge(1., 0., attempted))
                    call b_ptr%spproj_field%set(iptcl, 'cont_inpl_improved',  merge(1., 0., improved))
                endif
                call strategy2Dsrch(iptcl_batch)%ptr%kill
                ithr = omp_get_thread_num() + 1
                ts_thr(ithr) = ts_thr(ithr) + (omp_get_wtime() - t_ptcl)
            enddo
            !$omp end parallel do
            wall_search = wall_search + (omp_get_wtime() - t_region)
            if( ctrl%do_bench )then
                rt_align = rt_align + toc(t_align)
                t_cavg   = tic()
            endif
            call restore_class_averages_for_batch()
            if( ctrl%do_bench ) rt_cavg = rt_cavg + toc(t_cavg)
        enddo
        call reads%kill
        ! the particle buffer the shared-memory probability table filled for this iteration
        call release_pass_buffer
        ts_search = sum(ts_thr)
        deallocate(ts_thr)
        if( ctrl%do_bench ) reads_srch = ptcl_read_counts_since(reads0)
        call cleanup_search_state(strategy2Dsrch, pinds, batches, eulprob_obj_part, batchsz_max, orientation)
        rt_sigma2_write = 0.0
        if( p_ptr%cc_objfun == OBJFUN_EUCLID )then
            if( ctrl%do_bench ) t_sigma2_write = tic()
            call b_ptr%esig%write_sigma2
            if( ctrl%do_bench ) rt_sigma2_write = toc(t_sigma2_write)
        endif
        call write_orientations()
        if( l_distr_worker_glob )then
            if( ctrl%l_restore_cavgs ) call cavger_write_contribution(ctrl%l_frac_restore)
            call cavger_kill
        endif
        l_frac = ctrl%l_frac_restore
        call b_ptr%esig%kill
        call b_ptr%pftc%kill
        call qsys_declare_part_finished(p_ptr, string('simple_strategy2D_matcher :: refine2D_exec'))
        call maybe_write_bench(which_iter)

contains

        subroutine init_ctrl()
            ctrl%refine_flag       = trim(p_ptr%refine)
            ctrl%l_snhc            = str_has_substr(ctrl%refine_flag, 'snhc')
            ctrl%l_greedy          = str_has_substr(ctrl%refine_flag, 'greedy')
            ctrl%l_stream          = (trim(p_ptr%stream2d) == 'yes')
            ctrl%l_sample_updates  = p_ptr%l_update_frac
            ctrl%l_frac_restore    = refine2D_blends_carryover(p_ptr, which_iter)
            ctrl%l_prob_align      = p_ptr%l_prob_align_mode
            ctrl%l_restore_cavgs   = (trim(p_ptr%restore_cavgs) == 'yes')
            ctrl%l_require_full_assignment = refine2D_requires_full_assignment(p_ptr)
            ctrl%l_np_cls_defined  = cline%defined('nptcls_per_cls')
            ctrl%do_bench          = L_BENCH_GLOB
            if( p_ptr%extr_iter == 1 )then
                ctrl%l_greedy       = .true.
                ctrl%l_snhc         = .false.
            else if( p_ptr%extr_iter > p_ptr%extr_lim )then
                if( trim(ctrl%refine_flag) == 'snhc_smpl' ) ctrl%refine_flag = 'snhc'
            endif
            if( ctrl%l_stream )then
                if( ctrl%l_frac_restore )then
                    p_ptr%l_update_frac   = .true.
                    ctrl%l_sample_updates = .true.
                else
                    p_ptr%update_frac     = 1.0
                    p_ptr%l_update_frac   = .false.
                    ctrl%l_sample_updates = .false.
                endif
                if( trim(ctrl%refine_flag) == 'snhc' ) ctrl%refine_flag = 'snhc_smpl'
            endif
            ctrl%l_ctf = b_ptr%spproj%get_ctfflag('ptcl2D',iptcl=p_ptr%fromp) .ne. 'no'
            ! Serve both the alignment and the class-average restoration from the
            ! downscaled cache. Decided here because it sizes the image buffers and the
            ! class-averager work images below. Mandatory when requested: a rank that
            ! quietly fell back would contribute partial sums built by a different
            ! preprocessing path than its peers.
            call ptcl_cache_assert_ready(p_ptr, b_ptr)
            ctrl%l_cached = ptcl_cache_in_use(p_ptr, b_ptr)
            if( ctrl%l_cached .and. p_ptr%part == 1 ) &
                write(logfhandle,'(A)') '>>> REFINE2D: reading particles from the downscaled cache'
            ! The class averages live at box_crop, so their sums take the Fourier-cropped particles the cache would
            ! serve, made in memory from the full-size read; the stream's 2D pool keeps the full-size path
            ctrl%l_crop_sums = crop_particles2D(p_ptr, ctrl%l_cached)
        end subroutine init_ctrl

        subroutine sample_particles_for_update()
            if( allocated(pinds) ) deallocate(pinds)
            if( ctrl%l_prob_align )then
                ! prob_align2D owns the outer subset sampling in probabilistic mode;
                ! refine2D only reproduces that same subset for the downstream update
                ! (a distributed partition may hold none of it)
                call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top], nptcls2update, pinds, &
                    &allow_empty=p_ptr%l_distr_worker)
            else
                call sample_ptcls4update2D(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], ctrl%l_sample_updates, nptcls2update, pinds)
            endif
        end subroutine sample_particles_for_update

        subroutine compute_neigh_frac(neighfrac)
            real, intent(out) :: neighfrac
            neighfrac = 0.0
            if( p_ptr%extr_iter > p_ptr%extr_lim )then
                ! done
            else
                if( ctrl%l_snhc )then
                    neighfrac = extremal_decay2D( p_ptr%extr_iter, p_ptr%extr_lim )
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A,F8.2)') &
                        '>>> STOCHASTIC NEIGHBOURHOOD SIZE(%):', 100.0*(1.0-neighfrac)
                endif
            endif
        end subroutine compute_neigh_frac

        !> batches of the particles to update; a partition without any still runs the pass (no batch), which writes
        !! its orientations, sigma2 range and empty class sums with the references' centering offsets
        subroutine prepare_batches()
            if( nptcls2update < 1 )then
                batchsz_max = 1
                nbatches    = 0
                allocate(batches(0,2))
            else
                batchsz_max = min(nptcls2update, p_ptr%nthr * BATCHTHRSZ)
                nbatches    = ceiling(real(nptcls2update) / real(batchsz_max))
                batches     = split_nobjs_even(nptcls2update, nbatches)
                batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
            endif
            allocate(incr_shifts(2,batchsz_max), source=0.0)
        end subroutine prepare_batches

        !> a worker stops without halves; shared memory makes the project's partition
        subroutine ensure_even_odd_partition()
            call ensure_field_halves(b_ptr%spproj, p_ptr%oritype, p_ptr%projfile)
        end subroutine ensure_even_odd_partition

        subroutine prepare_class_averages_and_restoration()
            call cavger_new(p_ptr, b_ptr)
            if( .not. cline%defined('refs') )then
                THROW_HARD('need refs to be part of command line for refine2D execution')
            endif
            call cavger_read_all
            call cavger_init_online(batchsz_max, cropped_ptcls=ctrl%l_cached .or. ctrl%l_crop_sums)
        end subroutine prepare_class_averages_and_restoration

        subroutine prepare_alignment_references(batchsz_max)
            integer, intent(in) :: batchsz_max
            if( str_has_substr(ctrl%refine_flag, '_many') )then
                call prep_pftc4align2D(p_ptr, b_ptr, ptcl_match_imgs_pad, batchsz_max, which_iter, &
                                        &nmany_refs=s2D%snhc_nrefs_bound)
            else
                call prep_pftc4align2D(p_ptr, b_ptr, ptcl_match_imgs_pad, batchsz_max, which_iter)
            endif
        end subroutine prepare_alignment_references

        subroutine build_batch_particles_local()
            if( ctrl%do_bench ) t_build_batch_particles2D = tic()
            call build_batch_particles2D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), reads, &
                ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad, crop_sums=ctrl%l_crop_sums)
            if( ctrl%do_bench ) rt_build_batch_particles2D = rt_build_batch_particles2D + toc(t_build_batch_particles2D)
        end subroutine build_batch_particles_local

        subroutine allocate_strategy_for_particle(iptcl, iptcl_batch)
            integer, intent(in) :: iptcl, iptcl_batch
            logical :: first_or_unsearched, has_been_searched, l_fresh_start
            has_been_searched = b_ptr%spproj_field%has_been_searched(iptcl)
            l_fresh_start     = is_fresh_2D_start(p_ptr, p_ptr%which_iter)
            first_or_unsearched = l_fresh_start .or. (.not. has_been_searched)
            if( ctrl%l_prob_align )then
                allocate(strategy2D_prob :: strategy2Dsrch(iptcl_batch)%ptr)
            else if( ctrl%l_stream )then
                if( first_or_unsearched )then
                    allocate(strategy2D_greedy :: strategy2Dsrch(iptcl_batch)%ptr)
                else
                    select case(trim(ctrl%refine_flag))
                    case('greedy')
                        allocate(strategy2D_greedy         :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('greedy_smpl')
                        allocate(strategy2D_greedy_smpl    :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_smpl')
                        allocate(strategy2D_snhc_smpl      :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_smpl_many')
                        allocate(strategy2D_snhc_smpl_many :: strategy2Dsrch(iptcl_batch)%ptr)
                    case default
                        allocate(strategy2D_snhc           :: strategy2Dsrch(iptcl_batch)%ptr)
                    end select
                endif
            else
                if( str_has_substr(ctrl%refine_flag,'inpl') )then
                    if( ctrl%refine_flag == 'inpl' )then
                        allocate(strategy2D_inpl :: strategy2Dsrch(iptcl_batch)%ptr)
                    else if( ctrl%refine_flag == 'inpl_smpl' )then
                        allocate(strategy2D_inpl_smpl :: strategy2Dsrch(iptcl_batch)%ptr)
                    endif
                else if( ctrl%l_greedy .or. first_or_unsearched )then
                    if( trim(p_ptr%tseries) == 'yes' )then
                        if( ctrl%l_np_cls_defined )then
                            allocate(strategy2D_tseries :: strategy2Dsrch(iptcl_batch)%ptr)
                        else
                            allocate(strategy2D_greedy  :: strategy2Dsrch(iptcl_batch)%ptr)
                        endif
                    else
                        select case(trim(ctrl%refine_flag))
                        case('greedy_smpl')
                            allocate(strategy2D_greedy_smpl :: strategy2Dsrch(iptcl_batch)%ptr)
                        case default
                            allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                        end select
                    endif
                else
                    select case(trim(ctrl%refine_flag))
                    case('snhc_smpl')
                        allocate(strategy2D_snhc_smpl      :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_smpl_many')
                        allocate(strategy2D_snhc_smpl_many :: strategy2Dsrch(iptcl_batch)%ptr)
                    case default
                        allocate(strategy2D_snhc           :: strategy2Dsrch(iptcl_batch)%ptr)
                    end select
                endif
            endif
        end subroutine allocate_strategy_for_particle

        subroutine restore_class_averages_for_batch()
            integer(timer_int_kind) :: t_update
            call cavger_transf_oridat(batchsz, pinds(batch_start:batch_end), updated_only=.true.)
            if( ctrl%do_bench ) t_update = tic()
            call cavger_update_sums(batchsz, ptcl_imgs(1:batchsz), ts_sums, wall_sums)
            if( ctrl%do_bench ) rt_cavg_interp_splat = rt_cavg_interp_splat + toc(t_update)
        end subroutine restore_class_averages_for_batch

        subroutine cleanup_search_state(strategy2Dsrch, pinds, batches, eulprob_obj_part, batchsz_max, orientation)
            type(strategy2D_per_ptcl), allocatable, intent(inout) :: strategy2Dsrch(:)
            integer, allocatable, intent(inout) :: pinds(:), batches(:,:)
            type(eul_prob_tab2D), intent(inout) :: eulprob_obj_part
            integer,              intent(in)    :: batchsz_max
            type(ori),            intent(inout) :: orientation
            integer :: i
            call clean_strategy2D
            call orientation%kill
            do i = 1, batchsz_max
                nullify(strategy2Dsrch(i)%ptr)
            end do
            call clean_batch_particles2D(b_ptr, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
            deallocate(strategy2Dsrch, pinds, batches)
            if( ctrl%l_prob_align ) call eulprob_obj_part%kill
            call cavger_dealloc_online
        end subroutine cleanup_search_state

        subroutine write_orientations()
            if( p_ptr%top < p_ptr%fromp )then
                THROW_HARD('invalid output write range in refine2D_exec: TOP < FROMP')
            endif
            call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, b_ptr%spproj_field, &
                [p_ptr%fromp,p_ptr%top], isegment=PTCL2D_SEG)
            p_ptr%oritab = p_ptr%outfile
        end subroutine write_orientations

        subroutine maybe_write_bench(which_iter)
            use, intrinsic :: iso_fortran_env, only: int64, real64
            use simple_bench_io, only: open_bench_report
            integer, intent(in) :: which_iter
            integer(int64) :: peak_rss
            real(real64)    :: peak_rss_gib
            if( .not. ctrl%do_bench ) return
            rt_tot = toc(t_tot)
            peak_rss = get_peak_rss_bytes()
            peak_rss_gib = -1.0_real64
            if( peak_rss >= 0_int64 ) peak_rss_gib = real(peak_rss,real64) / real(1024_int64**3,real64)
            ! every partition writes its own report, so worker time and load imbalance can be aggregated
            call open_bench_report(fnr, 'REFINE2D', iter=which_iter, part=p_ptr%part, numlen=p_ptr%numlen)
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'match2D refine mode                 : ', trim(ctrl%refine_flag)
            write(fnr,'(a,l1)') 'match2D probabilistic alignment     : ', ctrl%l_prob_align
            write(fnr,'(a,l1)') 'match2D sample updates              : ', ctrl%l_sample_updates
            write(fnr,'(a,l1)') 'match2D restore class averages      : ', ctrl%l_restore_cavgs
            write(fnr,'(a,l1)') 'match2D require full assignment     : ', ctrl%l_require_full_assignment
            write(fnr,'(a,i0)') 'match2D nclasses                    : ', p_ptr%ncls
            write(fnr,'(a,i0)') 'match2D kfrom                       : ', p_ptr%kfromto(1)
            write(fnr,'(a,i0)') 'match2D kto                         : ', p_ptr%kfromto(2)
            write(fnr,'(a,i0)') 'match2D process partition           : ', p_ptr%part
            write(fnr,'(a,i0)') 'match2D process pid                 : ', p_ptr%pid
            write(fnr,'(a,i0)') 'match2D peak RSS (bytes)            : ', peak_rss
            write(fnr,'(a,f0.3)') 'match2D peak RSS (GiB)              : ', peak_rss_gib
            write(fnr,'(a,i0)') 'match2D node memory in use (bytes)  : ', get_node_mem_in_use_bytes()
            write(fnr,'(a,l1)') 'match2D particles from the 2D cache : ', ctrl%l_cached
            write(fnr,'(a,l1)') 'match2D class sums from cropped particles : ', ctrl%l_cached .or. ctrl%l_crop_sums
            call write_sigma2_bench_context(fnr, 'match2D', nptcls2update)
            call write_ptcl_read_bench_context(fnr, 'match2D search', reads_srch)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'match2D startup/setup               :', rt_startup
            write(fnr,'(a,1x,f0.2)') 'match2D particle allocation         :', rt_alloc_ptcl_imgs2D
            write(fnr,'(a,1x,f0.2)') 'match2D reference preparation       :', rt_prep_pftc_refs2D
            write(fnr,'(a,1x,f0.2)') 'match2D particle reading            :', reads_srch%rt_read + reads_srch%rt_cache
            write(fnr,'(a,1x,f0.2)') 'match2D particle preparation        :', &
                &rt_build_batch_particles2D - reads_srch%rt_read - reads_srch%rt_cache
            write(fnr,'(a,1x,f0.2)') 'match2D alignment search            :', rt_align
            write(fnr,'(a,1x,f0.2)') 'match2D class averaging             :', rt_cavg
            write(fnr,'(a,1x,f0.2)') 'match2D cavg FFT/CTF/interpolation  :', rt_cavg_interp_splat
            call write_sigma2_bench_timings(fnr, 'match2D', rt_sigma2_write)
            write(fnr,'(a,1x,f0.2)') 'match2D total time                  :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'match2D % accounted for             :', &
                ((rt_startup + rt_alloc_ptcl_imgs2D + rt_prep_pftc_refs2D + rt_build_batch_particles2D + &
                  rt_align + rt_cavg + rt_sigma2_write) / rt_tot) * 100.
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** THREAD-SECONDS (summed over threads; not wall time) ***'
            write(fnr,'(a,1x,f0.2)') 'match2D search                      :', ts_search
            write(fnr,'(a,1x,f0.2)') 'match2D search region wall time     :', wall_search
            write(fnr,'(a,1x,i0)')   'match2D search threads              :', nthr_team
            if( wall_search > 0.d0 )then
                write(fnr,'(a,1x,f0.3)') 'match2D search load balance         :', ts_search / (wall_search * real(nthr_team,dp))
            endif
            if( ctrl%l_restore_cavgs ) call cavger_write_sums_load(fnr, 'match2D', ts_sums, wall_sums)

            call fclose(fnr)
        end subroutine maybe_write_bench

    end subroutine refine2D_exec

end module simple_strategy2D_matcher
