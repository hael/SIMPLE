!@descr: high-level search routines for the cluster2D and abinitio2D applications
module simple_strategy2D_matcher
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api,          only: dp
use simple_pftc_srch_api
use simple_classaverager
use simple_binoris_io,               only: binwrite_oritab
use simple_progress,                 only: progressfile_update
use simple_strategy2D_alloc,         only: clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob, &
                                           s2D, set_strategy2D_stoch_bound, is_fresh_2D_start
use simple_builder,                  only: builder
use simple_qsys_funs,                only: qsys_job_finished
use simple_syslib,                   only: get_peak_rss_bytes
use simple_strategy2D,               only: strategy2D, strategy2D_per_ptcl
use simple_matcher_pftc_prep,        only: prep_pftc4align2D
use simple_matcher_smpl_and_lplims,  only: set_bp_range2d, sample_ptcls4update2D, cluster2D_requires_full_assignment, &
                                           all_active_ptcls_2D_assigned
use simple_matcher_ptcl_batch,       only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
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

public :: cluster2D_exec
public :: set_b_p_ptrs2D
public :: reset_sgd_stream_dispatch_count, get_sgd_stream_dispatch_count
public :: ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
private
#include "simple_local_flags.inc"

type(image),   allocatable :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
class(builder),    pointer :: b_ptr => null()
class(parameters), pointer :: p_ptr => null()
real(timer_int_kind)       :: rt_startup, rt_alloc_ptcl_imgs2D, rt_prep_pftc_refs2D
real(timer_int_kind)       :: rt_build_batch_particles2D, rt_align, rt_cavg, rt_tot
real(timer_int_kind)       :: rt_cavg_interp_splat
integer(timer_int_kind)    :: t_startup, t_alloc_ptcl_imgs2D, t_prep_pftc_refs2D
integer(timer_int_kind)    :: t_build_batch_particles2D, t_align, t_cavg, t_tot
type(string)               :: benchfname
integer                    :: sgd_stream_dispatch_count = 0

type :: cluster2D_ctrl
    character(len=:), allocatable :: refine_flag
    logical :: l_sample_updates
    logical :: l_frac_restore
    logical :: l_ctf
    logical :: l_snhc
    logical :: l_stream
    logical :: l_greedy
    logical :: l_sgd_stream
    logical :: l_np_cls_defined
    logical :: l_prob_align
    logical :: l_restore_cavgs
    logical :: l_require_full_assignment
    logical :: l_cached
    logical :: do_bench
  contains
    procedure :: display
end type cluster2D_ctrl

contains

    subroutine reset_sgd_stream_dispatch_count()
        sgd_stream_dispatch_count = 0
    end subroutine reset_sgd_stream_dispatch_count

    integer function get_sgd_stream_dispatch_count()
        get_sgd_stream_dispatch_count = sgd_stream_dispatch_count
    end function get_sgd_stream_dispatch_count

    subroutine set_b_p_ptrs2D( params, build )
        class(parameters), target, intent(in) :: params
        class(builder),    target, intent(in) :: build
        p_ptr => params
        b_ptr => build
    end subroutine set_b_p_ptrs2D

    !>  \brief  is the prime2D algorithm
    subroutine cluster2D_exec( params, build, cline, which_iter, converged )
        use simple_convergence, only: convergence
        use simple_decay_funs,  only: extremal_decay2D
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        class(cmdline),            intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(inout) :: converged
        type(strategy2D_per_ptcl), allocatable   :: strategy2Dsrch(:)
        real,                      allocatable   :: states(:), incr_shifts(:,:)
        integer,                   allocatable   :: pinds(:), batches(:,:)
        type(eul_prob_tab2D),      target        :: eulprob_obj_part
        type(cluster2D_ctrl)                     :: ctrl
        type(ori)             :: orientation
        type(convergence)     :: conv
        type(strategy2D_spec) :: strategy2Dspec
        logical :: attempted, improved, no_improvement, invalid
        real    :: frac_srch_space, neigh_frac
        integer :: iptcl, fnr, iptcl_map, iptcl_batch, ibatch, nptcls2update
        integer :: batchsz_max, batchsz, nbatches, batch_start, batch_end
        p_ptr => params
        b_ptr => build
        call init_ctrl()
        if( trim(p_ptr%inpl_cont) == 'yes' )then
            call b_ptr%spproj_field%set_all2single('cont_inpl_attempted', 0.)
            call b_ptr%spproj_field%set_all2single('cont_inpl_improved',  0.)
        endif
        if( p_ptr%sgd_diagnostic )then
            write(logfhandle,'(A,1X,A,I0,1X,A,L1,1X,A,L1,1X,A,1X,A)') &
                '>>> SEARCH DIAG: SGD matcher dispatch:', 'iteration=', p_ptr%which_iter, &
                'requested=', p_ptr%l_sgd, 'active=', ctrl%l_sgd_stream, &
                'refine=', trim(ctrl%refine_flag)
        endif
        if( ctrl%l_sgd_stream )then
            sgd_stream_dispatch_count = sgd_stream_dispatch_count + 1
            if( p_ptr%sgd_diagnostic )then
                write(logfhandle,'(A)') '>>> JOINT2D SGD STREAM ACTIVE: assignment=hard_class_angle objective=raw_euclidean'
                write(logfhandle,'(A)') '>>> JOINT2D SGD STREAM PATH: prob_align2D=off softmax=off topk=off optimizer=bounded_direct_gradient'
                write(logfhandle,'(A,I0,1X,A,F8.4,1X,A,I0)') &
                    '>>> JOINT2D SGD STREAM SETTINGS: iteration=', p_ptr%which_iter, &
                    'eta_shift=', p_ptr%sgd_eta_shift, 'shift_steps=', p_ptr%sgd_shift_its
            endif
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
        ! ptcl_imgs holds the raw images handed to class-average restoration; from the
        ! cache they arrive already Fourier-cropped, so both buffers shrink to box_crop
        if( ctrl%l_cached )then
            call alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max, &
                &imgbatch_box=p_ptr%box_crop, imgbatch_smpd=p_ptr%smpd_crop)
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
        rt_align = 0.0
        rt_cavg  = 0.0
        rt_cavg_interp_splat = 0.0
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles_local()
            call prep_strategy2D_batch( p_ptr, b_ptr%spproj, which_iter, batchsz, pinds(batch_start:batch_end) )
            if( ctrl%do_bench ) t_align = tic()
            !$omp parallel do private(iptcl,iptcl_batch,iptcl_map,orientation,strategy2Dspec,&
            !$omp attempted,improved,no_improvement,invalid)&
            !$omp default(shared) schedule(static) proc_bind(close)
            do iptcl_batch = 1, batchsz
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
            enddo
            !$omp end parallel do
            if( ctrl%do_bench )then
                rt_align = rt_align + toc(t_align)
                t_cavg   = tic()
            endif
            if( ctrl%l_sgd_stream .and. p_ptr%sgd_diagnostic ) &
                call report_stream_batch(ibatch, batch_start, batch_end, strategy2Dsrch)
            call restore_class_averages_for_batch()
            if( ctrl%do_bench ) rt_cavg = rt_cavg + toc(t_cavg)
        enddo
        call cleanup_search_state(strategy2Dsrch, pinds, batches, eulprob_obj_part, batchsz_max, orientation)
        if( p_ptr%cc_objfun == OBJFUN_EUCLID ) call b_ptr%esig%write_sigma2
        call write_orientations()
        call finalize_restoration_and_convergence(states, cline, conv, which_iter, converged)
        call b_ptr%esig%kill
        call b_ptr%pftc%kill
        call qsys_job_finished(p_ptr, string('simple_strategy2D_matcher :: cluster2D_exec'))
        call maybe_write_bench(which_iter)

contains

        subroutine init_ctrl()
            ctrl%refine_flag       = trim(p_ptr%refine)
            ctrl%l_snhc            = str_has_substr(ctrl%refine_flag, 'snhc')
            ctrl%l_greedy          = str_has_substr(ctrl%refine_flag, 'greedy')
            ctrl%l_sgd_stream      = p_ptr%l_sgd_streaming_active
            ctrl%l_stream          = (trim(p_ptr%stream2d) == 'yes')
            ctrl%l_sample_updates  = p_ptr%l_update_frac
            ctrl%l_frac_restore    = ctrl%l_sample_updates
            ctrl%l_prob_align      = p_ptr%l_prob_align_mode .and. .not. p_ptr%l_sgd_streaming_active
            if( p_ptr%l_sgd_streaming_active )then
                ! Design A keeps the existing discrete class/angle search but
                ! consumes each score immediately; no probabilistic table is
                ! built and no SoftMax/top-K assignment is available here.
                ctrl%refine_flag = 'greedy'
                ctrl%l_greedy    = .true.
            endif
            ctrl%l_restore_cavgs   = (trim(p_ptr%restore_cavgs) == 'yes')
            ctrl%l_require_full_assignment = cluster2D_requires_full_assignment(p_ptr)
            ctrl%l_np_cls_defined  = cline%defined('nptcls_per_cls')
            ctrl%do_bench          = L_BENCH_GLOB
            if( p_ptr%startit == 1 )then
                ctrl%l_frac_restore = .false.
            endif
            if( ctrl%l_sgd_stream )then
                ! SGD uses its own mini-batch fraction, including the first
                ! active iteration; do not inherit update_frac or sticky sampling.
                ctrl%l_sample_updates = .true.
                ctrl%l_frac_restore   = .true.
            endif
            if( p_ptr%extr_iter == 1 )then
                ctrl%l_greedy       = .true.
                ctrl%l_snhc         = .false.
            else if( p_ptr%extr_iter > p_ptr%extr_lim )then
                if( trim(ctrl%refine_flag) == 'snhc_smpl' ) ctrl%refine_flag = 'snhc'
            endif
            if( ctrl%l_stream )then
                ctrl%l_sample_updates = .false.
                ctrl%l_frac_restore   = .false.
                if( ctrl%l_sgd_stream )then
                    ! sample_ptcls4update2D applies sgd_update_frac and the
                    ! existing NSAMPLE_DEFAULT_2D cap on every active iteration.
                    ctrl%l_sample_updates = .true.
                    ctrl%l_frac_restore   = .true.
                else
                    ! Preserve legacy stream behavior for non-SGD runs.
                    if( (which_iter > 1) .and. (p_ptr%update_frac < 0.99) )then
                        p_ptr%l_update_frac   = .true.
                        ctrl%l_sample_updates = .true.
                        ctrl%l_frac_restore   = .true.
                    else
                        p_ptr%update_frac     = 1.0
                        p_ptr%l_update_frac   = .false.
                        ctrl%l_sample_updates = .false.
                        ctrl%l_frac_restore   = .false.
                    endif
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
                write(logfhandle,'(A)') '>>> CLUSTER2D: reading particles from the downscaled cache'
        end subroutine init_ctrl

        subroutine sample_particles_for_update()
            if( allocated(pinds) ) deallocate(pinds)
            if( ctrl%l_prob_align )then
                ! prob_align2D owns the outer subset sampling in probabilistic mode;
                ! cluster2D only reproduces that same subset for the downstream update.
                call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top], nptcls2update, pinds)
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

        subroutine prepare_batches()
            batchsz_max = min(nptcls2update, p_ptr%nthr * BATCHTHRSZ)
            nbatches    = ceiling(real(nptcls2update) / real(batchsz_max))
            batches     = split_nobjs_even(nptcls2update, nbatches)
            batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
            allocate(incr_shifts(2,batchsz_max), source=0.0)
        end subroutine prepare_batches

        subroutine ensure_even_odd_partition()
            if( b_ptr%spproj_field%get_nevenodd() == 0 )then
                if( l_distr_worker_glob ) THROW_HARD('no eo partitioning available; cluster2D_exec')
                call b_ptr%spproj_field%partition_eo
                call b_ptr%spproj%write_segment_inside(p_ptr%oritype)
            endif
        end subroutine ensure_even_odd_partition

        subroutine prepare_class_averages_and_restoration()
            call cavger_new(p_ptr, b_ptr)
            if( .not. cline%defined('refs') )then
                THROW_HARD('need refs to be part of command line for cluster2D execution')
            endif
            call cavger_read_all
            call cavger_init_online(batchsz_max, ctrl%l_frac_restore, cropped_ptcls=ctrl%l_cached)
        end subroutine prepare_class_averages_and_restoration

        subroutine prepare_alignment_references(batchsz_max)
            integer, intent(in) :: batchsz_max
            if( str_has_substr(ctrl%refine_flag, '_many') )then
                call prep_pftc4align2D(p_ptr, b_ptr, ptcl_match_imgs_pad, batchsz_max, which_iter, ctrl%l_stream,&
                                        &ctrl%l_frac_restore, nmany_refs=s2D%snhc_nrefs_bound)
            else
                call prep_pftc4align2D(p_ptr, b_ptr, ptcl_match_imgs_pad, batchsz_max, which_iter, ctrl%l_stream,&
                                        &ctrl%l_frac_restore)
            endif
        end subroutine prepare_alignment_references

        subroutine build_batch_particles_local()
            if( ctrl%do_bench ) t_build_batch_particles2D = tic()
            call build_batch_particles2D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
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
                if( ctrl%l_sgd_stream )then
                    ! SGD stream uses a deterministic top-1 scan followed by
                    ! bounded shift descent.  Preserve the legacy strategy
                    ! dispatch for all other stream modes.
                    allocate(strategy2D_greedy :: strategy2Dsrch(iptcl_batch)%ptr)
                else if( first_or_unsearched )then
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
            call cavger_update_sums(batchsz, ptcl_imgs(1:batchsz))
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
                THROW_HARD('invalid output write range in cluster2D_exec: TOP < FROMP')
            endif
            call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, b_ptr%spproj_field, &
                [p_ptr%fromp,p_ptr%top], isegment=PTCL2D_SEG)
            p_ptr%oritab = p_ptr%outfile
        end subroutine write_orientations

        subroutine finalize_restoration_and_convergence(states, cline, conv, which_iter, converged)
            real, allocatable, intent(inout) :: states(:)
            class(cmdline),    intent(inout) :: cline
            type(convergence), intent(inout) :: conv
            integer,           intent(in)    :: which_iter
            logical,           intent(inout) :: converged
            integer :: n_unassigned
            logical :: l_full_assignment
            if( l_distr_worker_glob )then
                if( ctrl%l_restore_cavgs )then
                    call cavger_readwrite_partial_sums('write')
                endif
                call cavger_kill
            else
                if( ctrl%l_restore_cavgs )then
                    if( cline%defined('which_iter') )then
                        p_ptr%refs      = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//MRC_EXT
                        p_ptr%refs_even = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//'_even'//MRC_EXT
                        p_ptr%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//'_odd'//MRC_EXT
                    else
                        THROW_HARD('which_iter expected to be part of command line in shared-memory execution')
                    endif
                    call cavger_readwrite_partial_sums('write')
                    call cavger_restore_cavgs( p_ptr%frcs )
                    call cavger_gen2Dclassdoc
                    call cavger_write_merged( p_ptr%refs )
                    call cavger_write_eo( p_ptr%refs_even, p_ptr%refs_odd )
                    call cavger_kill
                    call cline%set('refs', p_ptr%refs)
                    call b_ptr%spproj%os_cls3D%new(p_ptr%ncls, is_ptcl=.false.)
                    states = b_ptr%spproj%os_cls2D%get_all('state')
                    call b_ptr%spproj%os_cls3D%set_all('state',states)
                    call b_ptr%spproj%write_segment_inside('cls2D', p_ptr%projfile)
                    call b_ptr%spproj%write_segment_inside('cls3D', p_ptr%projfile)
                    deallocate(states)
                endif
                converged = conv%check_conv2D(p_ptr, cline, b_ptr%spproj_field, b_ptr%spproj_field%get_n('class'), p_ptr%msk)
                converged = converged .and. (p_ptr%which_iter >= p_ptr%minits)
                converged = converged .or.  (p_ptr%which_iter >= p_ptr%maxits)
                if( ctrl%l_require_full_assignment )then
                    l_full_assignment = all_active_ptcls_2D_assigned(b_ptr%spproj_field, [p_ptr%fromp,p_ptr%top], n_unassigned)
                    if( .not. l_full_assignment )then
                        write(logfhandle,'(A,I8)') &
                            '>>> CLUSTER2D FULL-ASSIGNMENT COVERAGE: UNASSIGNED ACTIVE PARTICLES =', n_unassigned
                    endif
                    converged = converged .and. l_full_assignment
                endif
                if(.not. ctrl%l_stream) call progressfile_update(conv%get('progress'))
            endif
        end subroutine finalize_restoration_and_convergence

        subroutine maybe_write_bench(which_iter)
            use, intrinsic :: iso_fortran_env, only: int64, real64
            integer, intent(in) :: which_iter
            integer(int64) :: peak_rss
            real(real64)    :: peak_rss_gib
            if( .not. ctrl%do_bench ) return
            if( p_ptr%part /= 1 ) return
            rt_tot = toc(t_tot)
            peak_rss = get_peak_rss_bytes()
            peak_rss_gib = -1.0_real64
            if( peak_rss >= 0_int64 ) peak_rss_gib = real(peak_rss,real64) / real(1024_int64**3,real64)
            benchfname = string('CLUSTER2D_BENCH_ITER')//int2str_pad(which_iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
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
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'match2D startup/setup               :', rt_startup
            write(fnr,'(a,1x,f0.2)') 'match2D particle allocation         :', rt_alloc_ptcl_imgs2D
            write(fnr,'(a,1x,f0.2)') 'match2D reference preparation       :', rt_prep_pftc_refs2D
            write(fnr,'(a,1x,f0.2)') 'match2D particle preparation        :', rt_build_batch_particles2D
            write(fnr,'(a,1x,f0.2)') 'match2D alignment search            :', rt_align
            write(fnr,'(a,1x,f0.2)') 'match2D class averaging             :', rt_cavg
            write(fnr,'(a,1x,f0.2)') 'match2D cavg FFT/CTF/interpolation  :', rt_cavg_interp_splat
            write(fnr,'(a,1x,f0.2)') 'match2D total time                  :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'match2D % accounted for             :', &
                ((rt_startup + rt_alloc_ptcl_imgs2D + rt_prep_pftc_refs2D + rt_build_batch_particles2D + &
                  rt_align + rt_cavg) / rt_tot) * 100.

            call fclose(fnr)
        end subroutine maybe_write_bench

    end subroutine cluster2D_exec

    subroutine display( self )
        class(cluster2D_ctrl), intent(in) :: self
        write(logfhandle,'(a)') '>>> CLUSTER2D CONTROL FLAGS:'
        write(logfhandle,'(a,a)') 'refine_flag           : ', trim(self%refine_flag)
        write(logfhandle,'(a,l1)') 'l_sample_updates     : ', self%l_sample_updates
        write(logfhandle,'(a,l1)') 'l_frac_restore       : ', self%l_frac_restore
        write(logfhandle,'(a,l1)') 'l_ctf                : ', self%l_ctf
        write(logfhandle,'(a,l1)') 'l_snhc               : ', self%l_snhc
        write(logfhandle,'(a,l1)') 'l_stream             : ', self%l_stream
        write(logfhandle,'(a,l1)') 'l_greedy             : ', self%l_greedy
        write(logfhandle,'(a,l1)') 'l_sgd_stream         : ', self%l_sgd_stream
        write(logfhandle,'(a,l1)') 'l_np_cls_defined     : ', self%l_np_cls_defined
        write(logfhandle,'(a,l1)') 'l_prob_align         : ', self%l_prob_align
        write(logfhandle,'(a,l1)') 'l_restore_cavgs      : ', self%l_restore_cavgs
        write(logfhandle,'(a,l1)') 'l_require_full_assignment : ', self%l_require_full_assignment
        write(logfhandle,'(a,l1)') 'do_bench             : ', self%do_bench
    end subroutine display

    subroutine report_stream_batch(ibatch, batch_start, batch_end, strategy2Dsrch)
        integer, intent(in) :: ibatch, batch_start, batch_end
        type(strategy2D_per_ptcl), intent(in) :: strategy2Dsrch(:)
        integer :: j, nused, naccepted, nimproved, nfinite, nbounded, nnonmonotonic
        real(dp) :: initial_sum, final_sum, max_abs_shift
        nused = 0
        naccepted = 0
        nimproved = 0
        nfinite = 0
        nbounded = 0
        nnonmonotonic = 0
        initial_sum = 0.
        final_sum = 0.
        max_abs_shift = 0.
        do j = 1, batch_end - batch_start + 1
            if( .not. associated(strategy2Dsrch(j)%ptr) ) cycle
            if( .not. strategy2Dsrch(j)%ptr%s%sgd_used ) cycle
            nused = nused + 1
            naccepted = naccepted + strategy2Dsrch(j)%ptr%s%sgd_accepted_steps
            initial_sum = initial_sum + strategy2Dsrch(j)%ptr%s%sgd_objective_initial
            final_sum = final_sum + strategy2Dsrch(j)%ptr%s%sgd_objective_final
            if( ieee_is_finite(strategy2Dsrch(j)%ptr%s%sgd_objective_initial) .and.&
                &ieee_is_finite(strategy2Dsrch(j)%ptr%s%sgd_objective_final) ) nfinite = nfinite + 1
            max_abs_shift = max(max_abs_shift, maxval(abs(real(strategy2Dsrch(j)%ptr%s%best_shvec,dp))))
            if( maxval(abs(strategy2Dsrch(j)%ptr%s%best_shvec)) <= strategy2Dsrch(j)%ptr%s%trs + 10.*epsilon(1.) )then
                nbounded = nbounded + 1
            endif
            if( strategy2Dsrch(j)%ptr%s%sgd_objective_final < &
                strategy2Dsrch(j)%ptr%s%sgd_objective_initial ) nimproved = nimproved + 1
            if( strategy2Dsrch(j)%ptr%s%sgd_objective_final > &
                strategy2Dsrch(j)%ptr%s%sgd_objective_initial + &
                100._dp*epsilon(1._dp)*max(1._dp,abs(strategy2Dsrch(j)%ptr%s%sgd_objective_initial)) )then
                nnonmonotonic = nnonmonotonic + 1
            endif
        enddo
        if( nused > 0 )then
            write(logfhandle,'(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,ES12.4,1X,A,ES12.4,1X,A,ES12.4)') &
                '>>> JOINT2D SGD STREAM BATCH: batch=', ibatch, 'particles=', batch_end-batch_start+1, &
                'updates=', nused, 'finite=', nfinite, 'bounded=', nbounded, 'accepted_steps=', naccepted, &
                'improved=', nimproved, 'nonmonotonic=', nnonmonotonic, &
                'mean_initial=', initial_sum/real(nused), 'mean_final=', final_sum/real(nused), &
                'max_abs_shift=', max_abs_shift
        else
            write(logfhandle,'(A,I0,1X,A,I0,1X,A)') &
                '>>> JOINT2D SGD STREAM BATCH: batch=', ibatch, 'particles=', batch_end-batch_start+1, &
                'updates=0'
        endif
    end subroutine report_stream_batch

end module simple_strategy2D_matcher
