!@descr: high-level search routines for the cluster2D and abinitio2D applications
module simple_strategy2D_matcher
use simple_pftc_srch_api
use simple_classaverager
use simple_binoris_io,              only: binwrite_oritab
use simple_progress,                only: progressfile_update
use simple_strategy2D_alloc,        only: clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob
use simple_builder,                 only: builder
use simple_qsys_funs,               only: qsys_job_finished
use simple_strategy2D,              only: strategy2D, strategy2D_per_ptcl
use simple_matcher_pftc_prep,       only: prep_pftc4align2D, prep_pftc4align2D_polar
use simple_matcher_smpl_and_lplims, only: set_bp_range2d, sample_ptcls4update2D
use simple_matcher_ptcl_batch,      only: prep_batch_particles2D, build_batch_particles2D, clean_batch_particles2D
use simple_strategy2D_greedy,       only: strategy2D_greedy
use simple_strategy2D_greedy_tree,  only: strategy2D_greedy_tree
use simple_strategy2D_greedy_smpl,  only: strategy2D_greedy_smpl
use simple_strategy2D_inpl,         only: strategy2D_inpl
use simple_strategy2D_inpl_smpl,    only: strategy2D_inpl_smpl
use simple_strategy2D_snhc,         only: strategy2D_snhc
use simple_strategy2D_snhc_ptree,   only: strategy2D_snhc_ptree
use simple_strategy2D_single_ptree, only: strategy2D_single_ptree
use simple_strategy2D_snhc_smpl,    only: strategy2D_snhc_smpl
use simple_strategy2D_prob,         only: strategy2D_prob
use simple_strategy2D_srch,         only: strategy2D_spec
use simple_strategy2D_tseries,      only: strategy2D_tseries
use simple_eul_prob_tab2D,          only: eul_prob_tab2D
implicit none

public :: cluster2D_exec
public :: set_b_p_ptrs2D
public :: ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
private
#include "simple_local_flags.inc"

type(image),   allocatable :: ptcl_imgs(:), ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
class(builder),    pointer :: b_ptr => null()
class(parameters), pointer :: p_ptr => null()
real(timer_int_kind)       :: rt_startup, rt_prep_batch_particles2D, rt_prep_pftc_refs2D
real(timer_int_kind)       :: rt_prep_strategy2D_batch
real(timer_int_kind)       :: rt_build_batch_particles2D, rt_align, rt_cavg, rt_projio, rt_tot
integer(timer_int_kind)    :: t, t_startup, t_prep_batch_particles2D, t_prep_pftc_refs2D
integer(timer_int_kind)    :: t_prep_strategy2D_batch
integer(timer_int_kind)    :: t_build_batch_particles2D, t_align, t_cavg, t_projio, t_tot
type(string)               :: benchfname

type :: cluster2D_ctrl
    character(len=:), allocatable :: refine_flag
    logical :: l_partial_sums
    logical :: l_update_frac
    logical :: l_ctf
    logical :: l_snhc
    logical :: l_stream
    logical :: l_greedy
    logical :: l_np_cls_defined
    logical :: l_alloc_read_cavgs
    logical :: l_prob_align
    logical :: l_polar
    logical :: l_restore_cavgs
    logical :: do_bench
end type cluster2D_ctrl

contains

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
        real    :: frac_srch_space, neigh_frac
        integer :: iptcl, fnr, updatecnt, iptcl_map, iptcl_batch, ibatch, nptcls2update
        integer :: batchsz_max, batchsz, nbatches, batch_start, batch_end
        p_ptr => params
        b_ptr => build
        call init_ctrl()
        if( ctrl%do_bench )then
            t_startup = tic()
            t_tot     = t_startup
        endif
        frac_srch_space = b_ptr%spproj_field%get_avg('frac')
        call sample_particles_for_update()
        call compute_neigh_frac()
        if( file_exists(p_ptr%frcs) ) call b_ptr%clsfrcs%read(p_ptr%frcs)
        call prepare_batches()
        call ensure_even_odd_partition()
        call prepare_class_averages_and_restoration()
        call set_bp_range2D(p_ptr, b_ptr, cline, which_iter, frac_srch_space)
        if( ctrl%do_bench )then
            rt_startup                 = toc(t_startup)
            rt_build_batch_particles2D = 0.0
            rt_prep_strategy2D_batch   = 0.0
            t_prep_batch_particles2D   = tic()
        endif
        call prep_batch_particles2D(p_ptr, b_ptr, batchsz_max, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        if( ctrl%do_bench )then
            rt_prep_batch_particles2D = toc(t_prep_batch_particles2D)
            t_prep_pftc_refs2D        = tic()
        endif
        call prepare_alignment_references(batchsz_max)
        if( ctrl%do_bench ) rt_prep_pftc_refs2D = toc(t_prep_pftc_refs2D)
        if( ctrl%l_polar )then
            if( which_iter == 1 ) call b_ptr%pftc%polar_cavger_new(.false.)
            call b_ptr%pftc%polar_cavger_zero_pft_refs
        endif
        call prep_strategy2D_glob(p_ptr, b_ptr%spproj, b_ptr%pftc%get_nrots(), neigh_frac)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> STRATEGY2D OBJECTS ALLOCATED'
        if( ctrl%l_prob_align )then
            call eulprob_obj_part%new(p_ptr, b_ptr, pinds)
            call eulprob_obj_part%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        endif
        call b_ptr%spproj_field%set_all2single('w', 1.0)
        allocate(strategy2Dsrch(batchsz_max))
        rt_align = 0.0
        rt_cavg  = 0.0
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles_local()
            if( ctrl%do_bench ) t_prep_strategy2D_batch = tic()
            call prep_strategy2D_batch( p_ptr, b_ptr%spproj, which_iter, batchsz, pinds(batch_start:batch_end) )
            if( ctrl%do_bench ) rt_prep_strategy2D_batch = rt_prep_strategy2D_batch + toc(t_prep_strategy2D_batch)
            if( ctrl%do_bench ) t_align = tic()
            !$omp parallel do private(iptcl,iptcl_batch,iptcl_map,updatecnt,orientation,strategy2Dspec)&
            !$omp default(shared) schedule(static) proc_bind(close)
            do iptcl_batch = 1, batchsz
                iptcl_map  = batch_start + iptcl_batch - 1
                iptcl      = pinds(iptcl_map)
                updatecnt  = b_ptr%spproj_field%get_updatecnt(iptcl)
                call allocate_strategy_for_particle(iptcl, iptcl_batch, updatecnt)
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
                call strategy2Dsrch(iptcl_batch)%ptr%kill
            enddo
            !$omp end parallel do
            if( ctrl%do_bench )then
                rt_align = rt_align + toc(t_align)
                t_cavg   = tic()
            endif
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
            ctrl%refine_flag     = trim(p_ptr%refine)
            ctrl%l_snhc          = str_has_substr(ctrl%refine_flag, 'snhc')
            ctrl%l_greedy        = str_has_substr(ctrl%refine_flag, 'greedy')
            ctrl%l_stream        = (trim(p_ptr%stream2d) == 'yes')
            ctrl%l_update_frac   = p_ptr%l_update_frac
            ctrl%l_partial_sums  = ctrl%l_update_frac
            ctrl%l_prob_align    = p_ptr%l_prob_align_mode
            ctrl%l_polar         = p_ptr%l_polar
            ctrl%l_restore_cavgs = (trim(p_ptr%restore_cavgs) == 'yes')
            ctrl%l_np_cls_defined= cline%defined('nptcls_per_cls')
            ctrl%do_bench        = L_BENCH_GLOB
            if( p_ptr%extr_iter == 1 )then
                ctrl%l_greedy       = .true.
                ctrl%l_snhc         = .false.
                ctrl%l_partial_sums = .false.
            else if( p_ptr%extr_iter > p_ptr%extr_lim )then
                if( trim(ctrl%refine_flag) == 'snhc_smpl' ) ctrl%refine_flag = 'snhc'
            endif
            if( ctrl%l_stream )then
                ctrl%l_update_frac = .false.
                if( (which_iter > 1) .and. (p_ptr%update_frac < 0.99) )then
                    p_ptr%l_update_frac = .true.
                    ctrl%l_partial_sums = .true.
                else
                    p_ptr%update_frac   = 1.0
                    p_ptr%l_update_frac = .false.
                    ctrl%l_partial_sums = .false.
                endif
                if( trim(ctrl%refine_flag) == 'snhc' ) ctrl%refine_flag = 'snhc_smpl'
            endif
            ctrl%l_ctf = b_ptr%spproj%get_ctfflag('ptcl2D',iptcl=p_ptr%fromp) .ne. 'no'
        end subroutine init_ctrl

        subroutine sample_particles_for_update()
            if( allocated(pinds) ) deallocate(pinds)
            if( ctrl%l_prob_align )then
                call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top], nptcls2update, pinds)
            else
                call sample_ptcls4update2D(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], ctrl%l_update_frac, nptcls2update, pinds)
            endif
        end subroutine sample_particles_for_update

        subroutine compute_neigh_frac()
            neigh_frac = 0.0
            if( p_ptr%extr_iter > p_ptr%extr_lim )then
                ! done
            else
                if( ctrl%l_snhc )then
                    neigh_frac = extremal_decay2D( p_ptr%extr_iter, p_ptr%extr_lim )
                    if( L_VERBOSE_GLOB ) write(logfhandle,'(A,F8.2)') &
                        '>>> STOCHASTIC NEIGHBOURHOOD SIZE(%):', 100.0*(1.0-neigh_frac)
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
            if( ctrl%l_polar .and. which_iter > 1 )then
                ! references are read in prep_pftc4align2D_polar below
            else
                ctrl%l_alloc_read_cavgs = l_distr_worker_glob .or. (which_iter == 1)
                call cavger_new(p_ptr, b_ptr, alloccavgs=ctrl%l_alloc_read_cavgs)
                if( ctrl%l_alloc_read_cavgs )then
                    if( .not. cline%defined('refs') )then
                        THROW_HARD('need refs to be part of command line for cluster2D execution')
                    endif
                    call cavger_read_all
                endif
            endif
            if( .not. ctrl%l_polar ) call cavger_init_online(batchsz_max, ctrl%l_partial_sums)
        end subroutine prepare_class_averages_and_restoration

        subroutine prepare_alignment_references(batchsz_max)
            integer, intent(in) :: batchsz_max
            if( ctrl%l_polar .and. which_iter > 1 )then
                call prep_pftc4align2D_polar(p_ptr, b_ptr, batchsz_max, which_iter, ctrl%l_stream)
            else
                call prep_pftc4align2D(p_ptr, b_ptr, ptcl_match_imgs_pad, batchsz_max, which_iter, ctrl%l_stream)
            endif
        end subroutine prepare_alignment_references

        subroutine build_batch_particles_local()
            if( ctrl%do_bench ) t_build_batch_particles2D = tic()
            call build_batch_particles2D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
            if( ctrl%do_bench ) rt_build_batch_particles2D = rt_build_batch_particles2D + toc(t_build_batch_particles2D)
        end subroutine build_batch_particles_local

        subroutine allocate_strategy_for_particle(iptcl, iptcl_batch, updatecnt)
            integer, intent(in) :: iptcl, iptcl_batch, updatecnt
            logical :: first_or_unsearched
            first_or_unsearched = (updatecnt == 1 .or. (.not. b_ptr%spproj_field%has_been_searched(iptcl)))
            if( ctrl%l_prob_align )then
                allocate(strategy2D_prob :: strategy2Dsrch(iptcl_batch)%ptr)
            else if( ctrl%l_stream )then
                if( first_or_unsearched )then
                    allocate(strategy2D_greedy :: strategy2Dsrch(iptcl_batch)%ptr)
                else
                    select case(trim(ctrl%refine_flag))
                    case('greedy')
                        allocate(strategy2D_greedy       :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('greedy_tree')
                        allocate(strategy2D_greedy_tree  :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('greedy_smpl')
                        allocate(strategy2D_greedy_smpl  :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_ptree')
                        allocate(strategy2D_snhc_ptree   :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('single_ptree')
                        allocate(strategy2D_single_ptree :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_smpl')
                        allocate(strategy2D_snhc_smpl    :: strategy2Dsrch(iptcl_batch)%ptr)
                    case default
                        allocate(strategy2D_snhc         :: strategy2Dsrch(iptcl_batch)%ptr)
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
                        case('greedy_tree')
                            allocate(strategy2D_greedy_tree :: strategy2Dsrch(iptcl_batch)%ptr)
                        case('greedy_smpl')
                            allocate(strategy2D_greedy_smpl :: strategy2Dsrch(iptcl_batch)%ptr)
                        case default
                            allocate(strategy2D_greedy      :: strategy2Dsrch(iptcl_batch)%ptr)
                        end select
                    endif
                else
                    select case(trim(ctrl%refine_flag))
                    case('snhc_ptree')
                        allocate(strategy2D_snhc_ptree   :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('single_ptree')
                        allocate(strategy2D_single_ptree :: strategy2Dsrch(iptcl_batch)%ptr)
                    case('snhc_smpl')
                        allocate(strategy2D_snhc_smpl    :: strategy2Dsrch(iptcl_batch)%ptr)
                    case default
                        allocate(strategy2D_snhc         :: strategy2Dsrch(iptcl_batch)%ptr)
                    end select
                endif
            endif
        end subroutine allocate_strategy_for_particle

        subroutine restore_class_averages_for_batch()
            if( ctrl%l_polar )then
                call b_ptr%pftc%polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end), &
                    b_ptr%spproj, incr_shifts(:,1:batchsz))
            else
                call cavger_transf_oridat(batchsz, pinds(batch_start:batch_end))
                call cavger_update_sums(batchsz, ptcl_imgs(1:batchsz))
            endif
        end subroutine restore_class_averages_for_batch

        subroutine cleanup_search_state(strategy2Dsrch, pinds, batches, eulprob_obj_part, batchsz_max, orientation)
            type(strategy2D_per_ptcl), allocatable, intent(inout) :: strategy2Dsrch(:)
            integer, allocatable, intent(inout) :: pinds(:), batches(:,:)
            type(eul_prob_tab2D), intent(inout) :: eulprob_obj_part
            integer, intent(in) :: batchsz_max
            type(ori), intent(inout) :: orientation
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
            if( ctrl%do_bench ) t_projio = tic()
            call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, b_ptr%spproj_field, &
                [p_ptr%fromp,p_ptr%top], isegment=PTCL2D_SEG)
            p_ptr%oritab = p_ptr%outfile
            if( ctrl%do_bench ) rt_projio = toc(t_projio)
        end subroutine write_orientations

        subroutine finalize_restoration_and_convergence(states, cline, conv, which_iter, converged)
            real, allocatable, intent(inout) :: states(:)
            class(cmdline), intent(inout) :: cline
            type(convergence), intent(inout) :: conv
            integer, intent(in) :: which_iter
            logical, intent(inout) :: converged
            if( l_distr_worker_glob )then
                if( ctrl%l_restore_cavgs )then
                    if( ctrl%l_polar )then
                        call b_ptr%pftc%polar_cavger_readwrite_partial_sums('write')
                    else
                        call cavger_readwrite_partial_sums('write')
                    endif
                endif
                call cavger_kill
                call b_ptr%pftc%polar_cavger_kill
            else
                converged = conv%check_conv2D(p_ptr, cline, b_ptr%spproj_field, b_ptr%spproj_field%get_n('class'), p_ptr%msk)
                converged = converged .and. (p_ptr%which_iter >= p_ptr%minits)
                converged = converged .or.  (p_ptr%which_iter >= p_ptr%maxits)
                if(.not. ctrl%l_stream) call progressfile_update(conv%get('progress'))
                if( ctrl%l_restore_cavgs )then
                    if( cline%defined('which_iter') )then
                        p_ptr%refs      = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//MRC_EXT
                        p_ptr%refs_even = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//'_even'//MRC_EXT
                        p_ptr%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//'_odd'//MRC_EXT
                    else
                        THROW_HARD('which_iter expected to be part of command line in shared-memory execution')
                    endif

                    if( ctrl%l_polar )then
                        if( which_iter == 1 ) call cavger_kill
                        call b_ptr%pftc%polar_cavger_merge_eos_and_norm2D(b_ptr%clsfrcs, string(FRCS_FILE))
                        call b_ptr%pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
                        call b_ptr%pftc%polar_cavger_gen2Dclassdoc(b_ptr%spproj, b_ptr%clsfrcs)
                        call b_ptr%pftc%polar_cavger_kill
                    else
                        call cavger_restore_cavgs( p_ptr%frcs )
                        call cavger_gen2Dclassdoc
                        call cavger_write_merged( p_ptr%refs )
                        if( ctrl%l_stream )then
                            call cavger_write_eo( p_ptr%refs_even, p_ptr%refs_odd )
                            call cavger_readwrite_partial_sums( 'write' )
                        endif
                        call cavger_kill(dealloccavgs=.false.)
                    endif
                    call cline%set('refs', p_ptr%refs)
                    call b_ptr%spproj%os_cls3D%new(p_ptr%ncls, is_ptcl=.false.)
                    states = b_ptr%spproj%os_cls2D%get_all('state')
                    call b_ptr%spproj%os_cls3D%set_all('state',states)
                    call b_ptr%spproj%write_segment_inside('cls2D', p_ptr%projfile)
                    call b_ptr%spproj%write_segment_inside('cls3D', p_ptr%projfile)
                    deallocate(states)
                endif
            endif
        end subroutine finalize_restoration_and_convergence

        subroutine maybe_write_bench(which_iter)
            integer, intent(in) :: which_iter
            if( .not. ctrl%do_bench ) return
            if( p_ptr%part /= 1 ) return
            rt_tot = toc(t_tot)
            benchfname = 'CLUSTER2D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'startup_overhead     : ', rt_startup
            write(fnr,'(a,1x,f9.2)') 'prep_batch_particles : ', rt_prep_batch_particles2D
            write(fnr,'(a,1x,f9.2)') 'prep_pftc_refs       : ', rt_prep_pftc_refs2D
            write(fnr,'(a,1x,f9.2)') 'build_batch_particles: ', rt_build_batch_particles2D
            write(fnr,'(a,1x,f9.2)') 'prep_strategy2D_batch: ', rt_prep_strategy2D_batch
            write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
            write(fnr,'(a,1x,f9.2)') 'class averaging      : ', rt_cavg
            write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', rt_projio
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'startup_overhead     : ', (rt_startup/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'prep_batch_particles : ', (rt_prep_batch_particles2D/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'prep_pftc_refs       : ', (rt_prep_pftc_refs2D/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'build_batch_particles: ', (rt_build_batch_particles2D/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'prep_strategy2D_batch: ', (rt_prep_strategy2D_batch/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'class averaging      : ', (rt_cavg/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'project file I/O     : ', (rt_projio/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for      : ', &
                ((rt_startup+rt_prep_batch_particles2D+rt_prep_pftc_refs2D+rt_build_batch_particles2D+ &
                  rt_prep_strategy2D_batch+rt_align+rt_cavg+rt_projio)/rt_tot) * 100.

            call fclose(fnr)
        end subroutine maybe_write_bench

    end subroutine cluster2D_exec

end module simple_strategy2D_matcher