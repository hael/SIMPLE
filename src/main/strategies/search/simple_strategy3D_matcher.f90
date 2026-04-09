module simple_strategy3D_matcher
use simple_pftc_srch_api
use simple_strategy3D_alloc,        only: clean_strategy3D, prep_strategy3D, s3D
use simple_binoris_io,              only: binwrite_oritab
use simple_builder,                 only: builder
use simple_convergence,             only: convergence
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_eul_prob_tab,            only: eul_prob_tab
use simple_matcher_2Dprep,          only: prepimg4align, prepimg4align_bench
use simple_matcher_3Drec,           only: init_rec, prep_imgs4rec, update_rec, finalize_rec
use simple_matcher_ptcl_batch,      only: prep_sigmas_alloc_ptcl_imgs, build_batch_particles3D
use simple_matcher_ptcl_io,         only: killimgbatch
use simple_matcher_pftc_prep,       only: prep_pftc4align3D_polar
use simple_matcher_refvol_utils,    only: read_mask_filter_reproject_refvols
use simple_matcher_smpl_and_lplims, only: set_bp_range3D, sample_ptcls4fillin, sample_ptcls4update3D
use simple_qsys_funs,               only: qsys_job_finished
use simple_strategy3D_eval,         only: strategy3D_eval
use simple_strategy3D_greedy_smpl,  only: strategy3D_greedy_smpl
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_prob,         only: strategy3D_prob
use simple_strategy3D_ptree_neigh_states, only: strategy3D_ptree_neigh_states
use simple_strategy3D_ptree_neigh,  only: strategy3D_ptree_neigh
use simple_strategy3D_ptree,        only: strategy3D_ptree
use simple_strategy3D_shc_ptree,    only: strategy3D_shc_ptree
use simple_strategy3D_shc_smpl,     only: strategy3D_shc_smpl
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_snhc_smpl,    only: strategy3D_snhc_smpl
use simple_strategy3D_srch,         only: strategy3D_spec
use simple_strategy3D,              only: strategy3D
implicit none

public :: refine3D_exec
private
#include "simple_local_flags.inc"

type :: refine3D_ctrl
    character(len=:), allocatable :: refine_mode
    character(len=:), allocatable :: polar_mode
    character(len=:), allocatable :: oritype
    logical :: do_restore
    logical :: do_polar_prepare
    logical :: do_prob_align
    logical :: do_polar
    logical :: do_sigma_mode
    logical :: do_write_oris
    logical :: do_bench
end type refine3D_ctrl

contains

    subroutine refine3D_exec( params, build, cline, which_iter, converged )
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(inout) :: converged
        class(parameters), pointer :: p_ptr => null()
        class(builder),    pointer :: b_ptr => null()
        type(eul_prob_tab), target :: eulprob_obj_part
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        type(strategy3D_spec),     allocatable :: strategy3Dspecs(:)
        type(image),               allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:), ptcl_rec_imgs(:)
        type(fplane_type),         allocatable :: fpls(:)
        type(class_sample),        allocatable :: clssmp(:)
        integer,                   allocatable :: batches(:,:), cnt_greedy(:), cnt_all(:), pinds(:)
        real,                      allocatable :: incr_shifts(:,:)
        type(convergence)   :: conv
        type(ori)           :: orientation
        type(refine3D_ctrl) :: ctrl
        real                :: frac_greedy
        integer             :: nbatches, batchsz_max, batch_start, batch_end, batchsz, nrefs
        integer             :: iptcl, fnr, ithr, iptcl_batch, iptcl_map, ibatch, nptcls2update
        logical             :: doprint, has_been_searched

        ! benchmarking
        type(string) :: benchfname
        integer(timer_int_kind) :: t_startup, t_build_batch_ptcls, t_prep_orisrch, t_align, t_rec, t_tot, t_projio
        integer(timer_int_kind) :: t_prep_pftc_polar_mode, t_prep_sigmas_alloc_ptcl_imgs
        integer(timer_int_kind) :: t_prep_reproj_refvols, t_memoize_refs
        real(timer_int_kind)    :: rt_startup, rt_build_batch_ptcls, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
        real(timer_int_kind)    :: rt_prep_pftc_polar_mode, rt_prep_sigmas_alloc_ptcl_imgs
        real(timer_int_kind)    :: rt_prep_reproj_refvols, rt_memoize_refs

        p_ptr => params
        b_ptr => build
        call init_ctrl()
        if( ctrl%do_bench )then
            t_startup = tic()
            t_tot     = t_startup
        endif
        call ensure_even_odd_partition()
        has_been_searched = .not. b_ptr%spproj%is_virgin_field(p_ptr%oritype)
        call set_bp_range3D(p_ptr, b_ptr, cline)
        call sample_particles_for_update()
        call prepare_batches()
        if( ctrl%do_bench )then
            rt_startup = toc(t_startup)
            rt_prep_pftc_polar_mode        = 0.0
            rt_prep_sigmas_alloc_ptcl_imgs = 0.0
            rt_prep_reproj_refvols         = 0.0
        endif
        call prepare_refs_sigmas_and_pftc()
        if( ctrl%do_bench ) t_memoize_refs = tic()
        if( .not. ctrl%do_prob_align ) call build%pftc%memoize_refs
        if( ctrl%do_bench )then
            rt_memoize_refs = toc(t_memoize_refs)
            t_prep_orisrch  = tic()
        endif
        call prep_strategy3D(p_ptr, b_ptr)
        allocate(strategy3Dspecs(batchsz_max), strategy3Dsrch(batchsz_max))
        if( ctrl%do_prob_align )then
            call eulprob_obj_part%new(p_ptr, b_ptr, pinds)
            call eulprob_obj_part%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        endif
        if( ctrl%do_bench )then
            rt_prep_orisrch     = toc(t_prep_orisrch)
            rt_build_batch_ptcls= 0.0
            rt_align            = 0.0
            rt_rec              = 0.0
        endif
        call maybe_init_reconstruction()
        allocate(cnt_greedy(p_ptr%nthr), cnt_all(p_ptr%nthr), source=0)
        allocate(incr_shifts(2,batchsz_max), source=0.0)
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles_local()
            if( ctrl%do_bench ) t_align = tic()
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map,ithr,orientation) &
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch = 1, batchsz
                iptcl_map = batch_start + iptcl_batch - 1
                iptcl     = pinds(iptcl_map)
                ithr      = omp_get_thread_num() + 1
                cnt_all(ithr) = cnt_all(ithr) + 1
                strategy3Dspecs(iptcl_batch)%iptcl     = iptcl
                strategy3Dspecs(iptcl_batch)%iptcl_map = iptcl_map
                if( ctrl%do_prob_align ) strategy3Dspecs(iptcl_batch)%eulprob_obj_part => eulprob_obj_part
                call choose_and_run_strategy(iptcl, iptcl_batch, iptcl_map, ithr, has_been_searched)
                if ( p_ptr%cc_objfun == OBJFUN_EUCLID ) then
                    call b_ptr%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch))
                    call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'proj')
                endif
            enddo
            !$omp end parallel do
            if( ctrl%do_bench ) rt_align = rt_align + toc(t_align)
            call maybe_restore_batch()
        enddo
        frac_greedy = 0.0
        if( any(cnt_greedy > 0) .and. any(cnt_all > 0) )then
            frac_greedy = real(sum(cnt_greedy)) / real(sum(cnt_all))
        endif
        call b_ptr%spproj_field%set_all2single('frac_greedy', frac_greedy)
        do iptcl_batch = 1, batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
        enddo
        deallocate(strategy3Dsrch, strategy3Dspecs, batches)
        call eulprob_obj_part%kill
        call deallocate_class_samples(clssmp)
        if( p_ptr%cc_objfun == OBJFUN_EUCLID ) call b_ptr%esig%write_sigma2
        if( trim(p_ptr%ptclw) == 'yes' )then
            ! not supported
        else
            if( trim(p_ptr%cavgw) == 'yes' )then
                call b_ptr%spproj_field%calc_cavg_soft_weights(p_ptr%frac)
            else
                call b_ptr%spproj_field%calc_hard_weights(p_ptr%frac)
            endif
        endif
        call clean_strategy3D
        call b_ptr%vol%kill
        call orientation%kill
        call dealloc_imgarr(ptcl_rec_imgs)
        call dealloc_imgarr(ptcl_match_imgs)
        call dealloc_imgarr(ptcl_match_imgs_pad)
        call killimgbatch(b_ptr)
        call maybe_write_orientations()

        if( ctrl%do_restore )then
            if( ctrl%do_bench ) t_rec = tic()
            if( ctrl%do_polar )then
                if( l_distr_worker_glob )then
                    call b_ptr%pftc%polar_cavger_readwrite_partial_sums('write')
                else
                    p_ptr%refs = CAVGS_ITER_FBODY//int2str_pad(p_ptr%which_iter,3)//MRC_EXT
                    select case(ctrl%polar_mode)
                    case('new')
                        call b_ptr%pftc%polar_cavger_merge_eos_and_norm_new(b_ptr%eulspace, &
                            b_ptr%pgrpsyms, cline, b_ptr%spproj_field%get_update_frac())
                    case default
                        call b_ptr%pftc%polar_cavger_merge_eos_and_norm(b_ptr%eulspace, &
                            b_ptr%pgrpsyms, cline, b_ptr%spproj_field%get_update_frac())
                    end select
                    call b_ptr%pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
                    call b_ptr%pftc%polar_cavger_kill
                endif
            else
                call finalize_rec(params, build, cline, fpls)
            endif
            if( ctrl%do_bench ) rt_rec = rt_rec + toc(t_rec)
        endif
        call b_ptr%pftc%kill
        call b_ptr%esig%kill
        call qsys_job_finished(p_ptr, string('simple_strategy3D_matcher :: refine3D_exec'))
        if( .not. p_ptr%l_distr_worker .and. ctrl%refine_mode /= 'sigma' )then
            converged = conv%check_conv3D(p_ptr, cline, b_ptr%spproj_field, p_ptr%msk)
        endif

        if( ctrl%do_bench )then
            rt_tot = toc(t_tot)
            doprint = .true.
            if( p_ptr%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'REFINE3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'startup_overhead : ',          rt_startup
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles3D : ',   rt_build_batch_ptcls
                write(fnr,'(a,1x,f9.2)') 'prep_pftc_polar_mode : ',      rt_prep_pftc_polar_mode
                write(fnr,'(a,1x,f9.2)') 'prep_sigmas_alloc_ptcls : ',   rt_prep_sigmas_alloc_ptcl_imgs
                write(fnr,'(a,1x,f9.2)') 'prepare_reproject_refvols : ', rt_prep_reproj_refvols
                write(fnr,'(a,1x,f9.2)') 'memoize_refs : ',              rt_memoize_refs
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation : ',     rt_prep_orisrch
                write(fnr,'(a,1x,f9.2)') '3D alignment : ',              rt_align
                write(fnr,'(a,1x,f9.2)') 'project file I/O : ',          rt_projio
                write(fnr,'(a,1x,f9.2)') 'reconstruction : ',            rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time : ',                rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'startup_overhead : ',          (rt_startup/rt_tot)                     * 100.
                write(fnr,'(a,1x,f9.2)') 'build_batch_particles3D : ',   (rt_build_batch_ptcls/rt_tot)           * 100.
                write(fnr,'(a,1x,f9.2)') 'prep_pftc_polar_mode : ',      (rt_prep_pftc_polar_mode/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') 'prep_sigmas_alloc_ptcls : ',   (rt_prep_sigmas_alloc_ptcl_imgs/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'prepare_reproject_refvols : ', (rt_prep_reproj_refvols/rt_tot)         * 100.
                write(fnr,'(a,1x,f9.2)') 'memoize_refs : ',              (rt_memoize_refs/rt_tot)                * 100.
                write(fnr,'(a,1x,f9.2)') 'orisrch3D preparation : ',     (rt_prep_orisrch/rt_tot)                * 100.
                write(fnr,'(a,1x,f9.2)') '3D alignment : ',              (rt_align/rt_tot)                       * 100.
                write(fnr,'(a,1x,f9.2)') 'project file I/O : ',          (rt_projio/rt_tot)                      * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction : ',            (rt_rec/rt_tot)                         * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for : ', &
                    ((rt_startup+rt_build_batch_ptcls+rt_prep_pftc_polar_mode+rt_prep_sigmas_alloc_ptcl_imgs+ &
                    rt_prep_reproj_refvols+rt_memoize_refs+rt_prep_orisrch+rt_align+rt_projio+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif

    contains

        subroutine init_ctrl()
            ctrl%refine_mode   = trim(p_ptr%refine)
            ctrl%polar_mode    = trim(p_ptr%polar)
            ctrl%oritype       = trim(p_ptr%oritype)
            ctrl%do_polar      = p_ptr%l_polar
            ctrl%do_prob_align = p_ptr%l_prob_align_mode
            ctrl%do_bench      = L_BENCH_GLOB
            ctrl%do_sigma_mode = (ctrl%refine_mode == 'sigma')
            ctrl%do_write_oris = .not. ctrl%do_sigma_mode
            select case(ctrl%refine_mode)
                case('eval','sigma')
                    ctrl%do_restore = .false.
                case default
                    if( ctrl%do_polar )then
                        ctrl%do_restore = .true.
                    else
                        ctrl%do_restore = (trim(p_ptr%volrec) == 'yes')
                    endif
            end select
            ctrl%do_polar_prepare = ctrl%do_polar .and. .not. cline%defined('vol1')
        end subroutine init_ctrl

        subroutine ensure_even_odd_partition()
            if( b_ptr%spproj_field%get_nevenodd() == 0 )then
                if( l_distr_worker_glob ) THROW_HARD('no eo partitioning available; refine3D_exec')
                call b_ptr%spproj_field%partition_eo
                call b_ptr%spproj%write_segment_inside(p_ptr%oritype)
            endif
        end subroutine ensure_even_odd_partition

        subroutine sample_particles_for_update()
            if( allocated(pinds) ) deallocate(pinds)
            if( ctrl%do_prob_align )then
                call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top], nptcls2update, pinds)
            else
                if( p_ptr%l_fillin .and. mod(which_iter,5) == 0 )then
                    call sample_ptcls4fillin(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls2update, pinds)
                else
                    call sample_ptcls4update3D(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls2update, pinds)
                endif
            endif
        end subroutine sample_particles_for_update

        subroutine prepare_batches()
            batchsz_max = min(nptcls2update, p_ptr%nthr * BATCHTHRSZ)
            nbatches    = ceiling(real(nptcls2update) / real(batchsz_max))
            batches     = split_nobjs_even(nptcls2update, nbatches)
            batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        end subroutine prepare_batches

        subroutine prepare_refs_sigmas_and_pftc()
            if( ctrl%do_prob_align .and. (.not. ctrl%do_polar_prepare) )then
                nrefs = p_ptr%nspace * p_ptr%nstates
                call b_ptr%pftc%new(p_ptr, nrefs, [1,batchsz_max], p_ptr%kfromto)
            endif
            if( ctrl%do_polar_prepare )then
                if( ctrl%do_bench ) t_prep_pftc_polar_mode = tic()
                call prep_pftc4align3D_polar(p_ptr, b_ptr, cline, batchsz_max)
                if( ctrl%do_bench ) rt_prep_pftc_polar_mode = toc(t_prep_pftc_polar_mode)
            endif
            if( ctrl%do_bench ) t_prep_sigmas_alloc_ptcl_imgs = tic()
            call prep_sigmas_alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
            if( ctrl%do_bench ) rt_prep_sigmas_alloc_ptcl_imgs = toc(t_prep_sigmas_alloc_ptcl_imgs)
            if( (.not. ctrl%do_prob_align) .and. (.not. ctrl%do_polar_prepare) )then
                if( ctrl%do_bench ) t_prep_reproj_refvols = tic()
                call read_mask_filter_reproject_refvols(p_ptr, b_ptr, cline, batchsz_max)
                if( ctrl%do_bench ) rt_prep_reproj_refvols = toc(t_prep_reproj_refvols)
            endif
            call build%vol%kill
            call build%vol_odd%kill
            call build%vol2%kill
            if( ctrl%do_polar .and. ctrl%do_restore )then
                if( cline%defined('vol1') )then
                    call b_ptr%pftc%polar_cavger_new(.true.)
                    if( p_ptr%l_trail_rec )then
                        call b_ptr%pftc%polar_cavger_write_eo_pftcrefs(string(POLAR_REFS_FBODY))
                    endif
                endif
                call b_ptr%pftc%polar_cavger_zero_pft_refs
                if( file_exists(p_ptr%frcs) )then
                    call b_ptr%clsfrcs%read(p_ptr%frcs)
                else
                    call b_ptr%clsfrcs%new(p_ptr%nspace, p_ptr%box_crop, p_ptr%smpd_crop, p_ptr%nstates)
                endif
            endif
        end subroutine prepare_refs_sigmas_and_pftc

        subroutine maybe_init_reconstruction()
            if( .not. ctrl%do_restore ) return
            if( ctrl%do_polar )then
                if( ctrl%polar_mode == 'new' )then
                    call init_rec(params, build, batchsz_max, fpls)
                    call alloc_imgarr(batchsz_max, [p_ptr%box,p_ptr%box,1], p_ptr%smpd, ptcl_rec_imgs)
                endif
            else
                call init_rec(params, build, batchsz_max, fpls)
                call alloc_imgarr(batchsz_max, [p_ptr%box,p_ptr%box,1], p_ptr%smpd, ptcl_rec_imgs)
            endif
        end subroutine maybe_init_reconstruction

        subroutine build_batch_particles_local()
            if( ctrl%do_bench ) t_build_batch_ptcls = tic()
            if( ctrl%do_polar )then
                select case(ctrl%polar_mode)
                case('new')
                    call build_batch_particles3D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                        ptcl_match_imgs, ptcl_match_imgs_pad, imgs4rec=ptcl_rec_imgs(1:batchsz))
                case default
                    call build_batch_particles3D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                        ptcl_match_imgs, ptcl_match_imgs_pad)
                end select
            else
                if( ctrl%do_restore )then
                    call build_batch_particles3D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                        ptcl_match_imgs, ptcl_match_imgs_pad, imgs4rec=ptcl_rec_imgs(1:batchsz))
                else
                    call build_batch_particles3D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                        ptcl_match_imgs, ptcl_match_imgs_pad)
                endif
            endif
            if( ctrl%do_bench ) rt_build_batch_ptcls = rt_build_batch_ptcls + toc(t_build_batch_ptcls)
        end subroutine build_batch_particles_local

        subroutine choose_and_run_strategy(iptcl, iptcl_batch, iptcl_map, ithr, has_been_searched)
            integer, intent(in) :: iptcl, iptcl_batch, iptcl_map, ithr
            logical, intent(in) :: has_been_searched
            select case(ctrl%refine_mode)
            case('shc')
                if( .not. has_been_searched )then
                    allocate(strategy3D_greedy :: strategy3Dsrch(iptcl_batch)%ptr)
                    cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                else
                    if( ran3() < GREEDY_FREQ )then
                        allocate(strategy3D_greedy :: strategy3Dsrch(iptcl_batch)%ptr)
                        cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                    else
                        allocate(strategy3D_shc :: strategy3Dsrch(iptcl_batch)%ptr)
                    endif
                endif
            case('shc_smpl')
                if( b_ptr%spproj_field%is_first_update(which_iter, iptcl) )then
                    allocate(strategy3D_greedy_smpl :: strategy3Dsrch(iptcl_batch)%ptr)
                    cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                else
                    allocate(strategy3D_shc_smpl :: strategy3Dsrch(iptcl_batch)%ptr)
                endif
            case('snhc_smpl')
                if( b_ptr%spproj_field%is_first_update(which_iter, iptcl) )then
                    allocate(strategy3D_greedy_smpl :: strategy3Dsrch(iptcl_batch)%ptr)
                    cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                else
                    allocate(strategy3D_snhc_smpl :: strategy3Dsrch(iptcl_batch)%ptr)
                endif
            case('eval')
                allocate(strategy3D_eval :: strategy3Dsrch(iptcl_batch)%ptr)
            case('neigh')
                allocate(strategy3D_greedy_sub :: strategy3Dsrch(iptcl_batch)%ptr)
            case('greedy')
                allocate(strategy3D_greedy :: strategy3Dsrch(iptcl_batch)%ptr)
            case('prob','prob_state','prob_neigh')
                allocate(strategy3D_prob :: strategy3Dsrch(iptcl_batch)%ptr)
            case('ptree')
                allocate(strategy3D_ptree :: strategy3Dsrch(iptcl_batch)%ptr)
            case('ptree_neigh')
                allocate(strategy3D_ptree_neigh :: strategy3Dsrch(iptcl_batch)%ptr)
            case('ptree_neigh_states')
                allocate(strategy3D_ptree_neigh_states :: strategy3Dsrch(iptcl_batch)%ptr)
            case('shc_ptree')
                allocate(strategy3D_shc_ptree :: strategy3Dsrch(iptcl_batch)%ptr)
            case('sigma')
                call b_ptr%spproj_field%get_ori(iptcl, orientation)
                call b_ptr%spproj_field%set(iptcl, 'proj', b_ptr%eulspace%find_closest_proj(orientation))
            case default
                THROW_HARD('refinement mode: '//trim(ctrl%refine_mode)//' unsupported')
            end select
            if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                call strategy3Dsrch(iptcl_batch)%ptr%new(p_ptr, strategy3Dspecs(iptcl_batch), b_ptr)
                call strategy3Dsrch(iptcl_batch)%ptr%srch(b_ptr%spproj_field, ithr)
                incr_shifts(:,iptcl_batch) = b_ptr%spproj_field%get_2Dshift(iptcl) - &
                    strategy3Dsrch(iptcl_batch)%ptr%s%prev_shvec
                call strategy3Dsrch(iptcl_batch)%ptr%kill
            endif
        end subroutine choose_and_run_strategy

        subroutine maybe_restore_batch()
            if( .not. ctrl%do_restore ) return
            if( ctrl%do_bench ) t_rec = tic()
            if( ctrl%do_polar )then
                select case(ctrl%polar_mode)
                case('new')
                    call prep_imgs4rec(params, b_ptr, batchsz, ptcl_rec_imgs(:batchsz), &
                        pinds(batch_start:batch_end), fpls(:batchsz))
                    call b_ptr%pftc%polar_cavger_insert_ptcls_oversamp(b_ptr%eulspace, b_ptr%spproj_field, &
                        b_ptr%pgrpsyms, batchsz, pinds(batch_start:batch_end), fpls(:batchsz))
                case default
                    call b_ptr%pftc%polar_cavger_update_sums(batchsz, pinds(batch_start:batch_end), &
                        b_ptr%spproj, incr_shifts(:,1:batchsz), is3D=.true.)
                end select
            else
                call prep_imgs4rec(params, b_ptr, batchsz, ptcl_rec_imgs(:batchsz), &
                    pinds(batch_start:batch_end), fpls(:batchsz))
                call update_rec(params, b_ptr, batchsz, pinds(batch_start:batch_end), fpls(:batchsz))
            endif
            if( ctrl%do_bench ) rt_rec = rt_rec + toc(t_rec)
        end subroutine maybe_restore_batch

        subroutine maybe_write_orientations()
            if( .not. ctrl%do_write_oris ) return
            if( ctrl%do_bench ) t_projio = tic()
            if( p_ptr%top < p_ptr%fromp )then
                THROW_HARD('invalid output write range in refine3D_exec: TOP < FROMP')
            endif
            select case(ctrl%oritype)
                case('ptcl3D')
                    call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, b_ptr%spproj_field, &
                        [p_ptr%fromp,p_ptr%top], isegment=PTCL3D_SEG)
                case('cls3D')
                    call binwrite_oritab(p_ptr%outfile, b_ptr%spproj, b_ptr%spproj_field, &
                        [p_ptr%fromp,p_ptr%top], isegment=CLS3D_SEG)
                case default
                    THROW_HARD('unsupported oritype: '//trim(ctrl%oritype)//'; refine3D_exec')
            end select
            p_ptr%oritab = p_ptr%outfile
            if( ctrl%do_bench ) rt_projio = toc(t_projio)
        end subroutine maybe_write_orientations

    end subroutine refine3D_exec

end module simple_strategy3D_matcher