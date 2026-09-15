!@descr: high-level particle matching and partial-reconstruction orchestration for refine3D workers
module simple_strategy3D_matcher
use, intrinsic :: iso_fortran_env, only: int64, real64
use simple_pftc_srch_api
use simple_matcher_refvol_utils
use simple_matcher_ptcl_batch
use simple_strategy3D_alloc,        only: clean_strategy3D, prep_strategy3D, s3D
use simple_binoris_io,              only: binwrite_oritab
use simple_builder,                 only: builder
use simple_euclid_sigma2,           only: euclid_sigma2
use simple_eul_prob_tab,            only: eul_prob_tab
use simple_matcher_2Dprep,          only: prepimg4align
use simple_matcher_3Drec,           only: calc_3Drec, calc_projdir3Drec
use simple_rec3D_pcg_strategy,      only: execute_rec3D_pcg_worker
use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4missing3D, sample_ptcls4update3D
use simple_qsys_funs,               only: qsys_job_finished
use simple_refine3D_fnames,         only: refine3D_bench_fname
use simple_syslib,                  only: get_peak_rss_bytes, get_current_rss_bytes
use simple_strategy3D_eval,         only: strategy3D_eval
use simple_strategy3D_greedy_smpl,  only: strategy3D_greedy_smpl
use simple_strategy3D_greedy_sub,   only: strategy3D_greedy_sub
use simple_strategy3D_greedy,       only: strategy3D_greedy
use simple_strategy3D_greedy_inpl,  only: strategy3D_greedy_inpl
use simple_strategy3D_prob,         only: strategy3D_prob
use simple_strategy3D_shc_smpl,     only: strategy3D_shc_smpl
use simple_strategy3D_shc,          only: strategy3D_shc
use simple_strategy3D_snhc_smpl,    only: strategy3D_snhc_smpl
use simple_strategy3D_srch,         only: strategy3D_spec
use simple_strategy3D,              only: strategy3D
use simple_ori_utils,               only: dm2euler
use simple_pose_cont_refine3D_adapter, only: pose_cont_reference_workspace, &
    &pose_cont_pose, pose_cont_config, pose_cont_limits, pose_cont_transaction_result, &
    &cartesian_pose_data, prepare_pose_cont_observation, shift_native_to_crop, &
    &shift_crop_to_native, LM_ACCEPTED_IMPROVEMENT, &
    &POSE_CONT_ROUTE_SHIFT_THEN_JOINT, POSE_CONT_ROUTE_JOINT
implicit none

public :: refine3D_exec
private
#include "simple_local_flags.inc"

type :: refine3D_ctrl
    character(len=:), allocatable :: refine_mode
    character(len=:), allocatable :: oritype
    logical :: do_write_partial_recs
    logical :: do_prob_align
    logical :: do_projrec
    logical :: do_sigma_mode
    logical :: do_emit_sigma
    logical :: do_write_oris
    logical :: do_bench
    logical :: do_pose_cont
  contains
    procedure :: print_flags
end type refine3D_ctrl

contains

    subroutine refine3D_exec( params, build, cline, which_iter, converged, l_write_partial_recs )
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(inout) :: converged
        logical,                   intent(in) :: l_write_partial_recs
        class(parameters), pointer :: p_ptr => null()
        class(builder),    pointer :: b_ptr => null()
        type(eul_prob_tab), target :: eulprob_obj_part
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        type(strategy3D_spec),     allocatable :: strategy3Dspecs(:)
        type(image),               allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:)
        integer,                   allocatable :: batches(:,:), cnt_greedy(:), cnt_all(:), pinds(:)
        real,                      allocatable :: incr_shifts(:,:)
        type(ori)           :: orientation
        type(refine3D_ctrl) :: ctrl
        type(pose_cont_reference_workspace) :: pose_cont_refs
        type(pose_cont_config) :: pose_config
        type(pose_cont_limits) :: pose_limits
        real                :: frac_greedy
        real(dp)            :: pose_cont_crop_scale
        integer             :: nbatches, batchsz_max, batch_start, batch_end, batchsz
        integer             :: iptcl, fnr, ithr, iptcl_batch, iptcl_map, ibatch, nptcls2update
        logical             :: has_been_searched
        ! benchmarking
        integer(int64) :: peak_rss, rss_after_teardown, rss_after_reconstruction
        real(real64)    :: peak_rss_gib, rss_after_teardown_gib, rss_after_reconstruction_gib
        type(string) :: benchfname
        integer(timer_int_kind) :: t_startup, t_build_batch_ptcls, t_prep_orisrch, t_align, t_rec, t_tot, t_projio
        integer(timer_int_kind) :: t_alloc_ptcl_imgs
        integer(timer_int_kind) :: t_prep_refs, t_memoize_refs
        real(timer_int_kind)    :: rt_startup, rt_build_batch_ptcls, rt_prep_orisrch, rt_align, rt_rec, rt_tot, rt_projio
        real(timer_int_kind)    :: rt_alloc_ptcl_imgs
        real(timer_int_kind)    :: rt_prep_refs, rt_memoize_refs, rt_rec_accum, rt_rec_write
        p_ptr => params
        b_ptr => build
        rss_after_teardown      = -1_int64
        rss_after_reconstruction = -1_int64
        call init_ctrl()
        converged = .false.
        if( ctrl%do_bench )then
            t_startup = tic()
            t_tot     = t_startup
        endif
        call ensure_even_odd_partition()
        has_been_searched = .not. b_ptr%spproj%is_virgin_field(p_ptr%oritype)
        call adopt_reprojection_model_range(p_ptr, b_ptr)
        call sample_particles_for_update( pinds, nptcls2update )
        if( nptcls2update < 1 )then
            if( p_ptr%l_update_missing )then
                write(logfhandle,'(A)') '>>> MATCH3D: no missing particles selected for update'
                ! Canonical consolidation still expects one range from every
                ! scheduled partition. Emit the unchanged committed slice so
                ! an empty update is a valid transaction rather than a
                ! missing-file failure.
                if( ctrl%do_emit_sigma )then
                    call prep_sigmas_objfun(p_ptr, b_ptr)
                    call b_ptr%esig%write_sigma2
                endif
                converged = .true.
                call qsys_job_finished(p_ptr, string('simple_strategy3D_matcher :: refine3D_exec'))
                return
            endif
            THROW_HARD('No particles selected for 3D update')
        endif
        call prepare_particles_batches( nptcls2update )
        if( ctrl%do_bench )then
            rt_startup = toc(t_startup)
            rt_build_batch_ptcls = 0.0
            rt_alloc_ptcl_imgs   = 0.0
            rt_prep_refs = 0.0
            rt_memoize_refs      = 0.0
            rt_prep_orisrch      = 0.0
            rt_align             = 0.0
            rt_projio            = 0.0
            rt_rec               = 0.0
            rt_rec_accum         = 0.0
            rt_rec_write         = 0.0
        endif
        call prepare_refs_sigmas_and_pftc()
        if( ctrl%do_bench ) t_memoize_refs = tic()
        if( .not. ctrl%do_prob_align ) call build%pftc%memoize_refs(eulspace=build%eulspace)
        if( ctrl%do_bench )then
            rt_memoize_refs = toc(t_memoize_refs)
            t_prep_orisrch  = tic()
        endif
        call prep_strategy3D(p_ptr, b_ptr)
        allocate(strategy3Dspecs(batchsz_max), strategy3Dsrch(batchsz_max))
        if( ctrl%do_prob_align )then
            call eulprob_obj_part%new_assignment(p_ptr, b_ptr, pinds)
            call eulprob_obj_part%read_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        endif
        if( ctrl%do_bench )then
            rt_prep_orisrch     = toc(t_prep_orisrch)
            rt_build_batch_ptcls= 0.0
            rt_align            = 0.0
        endif
        allocate(cnt_greedy(p_ptr%nthr), cnt_all(p_ptr%nthr), source=0)
        if( trim(p_ptr%inpl_cont) == 'yes' )then
            call b_ptr%spproj_field%set_all2single('cont_inpl_attempted', 0.)
            call b_ptr%spproj_field%set_all2single('cont_inpl_improved',  0.)
        endif
        allocate(incr_shifts(2,batchsz_max), source=0.0)
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            incr_shifts(:,1:batchsz) = 0.0
            call build_batch_particles_local()
            if( ctrl%do_bench ) t_align = tic()
            !$omp parallel do default(shared) private(iptcl,iptcl_batch,iptcl_map,ithr,orientation) &
            !$omp schedule(static) proc_bind(close)
            do iptcl_batch = 1, batchsz
                iptcl_map     = batch_start + iptcl_batch - 1
                iptcl         = pinds(iptcl_map)
                ithr          = omp_get_thread_num() + 1
                cnt_all(ithr) = cnt_all(ithr) + 1
                strategy3Dspecs(iptcl_batch)%iptcl     = iptcl
                strategy3Dspecs(iptcl_batch)%iptcl_map = iptcl_map
                if( ctrl%do_prob_align ) strategy3Dspecs(iptcl_batch)%eulprob_obj_part => eulprob_obj_part
                call choose_and_run_strategy(iptcl, iptcl_batch, ithr, has_been_searched)
                if( ctrl%do_emit_sigma )then
                    call b_ptr%spproj_field%get_ori(iptcl, orientation)
                    call orientation%set_shift(incr_shifts(:,iptcl_batch))
                    call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'proj')
                endif
            enddo
            !$omp end parallel do
            if( ctrl%do_bench ) rt_align = rt_align + toc(t_align)
        enddo
        frac_greedy = 0.0
        if( any(cnt_greedy > 0) .and. any(cnt_all > 0) )then
            frac_greedy = real(sum(cnt_greedy)) / real(sum(cnt_all))
        endif
        call b_ptr%spproj_field%set_all2single('frac_greedy', frac_greedy)
        if( ctrl%do_emit_sigma ) call b_ptr%esig%write_sigma2
        if( ctrl%do_projrec ) call b_ptr%spproj_field%set_projs(b_ptr%eulspace)
        call maybe_write_orientations()
        do iptcl_batch = 1, batchsz_max
            nullify(strategy3Dsrch(iptcl_batch)%ptr)
        enddo
        deallocate(strategy3Dsrch, strategy3Dspecs, batches)
        deallocate(cnt_greedy, cnt_all, incr_shifts)
        call eulprob_obj_part%kill
        call clean_strategy3D
        call b_ptr%kill_strategy3D_tbox
        call b_ptr%vol%kill
        call orientation%kill
        call clean_batch_particles3D(b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad)
        ! Registration is complete.  Release the all-state reprojection model,
        ! particle PFTs, memoized correlations, and PFTC thread workspaces
        ! before constructing the first state reconstruction.
        call b_ptr%pftc%kill
        if( b_ptr%pftc%exists() ) THROW_HARD('PFTC still allocated at reconstruction phase boundary')
        if( ctrl%do_bench ) rss_after_teardown = get_current_rss_bytes()
        if( ctrl%do_write_partial_recs )then
            if( ctrl%do_bench ) t_rec = tic()
            if( trim(params%rec_backend) == 'pcg' )then
                call execute_rec3D_pcg_worker(params, build, cline, pinds)
            else if( ctrl%do_projrec )then
                call calc_projdir3Drec(params, build, cline, nptcls2update, pinds)
            else
                call calc_3Drec(params, build, cline, nptcls2update, pinds)
            endif
            if( ctrl%do_bench ) rt_rec_write = rt_rec_write + toc(t_rec)
        endif
        call b_ptr%esig%kill
        call pose_cont_refs%kill
        if( ctrl%do_bench ) rss_after_reconstruction = get_current_rss_bytes()
        call qsys_job_finished(p_ptr, string('simple_strategy3D_matcher :: refine3D_exec'))
        if( ctrl%do_bench )then
            rt_rec = rt_rec_accum + rt_rec_write
            rt_tot = toc(t_tot)
            peak_rss = get_peak_rss_bytes()
            peak_rss_gib = -1.0_real64
            if( peak_rss >= 0_int64 ) peak_rss_gib = real(peak_rss,real64) / real(1024_int64**3,real64)
            rss_after_teardown_gib = -1.0_real64
            if( rss_after_teardown >= 0_int64 )then
                rss_after_teardown_gib = real(rss_after_teardown,real64) / real(1024_int64**3,real64)
            endif
            rss_after_reconstruction_gib = -1.0_real64
            if( rss_after_reconstruction >= 0_int64 )then
                rss_after_reconstruction_gib = real(rss_after_reconstruction,real64) / real(1024_int64**3,real64)
            endif
            ! every partition writes its own collision-free record so worker time,
            ! load imbalance and memory can be aggregated; partition 1 also keeps
            ! the legacy per-iteration file the existing parsers read
            call write_bench_file(refine3D_bench_fname(which_iter, p_ptr%part, p_ptr%numlen))
            if( p_ptr%part == 1 ) call write_bench_file(refine3D_bench_fname(which_iter))
        endif

    contains

        subroutine write_bench_file( fname )
            type(string), intent(in) :: fname
            benchfname = fname
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'match3D refine mode                 : ', trim(ctrl%refine_mode)
            write(fnr,'(a,l1)') 'match3D write partial outputs       : ', ctrl%do_write_partial_recs
            write(fnr,'(a,i0)') 'match3D nspace                      : ', p_ptr%nspace
            write(fnr,'(a,i0)') 'match3D nstates                     : ', p_ptr%nstates
            write(fnr,'(a,i0)') 'match3D kfrom                       : ', p_ptr%kfromto(1)
            write(fnr,'(a,i0)') 'match3D kto                         : ', p_ptr%kfromto(2)
            write(fnr,'(a,i0)') 'match3D process partition           : ', p_ptr%part
            write(fnr,'(a,i0)') 'match3D process pid                 : ', p_ptr%pid
            write(fnr,'(a,i0)') 'match3D nparts                      : ', p_ptr%nparts
            write(fnr,'(a,i0)') 'match3D worker threads              : ', p_ptr%nthr
            write(fnr,'(a,i0)') 'match3D box                         : ', p_ptr%box
            write(fnr,'(a,i0)') 'match3D box_crop                    : ', p_ptr%box_crop
            write(fnr,'(a,a)')  'match3D rec_backend                 : ', trim(p_ptr%rec_backend)
            write(fnr,'(a,i0)') 'match3D maxits_pcg                  : ', p_ptr%maxits_pcg
            write(fnr,'(a,es12.4)') 'match3D rtol                        : ', p_ptr%rtol
            write(fnr,'(a,i0)') 'match3D peak RSS (bytes)            : ', peak_rss
            write(fnr,'(a,f0.3)') 'match3D peak RSS (GiB)              : ', peak_rss_gib
            write(fnr,'(a,i0)') 'match3D RSS after align teardown (bytes): ', rss_after_teardown
            write(fnr,'(a,f0.3)') 'match3D RSS after align teardown (GiB)  : ', rss_after_teardown_gib
            write(fnr,'(a,i0)') 'match3D RSS after reconstruction (bytes): ', rss_after_reconstruction
            write(fnr,'(a,f0.3)') 'match3D RSS after reconstruction (GiB)  : ', rss_after_reconstruction_gib
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'match3D startup/setup              :', rt_startup
            write(fnr,'(a,1x,f0.2)') 'match3D particle preparation       :', rt_build_batch_ptcls + rt_alloc_ptcl_imgs
            write(fnr,'(a,1x,f0.2)') 'match3D reference preparation      :', rt_prep_refs + rt_memoize_refs
            write(fnr,'(a,1x,f0.2)') 'match3D orientation search         :', rt_prep_orisrch + rt_align
            write(fnr,'(a,1x,f0.2)') 'match3D project metadata I/O       :', rt_projio
            write(fnr,'(a,1x,f0.2)') 'match3D partial reconstruction     :', rt_rec
            write(fnr,'(a,1x,f0.2)') 'match3D partial reconstruction thread-s:', rt_rec * real(p_ptr%nthr, kind(rt_rec))
            write(fnr,'(a,1x,f0.2)') 'match3D total time                 :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'match3D total thread-s             :', rt_tot * real(p_ptr%nthr, kind(rt_tot))
            write(fnr,'(a,1x,f0.2)') 'match3D % accounted for            :', &
                &((rt_startup + rt_build_batch_ptcls + rt_alloc_ptcl_imgs + rt_prep_refs + &
                &  rt_memoize_refs + rt_prep_orisrch + rt_align + rt_projio + rt_rec) / rt_tot) * 100.
            call fclose(fnr)
        end subroutine write_bench_file

        subroutine init_ctrl()
            ctrl%refine_mode   = trim(p_ptr%refine)
            ctrl%oritype       = trim(p_ptr%oritype)
            ctrl%do_prob_align = p_ptr%l_prob_align_mode
            ctrl%do_projrec    = trim(p_ptr%projrec) == 'yes'
            ctrl%do_bench      = L_BENCH_GLOB
            ctrl%do_pose_cont  = trim(p_ptr%pose_cont) == 'yes'
            ctrl%do_sigma_mode = (ctrl%refine_mode == 'sigma')
            ctrl%do_emit_sigma = p_ptr%cc_objfun == OBJFUN_EUCLID .or. trim(p_ptr%cc_emit_sigma) == 'yes'
            ctrl%do_write_oris = .not. ctrl%do_sigma_mode
            if( ctrl%do_pose_cont )then
                if( p_ptr%cc_objfun /= OBJFUN_EUCLID ) &
                    &THROW_HARD('pose_cont requires objfun=euclid')
                if( ctrl%do_prob_align .or. ctrl%do_sigma_mode .or. ctrl%refine_mode == 'eval' ) &
                    &THROW_HARD('pose_cont requires an ordinary pose-search refinement mode')
                select case(trim(p_ptr%pose_cont_route))
                    case('shift_then_joint')
                        pose_config%route = POSE_CONT_ROUTE_SHIFT_THEN_JOINT
                    case('joint')
                        pose_config%route = POSE_CONT_ROUTE_JOINT
                    case DEFAULT
                        THROW_HARD('unsupported pose_cont_route')
                end select
                ! The adapter works in cropped-box pixels; express the one-native-
                ! pixel proposal and five-native-pixel capture bounds on that grid.
                pose_cont_crop_scale = real(p_ptr%box_crop,dp)/real(p_ptr%box,dp)
                pose_limits = pose_cont_limits(shift_step_bound=pose_cont_crop_scale, &
                    &max_total_shift=5._dp*pose_cont_crop_scale)
            endif
            select case(ctrl%refine_mode)
                case('eval','sigma')
                    ctrl%do_write_partial_recs = .false.
                case default
                    ctrl%do_write_partial_recs = l_write_partial_recs
            end select
        end subroutine init_ctrl

        subroutine ensure_even_odd_partition()
            if( b_ptr%spproj_field%get_nevenodd() == 0 )then
                if( l_distr_worker_glob ) THROW_HARD('no eo partitioning available; refine3D_exec')
                call b_ptr%spproj_field%partition_eo
                call b_ptr%spproj%write_segment_inside(p_ptr%oritype)
            endif
        end subroutine ensure_even_odd_partition

        subroutine sample_particles_for_update( pinds_local, nptcls )
            integer, allocatable, intent(out) :: pinds_local(:)
            integer,              intent(out) :: nptcls
            if( allocated(pinds_local) ) deallocate(pinds_local)
            if( ctrl%do_prob_align )then
                if( p_ptr%l_update_missing )then
                    THROW_HARD('update_missing requires matcher-owned assignment; use a non-probabilistic refine mode')
                endif
                call b_ptr%spproj_field%sample4update_reprod([p_ptr%fromp,p_ptr%top], nptcls, pinds_local)
            else
                if( p_ptr%l_update_missing )then
                    call sample_ptcls4missing3D(b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls, pinds_local)
                else if( p_ptr%l_fillin .and. mod(which_iter,5) == 0 )then
                    call sample_ptcls4fillin(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls, pinds_local)
                else
                    call sample_ptcls4update3D(p_ptr, b_ptr, [p_ptr%fromp,p_ptr%top], .true., nptcls, pinds_local)
                endif
            endif
        end subroutine sample_particles_for_update

        subroutine prepare_particles_batches( nptcls )
            integer, intent(in) :: nptcls
            batchsz_max = min(nptcls, p_ptr%nthr * BATCHTHRSZ)
            nbatches    = ceiling(real(nptcls) / real(batchsz_max))
            batches     = split_nobjs_even(nptcls, nbatches)
            batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        end subroutine prepare_particles_batches

        subroutine prepare_refs_sigmas_and_pftc()
            if( ctrl%do_bench ) t_prep_refs = tic()
            call read_reprojection_model(p_ptr, b_ptr, batchsz_max)
            ! Real-space artifacts bridge the reference-materialization process
            ! and every matcher worker; each worker loads them only once here.
            if( ctrl%do_pose_cont ) &
                &call pose_cont_refs%new_from_artifacts(p_ptr%nstates,p_ptr%box_crop,p_ptr%smpd_crop)
            call prep_sigmas_objfun(p_ptr, b_ptr)
            if( ctrl%do_bench ) rt_prep_refs = toc(t_prep_refs)
            if( ctrl%do_bench ) t_alloc_ptcl_imgs = tic()
            call alloc_ptcl_imgs(p_ptr, b_ptr, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
            if( ctrl%do_bench ) rt_alloc_ptcl_imgs = toc(t_alloc_ptcl_imgs)
            call build%vol%kill
            call build%vol_odd%kill
            call build%vol2%kill
        end subroutine prepare_refs_sigmas_and_pftc

        subroutine build_batch_particles_local()
            if( ctrl%do_bench ) t_build_batch_ptcls = tic()
            call build_batch_particles3D(p_ptr, b_ptr, batchsz, pinds(batch_start:batch_end), &
                ptcl_match_imgs, ptcl_match_imgs_pad)
            if( ctrl%do_bench ) rt_build_batch_ptcls = rt_build_batch_ptcls + toc(t_build_batch_ptcls)
        end subroutine build_batch_particles_local

        subroutine choose_and_run_strategy(iptcl, iptcl_batch, ithr, has_been_searched)
            integer, intent(in) :: iptcl, iptcl_batch, ithr
            logical, intent(in) :: has_been_searched
            type(ori) :: o_sigma ! procedure-local: thread-safe, unlike the host's orientation
            logical :: attempted, improved, no_improvement, invalid
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
                        allocate(strategy3D_greedy_smpl    :: strategy3Dsrch(iptcl_batch)%ptr)
                        cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                    else
                        allocate(strategy3D_shc_smpl       :: strategy3Dsrch(iptcl_batch)%ptr)
                    endif
                case('snhc_smpl')
                    if( b_ptr%spproj_field%is_first_update(which_iter, iptcl) )then
                        allocate(strategy3D_greedy_smpl    :: strategy3Dsrch(iptcl_batch)%ptr)
                        cnt_greedy(ithr) = cnt_greedy(ithr) + 1
                    else
                        allocate(strategy3D_snhc_smpl      :: strategy3Dsrch(iptcl_batch)%ptr)
                    endif
                case('eval')
                    allocate(strategy3D_eval               :: strategy3Dsrch(iptcl_batch)%ptr)
                case('neigh')
                    allocate(strategy3D_greedy_sub         :: strategy3Dsrch(iptcl_batch)%ptr)
                case('greedy')
                    allocate(strategy3D_greedy             :: strategy3Dsrch(iptcl_batch)%ptr)
                case('greedy_inpl')
                    allocate(strategy3D_greedy_inpl        :: strategy3Dsrch(iptcl_batch)%ptr)
                case('prob','prob_state','prob_neigh')
                    allocate(strategy3D_prob               :: strategy3Dsrch(iptcl_batch)%ptr)
                case('sigma')
                    ! residual-only pass: no search, the particle's projection
                    ! direction is the closest one to its stored orientation.
                    ! A contained procedure reaches the HOST's orientation, not
                    ! the caller's OpenMP-private copy, so a procedure-local
                    ! ori is mandatory here (shared-ori double free, 2026-09-07)
                    call b_ptr%spproj_field%get_ori(iptcl, o_sigma)
                    call b_ptr%spproj_field%set(iptcl, 'proj', b_ptr%eulspace%find_closest_proj(o_sigma))
                    call o_sigma%kill
                case default
                    THROW_HARD('refinement mode: '//trim(ctrl%refine_mode)//' unsupported')
            end select
            if( associated(strategy3Dsrch(iptcl_batch)%ptr) )then
                call strategy3Dsrch(iptcl_batch)%ptr%new(p_ptr, strategy3Dspecs(iptcl_batch), b_ptr)
                call strategy3Dsrch(iptcl_batch)%ptr%srch(b_ptr%spproj_field, ithr)
                if( ctrl%do_pose_cont ) call run_pose_cont_after_pftc(iptcl,iptcl_batch,ithr)
                if( trim(p_ptr%inpl_cont) == 'yes' )then
                    call strategy3Dsrch(iptcl_batch)%ptr%s%get_continuous_route_status( &
                        &attempted, improved, no_improvement, invalid)
                    call b_ptr%spproj_field%set(iptcl, 'cont_inpl_attempted', merge(1., 0., attempted))
                    call b_ptr%spproj_field%set(iptcl, 'cont_inpl_improved',  merge(1., 0., improved))
                endif
                incr_shifts(:,iptcl_batch) = b_ptr%spproj_field%get_2Dshift(iptcl) - &
                    strategy3Dsrch(iptcl_batch)%ptr%s%prev_shvec
                call strategy3Dsrch(iptcl_batch)%ptr%kill
                deallocate(strategy3Dsrch(iptcl_batch)%ptr)
                nullify(strategy3Dsrch(iptcl_batch)%ptr)
            endif
        end subroutine choose_and_run_strategy

        !> Run transactional Cartesian LM from the established matcher winner.
        !! With inpl_cont=yes this seed includes its accepted polish; otherwise
        !! it is the ordinary PFTC result. Only an accepted pose_cont result is
        !! committed, so rejected or invalid transactions preserve that seed.
        subroutine run_pose_cont_after_pftc(iptcl, iptcl_batch, ithr)
            integer, intent(in) :: iptcl, iptcl_batch, ithr
            type(ori) :: winner
            type(ctfparams) :: ctfparms, cropped_ctfparms
            type(cartesian_pose_data) :: data
            type(pose_cont_pose) :: seed
            type(pose_cont_transaction_result) :: result
            complex, allocatable :: observed(:,:)
            real, allocatable :: sigma2(:)
            real :: euler(3), shift_native(2)
            integer :: state, eo
            logical :: even

            ! Stage 1: capture the authoritative PFTC/inpl_cont winner.
            call b_ptr%spproj_field%get_ori(iptcl,winner)
            state = winner%get_state()
            eo = winner%get_eo()
            select case(eo)
                case(0)
                    even = .true.
                case(1)
                    even = .false.
                case default
                    THROW_HARD('pose_cont requires an even/odd half-set assignment')
            end select

            ! Stage 2: prepare the cropped Cartesian observation and its
            ! per-shell noise weights for the local objective.
            ctfparms = b_ptr%spproj%get_ctfparams(p_ptr%oritype,iptcl)
            call prepare_pose_cont_observation(b_ptr%imgbatch(iptcl_batch),b_ptr%lmsk, &
                &ptcl_match_imgs(ithr),p_ptr%msk_crop,p_ptr%smpd_crop,ctfparms, &
                &observed,cropped_ctfparms)
            if( .not. allocated(b_ptr%esig%sigma2_noise) ) &
                &THROW_HARD('pose_cont requires allocated sigma2 noise')
            if( p_ptr%kfromto(1) < lbound(b_ptr%esig%sigma2_noise,1) .or. &
                &p_ptr%kfromto(2) > ubound(b_ptr%esig%sigma2_noise,1) ) &
                &THROW_HARD('pose_cont shell range exceeds sigma2 noise bounds')
            if( iptcl < lbound(b_ptr%esig%sigma2_noise,2) .or. &
                &iptcl > ubound(b_ptr%esig%sigma2_noise,2) ) &
                &THROW_HARD('pose_cont particle index exceeds sigma2 noise bounds')
            allocate(sigma2(0:p_ptr%kfromto(2)),source=1.)
            sigma2(p_ptr%kfromto(1):p_ptr%kfromto(2)) = &
                &b_ptr%esig%sigma2_noise(p_ptr%kfromto(1):p_ptr%kfromto(2),iptcl)
            call pose_cont_refs%prepare_particle(state,even,observed,cropped_ctfparms, &
                &sigma2,p_ptr%kfromto,data)
            seed%rotmat = real(winner%get_mat(),dp)
            seed%shift = real(shift_native_to_crop(winner%get_2Dshift(), &
                &p_ptr%box,p_ptr%box_crop),dp)

            ! Stage 3: run the selected local LM route transactionally.
            call pose_cont_refs%refine_particle(state,even,seed,data,pose_config,pose_limits,result)

            ! Stage 4: commit only an accepted Cartesian improvement. Rejected
            ! or invalid transactions leave the established pose unchanged.
            if( result%status == LM_ACCEPTED_IMPROVEMENT )then
                ! Convert the accepted Cartesian pose back to SIMPLE's native
                ! project coordinates before sigma evaluation/reconstruction.
                euler = real(dm2euler(result%pose%rotmat))
                shift_native = shift_crop_to_native(real(result%pose%shift), &
                    &p_ptr%box,p_ptr%box_crop)
                call winner%set_euler(euler)
                call winner%set_shift(shift_native)

                ! Keep the nearest discrete companions consistent for legacy
                ! PFTC/sigma consumers. All other fields, including corr,
                ! state, and half-set identity, remain those of the winner.
                call winner%set('proj',real(b_ptr%eulspace%find_closest_proj(winner)))
                call winner%set('inpl',real(b_ptr%pftc%get_roind(360.-winner%e3get())))
                call b_ptr%spproj_field%set_ori(iptcl,winner)
            endif
            call winner%kill
        end subroutine run_pose_cont_after_pftc

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

    ! debugging convenience function
    subroutine print_flags( ctrl )
        class(refine3D_ctrl), intent(in) :: ctrl
        write(logfhandle,*) 'refine_mode           : ', ctrl%refine_mode
        write(logfhandle,*) 'oritype               : ', ctrl%oritype
        write(logfhandle,*) 'do_write_partial_recs : ', ctrl%do_write_partial_recs
        write(logfhandle,*) 'do_prob_align         : ', ctrl%do_prob_align
        write(logfhandle,*) 'do_projrec            : ', ctrl%do_projrec
        write(logfhandle,*) 'do_sigma_mode         : ', ctrl%do_sigma_mode
        write(logfhandle,*) 'do_write_oris         : ', ctrl%do_write_oris
        write(logfhandle,*) 'do_bench              : ', ctrl%do_bench
        write(logfhandle,*) 'do_pose_cont          : ', ctrl%do_pose_cont
    end subroutine print_flags

end module simple_strategy3D_matcher
