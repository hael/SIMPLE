!@descr: the pose student's strategy routines: the teacher's rows of the selected particles at the student grid and band, their resolvability and their features, the student's training (learn_pose_student), and its prediction of the mode, the cone and the confidence (predict_poses)
! The teacher is the probabilistic machinery's posterior (doc/implementation_notes/completed/learned_pose_state_assignment.md,
! section 4): the dense table's rows as the production fill computes them against the references prepared as
! refine3D's matcher prepares them, batch by batch, each row read before the batch's assignment releases it.
module simple_strategy_pose_student
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int32
use simple_core_module_api
use simple_builder,            only: builder
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_image,              only: image
use simple_ori,                only: ori
use simple_eo_partition,       only: ensure_field_halves
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_utils, only: prob_row_record, prob_row_file, prob_fill_profile, write_fill_profile_bench, &
    &temperature_or_argmin
use simple_matcher_refvol_utils,    only: materialize_reprojection_model_from_volumes, read_reprojection_model, &
    &complete_volume_source_available
use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D, &
    &check_pass_team
use simple_matcher_ptcl_io,         only: ptcl_read_session
use simple_pftc_shsrch_grad,        only: pftc_shsrch_grad
use simple_matcher_smpl_and_lplims, only: set_bp_range3D, adopt_run_seed
use simple_sigma2_bootstrap,        only: ensure_sigma2_for_iteration
use simple_qsys_funs,               only: qsys_declare_part_finished
use simple_polar_feats,             only: polar_feats, polar_feat, polar_feats_cache, polar_feats_contract, &
    &polar_feats_hash_ints, POLAR_FEATS_NHEAD
use simple_pose_student,  only: pose_student, student_pred
use simple_polar_student, only: student_train_profile
use simple_strategy_student_common, only: student_holdout_split, student_conf_threshold, student_geometry_mismatch, &
    &student_part_pinds, file_hash, choose_student_design, fname_base, student_predict_parts, student_sigma2_for_pass
use simple_defs_hhcnn,    only: NSPACE_STUDENT, KMAX_STUDENT, STUDENT_CONF_FRACTION, STUDENT_DESIGN_QUANTILE, STUDENT_KMAX_ROW, &
    &STUDENT_WIDTHS, STUDENT_WIDEST_MIN_ROWS, STUDENT_MAX_EPOCHS, STUDENT_PATIENCE_EPOCHS, STUDENT_SHIFT_AUG_NSIG, &
    &STUDENT_SHIFT_REG_HALFWIDTH, STUDENT_SHIFT_REG_MAXITS, POSE_SHIFT_CHECK_NPTCLS, POSE_SHIFT_CHECK_MAGS, &
    &POSE_SHIFT_CHECK_NDIR, POSE_SHIFT_CHECK_NGROUPS
use simple_rnd,          only: rnd_run_seed
use simple_math,         only: hpsort
use simple_type_defs,    only: OBJFUN_EUCLID, PTCL3D_SEG
use simple_binoris_io,   only: binwrite_oritab
implicit none

public :: run_learn_pose_student, learn_pose_student_exec, run_predict_poses, predict_poses_exec
private
#include "simple_local_flags.inc"



contains

    ! ===== LEARN_POSE_STUDENT =====

    !> learn_pose_student on its own parameters and builder, built from cline
    subroutine run_learn_pose_student( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target :: params
        type(builder)            :: build
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        call learn_pose_student_exec(params, build, cline)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_LEARN_POSE_STUDENT NORMAL STOP ****', print_simple=.false.)
    end subroutine run_learn_pose_student

    !> The teacher's rows of the project's selected particles (state 1 of ptcl3D) at the student grid (nspace) and
    !! band (lp) against the references vol1, prepared as refine3D's matcher prepares them, filled batch by batch
    !! with the production fill; each row's record with its resolvability, its whole row, its error against the
    !! particle's stored direction and its feature tensor are written beside fname. Then the student is trained on
    !! the rows and the features (every tenth particle held out for the early stopping) and written as the model
    !! fname with its weights beside it. Project: reads the particle field; writes what refine3D's setup writes
    !! (the halves, the sigma2 bootstrap, the run seed), nothing else
    subroutine learn_pose_student_exec( params, build, cline )
        type(parameters), target, intent(inout) :: params
        type(builder),    target, intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        class(eul_prob_tab), allocatable :: tab
        type(ptcl_read_session)  :: reads
        type(image), allocatable :: ptcl_imgs(:), ptcl_imgs_pad(:)
        type(polar_feats)        :: feats
        type(polar_feat)         :: feat
        type(polar_feats_cache)  :: cache
        type(prob_row_file)      :: rows
        type(prob_row_record)    :: rec
        type(prob_fill_profile)  :: profile
        type(pose_student)       :: student
        type(ori)                :: o
        type(string)             :: base, fname_rows, fname_feats, fname_stats, fname_header, fname_grid, projpath
        real,    allocatable     :: dists(:), wk(:,:)
        real,    allocatable     :: x_all(:,:,:,:), ctf2_all(:,:), side_all(:), dists_all(:,:), temps_all(:), xrefs(:,:,:,:), ctf_one(:,:)
        complex, allocatable     :: pft_ref(:,:)
        real,    allocatable     :: prob(:,:), conf_held(:), radius_held(:)
        logical, allocatable     :: held(:), recov(:,:)
        integer, allocatable     :: pinds(:), batches(:,:), statbuf(:), mode_all(:), held_inds(:), group_held(:)
        integer, allocatable     :: k99_held(:), mode_held_s(:), group_of(:), stkinds(:)
        type(student_pred)       :: pred
        integer(int64)           :: key(POLAR_FEATS_NHEAD), mtime, vol_hash
        integer(timer_int_kind)  :: t0, t_all
        real(timer_int_kind)     :: rt_prep, rt_batch, rt_fill, rt_rows, rt_all, rt_train, rt_cand(size(STUDENT_WIDTHS))
        real(timer_int_kind)     :: rt_heldpred(size(STUDENT_WIDTHS)), rt_boot(size(STUDENT_WIDTHS)), rt_write, rt_kept
        type(student_train_profile) :: prof_cand(size(STUDENT_WIDTHS))
        integer(timer_int_kind)  :: tb
        integer                  :: ntried, ntry
        real     :: euls_ref(3), euls_mode(3), euls_sym(3), euldist, inpldist, temp, dmin, w_ref, invalid, shift_aug
        real     :: recov_cand(size(STUDENT_WIDTHS)), loss_cand(size(STUDENT_WIDTHS)), gain_q(size(STUDENT_WIDTHS))
        real     :: conf_threshold, heldout_seeded, heldout_recovered, cone_bound
        real(dp) :: z
        integer  :: nptcls, batchsz_max, nbatches, ibatch, b0, b1, batchsz, i, iptcl, iproj_ref, iproj_mode, rank_ref
        integer  :: stat_unit, grid_unit, head_unit, io_stat, status, stat_status, j, nvalid_ref, nk, nheld
        integer  :: nptcls_train, nheldout, epochs_run, best_epoch, nweights, iw, iw_best, nthr_here
        integer  :: ithr, epoch_cand(size(STUDENT_WIDTHS)), held_unit
        real     :: loss_train, loss_held, mode_held, time_train_s
        logical  :: l_rebuilt, l_cache, l_by_stack
        character(len=:), allocatable :: msg
        type(pose_student)       :: reread
        integer                  :: nsel
        t_all = tic()
        call check_pass_team('learn_pose_student')
        ! the contract
        if( params%cc_objfun /= OBJFUN_EUCLID ) THROW_HARD('the teacher is the Euclidean posterior: objfun=euclid; learn_pose_student')
        if( .not. params%l_lpset ) THROW_HARD('the band is lp; learn_pose_student')
        if( params%nstates /= 1 ) THROW_HARD('the pose student is single-state; learn_pose_student')
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('the rows are ptcl3D''s; learn_pose_student')
        if( .not. complete_volume_source_available(params) ) THROW_HARD('the references are vol1; learn_pose_student')
        base = fname_base(params%fname)
        fname_rows   = base%to_char()//'_rows.bin'
        fname_feats  = base%to_char()//'_feats.bin'
        fname_stats  = base%to_char()//'_rows.txt'
        fname_header = base%to_char()//'_rows_header.txt'
        fname_grid   = base%to_char()//'_grid.txt'
        ! the selected particles
        call build%spproj_field%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no selected particle; learn_pose_student')
        nptcls = size(pinds)
        if( nptcls < 1 ) THROW_HARD('no selected particle; learn_pose_student')
        write(logfhandle,'(A,I0,A,I0,A,F7.2,A)') '>>> POSE STUDENT ROWS: ', nptcls, ' particles at nspace ', params%nspace, &
            &' and lp ', params%lp, ' A'
        ! the held-out micrographs: every tenth of the selected particles' stacks, every tenth selected particle when they
        ! lie on fewer stacks; a selection that leaves a cohort empty is refused before the rows are made
        allocate(held(nptcls), source=.false.)
        allocate(group_of(nptcls), stkinds(nptcls), source=0)
        do i = 1,nptcls
            stkinds(i) = build%spproj_field%get_int(pinds(i), 'stkind')
        enddo
        call student_holdout_split(nptcls, stkinds, held, group_of, l_by_stack, msg)
        if( len(msg) > 0 ) THROW_HARD(msg)
        deallocate(stkinds)
        ! refine3D's setup: the halves (their one maker), the sigma2 state, the run seed of the keyed draws
        t0 = tic()
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, 1, params%box, params%smpd, &
            &'POSE STUDENT SIGMA2', l_rebuilt)
        if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
        call adopt_run_seed(params, build)
        ! the band and the references at the student grid: the matcher's reprojection model from the volume
        call set_bp_range3D(params, build, cline)
        call materialize_reprojection_model_from_volumes(params, build, cline, cleanup_pftc=.true.)
        ! the batches and the calculator, as refine3D_exec prepares them
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = 1 + (nptcls - 1) / batchsz_max
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call read_reprojection_model(params, build, batchsz_max)
        call prep_sigmas_objfun(params, build)
        call alloc_ptcl_imgs(params, build, ptcl_imgs, ptcl_imgs_pad, batchsz_max)
        call build%vol%kill
        call build%vol_odd%kill
        call build%vol2%kill
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        ! the features' geometry and the frozen weights: the calculator's band and the spectra the rows are whitened with
        call build%pftc%get_sigma2_weights(wk)
        call feats%new(build%pftc%get_pftsz(), params%kfromto, params%box, wk)
        nk = params%kfromto(2) - params%kfromto(1) + 1
        ! the grid's reprojections through the features (no CTF): the readout's templates start from them (ruling R3)
        allocate(xrefs(feats%get_pftsz(), nk, 2, params%nspace), source=0.)
        allocate(pft_ref(feats%get_pftsz_in(), params%kfromto(1):params%kfromto(2)))
        allocate(ctf_one(feats%get_pftsz_in(), params%kfromto(1):params%kfromto(2)), source=1.)
        do i = 1,params%nspace
            call build%pftc%get_ref_pft(i, .true., pft_ref)
            call feats%build_from_pft(pft_ref, ctf_one, .true., feat)
            xrefs(:,:,:,i) = feat%x
        enddo
        deallocate(pft_ref, ctf_one)
        projpath = simple_abspath(params%projfile)
        call simple_file_stat(projpath, stat_status, statbuf)
        if( stat_status /= 0 ) THROW_HARD('could not stat the project file; learn_pose_student')
        mtime = int(statbuf(10),int64)
        key   = polar_feats_contract(projpath%to_char(), mtime, build%spproj_field%get_noris(), params%box, params%box_crop, &
            &params%smpd, pinds, build%pftc%get_pftsz(), params%kfromto, wk, trim(params%ctf) /= 'no')
        call cache%new(fname_feats)
        call cache%create(feats, maxval(pinds), status)
        l_cache   = status == 0
        vol_hash  = file_hash(params%vols(1))
        ! the table of the pass: the iteration's model from the complete field, the assignment opened
        allocate(eul_prob_tab :: tab)
        call tab%new_pipeline(params, build, pinds)
        call tab%set_model_from_field
        call tab%open_assignment
        call rows%open_write(fname_rows, nptcls, params%nspace)
        call fopen(stat_unit, fname_stats, status='replace', action='write', iostat=io_stat)
        call fileiochk('learn_pose_student; open '//fname_stats%to_char(), io_stat)
        write(stat_unit,'(A)') '# iptcl eo temp nvalid iproj_mode pmax neff k99 nkept wtail lmin inpl_mode x_mode y_mode '//&
            &'dist_mode_deg iproj_ref w_ref rank_ref in_cone in_k99 wsq logpow dfx tazim'
        call write_grid
        ! the training set in memory: the tensors, the rows' losses and temperatures
        allocate(x_all(feats%get_pftsz(), nk, 2, nptcls), ctf2_all(nk, nptcls), side_all(nptcls), source=0.)
        allocate(dists_all(params%nspace, nptcls), temps_all(nptcls), source=0.)
        allocate(mode_all(nptcls), source=0)
        rt_prep  = toc(t0)
        rt_batch = 0.
        rt_fill  = 0.
        rt_rows  = 0.
        allocate(dists(params%nspace))
        call reads%new(params)
        call reads%expect(pinds)
        invalid = 0.1 * huge(invalid)
        do ibatch = 1,nbatches
            b0      = batches(ibatch,1)
            b1      = batches(ibatch,2)
            batchsz = b1 - b0 + 1
            t0 = tic()
            call build_batch_particles3D(params, build, batchsz, pinds(b0:b1), ptcl_imgs, ptcl_imgs_pad, reads)
            rt_batch = rt_batch + toc(t0)
            t0 = tic()
            call tab%fill_tab_range(b0, b1)
            rt_fill = rt_fill + toc(t0)
            t0 = tic()
            do i = b0,b1
                iptcl = pinds(i)
                call tab%get_posterior_row(i, params%prob_beta, STUDENT_KMAX_ROW, rec, dists)
                call rows%put(rec, dists)
                call feats%build(build%pftc, iptcl, feat)
                if( l_cache ) call cache%put(iptcl, feat)
                x_all(:,:,:,i)  = feat%x
                ctf2_all(:,i)   = feat%ctf2
                side_all(i)     = feat%logpow
                dists_all(:,i)  = dists(1:params%nspace)
                temps_all(i)    = rec%temp
                if( rec%imin > 0 ) mode_all(i) = tab%ref_proj(rec%imin)
                call write_stats_line
            enddo
            rt_rows = rt_rows + toc(t0)
            call tab%assign_batch(b0, b1)
        enddo
        call tab%close_assignment
        profile = tab%get_fill_profile()
        call rows%close
        call fclose(stat_unit)
        if( l_cache ) call cache%seal(key)
        call reads%kill
        call feat%kill
        call rec%kill
        call clean_batch_particles3D(build, ptcl_imgs, ptcl_imgs_pad)
        call tab%kill
        deallocate(tab)
        call build%pftc%kill
        call build%esig%kill
        write(logfhandle,'(A,F9.2,A,F9.2,A)') '>>> POSE STUDENT ROWS: fill ', rt_fill, ' s wall, ', &
            &rt_fill / real(nptcls,timer_int_kind) * 1000., ' ms per particle'
        ! the students, cheapest first (the widest only on a large set), trained on the rows with the shift augmentation of
        ! the centring's disk; a wider design replaces the kept one when it recovers clearly more held-out seeds (a seed recovered, ruling R1: the
        ! teacher's direction or its mirror partner inside a folded cone within the cap and within the route's bound,
        ! the dense search's cell prob_athres or the grid's spacing), the paired difference over the micrographs
        ! resampled by micrograph
        t0 = tic()
        shift_aug = 0.
        if( params%l_shift_map ) shift_aug = STUDENT_SHIFT_AUG_NSIG * params%shift_sigma / params%smpd
        nheld = count(held)
        allocate(held_inds(nheld), group_held(nheld), k99_held(nheld), mode_held_s(nheld), source=0)
        allocate(conf_held(nheld), radius_held(nheld), source=0.)
        allocate(recov(nheld, size(STUDENT_WIDTHS)), source=.false.)
        cone_bound = 0.
        j = 0
        do i = 1,nptcls
            if( .not. held(i) ) cycle
            j = j + 1
            held_inds(j)  = i
            group_held(j) = group_of(i)
        enddo
        nthr_here = params%nthr
        !$ nthr_here = omp_get_max_threads()
        allocate(prob(params%nspace, nthr_here), source=0.)
        iw_best    = 0
        recov_cand = 0.
        loss_cand  = 0.
        gain_q     = 0.
        rt_cand    = 0.
        epoch_cand = 0
        rt_heldpred = 0.
        rt_boot     = 0.
        rt_write    = 0.
        ntry = size(STUDENT_WIDTHS)
        if( nptcls - nheld < STUDENT_WIDEST_MIN_ROWS ) ntry = size(STUDENT_WIDTHS) - 1
        call choose_student_design('POSE STUDENT', STUDENT_WIDTHS, ntry, nheld, group_held, rnd_run_seed(), try_design, &
            &keep_design, recov, iw_best, ntried, gain_q, rt_boot)
        ! the kept design back from its file (the loop's last candidate may have lost), its threshold on the held-out
        ! rows and the held-out table
        tb = tic()
        call student%kill
        call student%read(params%fname)
        call student%prepare_predict(nthr_here)
        !$omp parallel do default(shared) private(j,i,ithr,pred) schedule(dynamic) proc_bind(close)
        do j = 1,nheld
            i    = held_inds(j)
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call student%predict_tensor(ithr, x_all(:,:,:,i), ctf2_all(:,i), side_all(i), prob(:,ithr))
            call student%cone_of(prob(:,ithr), pred)
            conf_held(j)      = pred%conf
            k99_held(j)       = pred%k99
            radius_held(j)    = pred%radius
            mode_held_s(j)    = pred%mode
            recov(j,iw_best)  = l_recovered(pred, mode_all(i))
        enddo
        !$omp end parallel do
        call student%set_cone_bound(cone_bound)
        call set_conf_threshold
        call student%set_threshold(conf_threshold, STUDENT_CONF_FRACTION, heldout_seeded, heldout_recovered)
        call student%write(params%fname)
        call check_reread_population
        call write_heldout_table
        call student%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
            &time_train_s, nweights)
        rt_kept  = toc(tb)
        rt_train = toc(t0)
        write(logfhandle,'(A,I0,A,F8.5,A,F8.4,A,F8.4,A,F6.2,A,F9.1,A)') '>>> POSE STUDENT KEPT: width ', STUDENT_WIDTHS(iw_best), &
            &', confidence threshold ', conf_threshold, ' seeding ', heldout_seeded, ' of the held-out rows with ', &
            &heldout_recovered, ' of their teacher directions in the folded cone (its radius within ', cone_bound, &
            &' deg); the designs ', rt_train, ' s'
        rt_all = toc(t_all)
        call write_header
        if( L_BENCH_GLOB ) call write_learn_bench
        call student%kill
        call feats%kill
        call cache%kill
        call o%kill
        deallocate(dists, pinds, batches, wk, x_all, ctf2_all, side_all, dists_all, temps_all, held, mode_all, group_of, xrefs)
        deallocate(held_inds, group_held, conf_held, radius_held, k99_held, mode_held_s, recov, prob)
        if( allocated(statbuf) ) deallocate(statbuf)

    contains

        !> BENCH/POSE_STUDENT_LEARN_BENCH.txt: the context, the rows' costs, every design's training broken into the
        !! batch assembly, the gradients (the parallel region), the optimiser steps and the held-out evaluation, the
        !! design choice and the kept model, with the share accounted for
        subroutine write_learn_bench
            use simple_bench_io, only: open_bench_report
            integer  :: fnr, jw
            real(dp) :: rt_designs, acc
            call open_bench_report(fnr, 'POSE_STUDENT_LEARN')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,i0)')   'student particles (rows)            : ', nptcls
            write(fnr,'(a,i0)')   'student held-out particles          : ', nheld
            write(fnr,'(a,i0)')   'student grid directions (nspace)    : ', params%nspace
            write(fnr,'(a,i0)')   'student classes (mirror pairs)      : ', params%nspace / 2
            write(fnr,'(a,f0.2)') 'student band (lp, A)                : ', params%lp
            write(fnr,'(a,i0)')   'student rings                       : ', params%kfromto(2) - params%kfromto(1) + 1
            write(fnr,'(a,i0,a,i0)') 'student half-ring angles (calc/net) : ', feats%get_pftsz_in(), ' / ', feats%get_pftsz()
            write(fnr,'(a,i0,a,i0)') 'student box / box_crop              : ', params%box, ' / ', params%box_crop
            write(fnr,'(a,i0)')   'student threads                     : ', nthr_here
            write(fnr,'(a,l1)')   'student feature cache used          : ', l_cache
            write(fnr,'(a,f0.2)') 'student shift augmentation (px)     : ', shift_aug
            write(fnr,'(a,i0)')   'student designs trained             : ', ntried
            write(fnr,'(a,i0)')   'student kept width                  : ', STUDENT_WIDTHS(iw_best)
            write(fnr,'(a,i0)')   'student kept weights                : ', nweights
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'rows setup (halves, sigma2, references, calculator):', rt_prep
            write(fnr,'(a,1x,f0.2)') 'rows particle reading and preparation          :', rt_batch
            write(fnr,'(a,1x,f0.2)') 'rows table fill (the teacher)                  :', rt_fill
            write(fnr,'(a,1x,f0.3)') 'rows table fill per particle (ms)              :', 1000. * rt_fill / real(max(1,nptcls))
            write(fnr,'(a,1x,f0.2)') 'rows records and features                      :', rt_rows
            rt_designs = 0.d0
            do jw = 1,ntried
                acc = prof_cand(jw)%rt_assemble + prof_cand(jw)%net%rt_grad + prof_cand(jw)%net%rt_step + &
                    &prof_cand(jw)%net%rt_eval
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), ' training total               :', rt_cand(jw)
                write(fnr,'(a,i0,a,1x,i0,1x,i0,1x,i0)') 'design width ', STUDENT_WIDTHS(jw), &
                    &' epochs / batch / batches per epoch:', prof_cand(jw)%epochs, prof_cand(jw)%batch_size, &
                    &prof_cand(jw)%nbatches_epoch
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), '   batch assembly             :', &
                    &prof_cand(jw)%rt_assemble
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), '   gradients (parallel)       :', &
                    &prof_cand(jw)%net%rt_grad
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), '   reduction and Adam steps   :', &
                    &prof_cand(jw)%net%rt_step
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), '   held-out evaluation        :', &
                    &prof_cand(jw)%net%rt_eval
                write(fnr,'(a,i0,a,1x,i0)')   'design width ', STUDENT_WIDTHS(jw), '   optimiser steps            :', &
                    &prof_cand(jw)%net%nsteps
                if( prof_cand(jw)%net%nsteps > 0 ) write(fnr,'(a,i0,a,1x,f0.1,1x,i0,1x,f0.3)') 'design width ', &
                    &STUDENT_WIDTHS(jw), '   chunks per step / team / occupancy:', &
                    &real(prof_cand(jw)%net%nchunks,dp) / real(prof_cand(jw)%net%nsteps,dp), prof_cand(jw)%net%nthr, &
                    &real(prof_cand(jw)%net%nchunks,dp) / max(1.d0, real(prof_cand(jw)%net%nslots,dp))
                if( prof_cand(jw)%net%nsamples_train > 0 ) write(fnr,'(a,i0,a,1x,f0.1)') 'design width ', &
                    &STUDENT_WIDTHS(jw), '   training samples per second:', &
                    &real(prof_cand(jw)%net%nsamples_train,dp) / max(prof_cand(jw)%net%rt_grad + prof_cand(jw)%net%rt_step, 1.d-9)
                if( rt_cand(jw) > 0. ) write(fnr,'(a,i0,a,1x,f0.1)') 'design width ', STUDENT_WIDTHS(jw), &
                    &'   % of its training accounted  :', 100.d0 * acc / real(rt_cand(jw),dp)
                write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), '   held-out prediction        :', rt_heldpred(jw)
                if( jw > 1 ) write(fnr,'(a,i0,a,1x,f0.2)') 'design width ', STUDENT_WIDTHS(jw), &
                    &'   resampled design choice    :', rt_boot(jw)
                rt_designs = rt_designs + real(rt_cand(jw) + rt_heldpred(jw) + rt_boot(jw),dp)
            enddo
            write(fnr,'(a,1x,f0.2)') 'model writes                                   :', rt_write
            write(fnr,'(a,1x,f0.2)') 'kept model: reread, held-out, threshold, table :', rt_kept
            write(fnr,'(a,1x,f0.2)') 'designs, choice and kept model                 :', rt_train
            write(fnr,'(a,1x,f0.2)') 'total time                                     :', rt_all
            write(fnr,'(a,1x,f0.2)') '% accounted for                                :', 100. * real(rt_prep + rt_batch + rt_fill + &
                &rt_rows + real(rt_designs) + rt_write + rt_kept) / max(rt_all, real(1.e-9,timer_int_kind))
            write(fnr,'(a)') ''
            call write_fill_profile_bench(fnr, 'rows', profile, params%nstates)
            call fclose(fnr)
        end subroutine write_learn_bench

        !> design iw (choose_student_design's trial): built with the grid's reprojections as its templates, trained on
        !! the rows, and judged on the held-out rows, recov_iw the seeds it recovers
        subroutine try_design( iw_here, recov_iw )
            integer, intent(in)    :: iw_here
            logical, intent(inout) :: recov_iw(:)
            integer :: jj, ii, ithr_loc
            type(student_pred) :: pred_loc
            call student%new(feats%get_pftsz_in(), params%kfromto, params%box, params%box_crop, params%smpd, params%mskdiam, &
                &trim(params%pgrp), params%nspace, params%lp, wk, trim(params%ctf) /= 'no', STUDENT_WIDTHS(iw_here), rnd_run_seed())
            call student%init_readout(xrefs)
            if( l_by_stack )then
                call student%set_provenance(projpath%to_char(), polar_feats_hash_ints(pinds), vol_hash, 'every_10th_selected_stack')
            else
                call student%set_provenance(projpath%to_char(), polar_feats_hash_ints(pinds), vol_hash, 'every_10th_selected_particle')
            endif
            if( iw_here == 1 ) cone_bound = max(params%prob_athres, student%grid_spacing())
            call student%set_cone_bound(cone_bound)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.2,A)') '>>> POSE STUDENT TRAINING: width ', STUDENT_WIDTHS(iw_here), ', ', &
                &nptcls - nheld, ' rows, ', nheld, ' held out, shifts augmented within ', shift_aug, ' px'
            call student%train(nptcls, x_all, ctf2_all, side_all, dists_all, temps_all, held, params%prob_beta, &
                &shift_aug, STUDENT_MAX_EPOCHS, STUDENT_PATIENCE_EPOCHS, params%nthr, logfhandle)
            call student%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
                &time_train_s, nweights)
            rt_cand(iw_here)    = time_train_s
            loss_cand(iw_here)  = loss_held
            epoch_cand(iw_here) = best_epoch
            prof_cand(iw_here)  = student%get_train_profile()
            ! the held-out seeds this design recovers
            tb = tic()
            call student%prepare_predict(nthr_here)
            !$omp parallel do default(shared) private(jj,ii,ithr_loc,pred_loc) schedule(dynamic) proc_bind(close)
            do jj = 1,nheld
                ii       = held_inds(jj)
                ithr_loc = 1
                !$ ithr_loc = omp_get_thread_num() + 1
                call student%predict_tensor(ithr_loc, x_all(:,:,:,ii), ctf2_all(:,ii), side_all(ii), prob(:,ithr_loc))
                call student%cone_of(prob(:,ithr_loc), pred_loc)
                conf_held(jj)   = pred_loc%conf
                k99_held(jj)    = pred_loc%k99
                radius_held(jj) = pred_loc%radius
                mode_held_s(jj) = pred_loc%mode
                recov_iw(jj)    = l_recovered(pred_loc, mode_all(ii))
            enddo
            !$omp end parallel do
            recov_cand(iw_here)  = real(count(recov_iw(1:nheld))) / real(nheld)
            rt_heldpred(iw_here) = toc(tb)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F9.5,A,F7.4,A,F7.4,A,F9.1,A)') '>>> POSE STUDENT TRAINED: width ', &
                &STUDENT_WIDTHS(iw_here), ', best epoch ', best_epoch, ' of ', epochs_run, ', held-out loss ', loss_held, &
                &', held-out mode agreement ', mode_held, ', held-out seeds recovered ', recov_cand(iw_here), ', ', time_train_s, ' s'
        end subroutine try_design

        !> design iw kept (choose_student_design): its model written
        subroutine keep_design( iw_here )
            integer, intent(in) :: iw_here
            integer(timer_int_kind) :: tw
            tw = tic()
            call student%write(params%fname)
            rt_write = rt_write + toc(tw)
        end subroutine keep_design

        !> a held-out seed recovered (ruling R1): the folded cone within the cap and the route's bound, and the
        !! teacher's direction or its mirror partner among the cone's folded directions
        logical function l_recovered( pred_loc, imode_loc )
            type(student_pred), intent(in) :: pred_loc
            integer,            intent(in) :: imode_loc
            integer :: jj, ip
            l_recovered = .false.
            if( pred_loc%k99 > KMAX_STUDENT .or. pred_loc%radius > cone_bound ) return
            ip = student%fold_index(imode_loc)
            do jj = 1,pred_loc%ncone
                if( student%fold_index(pred_loc%cone(jj)) == ip ) l_recovered = .true.
            enddo
        end function l_recovered

        !> the held-out rows the route would seed: a cone within the cap and within the bound bound_here
        subroutine eligible_held( bound_here, elig )
            real,    intent(in)  :: bound_here
            logical, intent(out) :: elig(nheld)
            integer :: jj
            do jj = 1,nheld
                elig(jj) = k99_held(jj) <= KMAX_STUDENT .and. radius_held(jj) <= bound_here
            enddo
        end subroutine eligible_held

        !> the confidence threshold on the held-out rows the route would seed (student_conf_threshold), the seeded
        !! share of the held-out rows and the recovered share among them
        subroutine set_conf_threshold
            logical :: elig(nheld)
            integer :: nrec
            call eligible_held(cone_bound, elig)
            call student_conf_threshold(nheld, conf_held, elig, recov(:,iw_best), STUDENT_CONF_FRACTION, conf_threshold, nsel, nrec)
            heldout_seeded    = real(nsel) / real(nheld)
            heldout_recovered = 0.
            if( nsel > 0 ) heldout_recovered = real(nrec) / real(nsel)
        end subroutine set_conf_threshold

        !> the model read back from its file seeds exactly the calibrated held-out rows (its threshold and bound
        !! survive the file)
        subroutine check_reread_population
            logical :: elig(nheld)
            real    :: thr, cf, hs, hr, bound_read
            call reread%read(params%fname)
            call reread%get_threshold(thr, cf, hs, hr)
            bound_read = reread%get_cone_bound()
            call reread%kill
            call eligible_held(bound_read, elig)
            if( count(elig .and. conf_held >= thr) /= nsel ) &
                &THROW_HARD('the model read back seeds other held-out rows than the calibration; learn_pose_student')
        end subroutine check_reread_population

        !> the held-out table: per held-out row its micrograph, confidence, 99% set size, the teacher's and the
        !! student's modes and whether the teacher's lies in the cone, for the kept design
        subroutine write_heldout_table
            type(string) :: fname_held
            integer :: jj, i_here
            fname_held = base%to_char()//'_heldout.txt'
            call fopen(held_unit, fname_held, status='replace', action='write', iostat=io_stat)
            call fileiochk('learn_pose_student; open '//fname_held%to_char(), io_stat)
            write(held_unit,'(A)') '# iptcl group conf k99 radius_deg teacher_mode student_mode recovered'
            do jj = 1,nheld
                i_here = held_inds(jj)
                write(held_unit,'(I0,1X,I0,1X,F8.5,1X,I0,1X,F8.3,1X,I0,1X,I0,1X,I0)') pinds(i_here), group_held(jj), conf_held(jj), &
                    &k99_held(jj), radius_held(jj), mode_all(i_here), mode_held_s(jj), merge(1, 0, recov(jj,iw_best))
            enddo
            call fclose(held_unit)
            call fname_held%kill
        end subroutine write_heldout_table

        !> the row's statistics and its error against the particle's stored direction (the truth on a truth set, the
        !! reference alignment on a real one): the angle between the mode's direction and the stored one across the
        !! point group, the stored direction's nearest grid point with its mass and rank in the row, whether it lies in
        !! the cone (the KMAX_STUDENT entries of largest mass) and in the 99% set
        subroutine write_stats_line
            call build%spproj_field%get_ori(iptcl, o)
            euls_ref   = o%get_euler()
            iproj_ref  = build%eulspace%find_closest_proj(o)
            iproj_mode = 0
            euldist    = -1.
            if( rec%imin > 0 )then
                iproj_mode = tab%ref_proj(rec%imin)
                euls_mode  = build%eulspace%get_euler(iproj_mode)
                call build%pgrpsyms%sym_dists(euls_mode, euls_ref, euls_sym, euldist, inpldist)
            endif
            ! the stored direction's mass and rank in the row
            temp       = temperature_or_argmin(rec%temp)
            w_ref      = 0.
            rank_ref   = 0
            nvalid_ref = 0
            if( rec%imin > 0 .and. valid_dist(dists(iproj_ref)) )then
                dmin = dists(rec%imin)
                z    = 0.d0
                do j = 1,params%nspace
                    if( .not. valid_dist(dists(j)) ) cycle
                    nvalid_ref = nvalid_ref + 1
                    if( temp * (dists(j) - dmin) < 80. ) z = z + exp(-real(temp,dp) * real(dists(j) - dmin,dp))
                    if( dists(j) < dists(iproj_ref) ) rank_ref = rank_ref + 1
                enddo
                rank_ref = rank_ref + 1
                if( z > 0.d0 .and. temp * (dists(iproj_ref) - dmin) < 80. )&
                    &w_ref = real(exp(-real(temp,dp) * real(dists(iproj_ref) - dmin,dp)) / z)
            endif
            write(stat_unit,'(I0,1X,I0,1X,ES13.6,1X,I0,1X,I0,1X,F8.5,1X,F10.3,1X,I0,1X,I0,1X,ES11.4,1X,ES13.6,1X,I0,1X,F8.3,1X,F8.3,'//&
                &'1X,F8.3,1X,I0,1X,ES11.4,1X,I0,1X,I0,1X,I0,1X,ES13.6,1X,F9.4,1X,F8.4,1X,F8.5)') &
                &iptcl, build%spproj_field%get_eo(iptcl), rec%temp, rec%nvalid, iproj_mode, rec%pmax, rec%neff, rec%k99, &
                &rec%nkept, rec%wtail, rec%lmin, rec%inpl_mode, rec%x_mode, rec%y_mode, euldist, iproj_ref, w_ref, rank_ref, &
                &merge(1, 0, rank_ref >= 1 .and. rank_ref <= KMAX_STUDENT), merge(1, 0, rank_ref >= 1 .and. rank_ref <= rec%k99), &
                &feat%wsq, feat%logpow, build%spproj_field%get(iptcl, 'dfx'), feat%tazim
        end subroutine write_stats_line

        logical function valid_dist( d )
            real, intent(in) :: d
            valid_dist = ieee_is_finite(d) .and. d < invalid
        end function valid_dist

        !> the student grid: one line per direction, its index and Euler angles
        subroutine write_grid
            integer :: iproj
            real    :: euls(3)
            call fopen(grid_unit, fname_grid, status='replace', action='write', iostat=io_stat)
            call fileiochk('learn_pose_student; open '//fname_grid%to_char(), io_stat)
            write(grid_unit,'(A)') '# iproj e1 e2 e3'
            do iproj = 1,build%eulspace%get_noris()
                euls = build%eulspace%get_euler(iproj)
                write(grid_unit,'(I0,3(1X,F10.4))') iproj, euls
            enddo
            call fclose(grid_unit)
        end subroutine write_grid

        !> the rows' provenance and cost: the project, the selection, the grid, the band, the references and their
        !! hash, the objective's temperatures and shell weights, the fill's budget and its cost per particle, and the
        !! training's record
        subroutine write_header
            integer :: k
            call fopen(head_unit, fname_header, status='replace', action='write', iostat=io_stat)
            call fileiochk('learn_pose_student; open '//fname_header%to_char(), io_stat)
            write(head_unit,'(A)')        '# pose student teacher rows'
            write(head_unit,'(A)')        'rows_version=1'
            write(head_unit,'(A,A)')      'projfile=', projpath%to_char()
            write(head_unit,'(A,I0)')     'project_mtime=', mtime
            write(head_unit,'(A,I0)')     'project_rows=', build%spproj_field%get_noris()
            write(head_unit,'(A,I0)')     'nptcls=', nptcls
            write(head_unit,'(A,I0)')     'selection_hash=', polar_feats_hash_ints(pinds)
            write(head_unit,'(A,I0)')     'nspace=', params%nspace
            write(head_unit,'(A,A)')      'pgrp=', trim(params%pgrp)
            write(head_unit,'(A,F10.4)')  'lp=', params%lp
            write(head_unit,'(A,I0,A,I0)') 'kfromto=', params%kfromto(1), ',', params%kfromto(2)
            write(head_unit,'(A,I0)')     'pftsz=', feats%get_pftsz_in()
            write(head_unit,'(A,I0)')     'nrots=', feats%get_nrots()
            write(head_unit,'(A,I0)')     'box=', params%box
            write(head_unit,'(A,I0)')     'box_crop=', params%box_crop
            write(head_unit,'(A,F10.5)')  'smpd=', params%smpd
            write(head_unit,'(A,F10.3)')  'mskdiam=', params%mskdiam
            write(head_unit,'(A,A)')      'ctf=', trim(params%ctf)
            write(head_unit,'(A,A)')      'objfun=', trim(params%objfun)
            write(head_unit,'(A,F8.4)')   'prob_beta=', params%prob_beta
            write(head_unit,'(A,F8.3)')   'trs=', params%trs
            write(head_unit,'(A,F8.3)')   'shift_sigma=', params%shift_sigma
            write(head_unit,'(A,F8.3)')   'prob_athres=', params%prob_athres
            write(head_unit,'(A,I0)')     'projs_ns=', profile%projs_ns
            write(head_unit,'(A,I0)')     'kmax_row=', STUDENT_KMAX_ROW
            write(head_unit,'(A,I0)')     'kmax_student=', KMAX_STUDENT
            write(head_unit,'(A,A)')      'vol1=', params%vols(1)%to_char()
            write(head_unit,'(A,I0)')     'vol1_hash=', vol_hash
            write(head_unit,'(A,L1)')     'feats_cached=', l_cache
            write(head_unit,'(A,I0)')     'nthr=', params%nthr
            write(head_unit,'(A,F10.3)')  'time_prep_s=', rt_prep
            write(head_unit,'(A,F10.3)')  'time_batch_prep_s=', rt_batch
            write(head_unit,'(A,F10.3)')  'time_fill_s=', rt_fill
            write(head_unit,'(A,F10.3)')  'time_rows_feats_s=', rt_rows
            write(head_unit,'(A,F10.3)')  'time_train_s=', rt_train
            write(head_unit,'(A,F10.3)')  'time_total_s=', rt_all
            write(head_unit,'(A,F12.5)')  'fill_s_per_ptcl=', rt_fill / real(nptcls,timer_int_kind)
            write(head_unit,'(A,F12.5)')  'total_s_per_ptcl=', rt_all / real(nptcls,timer_int_kind)
            write(head_unit,'(A,A)')      'model=', params%fname%to_char()
            write(head_unit,'(A,I0)')     'student_width=', STUDENT_WIDTHS(iw_best)
            write(head_unit,'(A,I0)')     'student_nweights=', nweights
            do iw = 1,size(STUDENT_WIDTHS)
                if( rt_cand(iw) <= 0. ) cycle
                write(head_unit,'(A,I0,A,I0,A,F8.5,A,F9.5,A,I0,A,F9.1,A,F8.5)') 'candidate_width_', STUDENT_WIDTHS(iw), '=', &
                    &STUDENT_WIDTHS(iw), ' recovered=', recov_cand(iw), ' loss_heldout=', loss_cand(iw), ' best_epoch=', &
                    &epoch_cand(iw), ' time_s=', rt_cand(iw), ' gain_quantile=', gain_q(iw)
            enddo
            write(head_unit,'(A,A)')      'holdout_rule=', trim(merge('every_10th_stack            ', 'every_10th_selected_particle', l_by_stack))
            write(head_unit,'(A,F8.5)')   'conf_threshold=', conf_threshold
            write(head_unit,'(A,F8.5)')   'conf_fraction=', STUDENT_CONF_FRACTION
            write(head_unit,'(A,F8.5)')   'heldout_seeded=', heldout_seeded
            write(head_unit,'(A,F8.5)')   'heldout_recovered=', heldout_recovered
            write(head_unit,'(A,F8.3)')   'cone_bound_deg=', cone_bound
            write(head_unit,'(A)')        'recovery=folded_sphere_within_cap_and_bound'
            write(head_unit,'(A,I0)')     'train_rows=', nptcls_train
            write(head_unit,'(A,I0)')     'heldout_rows=', nheldout
            write(head_unit,'(A,F8.3)')   'shift_aug_px=', shift_aug
            write(head_unit,'(A,I0)')     'max_epochs=', STUDENT_MAX_EPOCHS
            write(head_unit,'(A,I0)')     'patience_epochs=', STUDENT_PATIENCE_EPOCHS
            write(head_unit,'(A,I0)')     'epochs_run=', epochs_run
            write(head_unit,'(A,I0)')     'best_epoch=', best_epoch
            write(head_unit,'(A,ES14.6)') 'loss_train=', loss_train
            write(head_unit,'(A,ES14.6)') 'loss_heldout=', loss_held
            write(head_unit,'(A,F8.5)')   'mode_heldout=', mode_held
            write(head_unit,'(A)', advance='no') 'wk_even='
            do k = params%kfromto(1),params%kfromto(2)
                if( k > params%kfromto(1) ) write(head_unit,'(A)', advance='no') ','
                write(head_unit,'(ES14.6)', advance='no') wk(k,1)
            enddo
            write(head_unit,*)
            write(head_unit,'(A)', advance='no') 'wk_odd='
            do k = params%kfromto(1),params%kfromto(2)
                if( k > params%kfromto(1) ) write(head_unit,'(A)', advance='no') ','
                write(head_unit,'(ES14.6)', advance='no') wk(k,2)
            enddo
            write(head_unit,*)
            call write_fill_profile_bench(head_unit, 'rows', profile, params%nstates)
            call fclose(head_unit)
        end subroutine write_header

    end subroutine learn_pose_student_exec

    ! ===== PREDICT_POSES =====

    !> predict_poses on its own parameters and builder: the model's grid, point group and band are set on the
    !! command line from its header before the parameters are made
    subroutine run_predict_poses( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target :: params
        type(builder)            :: build
        type(pose_student)       :: student
        type(string)             :: infile
        infile = cline%get_carg('infile')
        call student%read(infile)
        call cline%set('nspace', student%get_nspace())
        call cline%set('pgrp',   student%get_pgrp())
        call cline%set('lp',     student%get_lp())
        call student%kill
        call infile%kill
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            if( cline%get_iarg('nparts') > 1 )then
                call predict_poses_master(cline)
                return
            endif
            call cline%delete('nparts')
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=cline%defined('vol1'))
        call predict_poses_exec(params, build, cline)
        call build%kill_general_tbox
        if( params%l_distr_worker ) call qsys_declare_part_finished(params, string('predict_poses'))
        call simple_end('**** SIMPLE_PREDICT_POSES NORMAL STOP ****', print_simple=.false.)
    end subroutine run_predict_poses

    !> predict_poses over nparts parts: given vol1, the halves, sigma2 and the reprojection model made first; the parts'
    !! records merged and written once, their tables joined (doc/algorithms/pose_student.md)
    subroutine predict_poses_master( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target  :: params
        type(builder)             :: build
        type(string), allocatable :: part_bases(:)
        logical :: l_rebuilt
        call build%init_params_and_build_general_tbox(cline, params, do3d=cline%defined('vol1'))
        if( complete_volume_source_available(params) )then
            ! as predict_poses_exec prepares the registration, before any part reads it
            call ensure_field_halves(build%spproj, params%oritype, params%projfile)
            call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, 1, params%box, params%smpd, &
                &'POSE STUDENT SIGMA2', l_rebuilt)
            if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
            call set_bp_range3D(params, build, cline)
            call materialize_reprojection_model_from_volumes(params, build, cline, cleanup_pftc=.true.)
        endif
        call student_predict_parts(params, build, cline, fname_base(params%fname), part_bases)
        call build%spproj%write_segment_inside(params%oritype)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_PREDICT_POSES NORMAL STOP ****', print_simple=.false.)
    end subroutine predict_poses_master

    !> The student's prediction for the project's selected particles: per particle the features under the model's
    !! band and weights, the forward pass, the mode, the cone holding 99% of the posterior (at most KMAX_STUDENT
    !! directions) and the confidence, with the error against the particle's stored direction, written to the text
    !! table <fname>_pred.txt. Project: reads the particle field and the stacks; writes nothing
    subroutine predict_poses_exec( params, build, cline )
        type(parameters), target, intent(inout) :: params
        type(builder),    target, intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        integer, parameter       :: TH_NSTEPS = 6   !< the timed per-particle steps
        type(pose_student)       :: student
        type(ptcl_read_session)  :: reads
        type(image), allocatable :: ptcl_imgs(:), ptcl_imgs_pad(:)
        type(polar_feats)        :: feats
        type(string)             :: base, fname_pred
        real,    allocatable     :: wk(:,:), prob(:,:)
        integer, allocatable     :: pinds(:), batches(:,:), pinds_part(:)
        type(student_pred), allocatable :: preds(:)
        character(len=:), allocatable :: msg
        real,    allocatable     :: dist_mode(:), dist_cone(:), euls_pose(:,:)
        integer, allocatable     :: iproj_ref(:), rank_ref(:), in_cone(:)
        integer(timer_int_kind)  :: t0, t_all
        real(timer_int_kind)     :: rt_prep, rt_batch, rt_pred, rt_all, rt_table, rt_write, rt_chkout
        real(dp), allocatable    :: th_s(:,:)       ! per thread: features, forward, cone, reference error, shift check, registration
        integer,  allocatable    :: th_n(:)
        type(pftc_shsrch_grad), allocatable :: shsrch(:)
        real,    allocatable     :: sh_reg(:,:), corr_reg(:)
        real     :: lims_reg(2,2)
        logical  :: l_reg, l_rebuilt
        integer(int64)           :: clock_rate
        integer  :: nptcls, batchsz_max, nbatches, ibatch, b0, b1, batchsz, i, iptcl, nthr, pred_unit, io_stat, kfromto(2)
        integer  :: n_in_cone, n_within_two, n_capped, n_seeded, n_in_cone_seeded, route, n_wide, nchk, nfeat, chk_unit
        real(dp) :: ncone_sum, conf_sum, radius_sum
        real     :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered, spacing, cone, euls(3), cone_bound
        real,    allocatable :: chk_tv(:,:), chk_jacc(:,:), chk_fchg(:,:,:), chk_sorted(:)
        integer, allocatable :: chk_mode(:,:), chk_cone(:,:)
        type(string)         :: fname_chk
        logical  :: l_part
        t_all = tic()
        call check_pass_team('predict_poses')
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('the poses are ptcl3D''s; predict_poses')
        ! a part of a distributed caller (part, nparts, fromp, top, outfile: the commander checks they are given)
        ! writes its range to outfile; otherwise the run's own project copy is written
        l_part = cline%defined('part')
        t0 = tic()
        call student%read(params%infile)
        kfromto = student%get_kfromto()
        msg = student_geometry_mismatch(student%get_pftsz(), student%get_box(), student%get_box_crop(), student%get_smpd(), &
            &student%get_mskdiam(), student%get_with_ctf(), params%pftsz, params%box, params%box_crop, params%smpd, &
            &params%mskdiam, trim(params%ctf) /= 'no')
        if( len(msg) > 0 ) THROW_HARD(msg)
        params%kfromto = kfromto
        params%lp      = student%get_lp()
        call student%get_weights(wk)
        call feats%new(student%get_pftsz(), kfromto, params%box, wk)
        base = fname_base(params%fname)
        ! a part's tables carry its number: the parts of one caller share fname
        if( l_part ) base = base%to_char()//'_part'//int2str(params%part)
        fname_pred = base%to_char()//'_pred.txt'
        ! the selected particles (a part's range of them; a part may hold none and writes its range unchanged)
        l_reg = complete_volume_source_available(params)
        call build%spproj_field%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no selected particle; predict_poses')
        if( size(pinds) < 1 ) THROW_HARD('no selected particle; predict_poses')
        if( l_part )then
            call student_part_pinds(pinds, params%part, params%nparts, params%fromp, params%top, &
                &build%spproj_field%get_noris(), pinds_part, msg)
            if( len(msg) > 0 ) THROW_HARD(msg)
            call move_alloc(pinds_part, pinds)
            if( size(pinds) == 0 )then
                call write_empty_part
                return
            endif
        endif
        nptcls = size(pinds)
        if( nptcls < 1 ) THROW_HARD('no selected particle; predict_poses')
        call student%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
        spacing    = student%grid_spacing()
        cone_bound = student%get_cone_bound()
        if( .not. (cone_bound > 0.) )then
            ! a model trained before ruling R1 carries no bound: the parameter's, with the threshold calibrated without it
            cone_bound = max(params%prob_athres, spacing)
            THROW_WARN('the model carries no cone bound (trained before ruling R1): the route takes prob_athres or the grid spacing; its threshold was calibrated without the bound')
        endif
        fname_chk = base%to_char()//'_shiftcheck.txt'
        write(logfhandle,'(A,I0,A,I0,A,F7.2,A)') '>>> POSE STUDENT PREDICTION: ', nptcls, ' particles at nspace ', &
            &student%get_nspace(), ' and lp ', params%lp, ' A'
        ! the batches and the calculator: one reference slot (the features need the particles alone), or, given vol1,
        ! the map's reprojections at the student's grid and band with sigma2, as learn_pose_student prepares them,
        ! for the shift registration against the read-out pose
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = 1 + (nptcls - 1) / batchsz_max
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        if( l_reg )then
            call ensure_field_halves(build%spproj, params%oritype, params%projfile)
            call student_sigma2_for_pass(cline, params%projfile, params%oritype, 1, params%box, params%smpd, &
                &'POSE STUDENT SIGMA2', l_rebuilt)
            if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
            call set_bp_range3D(params, build, cline)
            if( any(params%kfromto /= kfromto) ) THROW_HARD('the map''s band differs from the model''s; predict_poses')
            ! a part reads the reprojection model its master made
            if( .not. l_part ) call materialize_reprojection_model_from_volumes(params, build, cline, cleanup_pftc=.true.)
            call read_reprojection_model(params, build, batchsz_max)
            call prep_sigmas_objfun(params, build)
            call alloc_ptcl_imgs(params, build, ptcl_imgs, ptcl_imgs_pad, batchsz_max)
            call build%vol%kill
            call build%vol_odd%kill
            call build%vol2%kill
            call build%pftc%memoize_refs(eulspace=build%eulspace)
        else
            call build%pftc%new(params, 1, [1,batchsz_max], kfromto)
            call alloc_ptcl_imgs(params, build, ptcl_imgs, ptcl_imgs_pad, batchsz_max)
        endif
        nthr = params%nthr
        !$ nthr = omp_get_max_threads()
        allocate(sh_reg(2,nptcls), corr_reg(nptcls), source=0.)
        if( l_reg )then
            lims_reg(:,1) = -STUDENT_SHIFT_REG_HALFWIDTH / params%smpd
            lims_reg(:,2) =  STUDENT_SHIFT_REG_HALFWIDTH / params%smpd
            allocate(shsrch(nthr))
            do i = 1,nthr
                call shsrch(i)%new_fixed(build, lims_reg, maxits=STUDENT_SHIFT_REG_MAXITS, coarse_init=.true.)
            enddo
            write(logfhandle,'(A,F6.1,A)') '>>> POSE STUDENT PREDICTION: shifts registered against the map''s reprojection at the read-out pose, within ', &
                &STUDENT_SHIFT_REG_HALFWIDTH, ' A of the stored shift'
        endif
        call student%prepare_predict(nthr)
        allocate(th_s(TH_NSTEPS, nthr), source=0.d0)
        allocate(th_n(nthr), source=0)
        call system_clock(count_rate=clock_rate)
        allocate(preds(nptcls), prob(student%get_nspace(), batchsz_max))
        allocate(dist_mode(nptcls), dist_cone(nptcls), source=-1.)
        allocate(euls_pose(3,nptcls), source=0.)
        allocate(iproj_ref(nptcls), rank_ref(nptcls), in_cone(nptcls), source=0)
        ! ruling R2's shift check on the first selected particles: the posterior and the pooled features of each
        ! particle shifted through the calculator's phase ramps against its own
        nchk  = min(nptcls, POSE_SHIFT_CHECK_NPTCLS)
        nfeat = student%get_nfeat()
        allocate(chk_tv(nchk, size(POSE_SHIFT_CHECK_MAGS)), chk_jacc(nchk, size(POSE_SHIFT_CHECK_MAGS)), &
            &chk_fchg(nchk, size(POSE_SHIFT_CHECK_MAGS), POSE_SHIFT_CHECK_NGROUPS), source=0.)
        allocate(chk_mode(nchk, size(POSE_SHIFT_CHECK_MAGS)), chk_cone(nchk, size(POSE_SHIFT_CHECK_MAGS)), source=0)
        rt_prep  = toc(t0)
        rt_batch = 0.
        rt_pred  = 0.
        call reads%new(params)
        call reads%expect(pinds)
        do ibatch = 1,nbatches
            b0      = batches(ibatch,1)
            b1      = batches(ibatch,2)
            batchsz = b1 - b0 + 1
            t0 = tic()
            call build_batch_particles3D(params, build, batchsz, pinds(b0:b1), ptcl_imgs, ptcl_imgs_pad, reads)
            rt_batch = rt_batch + toc(t0)
            t0 = tic()
            !$omp parallel do default(shared) private(i) schedule(dynamic) proc_bind(close)
            do i = b0,b1
                call predict_one(i, pinds(i), i - b0 + 1)
            enddo
            !$omp end parallel do
            rt_pred = rt_pred + toc(t0)
        enddo
        call reads%kill
        call clean_batch_particles3D(build, ptcl_imgs, ptcl_imgs_pad)
        call build%pftc%kill
        ! the table and its summary
        t0 = tic()
        call fopen(pred_unit, fname_pred, status='replace', action='write', iostat=io_stat)
        call fileiochk('predict_poses; open '//fname_pred%to_char(), io_stat)
        write(pred_unit,'(A)') '# iptcl eo mode pmax k99 ncone conf radius_deg iproj_ref dist_mode_deg dist_cone_deg rank_ref in_cone route e1 e2 e3 x y'
        n_in_cone    = 0
        n_within_two = 0
        n_capped     = 0
        n_seeded     = 0
        n_in_cone_seeded = 0
        n_wide       = 0
        ncone_sum    = 0.d0
        conf_sum     = 0.d0
        radius_sum   = 0.d0
        do i = 1,nptcls
            iptcl = pinds(i)
            ! the route (ruling R1): the neighbourhood search from the student's folded cone when the cone is within
            ! the cap, its radius within the model's bound (the dense search's cell) and the confidence at or above
            ! the model's threshold; the dense search otherwise. A seeded particle takes the mode (the heavier of the
            ! pair) as its direction (the in-plane rotation and the shift stay the search's; the search tries the
            ! antipode too), the cone's radius at least the grid's spacing; the confidence and the route in their slots
            route = 2
            if( preds(i)%k99 <= KMAX_STUDENT .and. preds(i)%radius <= cone_bound .and. preds(i)%conf >= conf_threshold ) route = 1
            cone = 0.
            if( route == 1 )then
                euls = student%grid_euler(preds(i)%mode)
                euls(3) = build%spproj_field%e3get(iptcl)
                call build%spproj_field%set_euler(iptcl, euls)
                cone = max(preds(i)%radius, spacing)
                n_seeded = n_seeded + 1
                n_in_cone_seeded = n_in_cone_seeded + in_cone(i)
            endif
            call build%spproj_field%set(iptcl, 'cone',  cone)
            call build%spproj_field%set(iptcl, 'pconf', preds(i)%conf)
            call build%spproj_field%set(iptcl, 'route', real(route))
            write(pred_unit,'(I0,1X,I0,1X,I0,1X,F8.5,1X,I0,1X,I0,1X,F8.5,1X,F8.3,1X,I0,1X,F8.3,1X,F8.3,1X,I0,1X,I0,1X,I0,5(1X,F8.3))') &
                &iptcl, build%spproj_field%get_eo(iptcl), preds(i)%mode, preds(i)%pmax, preds(i)%k99, preds(i)%ncone, &
                &preds(i)%conf, preds(i)%radius, iproj_ref(i), dist_mode(i), dist_cone(i), rank_ref(i), in_cone(i), route, &
                &euls_pose(:,i), sh_reg(:,i)
            n_in_cone = n_in_cone + in_cone(i)
            if( dist_mode(i) >= 0. .and. dist_mode(i) <= 12. ) n_within_two = n_within_two + 1
            if( preds(i)%k99 > KMAX_STUDENT ) n_capped = n_capped + 1
            if( preds(i)%radius > cone_bound ) n_wide = n_wide + 1
            ncone_sum  = ncone_sum  + real(preds(i)%ncone,dp)
            conf_sum   = conf_sum   + real(preds(i)%conf,dp)
            radius_sum = radius_sum + real(preds(i)%radius,dp)
        enddo
        call fclose(pred_unit)
        rt_table = toc(t0)
        ! the records: a part's range to its document, the run's own project copy otherwise
        t0 = tic()
        if( l_part )then
            call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [params%fromp, params%top], isegment=PTCL3D_SEG)
        else
            call build%spproj%write_segment_inside(params%oritype)
        endif
        rt_write = toc(t0)
        write(logfhandle,'(A,F8.4,A,F8.5,A,F8.4,A)') '>>> POSE STUDENT PREDICTION: seeded (the neighbourhood search from the cone) ', &
            &real(n_seeded) / real(nptcls), ' at the confidence threshold ', conf_threshold, ', the stored direction in the cone of ', &
            &real(n_in_cone_seeded) / real(max(1,n_seeded)), ' of the seeded'
        write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,F6.2,A)') '>>> POSE STUDENT PREDICTION: stored direction in the folded cone ', &
            &real(n_in_cone) / real(nptcls), ', folded mode within two grid steps (12 deg) ', real(n_within_two) / real(nptcls), &
            &', cone capped at KMAX_STUDENT ', real(n_capped) / real(nptcls), ', cone wider than the bound ', &
            &real(n_wide) / real(nptcls), ' (', cone_bound, ' deg)'
        write(logfhandle,'(A,F8.2,A,F8.4,A,F8.2,A)') '>>> POSE STUDENT PREDICTION: mean cone size ', ncone_sum / real(nptcls,dp), &
            &', mean confidence ', conf_sum / real(nptcls,dp), ', mean cone radius ', radius_sum / real(nptcls,dp), ' deg'
        write(logfhandle,'(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A)') '>>> POSE STUDENT PREDICTION: prep ', rt_prep, ' s, particles ', &
            &rt_batch, ' s, features and forward ', rt_pred, ' s, total ', rt_all, ' s'
        write(logfhandle,'(A,F9.3,A)') '>>> POSE STUDENT PREDICTION: ', rt_pred / real(nptcls,timer_int_kind) * 1000., &
            &' ms per particle (features and forward)'
        t0 = tic()
        call write_shift_check
        rt_chkout = toc(t0)
        rt_all = toc(t_all)
        if( L_BENCH_GLOB ) call write_predict_bench
        call student%kill
        call feats%kill
        if( allocated(shsrch) )then
            do i = 1,size(shsrch)
                call shsrch(i)%kill
            enddo
            deallocate(shsrch)
        endif
        deallocate(preds, prob, dist_mode, dist_cone, euls_pose, sh_reg, corr_reg, iproj_ref, rank_ref, in_cone, pinds, batches, wk)
        deallocate(chk_tv, chk_jacc, chk_fchg, chk_mode, chk_cone, th_s, th_n)
        call fname_chk%kill

    contains

        !> a part whose range holds no selected particle: the table's header and the range's records unchanged, the
        !! contribution the merge expects
        subroutine write_empty_part
            call fopen(pred_unit, fname_pred, status='replace', action='write', iostat=io_stat)
            call fileiochk('predict_poses; open '//fname_pred%to_char(), io_stat)
            write(pred_unit,'(A)') '# iptcl eo mode pmax k99 ncone conf radius_deg iproj_ref dist_mode_deg dist_cone_deg rank_ref in_cone route e1 e2 e3 x y'
            call fclose(pred_unit)
            call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [params%fromp, params%top], isegment=PTCL3D_SEG)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> POSE STUDENT PREDICTION: no selected particle in the part''s range ', &
                &params%fromp, '-', params%top, '; its records written unchanged'
            call student%kill
            call feats%kill
            deallocate(pinds, wk)
        end subroutine write_empty_part

        !> BENCH/POSE_STUDENT_PREDICT_BENCH.txt: the context, the phases, the prediction's per-particle steps summed
        !! over the threads (thread-seconds) and the team's load balance, with the share accounted for
        subroutine write_predict_bench
            use simple_bench_io, only: open_bench_report
            character(len=*), parameter :: STEP_NAMES(TH_NSTEPS) = [character(len=32) :: 'features (pooled tensor)', &
                &'forward pass (HHCNN)', 'folded cone and readout pose', 'reference direction error', 'shift check (R2)', &
                &'shift registration (vol1)']
            integer  :: fnr, k
            real(dp) :: th_tot(nthr), busy
            call open_bench_report(fnr, 'POSE_STUDENT_PREDICT')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,i0)')   'student particles predicted         : ', nptcls
            write(fnr,'(a,i0)')   'student grid directions (nspace)    : ', student%get_nspace()
            write(fnr,'(a,f0.2)') 'student band (lp, A)                : ', params%lp
            write(fnr,'(a,i0)')   'student rings                       : ', kfromto(2) - kfromto(1) + 1
            write(fnr,'(a,i0,a,i0)') 'student half-ring angles (calc/net) : ', feats%get_pftsz_in(), ' / ', feats%get_pftsz()
            write(fnr,'(a,i0,a,i0)') 'student box / box_crop              : ', params%box, ' / ', params%box_crop
            write(fnr,'(a,i0)')   'student threads                     : ', nthr
            write(fnr,'(a,i0)')   'student batches                     : ', nbatches
            write(fnr,'(a,i0)')   'student shift-checked particles     : ', nchk
            write(fnr,'(a,l1)')   'student distributed part            : ', l_part
            write(fnr,'(a,l1)')   'student shift registration (vol1)   : ', l_reg
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'model read and calculator setup      :', rt_prep
            write(fnr,'(a,1x,f0.2)') 'particle reading and preparation     :', rt_batch
            write(fnr,'(a,1x,f0.2)') 'prediction (parallel, wall)          :', rt_pred
            write(fnr,'(a,1x,f0.3)') 'prediction per particle (ms, wall)   :', 1000. * rt_pred / real(max(1,nptcls))
            write(fnr,'(a,1x,f0.2)') 'routes and prediction table          :', rt_table
            write(fnr,'(a,1x,f0.2)') 'records write                        :', rt_write
            write(fnr,'(a,1x,f0.2)') 'shift check table                    :', rt_chkout
            write(fnr,'(a,1x,f0.2)') 'total time                           :', rt_all
            write(fnr,'(a,1x,f0.2)') '% accounted for                      :', 100. * (rt_prep + rt_batch + rt_pred + rt_table + &
                &rt_write + rt_chkout) / max(rt_all, real(1.e-9,timer_int_kind))
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** PREDICTION STEPS (thread-seconds, ms per particle) ***'
            do k = 1,TH_NSTEPS
                write(fnr,'(a32,a,1x,f0.3,1x,f0.4)') STEP_NAMES(k), ':', sum(th_s(k,:)), &
                    &1000.d0 * sum(th_s(k,:)) / real(max(1,merge(nchk, nptcls, k == 5)),dp)
            enddo
            th_tot = sum(th_s, dim=1)
            busy   = sum(th_tot)
            write(fnr,'(a32,a,1x,f0.3)') 'all steps', ':', busy
            if( rt_pred > 0. ) write(fnr,'(a32,a,1x,f0.1)') 'team utilisation (%)', ':', &
                &100.d0 * busy / (real(nthr,dp) * real(rt_pred,dp))
            if( busy > 0.d0 ) write(fnr,'(a32,a,1x,f0.3)') 'load imbalance (max/mean thread)', ':', &
                &maxval(th_tot) / (busy / real(nthr,dp))
            write(fnr,'(a32,a,1x,i0,1x,i0)') 'particles per thread (min, max)', ':', minval(th_n), maxval(th_n)
            call fclose(fnr)
        end subroutine write_predict_bench

        !> one particle's features, forward pass, folded cone, the readout's pose (direction, mirror and in-plane angle,
        !! no search) and error against its stored direction (the smaller of a direction's and its mirror partner's
        !! angle), on the calling thread; the shift check on the first particles
        subroutine predict_one( i_here, iptcl_here, islot )
            integer, intent(in) :: i_here, iptcl_here, islot
            type(polar_feat) :: feat_here
            type(ori)        :: o_here
            real    :: euls_ref(3), d, cxy(3)
            integer :: ithr_here, j, ip_ref, ip_pose, irot
            integer(int64) :: c0
            ithr_here = 1
            !$ ithr_here = omp_get_thread_num() + 1
            call system_clock(count=c0)
            call feats%build(build%pftc, iptcl_here, feat_here)
            call lap(1, ithr_here, c0)
            call student%predict_feat(ithr_here, feat_here, prob(:,islot))
            call lap(2, ithr_here, c0)
            call student%cone_of(prob(:,islot), preds(i_here))
            call student%pose_of(ithr_here, preds(i_here)%mode, euls_pose(:,i_here), ip_pose)   ! before the shift check reuses the readout
            call lap(3, ithr_here, c0)
            ! the shift: the stored one, plus, given the map, the increment registered against the reference at the read-out
            ! pose (its direction and rotation fixed; the search's convention e3 = 360 - rot(irot))
            sh_reg(:,i_here) = build%spproj_field%get_2Dshift(iptcl_here)
            if( l_reg )then
                irot = build%pftc%get_roind(360. - euls_pose(3,i_here))
                call shsrch(ithr_here)%set_indices(ip_pose, iptcl_here)
                cxy = shsrch(ithr_here)%minimize(irot)
                if( irot > 0 )then
                    sh_reg(:,i_here) = sh_reg(:,i_here) + cxy(2:3)
                    corr_reg(i_here) = cxy(1)
                endif
                call lap(6, ithr_here, c0)
            endif
            if( i_here <= nchk )then
                call shift_check(i_here, ithr_here, islot, feat_here)
                call lap(5, ithr_here, c0)
            endif
            call build%spproj_field%get_ori(iptcl_here, o_here)
            euls_ref          = o_here%get_euler()
            iproj_ref(i_here) = student%grid_closest(o_here)
            ip_ref            = student%fold_index(iproj_ref(i_here))
            dist_mode(i_here) = min(student%angle_between(preds(i_here)%mode, euls_ref), &
                &student%angle_between(student%mirror_of(preds(i_here)%mode), euls_ref))
            rank_ref(i_here)  = 1 + count(prob(:,islot) > prob(ip_ref,islot))   ! the pair's mass sits on its unmirrored direction
            in_cone(i_here)   = 0
            dist_cone(i_here) = huge(d)
            do j = 1,preds(i_here)%ncone
                if( student%fold_index(preds(i_here)%cone(j)) == ip_ref ) in_cone(i_here) = 1
                d = min(student%angle_between(preds(i_here)%cone(j), euls_ref), &
                    &student%angle_between(student%mirror_of(preds(i_here)%cone(j)), euls_ref))
                dist_cone(i_here) = min(dist_cone(i_here), d)
            enddo
            call feat_here%kill
            call o_here%kill
            call lap(4, ithr_here, c0)
            th_n(ithr_here) = th_n(ithr_here) + 1
        end subroutine predict_one

        !> the time since c0 to step k of thread ithr_here, c0 moved on (system_clock: tic/toc share module state)
        subroutine lap( k, ithr_here, c0 )
            integer,        intent(in)    :: k, ithr_here
            integer(int64), intent(inout) :: c0
            integer(int64) :: c1
            call system_clock(count=c1)
            th_s(k,ithr_here) = th_s(k,ithr_here) + real(c1 - c0,dp) / real(clock_rate,dp)
            c0 = c1
        end subroutine lap

        !> ruling R2's shift check of one particle: its tensor shifted through the calculator's phase ramps by each
        !! magnitude in each direction, predicted again; per magnitude the largest total variation of the posterior,
        !! how many directions keep the folded mode and keep the unshifted mode inside the shifted cone, the mean
        !! Jaccard index of the cones' folded sets, and the relative change of each pooled feature group
        subroutine shift_check( i_here, ithr_here, islot, feat_here )
            integer,          intent(in) :: i_here, ithr_here, islot
            type(polar_feat), intent(in) :: feat_here
            type(student_pred) :: pred_sh
            real    :: x_sh(size(feat_here%x,1), size(feat_here%x,2), 2), prob_sh(size(prob,1)), f0(nfeat), f1(nfeat)
            real    :: shift(2), ang, tv, num, den
            integer :: im, id, g, lo, hi, nint_set, nuni, j, k, ip0
            integer :: set0(KMAX_STUDENT), set1(KMAX_STUDENT), n0, n1
            logical :: l_in
            call student%get_features(ithr_here, f0)
            n0 = 0
            do j = 1,preds(i_here)%ncone
                n0 = n0 + 1
                set0(n0) = student%fold_index(preds(i_here)%cone(j))
            enddo
            ip0 = student%fold_index(preds(i_here)%mode)
            do im = 1,size(POSE_SHIFT_CHECK_MAGS)
                chk_tv(i_here,im)   = 0.
                chk_jacc(i_here,im) = 0.
                chk_mode(i_here,im) = 0
                chk_cone(i_here,im) = 0
                chk_fchg(i_here,im,:) = 0.
                do id = 1,POSE_SHIFT_CHECK_NDIR
                    ang   = real(id - 1) * 360. / real(POSE_SHIFT_CHECK_NDIR) * (PI / 180.)
                    shift = POSE_SHIFT_CHECK_MAGS(im) * [cos(ang), sin(ang)]
                    call feats%shift_feat(feat_here%x, shift, x_sh)
                    call student%predict_tensor(ithr_here, x_sh, feat_here%ctf2, feat_here%logpow, prob_sh)
                    call student%get_features(ithr_here, f1)
                    call student%cone_of(prob_sh, pred_sh)
                    tv = 0.5 * sum(abs(prob_sh - prob(:,islot)))
                    chk_tv(i_here,im) = max(chk_tv(i_here,im), tv)
                    if( student%fold_index(pred_sh%mode) == ip0 ) chk_mode(i_here,im) = chk_mode(i_here,im) + 1
                    n1   = 0
                    l_in = .false.
                    do j = 1,pred_sh%ncone
                        n1 = n1 + 1
                        set1(n1) = student%fold_index(pred_sh%cone(j))
                        if( set1(n1) == ip0 ) l_in = .true.
                    enddo
                    if( l_in ) chk_cone(i_here,im) = chk_cone(i_here,im) + 1
                    nint_set = 0
                    do j = 1,n0
                        do k = 1,n1
                            if( set0(j) == set1(k) ) nint_set = nint_set + 1
                        enddo
                    enddo
                    nuni = n0 + n1 - nint_set
                    if( nuni > 0 ) chk_jacc(i_here,im) = chk_jacc(i_here,im) + real(nint_set) / real(nuni)
                    do g = 1,POSE_SHIFT_CHECK_NGROUPS
                        lo  = 1
                        hi  = nfeat
                        num = sum(abs(f1(lo:hi) - f0(lo:hi)))
                        den = sum(abs(f0(lo:hi)))
                        if( den > 0. ) chk_fchg(i_here,im,g) = chk_fchg(i_here,im,g) + num / den
                    enddo
                enddo
                chk_jacc(i_here,im)   = chk_jacc(i_here,im)   / real(POSE_SHIFT_CHECK_NDIR)
                chk_fchg(i_here,im,:) = chk_fchg(i_here,im,:) / real(POSE_SHIFT_CHECK_NDIR)
            enddo
        end subroutine shift_check

        !> the shift check's table and summary: per shift magnitude over the checked particles
        subroutine write_shift_check
            integer :: im, io_here, g
            real    :: tv_mean, tv_q90, mode_frac, cone_frac, jacc_mean, fchg(POSE_SHIFT_CHECK_NGROUPS)
            if( nchk < 1 ) return
            call fopen(chk_unit, fname_chk, status='replace', action='write', iostat=io_here)
            call fileiochk('predict_poses; open '//fname_chk%to_char(), io_here)
            write(chk_unit,'(A,I0,A,I0,A)') '# ruling R2 shift check: ', nchk, ' particles shifted through the calculator''s phase ramps in ', &
                &POSE_SHIFT_CHECK_NDIR, ' directions; per magnitude the posterior''s largest total variation (mean, 90th percentile), the share of'
            write(chk_unit,'(A)') '# shifted copies keeping the folded mode and holding the unshifted mode in their cone, the cones'' mean Jaccard index,'
            write(chk_unit,'(A)') '# and the relative change of the readout map (ruling R3''s template readout)'
            write(chk_unit,'(A)') '# shift_px tv_mean tv_q90 mode_kept cone_holds_mode jaccard fchg_readout_map'
            allocate(chk_sorted(nchk))
            do im = 1,size(POSE_SHIFT_CHECK_MAGS)
                chk_sorted = chk_tv(:,im)
                call hpsort(chk_sorted)
                tv_mean   = sum(chk_tv(:,im)) / real(nchk)
                tv_q90    = chk_sorted(max(1, nint(0.9 * real(nchk))))
                mode_frac = real(sum(chk_mode(:,im))) / real(nchk * POSE_SHIFT_CHECK_NDIR)
                cone_frac = real(sum(chk_cone(:,im))) / real(nchk * POSE_SHIFT_CHECK_NDIR)
                jacc_mean = sum(chk_jacc(:,im)) / real(nchk)
                do g = 1,POSE_SHIFT_CHECK_NGROUPS
                    fchg(g) = sum(chk_fchg(:,im,g)) / real(nchk)
                enddo
                write(chk_unit,'(F6.2,1X,F8.5,1X,F8.5,1X,F8.5,1X,F8.5,1X,F8.5,1X,F9.5)') POSE_SHIFT_CHECK_MAGS(im), &
                    &tv_mean, tv_q90, mode_frac, cone_frac, jacc_mean, fchg(1)
                write(logfhandle,'(A,F5.1,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F8.5)') &
                    &'>>> POSE STUDENT SHIFT CHECK: ', POSE_SHIFT_CHECK_MAGS(im), ' px: posterior TV mean ', tv_mean, ' q90 ', tv_q90, &
                    &', folded mode kept ', mode_frac, ', cone holds the mode ', cone_frac, ', cones'' Jaccard ', jacc_mean, &
                    &'; readout map change ', fchg(1)
            enddo
            call fclose(chk_unit)
            deallocate(chk_sorted)
        end subroutine write_shift_check

    end subroutine predict_poses_exec

end module simple_strategy_pose_student
