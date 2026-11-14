!@descr: the class student's strategy routines: the teacher's rows of the selected particles against frozen class averages and their features, the student's training (learn_class_student), and its prediction of the class, the in-plane angle, the candidates and the confidence (predict_classes)
! The teacher is the 2D probabilistic table's posterior against the class averages refs, prepared as refine2D's
! matcher prepares them, batch by batch, each row read before the batch's assignment releases it
! (doc/algorithms/class_student.md).
module simple_strategy_class_student
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_builder,                 only: builder
use simple_cmdline,                 only: cmdline
use simple_parameters,              only: parameters
use simple_image,                   only: image
use simple_eo_partition,            only: ensure_field_halves
use simple_classaverager,           only: cavger_new, cavger_read_all, cavger_kill
use simple_eul_prob_tab2D,          only: eul_prob_tab2D
use simple_eul_prob_tab_utils,      only: prob_row_record, prob_row_file
use simple_matcher_ptcl_batch,      only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D, check_pass_team
use simple_matcher_ptcl_io,         only: ptcl_read_session
use simple_matcher_pftc_prep,       only: prep_pftc4align2D
use simple_matcher_smpl_and_lplims, only: set_bp_range2D, adopt_run_seed
use simple_sigma2_bootstrap,        only: ensure_sigma2_for_iteration
use simple_qsys_funs,               only: qsys_declare_part_finished
use simple_pftc_shsrch_grad,        only: pftc_shsrch_grad
use simple_polar_feats,             only: polar_feats, polar_feat, polar_feats_hash_ints
use simple_class_student,           only: class_student, class_pred
use simple_strategy_student_common, only: student_holdout_split, student_conf_threshold, student_geometry_mismatch, &
    &file_hash, choose_student_design, fname_base, student_part_pinds, student_predict_parts, student_sigma2_for_pass
use simple_defs_hhcnn,              only: KMAX_CLASS, STUDENT_CONF_FRACTION, STUDENT_KMAX_ROW, STUDENT_WIDTHS, &
    &STUDENT_WIDEST_MIN_ROWS, STUDENT_MAX_EPOCHS, STUDENT_PATIENCE_EPOCHS, STUDENT_SHIFT_AUG_NSIG, &
    &STUDENT_SHIFT_REG_HALFWIDTH, STUDENT_SHIFT_REG_MAXITS, CLASS_STUDENT_TEACHER_ITER, CLASS_STUDENT_NDESIGNS
use simple_rnd,                     only: rnd_run_seed
use simple_type_defs,               only: OBJFUN_EUCLID, PTCL2D_SEG
use simple_binoris_io,              only: binwrite_oritab
use simple_class_candidates,        only: class_candidates
implicit none

public :: run_learn_class_student, learn_class_student_exec, run_predict_classes, predict_classes_exec
private
#include "simple_local_flags.inc"


contains

    ! ===== LEARN_CLASS_STUDENT =====

    !> learn_class_student on its own parameters and builder, built from cline
    subroutine run_learn_class_student( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target :: params
        type(builder)            :: build
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        call learn_class_student_exec(params, build, cline)
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_LEARN_CLASS_STUDENT NORMAL STOP ****', print_simple=.false.)
    end subroutine run_learn_class_student

    !> The dense 2D table's rows of the selected particles against the class averages refs at the band lp (refine2D's
    !! preparation and fill), written beside fname, and the student trained on them and written as fname. Project:
    !! reads the particle field; writes what refine2D's setup writes (the halves, the sigma2 bootstrap)
    subroutine learn_class_student_exec( params, build, cline )
        type(parameters), target, intent(inout) :: params
        type(builder),    target, intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        type(eul_prob_tab2D)     :: tab
        type(ptcl_read_session)  :: reads
        type(image), allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:), ptcl_imgs(:)
        type(polar_feats)        :: feats
        type(polar_feat)         :: feat
        type(prob_row_file)      :: rows
        type(prob_row_record)    :: rec
        type(class_student)      :: student, reread
        type(class_pred)         :: pred
        type(string)             :: base, fname_rows, fname_stats, fname_header, projpath
        real,    allocatable     :: dists(:), wk(:,:), x_all(:,:,:,:), ctf2_all(:,:), side_all(:), dists_all(:,:)
        real,    allocatable     :: temps_all(:), e3_all(:), xcavgs(:,:,:,:), ctf_one(:,:), prob(:,:)
        real,    allocatable     :: conf_held(:), e3_held(:)
        complex, allocatable     :: pft_ref(:,:)
        logical, allocatable     :: held(:), recov(:,:)
        integer, allocatable     :: pinds(:), batches(:,:), mode_all(:), held_inds(:), group_held(:), group_of(:), stkinds(:)
        integer, allocatable     :: k99_held(:), ncand_held(:), mode_held_s(:), states(:)
        integer(int64)           :: cavgs_hash
        integer(timer_int_kind)  :: t0, t_all
        real(timer_int_kind)     :: rt_prep, rt_batch, rt_fill, rt_rows, rt_train, rt_all
        real(timer_int_kind)     :: rt_cand(size(STUDENT_WIDTHS)), rt_boot(size(STUDENT_WIDTHS))
        real     :: recov_cand(size(STUDENT_WIDTHS)), loss_cand(size(STUDENT_WIDTHS)), gain_q(size(STUDENT_WIDTHS))
        integer  :: epoch_cand(size(STUDENT_WIDTHS))
        real     :: shift_aug, conf_threshold, heldout_seeded, heldout_recovered
        real     :: loss_train, loss_held, mode_held, time_train_s
        integer  :: nptcls, batchsz_max, nbatches, ibatch, b0, b1, batchsz, i, iptcl, nk, nheld, j, icls, ntry
        integer  :: stat_unit, head_unit, io_stat, nthr_here, iw_best, ntried, nsel, ncls
        integer  :: nptcls_train, nheldout, epochs_run, best_epoch, nweights
        logical  :: l_rebuilt, l_by_stack
        character(len=:), allocatable :: msg
        t_all = tic()
        call check_pass_team('learn_class_student')
        ! the contract
        if( params%cc_objfun /= OBJFUN_EUCLID ) THROW_HARD('the teacher is the Euclidean posterior: objfun=euclid; learn_class_student')
        if( .not. params%l_lpset ) THROW_HARD('the band is lp; learn_class_student')
        if( trim(params%oritype) /= 'ptcl2D' ) THROW_HARD('the rows are ptcl2D''s; learn_class_student')
        if( trim(params%refine) /= 'prob' ) THROW_HARD('the teacher is the dense table: refine=prob; learn_class_student')
        if( .not. cline%defined('refs') ) THROW_HARD('the references are the class averages refs; learn_class_student')
        ncls         = params%ncls
        base         = fname_base(params%fname)
        fname_rows   = base%to_char()//'_rows.bin'
        fname_stats  = base%to_char()//'_rows.txt'
        fname_header = base%to_char()//'_rows_header.txt'
        ! the selected particles and the held-out micrographs
        call build%spproj_field%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no selected particle; learn_class_student')
        nptcls = size(pinds)
        if( nptcls < 1 ) THROW_HARD('no selected particle; learn_class_student')
        write(logfhandle,'(A,I0,A,I0,A,F7.2,A)') '>>> CLASS STUDENT ROWS: ', nptcls, ' particles against ', ncls, &
            &' classes at lp ', params%lp, ' A'
        allocate(held(nptcls), source=.false.)
        allocate(group_of(nptcls), stkinds(nptcls), source=0)
        do i = 1,nptcls
            stkinds(i) = build%spproj_field%get_int(pinds(i), 'stkind')
        enddo
        call student_holdout_split(nptcls, stkinds, held, group_of, l_by_stack, msg)
        if( len(msg) > 0 ) THROW_HARD(msg)
        deallocate(stkinds)
        ! refine2D's setup: the halves, the sigma2 state, the run seed of the keyed draws
        t0 = tic()
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, CLASS_STUDENT_TEACHER_ITER, params%box, params%smpd, &
            &'CLASS STUDENT SIGMA2', l_rebuilt)
        if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
        call adopt_run_seed(params, build)
        params%which_iter = CLASS_STUDENT_TEACHER_ITER
        ! the band, the batches and the references, as refine2D_exec prepares them
        call cavger_new(params, build)
        call cavger_read_all
        call set_bp_range2D(params, build, cline, CLASS_STUDENT_TEACHER_ITER, 0.)
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = 1 + (nptcls - 1) / batchsz_max
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, CLASS_STUDENT_TEACHER_ITER)
        call cavger_kill
        ! the features' geometry and the frozen weights: the calculator's band and the spectra the rows are whitened with
        call build%pftc%get_sigma2_weights(wk)
        call feats%new(build%pftc%get_pftsz(), params%kfromto, params%box, wk)
        nk = params%kfromto(2) - params%kfromto(1) + 1
        ! the class averages through the features (no CTF): the readout's templates start from them
        allocate(xcavgs(feats%get_pftsz(), nk, 2, ncls), source=0.)
        allocate(pft_ref(feats%get_pftsz_in(), params%kfromto(1):params%kfromto(2)))
        allocate(ctf_one(feats%get_pftsz_in(), params%kfromto(1):params%kfromto(2)), source=1.)
        do icls = 1,ncls
            call build%pftc%get_ref_pft(icls, .true., pft_ref)
            call feats%build_from_pft(pft_ref, ctf_one, .true., feat)
            xcavgs(:,:,:,icls) = feat%x
        enddo
        deallocate(pft_ref, ctf_one)
        ! the classes' states in the teacher's project (all selected without a class document of this size)
        allocate(states(ncls), source=1)
        if( build%spproj%os_cls2D%get_noris() == ncls )then
            do icls = 1,ncls
                states(icls) = build%spproj%os_cls2D%get_state(icls)
            enddo
        endif
        projpath   = simple_abspath(params%projfile)
        cavgs_hash = file_hash(params%refs)
        ! the table of the pass: the iteration's model from the complete field, the assignment opened
        call tab%new_pipeline(params, build, pinds)
        call tab%set_model_from_field
        call tab%open_assignment
        call rows%open_write(fname_rows, nptcls, ncls)
        call fopen(stat_unit, fname_stats, status='replace', action='write', iostat=io_stat)
        call fileiochk('learn_class_student; open '//fname_stats%to_char(), io_stat)
        write(stat_unit,'(A)') '# iptcl eo temp nvalid mode pmax neff k99 nkept wtail lmin inpl_mode e3_mode x_mode y_mode '//&
            &'class_stored wsq logpow'
        ! the training set in memory: the tensors, the rows' losses and temperatures
        allocate(x_all(feats%get_pftsz(), nk, 2, nptcls), ctf2_all(nk, nptcls), side_all(nptcls), source=0.)
        allocate(dists_all(ncls, nptcls), temps_all(nptcls), source=0.)
        allocate(e3_all(nptcls), source=-1.)
        allocate(mode_all(nptcls), source=0)
        allocate(dists(ncls))
        rt_prep  = toc(t0)
        rt_batch = 0.
        rt_fill  = 0.
        rt_rows  = 0.
        call reads%new(params)
        call reads%expect(pinds)
        do ibatch = 1,nbatches
            b0      = batches(ibatch,1)
            b1      = batches(ibatch,2)
            batchsz = b1 - b0 + 1
            t0 = tic()
            call build_batch_particles2D(params, build, batchsz, pinds(b0:b1), reads, ptcl_match_imgs=ptcl_match_imgs, &
                &ptcl_match_imgs_pad=ptcl_match_imgs_pad)
            rt_batch = rt_batch + toc(t0)
            t0 = tic()
            call tab%fill_tab_range(b0, b1)
            rt_fill = rt_fill + toc(t0)
            t0 = tic()
            call tab%keep_batch_rows(b0, b1)
            do i = b0,b1
                iptcl = pinds(i)
                call tab%get_posterior_row(i, params%prob_beta, STUDENT_KMAX_ROW, rec, dists)
                call rows%put(rec, dists)
                call feats%build(build%pftc, iptcl, feat)
                x_all(:,:,:,i) = feat%x
                ctf2_all(:,i)  = feat%ctf2
                side_all(i)    = feat%logpow
                dists_all(:,i) = dists
                temps_all(i)   = rec%temp
                mode_all(i)    = rec%imin
                ! the search's convention: e3 = 360 - rot(irot)
                if( rec%inpl_mode > 0 ) e3_all(i) = 360. - build%pftc%get_rot(rec%inpl_mode)
                write(stat_unit,'(I0,1X,I0,1X,ES13.6,1X,I0,1X,I0,1X,F8.5,1X,F10.3,1X,I0,1X,I0,1X,ES11.4,1X,ES13.6,1X,I0,'//&
                    &'1X,F8.3,1X,F8.3,1X,F8.3,1X,I0,1X,ES13.6,1X,F9.4)') iptcl, build%spproj_field%get_eo(iptcl), rec%temp, &
                    &rec%nvalid, rec%imin, rec%pmax, rec%neff, rec%k99, rec%nkept, rec%wtail, rec%lmin, rec%inpl_mode, &
                    &e3_all(i), rec%x_mode, rec%y_mode, build%spproj_field%get_class(iptcl), feat%wsq, feat%logpow
            enddo
            rt_rows = rt_rows + toc(t0)
            call tab%assign_batch(b0, b1)
        enddo
        call tab%close_assignment
        call rows%close
        call fclose(stat_unit)
        call reads%kill
        call feat%kill
        call rec%kill
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        call tab%kill
        call build%pftc%kill
        call build%esig%kill
        write(logfhandle,'(A,F9.2,A,F9.2,A)') '>>> CLASS STUDENT ROWS: fill ', rt_fill, ' s wall, ', &
            &rt_fill / real(nptcls,timer_int_kind) * 1000., ' ms per particle'
        ! the students, cheapest first (the widest only on a large set), trained with the centring disk's shift
        ! augmentation; a wider design replaces the kept one when it recovers clearly more held-out particles (the
        ! teacher's class among the candidates of a particle whose candidates fit the cap)
        t0 = tic()
        shift_aug = 0.
        if( params%l_shift_map ) shift_aug = STUDENT_SHIFT_AUG_NSIG * params%shift_sigma / params%smpd
        nheld = count(held)
        allocate(held_inds(nheld), group_held(nheld), k99_held(nheld), ncand_held(nheld), mode_held_s(nheld), source=0)
        allocate(conf_held(nheld), e3_held(nheld), source=0.)
        allocate(recov(nheld, size(STUDENT_WIDTHS)), source=.false.)
        j = 0
        do i = 1,nptcls
            if( .not. held(i) ) cycle
            j = j + 1
            held_inds(j)  = i
            group_held(j) = group_of(i)
        enddo
        nthr_here = params%nthr
        !$ nthr_here = omp_get_max_threads()
        allocate(prob(ncls, nthr_here), source=0.)
        recov_cand = 0.
        loss_cand  = 0.
        gain_q     = 0.
        rt_cand    = 0.
        rt_boot    = 0.
        epoch_cand = 0
        ntry = min(CLASS_STUDENT_NDESIGNS, size(STUDENT_WIDTHS))
        if( nptcls - nheld < STUDENT_WIDEST_MIN_ROWS ) ntry = min(ntry, size(STUDENT_WIDTHS) - 1)
        call choose_student_design('CLASS STUDENT', STUDENT_WIDTHS, ntry, nheld, group_held, rnd_run_seed(), try_design, &
            &keep_design, recov, iw_best, ntried, gain_q, rt_boot)
        ! the kept design back from its file, its threshold on the held-out rows and the held-out table
        call student%kill
        call student%read(params%fname)
        call predict_held(recov(:,iw_best))
        call set_conf_threshold
        call student%set_threshold(conf_threshold, STUDENT_CONF_FRACTION, heldout_seeded, heldout_recovered)
        call student%write(params%fname)
        call check_reread_population
        call write_heldout_table
        call student%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
            &time_train_s, nweights)
        rt_train = toc(t0)
        write(logfhandle,'(A,I0,A,F8.5,A,F8.4,A,F8.4,A,F9.1,A)') '>>> CLASS STUDENT KEPT: width ', STUDENT_WIDTHS(iw_best), &
            &', confidence threshold ', conf_threshold, ' seeding ', heldout_seeded, ' of the held-out rows with ', &
            &heldout_recovered, ' of their teacher classes among the candidates; the designs ', rt_train, ' s'
        call log_inpl_error
        rt_all = toc(t_all)
        call write_header
        call student%kill
        call feats%kill
        deallocate(dists, pinds, batches, wk, x_all, ctf2_all, side_all, dists_all, temps_all, e3_all, held, mode_all, &
            &group_of, xcavgs, states, held_inds, group_held, conf_held, e3_held, k99_held, ncand_held, mode_held_s, recov, prob)

    contains

        !> design iw (choose_student_design's trial): built with the class averages as its templates, trained on the
        !! rows, and judged on the held-out rows, recov_iw the particles it recovers
        subroutine try_design( iw_here, recov_iw )
            integer, intent(in)    :: iw_here
            logical, intent(inout) :: recov_iw(:)
            call student%new(feats%get_pftsz_in(), params%kfromto, params%box, params%box_crop, params%smpd, &
                &params%mskdiam, params%lp, wk, trim(params%ctf) /= 'no', ncls, states, STUDENT_WIDTHS(iw_here), rnd_run_seed())
            call student%init_readout(xcavgs)
            if( l_by_stack )then
                call student%set_provenance(projpath%to_char(), polar_feats_hash_ints(pinds), cavgs_hash, 'every_10th_selected_stack')
            else
                call student%set_provenance(projpath%to_char(), polar_feats_hash_ints(pinds), cavgs_hash, 'every_10th_selected_particle')
            endif
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.2,A)') '>>> CLASS STUDENT TRAINING: width ', STUDENT_WIDTHS(iw_here), ', ', &
                &nptcls - nheld, ' rows, ', nheld, ' held out, shifts augmented within ', shift_aug, ' px'
            call student%train(nptcls, x_all, ctf2_all, side_all, dists_all, temps_all, held, params%prob_beta, &
                &shift_aug, STUDENT_MAX_EPOCHS, STUDENT_PATIENCE_EPOCHS, params%nthr, logfhandle)
            call student%get_record(nptcls_train, nheldout, epochs_run, best_epoch, loss_train, loss_held, mode_held, &
                &time_train_s, nweights)
            rt_cand(iw_here)    = time_train_s
            loss_cand(iw_here)  = loss_held
            epoch_cand(iw_here) = best_epoch
            call predict_held(recov_iw)
            recov_cand(iw_here) = real(count(recov_iw(1:nheld))) / real(nheld)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F9.5,A,F7.4,A,F7.4,A,F9.1,A)') '>>> CLASS STUDENT TRAINED: width ', &
                &STUDENT_WIDTHS(iw_here), ', best epoch ', best_epoch, ' of ', epochs_run, ', held-out loss ', loss_held, &
                &', held-out mode agreement ', mode_held, ', held-out classes recovered ', recov_cand(iw_here), ', ', &
                &time_train_s, ' s'
        end subroutine try_design

        !> design iw kept (choose_student_design): its model written
        subroutine keep_design( iw_here )
            integer, intent(in) :: iw_here
            call student%write(params%fname)
        end subroutine keep_design

        !> the student's prediction of every held-out row: the candidates, the confidence, the mode and its in-plane angle;
        !! recov_here whether the teacher's class lies among the candidates of a particle whose candidates fit the cap
        subroutine predict_held( recov_here )
            logical, intent(inout) :: recov_here(:)
            integer :: jj, ii, ithr_loc, c
            call student%prepare_predict(nthr_here)
            !$omp parallel do default(shared) private(jj,ii,ithr_loc,pred,c) schedule(dynamic) proc_bind(close)
            do jj = 1,nheld
                ii       = held_inds(jj)
                ithr_loc = 1
                !$ ithr_loc = omp_get_thread_num() + 1
                call student%predict_tensor(ithr_loc, x_all(:,:,:,ii), ctf2_all(:,ii), side_all(ii), prob(:,ithr_loc))
                call student%cands_of(prob(:,ithr_loc), pred)
                conf_held(jj)   = pred%conf
                k99_held(jj)    = pred%k99
                ncand_held(jj)  = pred%ncand
                mode_held_s(jj) = pred%mode
                e3_held(jj)     = student%rotation_of(ithr_loc, pred%mode)
                recov_here(jj)  = .false.
                if( pred%k99 <= KMAX_CLASS )then
                    do c = 1,pred%ncand
                        if( pred%cands(c) == mode_all(ii) ) recov_here(jj) = .true.
                    enddo
                endif
            enddo
            !$omp end parallel do
        end subroutine predict_held

        !> the confidence threshold on the held-out rows whose candidates fit the cap (student_conf_threshold), the seeded
        !! share of the held-out rows and the recovered share among them
        subroutine set_conf_threshold
            integer :: nrec
            call student_conf_threshold(nheld, conf_held, k99_held <= KMAX_CLASS, recov(:,iw_best), STUDENT_CONF_FRACTION, &
                &conf_threshold, nsel, nrec)
            heldout_seeded    = real(nsel) / real(nheld)
            heldout_recovered = 0.
            if( nsel > 0 ) heldout_recovered = real(nrec) / real(nsel)
        end subroutine set_conf_threshold

        !> the model read back from its file seeds exactly the calibrated held-out rows
        subroutine check_reread_population
            real :: thr, cf, hs, hr
            call reread%read(params%fname)
            call reread%get_threshold(thr, cf, hs, hr)
            call reread%kill
            if( count(k99_held <= KMAX_CLASS .and. conf_held >= thr) /= nsel ) &
                &THROW_HARD('the model read back seeds other held-out rows than the calibration; learn_class_student')
        end subroutine check_reread_population

        !> the held-out table of the kept design: per held-out row its micrograph, confidence, 99% set and candidates,
        !! the teacher's and the student's classes and in-plane angles, and whether it is recovered
        subroutine write_heldout_table
            type(string) :: fname_held
            integer :: jj, ii, held_unit
            fname_held = base%to_char()//'_heldout.txt'
            call fopen(held_unit, fname_held, status='replace', action='write', iostat=io_stat)
            call fileiochk('learn_class_student; open '//fname_held%to_char(), io_stat)
            write(held_unit,'(A)') '# iptcl group conf k99 ncand teacher_class student_class recovered teacher_e3 student_e3'
            do jj = 1,nheld
                ii = held_inds(jj)
                write(held_unit,'(I0,1X,I0,1X,F8.5,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,F8.3,1X,F8.3)') pinds(ii), group_held(jj), &
                    &conf_held(jj), k99_held(jj), ncand_held(jj), mode_all(ii), mode_held_s(jj), &
                    &merge(1, 0, recov(jj,iw_best)), e3_all(ii), e3_held(jj)
            enddo
            call fclose(held_unit)
            call fname_held%kill
        end subroutine write_heldout_table

        !> the student's in-plane angle against the teacher's on the held-out rows where their classes agree
        subroutine log_inpl_error
            real, allocatable :: d(:)
            integer :: jj, ii, n
            allocate(d(nheld), source=0.)
            n = 0
            do jj = 1,nheld
                ii = held_inds(jj)
                if( mode_held_s(jj) /= mode_all(ii) .or. e3_all(ii) < 0. ) cycle
                n    = n + 1
                d(n) = abs(e3_held(jj) - e3_all(ii))
                d(n) = min(d(n), 360. - d(n))
            enddo
            if( n > 0 ) write(logfhandle,'(A,I0,A,F8.3,A)') '>>> CLASS STUDENT KEPT: in-plane angle against the teacher''s on ', &
                &n, ' held-out rows of agreeing class, median ', median_nocopy(d(1:n)), ' deg'
            deallocate(d)
        end subroutine log_inpl_error

        !> the rows' provenance and cost and the training's record
        subroutine write_header
            integer :: iw
            call fopen(head_unit, fname_header, status='replace', action='write', iostat=io_stat)
            call fileiochk('learn_class_student; open '//fname_header%to_char(), io_stat)
            write(head_unit,'(A)')        '# class student teacher rows'
            write(head_unit,'(A)')        'rows_version=1'
            write(head_unit,'(A,A)')      'projfile=', projpath%to_char()
            write(head_unit,'(A,I0)')     'project_rows=', build%spproj_field%get_noris()
            write(head_unit,'(A,I0)')     'nptcls=', nptcls
            write(head_unit,'(A,I0)')     'selection_hash=', polar_feats_hash_ints(pinds)
            write(head_unit,'(A,I0)')     'ncls=', ncls
            write(head_unit,'(A,A)')      'refs=', params%refs%to_char()
            write(head_unit,'(A,I0)')     'cavgs_hash=', cavgs_hash
            write(head_unit,'(A,F10.4)')  'lp=', params%lp
            write(head_unit,'(A,I0,A,I0)') 'kfromto=', params%kfromto(1), ',', params%kfromto(2)
            write(head_unit,'(A,I0)')     'pftsz=', feats%get_pftsz_in()
            write(head_unit,'(A,I0)')     'box=', params%box
            write(head_unit,'(A,I0)')     'box_crop=', params%box_crop
            write(head_unit,'(A,F10.5)')  'smpd=', params%smpd
            write(head_unit,'(A,F10.3)')  'mskdiam=', params%mskdiam
            write(head_unit,'(A,A)')      'ctf=', trim(params%ctf)
            write(head_unit,'(A,F8.4)')   'prob_beta=', params%prob_beta
            write(head_unit,'(A,F8.3)')   'trs=', params%trs
            write(head_unit,'(A,F8.3)')   'shift_sigma=', params%shift_sigma
            write(head_unit,'(A,I0)')     'kmax_row=', STUDENT_KMAX_ROW
            write(head_unit,'(A,I0)')     'kmax_class=', KMAX_CLASS
            write(head_unit,'(A,I0)')     'nthr=', params%nthr
            write(head_unit,'(A,F10.3)')  'time_prep_s=', rt_prep
            write(head_unit,'(A,F10.3)')  'time_batch_prep_s=', rt_batch
            write(head_unit,'(A,F10.3)')  'time_fill_s=', rt_fill
            write(head_unit,'(A,F10.3)')  'time_rows_feats_s=', rt_rows
            write(head_unit,'(A,F10.3)')  'time_train_s=', rt_train
            write(head_unit,'(A,F10.3)')  'time_total_s=', rt_all
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
            write(head_unit,'(A)')        'recovery=teacher_class_among_candidates_within_cap'
            write(head_unit,'(A,I0)')     'train_rows=', nptcls_train
            write(head_unit,'(A,I0)')     'heldout_rows=', nheldout
            write(head_unit,'(A,F8.3)')   'shift_aug_px=', shift_aug
            write(head_unit,'(A,I0)')     'epochs_run=', epochs_run
            write(head_unit,'(A,I0)')     'best_epoch=', best_epoch
            write(head_unit,'(A,ES14.6)') 'loss_train=', loss_train
            write(head_unit,'(A,ES14.6)') 'loss_heldout=', loss_held
            write(head_unit,'(A,F8.5)')   'mode_heldout=', mode_held
            call fclose(head_unit)
        end subroutine write_header

    end subroutine learn_class_student_exec

    ! ===== PREDICT_CLASSES =====

    !> predict_classes on its own parameters and builder: the model's classes and band are set on the command line from
    !! its header before the parameters are made
    subroutine run_predict_classes( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target :: params
        type(builder)            :: build
        type(class_student)      :: student
        type(string)             :: infile
        infile = cline%get_carg('infile')
        call student%read(infile)
        call cline%set('ncls', student%get_ncls())
        call cline%set('lp',   student%get_lp())
        call student%kill
        call infile%kill
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            if( cline%get_iarg('nparts') > 1 )then
                call predict_classes_master(cline)
                return
            endif
            call cline%delete('nparts')
        endif
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        call predict_classes_exec(params, build, cline)
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        if( params%l_distr_worker ) call qsys_declare_part_finished(params, string('predict_classes'))
        call simple_end('**** SIMPLE_PREDICT_CLASSES NORMAL STOP ****', print_simple=.false.)
    end subroutine run_predict_classes

    !> predict_classes over nparts parts: given refs, the halves and sigma2 made first; the parts' records merged and
    !! written once, their tables joined, their candidates published as one artifact (doc/algorithms/class_student.md)
    subroutine predict_classes_master( cline )
        class(cmdline), intent(inout) :: cline
        type(parameters), target  :: params
        type(builder)             :: build
        type(class_candidates)    :: cands_art
        type(string)              :: base
        type(string), allocatable :: part_bases(:), part_cands(:)
        logical :: l_rebuilt
        integer :: ipart
        call build%init_params_and_build_spproj(cline, params)
        if( cline%defined('refs') )then
            call ensure_field_halves(build%spproj, params%oritype, params%projfile)
            call ensure_sigma2_for_iteration(cline, params%projfile, params%oritype, CLASS_STUDENT_TEACHER_ITER, &
                &params%box, params%smpd, 'CLASS STUDENT SIGMA2', l_rebuilt)
            if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
        endif
        base = fname_base(params%fname)
        call student_predict_parts(params, build, cline, base, part_bases)
        call build%spproj%write_segment_inside(params%oritype)
        allocate(part_cands(params%nparts))
        do ipart = 1,params%nparts
            part_cands(ipart) = part_bases(ipart)%to_char()//'_cands.bin'
        enddo
        call cands_art%publish_parts(build%spproj, build%spproj_field, params%projfile, &
            &string(base%to_char()//'_cands.bin'), part_cands, params%ncls)
        call cands_art%kill
        do ipart = 1,params%nparts
            call del_file(part_cands(ipart))
        enddo
        call build%spproj%kill
        call simple_end('**** SIMPLE_PREDICT_CLASSES NORMAL STOP ****', print_simple=.false.)
    end subroutine predict_classes_master

    !> The student's class, in-plane angle, candidates and confidence for the selected particles, the shift registered
    !! against the predicted class given refs; to <fname>_pred.txt, the records (class, e3, shift, pconf, route) and the
    !! project's class_candidates artifact
    subroutine predict_classes_exec( params, build, cline )
        type(parameters), target, intent(inout) :: params
        type(builder),    target, intent(inout) :: build
        class(cmdline),           intent(inout) :: cline
        type(class_student)      :: student
        type(ptcl_read_session)  :: reads
        type(image), allocatable :: ptcl_match_imgs(:), ptcl_match_imgs_pad(:), ptcl_imgs(:)
        type(polar_feats)        :: feats
        type(string)             :: base, fname_pred
        type(class_pred), allocatable :: preds(:)
        type(pftc_shsrch_grad), allocatable :: shsrch(:)
        real,    allocatable     :: wk(:,:), prob(:,:), e3_pred(:), sh_reg(:,:)
        integer, allocatable     :: pinds(:), batches(:,:), class_stored(:), routes(:), ncands(:), cands_all(:,:)
        type(class_candidates)   :: cands_art
        character(len=:), allocatable :: msg
        character(len=*), parameter :: PRED_HEADER = '# iptcl eo class pmax k99 ncand conf route e3 x y class_stored cands'
        integer(timer_int_kind)  :: t0, t_all
        real(timer_int_kind)     :: rt_prep, rt_batch, rt_pred, rt_write, rt_all
        real     :: conf_threshold, conf_fraction, heldout_seeded, heldout_recovered, lims_reg(2,2)
        integer  :: nptcls, batchsz_max, nbatches, ibatch, b0, b1, batchsz, i, iptcl, nthr, pred_unit, io_stat, kfromto(2)
        integer, allocatable     :: pinds_part(:)
        integer  :: route, n_seeded, n_capped, n_agree, j
        real(dp) :: ncand_sum, conf_sum
        logical  :: l_reg, l_rebuilt, l_part
        t_all = tic()
        call check_pass_team('predict_classes')
        ! a part of a distributed master (part, nparts, fromp, top, outfile: the commander checks they are given)
        ! writes its range's records to outfile and its candidates to a part file; otherwise the project is written
        l_part = cline%defined('part')
        if( trim(params%oritype) /= 'ptcl2D' ) THROW_HARD('the classes are ptcl2D''s; predict_classes')
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
        if( l_part ) base = base%to_char()//'_part'//int2str(params%part)
        fname_pred = base%to_char()//'_pred.txt'
        call build%spproj_field%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no selected particle; predict_classes')
        if( size(pinds) < 1 ) THROW_HARD('no selected particle; predict_classes')
        if( l_part )then
            call student_part_pinds(pinds, params%part, params%nparts, params%fromp, params%top, &
                &build%spproj_field%get_noris(), pinds_part, msg)
            if( len(msg) > 0 ) THROW_HARD(msg)
            call move_alloc(pinds_part, pinds)
        endif
        nptcls = size(pinds)
        if( nptcls < 1 )then
            ! a part whose range holds no selected particle: the table's header, the range's records and candidates
            call write_empty_part
            return
        endif
        call student%get_threshold(conf_threshold, conf_fraction, heldout_seeded, heldout_recovered)
        write(logfhandle,'(A,I0,A,I0,A,F7.2,A)') '>>> CLASS STUDENT PREDICTION: ', nptcls, ' particles over ', &
            &student%get_ncls(), ' classes at lp ', params%lp, ' A'
        ! the batches and the calculator: one reference slot (the features need the particles alone), or, given refs, the
        ! class averages with sigma2, as learn_class_student prepares them, for the shift registration
        l_reg       = cline%defined('refs')
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = 1 + (nptcls - 1) / batchsz_max
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        if( l_reg )then
            call ensure_field_halves(build%spproj, params%oritype, params%projfile)
            call student_sigma2_for_pass(cline, params%projfile, params%oritype, CLASS_STUDENT_TEACHER_ITER, params%box, &
                &params%smpd, 'CLASS STUDENT SIGMA2', l_rebuilt)
            if( l_rebuilt ) call build%spproj%read_segment('projinfo', params%projfile)
            params%which_iter = CLASS_STUDENT_TEACHER_ITER
            call cavger_new(params, build)
            call cavger_read_all
            call set_bp_range2D(params, build, cline, CLASS_STUDENT_TEACHER_ITER, 0.)
            if( any(params%kfromto /= kfromto) ) THROW_HARD('the class averages'' band differs from the model''s; predict_classes')
            call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, CLASS_STUDENT_TEACHER_ITER)
            call cavger_kill
        else
            call build%pftc%new(params, 1, [1,batchsz_max], kfromto)
        endif
        nthr = params%nthr
        !$ nthr = omp_get_max_threads()
        if( l_reg )then
            lims_reg(:,1) = -STUDENT_SHIFT_REG_HALFWIDTH / params%smpd
            lims_reg(:,2) =  STUDENT_SHIFT_REG_HALFWIDTH / params%smpd
            allocate(shsrch(nthr))
            do i = 1,nthr
                call shsrch(i)%new_fixed(build, lims_reg, maxits=STUDENT_SHIFT_REG_MAXITS, coarse_init=.true.)
            enddo
        endif
        call student%prepare_predict(nthr)
        allocate(preds(nptcls), prob(student%get_ncls(), batchsz_max))
        allocate(e3_pred(nptcls), source=0.)
        allocate(sh_reg(2,nptcls), source=0.)
        allocate(class_stored(nptcls), source=0)
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
            call build_batch_particles2D(params, build, batchsz, pinds(b0:b1), reads, ptcl_match_imgs=ptcl_match_imgs, &
                &ptcl_match_imgs_pad=ptcl_match_imgs_pad)
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
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        call build%pftc%kill
        call build%esig%kill
        ! the table and the records: the route seeds a particle whose candidates fit the cap at or above the model's
        ! confidence threshold; every particle takes the predicted class, angle and shift
        t0 = tic()
        call fopen(pred_unit, fname_pred, status='replace', action='write', iostat=io_stat)
        call fileiochk('predict_classes; open '//fname_pred%to_char(), io_stat)
        write(pred_unit,'(A)') PRED_HEADER
        allocate(routes(nptcls), ncands(nptcls), source=0)
        allocate(cands_all(KMAX_CLASS,nptcls), source=0)
        n_seeded  = 0
        n_capped  = 0
        n_agree   = 0
        ncand_sum = 0.d0
        conf_sum  = 0.d0
        do i = 1,nptcls
            iptcl = pinds(i)
            route = 2
            if( preds(i)%k99 <= KMAX_CLASS .and. preds(i)%conf >= conf_threshold ) route = 1
            if( route == 1 ) n_seeded = n_seeded + 1
            routes(i)      = route
            ncands(i)      = preds(i)%ncand
            cands_all(:,i) = preds(i)%cands
            if( preds(i)%k99 > KMAX_CLASS ) n_capped = n_capped + 1
            if( preds(i)%mode == class_stored(i) ) n_agree = n_agree + 1
            ncand_sum = ncand_sum + real(preds(i)%ncand,dp)
            conf_sum  = conf_sum  + real(preds(i)%conf,dp)
            call build%spproj_field%set(iptcl, 'class', preds(i)%mode)
            call build%spproj_field%e3set(iptcl, e3_pred(i))
            call build%spproj_field%set_shift(iptcl, sh_reg(:,i))
            call build%spproj_field%set(iptcl, 'pconf', preds(i)%conf)
            call build%spproj_field%set(iptcl, 'route', route)
            write(pred_unit,'(I0,1X,I0,1X,I0,1X,F8.5,1X,I0,1X,I0,1X,F8.5,1X,I0,3(1X,F8.3),1X,I0)', advance='no') iptcl, &
                &build%spproj_field%get_eo(iptcl), preds(i)%mode, preds(i)%pmax, preds(i)%k99, preds(i)%ncand, &
                &preds(i)%conf, route, e3_pred(i), sh_reg(:,i), class_stored(i)
            do j = 1,preds(i)%ncand
                write(pred_unit,'(1X,I0)', advance='no') preds(i)%cands(j)
            enddo
            write(pred_unit,*)
        enddo
        call fclose(pred_unit)
        if( l_part )then
            ! the range's records and candidates for the master, which owns the project
            call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [params%fromp, params%top], &
                &isegment=PTCL2D_SEG)
            call cands_art%write_part(build%spproj_field%get_noris(), params%fromp, params%top, &
                &string(base%to_char()//'_cands.bin'), pinds, routes, ncands, cands_all, student%get_ncls())
        else
            call build%spproj%write_segment_inside(params%oritype)
            ! the candidate artifact for the sparse 2D table, registered under the particles' layout digest
            call cands_art%publish(build%spproj, build%spproj_field, params%projfile, string(base%to_char()//'_cands.bin'), &
                &pinds, routes, ncands, cands_all, student%get_ncls())
        endif
        call cands_art%kill
        rt_write = toc(t0)
        rt_all   = toc(t_all)
        write(logfhandle,'(A,F8.4,A,F8.5,A,F8.4,A,F8.2,A,F8.4)') '>>> CLASS STUDENT PREDICTION: seeded ', &
            &real(n_seeded) / real(nptcls), ' at the confidence threshold ', conf_threshold, ', capped ', &
            &real(n_capped) / real(nptcls), ', mean candidates ', ncand_sum / real(nptcls,dp), ', mean confidence ', &
            &conf_sum / real(nptcls,dp)
        write(logfhandle,'(A,F8.4)') '>>> CLASS STUDENT PREDICTION: the stored class is the predicted one for ', &
            &real(n_agree) / real(nptcls)
        write(logfhandle,'(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A,F9.2,A)') '>>> CLASS STUDENT PREDICTION: prep ', rt_prep, &
            &' s, particles ', rt_batch, ' s, features, forward and registration ', rt_pred, ' s, records ', rt_write, &
            &' s, total ', rt_all, ' s'
        write(logfhandle,'(A,F9.3,A)') '>>> CLASS STUDENT PREDICTION: ', rt_pred / real(nptcls,timer_int_kind) * 1000., &
            &' ms per particle'
        call student%kill
        call feats%kill
        if( allocated(shsrch) )then
            do i = 1,size(shsrch)
                call shsrch(i)%kill
            enddo
            deallocate(shsrch)
        endif
        deallocate(preds, prob, e3_pred, sh_reg, class_stored, pinds, batches, wk, routes, ncands, cands_all)

    contains

        !> a part whose range holds no selected particle: the table's header, the range's records unchanged and its
        !! empty candidates, the contributions the master's merges expect
        subroutine write_empty_part
            integer :: none(0), none2(KMAX_CLASS,0)
            call fopen(pred_unit, fname_pred, status='replace', action='write', iostat=io_stat)
            call fileiochk('predict_classes; open '//fname_pred%to_char(), io_stat)
            write(pred_unit,'(A)') PRED_HEADER
            call fclose(pred_unit)
            call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [params%fromp, params%top], &
                &isegment=PTCL2D_SEG)
            call cands_art%write_part(build%spproj_field%get_noris(), params%fromp, params%top, &
                &string(base%to_char()//'_cands.bin'), none, none, none, none2, student%get_ncls())
            write(logfhandle,'(A,I0,A,I0,A)') '>>> CLASS STUDENT PREDICTION: no selected particle in the part''s range ', &
                &params%fromp, '-', params%top, '; its records written unchanged'
            call student%kill
            call feats%kill
            deallocate(pinds, wk)
        end subroutine write_empty_part

        !> one particle's features, forward pass, candidates and read-out angle, and, given the class averages, the
        !! shift increment registered against the predicted class at that angle, on the calling thread
        subroutine predict_one( i_here, iptcl_here, islot )
            integer, intent(in) :: i_here, iptcl_here, islot
            type(polar_feat) :: feat_here
            real    :: cxy(3)
            integer :: ithr_here, irot
            ithr_here = 1
            !$ ithr_here = omp_get_thread_num() + 1
            call feats%build(build%pftc, iptcl_here, feat_here)
            call student%predict_feat(ithr_here, feat_here, prob(:,islot))
            call student%cands_of(prob(:,islot), preds(i_here))
            e3_pred(i_here)      = student%rotation_of(ithr_here, preds(i_here)%mode)
            class_stored(i_here) = build%spproj_field%get_class(iptcl_here)
            sh_reg(:,i_here)     = build%spproj_field%get_2Dshift(iptcl_here)
            if( l_reg )then
                irot = build%pftc%get_roind(360. - e3_pred(i_here))
                call shsrch(ithr_here)%set_indices(preds(i_here)%mode, iptcl_here)
                cxy = shsrch(ithr_here)%minimize(irot)
                if( irot > 0 ) sh_reg(:,i_here) = sh_reg(:,i_here) + cxy(2:3)
            endif
            call feat_here%kill
        end subroutine predict_one

    end subroutine predict_classes_exec

end module simple_strategy_class_student
