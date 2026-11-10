!@descr: the pose student's split-half protocol in miniature on the blob phantom's simulated particles, a gate of simple_commanders_test_highlevel, and the HHCNN kernel benchmark
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_pose
implicit none
#include "simple_local_flags.inc"

contains

    !> The pose student's nightly route (doc/implementation_notes/completed/learned_pose_state_assignment.md,
    !! section 10 tier 6 and ruling R3 item 7): the split-half protocol in miniature on 2 000 simulated particles
    !! of the blob phantom of the multi-state test fixture (one state, asymmetric, SNR 1, no shift) with their truth
    !! poses. The training half trains the student against the truth map, the test
    !! half is predicted from blank poses, the neighbourhood search (prob_neigh geom) aligns every test particle
    !! from the student's cone where it has one and from its coarse cell where it has none, the same search runs
    !! on a copy without student records, the dense search (prob) runs on every row as the protocol's fallback,
    !! and the half is reconstructed from the protocol's merge (the seeded rows' poses from the cone search, the
    !! others' from the dense search) and from the dense search alone. The gate asserts the machinery (the
    !! route's contract per particle, the fold's bound on the cone, the rows without a cone identical between
    !! the search with and without student records, the seeded rows changed by their cones, finite poses) and
    !! reports the reads (the held-out recovery, the seeded fraction, the cone's coverage of the truth, the
    !! agreement with the truth per route beside the dense search's on the same rows, the maps' FSC against the
    !! truth, the costs): no floor is taken from a student that seeds nobody; a profile sets floors later.
    !! Every row goes to metrics.tsv. The classic generator is reset to the gate's seed before every search, so
    !! that in-process runs draw alike (the keyed draws share the run's seed).
    module subroutine exec_test_pose_student( self, cline )
        class(commander_test_pose_student), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        class(parameters), allocatable :: params
        type(string) :: cwd_saved, fixture_root
        integer      :: status
        logical      :: all_ok
        allocate(params)
        call simple_getcwd(cwd_saved)
        fixture_root = filepath(cwd_saved, 'test_pose_student_'//int2str(get_process_id()))
        if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
        call simple_mkdir(fixture_root)
        call simple_chdir(fixture_root, status)
        if( status /= 0 ) THROW_HARD('TEST_POSE_STUDENT FAILED: could not enter fixture directory')
        call params%new(cline)
        all_ok = .true.
        call run_pose_student_gate(params%nthr, filepath(cwd_saved, 'metrics.tsv'), all_ok)
        call simple_chdir(cwd_saved, status)
        if( status /= 0 ) THROW_HARD('TEST_POSE_STUDENT FAILED: could not restore original directory')
        if( all_ok )then
            call simple_rmdir(fixture_root)
            write(logfhandle,'(A)') 'PASS: the pose student''s split-half protocol in miniature'
            call simple_end('**** SIMPLE_TEST_POSE_STUDENT NORMAL STOP ****')
        else
            THROW_HARD('TEST_POSE_STUDENT FAILED')
        endif
    end subroutine exec_test_pose_student

    !> The HHCNN kernel benchmark on synthetic tensors at a production half ring: per ring count and width the forward
    !! and backward passes, the training step at one thread and at the team, the team's prediction throughput and the
    !! step's thread identity, to BENCH/POSE_STUDENT_KERNELS_BENCH.txt and the log; a measurement, no thresholds
    module subroutine exec_test_pose_student_bench( self, cline )
        use simple_hhcnn,    only: run_hhcnn_bench
        use simple_bench_io, only: open_bench_report
        !$ use omp_lib,      only: omp_get_max_threads, omp_set_num_threads
        class(commander_test_pose_student_bench), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        integer, parameter :: PFTSZ_BENCH     = 160              !< the half ring of a 270 A mask at 2 A
        integer, parameter :: NKS_BENCH(4)    = [12, 18, 25, 37] !< rings of the band at 20, 15, 12 and 8 A of a 300 A box
        integer, parameter :: WIDTHS_BENCH(3) = [16, 32, 64]     !< the designs learn_pose_student tries
        integer, parameter :: NREP_BENCH      = 5
        integer :: fnr, nthr, nspace
        logical :: ok_identity
        nthr = 1
        !$ nthr = omp_get_max_threads()
        if( cline%defined('nthr') ) nthr = cline%get_iarg('nthr')
        nspace = 1000
        if( cline%defined('nspace') ) nspace = cline%get_iarg('nspace')
        if( nspace < 2 .or. mod(nspace, 2) /= 0 ) THROW_HARD('nspace must be even: the classes are its mirror pairs')
        !$ call omp_set_num_threads(nthr)
        call open_bench_report(fnr, 'POSE_STUDENT_KERNELS')
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,i0)') 'student half-ring angles (pftsz)    : ', PFTSZ_BENCH
        write(fnr,'(a,i0)') 'student classes (mirror pairs)      : ', nspace / 2
        write(fnr,'(a,i0)') 'student team threads                : ', nthr
        write(fnr,'(a,i0)') 'student repetitions per measurement : ', NREP_BENCH
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** KERNELS (us per sample, ms per step of 256 samples) ***'
        call run_hhcnn_bench(PFTSZ_BENCH, NKS_BENCH, WIDTHS_BENCH, nspace / 2, nthr, NREP_BENCH, fnr, ok_identity)
        write(fnr,'(a)') ''
        write(fnr,'(a,l1)') 'training step identical at one thread and at the team: ', ok_identity
        call fclose(fnr)
        write(logfhandle,'(a,l1)') '>>> POSE STUDENT BENCH: training step identical at one thread and at the team: ', ok_identity
        call simple_end('**** SIMPLE_TEST_POSE_STUDENT_BENCH NORMAL STOP ****')
    end subroutine exec_test_pose_student_bench

    subroutine run_pose_student_gate( nthr, metrics_file, all_ok )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use simple_test_state_fixture,      only: state_fixture, FIX_KNOWN
        use simple_ui,                      only: make_ui
        use simple_image,                   only: image
        use simple_sym,                     only: sym
        use simple_ori_utils,               only: euler_mirror
        use simple_timer,                   only: tic, toc, timer_int_kind
        use simple_estimate_ssnr,           only: get_resolution
        use simple_commanders_pose_student, only: commander_learn_pose_student, commander_predict_poses
        use simple_commanders_refine3D,     only: commander_refine3D
        use simple_commanders_rec,          only: commander_rec3D
        use simple_test_gate,               only: test_gate
        integer,       intent(in)    :: nthr
        class(string), intent(in)    :: metrics_file
        logical,       intent(inout) :: all_ok
        character(len=*), parameter :: GATE_DIR   = 'pose_student_gate'
        character(len=*), parameter :: TRUTH_VOL  = 'truth.mrc'
        character(len=*), parameter :: PROJNAME   = 'pose_gate'
        character(len=*), parameter :: MODEL      = 'student.txt'
        character(len=*), parameter :: PRED_TABLE = 'student_pred.txt'
        integer, parameter :: NPTCLS  = 2000       !< 1 000 training rows (100 held out), 1 000 test rows
        integer, parameter :: NSPACE  = 1000       !< the student's grid (NSPACE_STUDENT)
        integer, parameter :: GATE_SEED = 20261108 !< the classic generator's seed before every search (in-process identity)
        real,    parameter :: SNR     = 1.         !< the fixture's phantom SNR
        real,    parameter :: LP      = 8.         !< the band: 16 rings at the blob box (64 at 2 A), two ring bands in the readout
        real,    parameter :: SHIFT_SIGMA = 2.     !< the centring width of the training's shift disk (3 widths): the 2 px of the measurements, not the 10 px default
        real,    parameter :: AGREE_DEG = 6., TWO_STEPS_DEG = 12.
        type(state_fixture)                 :: fix
        type(commander_learn_pose_student)  :: xlearn
        type(commander_predict_poses)       :: xpredict
        type(commander_refine3D)            :: xrefine3D
        type(commander_rec3D)               :: xrec3D
        class(cmdline),    allocatable :: cl
        class(sp_project), allocatable :: proj, pred, searched, control, dense
        type(image)         :: vol, vol_rec
        type(oris)          :: truth
        type(ori)           :: o_truth
        type(test_gate)     :: gate
        type(sym)           :: symop
        type(string)        :: root, gate_root, truth_abs, model_abs, proj_abs
        real, allocatable   :: corrs(:), res(:)
        integer :: box
        real    :: smpd, mskdiam, mskrad
        real    :: conf_threshold, cone_bound, heldout_seeded, heldout_recovered, time_train, spacing
        real    :: radius(NPTCLS), euls(3), euls_m(3), euls_t(3), euls_sym(3), d1, d2, inpld, fsc05, fsc0143
        real    :: wall_learn, wall_predict, wall_search, wall_control, wall_dense, wall_rec, fsc05_d, fsc0143_d
        integer :: i, j, k, status, route(NPTCLS), nsel, nseeded, ndense, n_in_cone, n_two, nshell
        integer :: n6(0:1), n12(0:1), nroute(0:1), nk12, n6d(0:1), n12d(0:1)
        logical :: l_route_ok, l_cone_ok, l_dense_ok, l_fold_ok, l_fallback_same, l_seeded_changed, l_finite, l_sel(NPTCLS)
        integer(timer_int_kind) :: t0
        allocate(cl, proj, pred, searched, control, dense)
        call make_ui
        write(logfhandle,'(a)') '>>> TEST_POSE_STUDENT: the split-half protocol in miniature'
        call simple_getcwd(root)
        if( file_exists(GATE_DIR) ) call simple_rmdir(GATE_DIR, status)
        call simple_mkdir(GATE_DIR)
        call simple_chdir(GATE_DIR, status)
        if( status /= 0 ) THROW_HARD('Could not enter '//GATE_DIR)
        call simple_getcwd(gate_root)
        call gate%new(metrics_file)
        call symop%new('c1')
        ! ---- the fixture: the blob phantom's particles at their truth poses (FIX_KNOWN), its truth map soft-masked
        !      the references of the teacher and of the search ----
        call fix%new('blob', FIX_KNOWN, 1, NPTCLS, SNR, nthr, sherr=0.)
        call fix%build
        box     = fix%box()
        smpd    = fix%smpd()
        mskdiam = fix%mskdiam()
        mskrad  = 0.5 * mskdiam / smpd
        call vol%new([box,box,box], smpd)
        call vol%read(fix%truth_fname(1))
        call vol%mask3D_soft(mskrad, backgr=0.)
        call vol%write(string(TRUTH_VOL))
        truth_abs = simple_abspath(string(TRUTH_VOL))
        call vol%fft                               ! the FSC reads the transforms
        call proj%read(fix%project_fname())
        if( proj%os_ptcl3D%get_noris() /= NPTCLS ) THROW_HARD('the fixture holds another particle count; test_pose_student')
        call truth%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call proj%os_ptcl3D%get_ori(i, o_truth)
            call truth%set_ori(i, o_truth)
            call proj%os_ptcl3D%set_shift(i, [0.,0.])
            call proj%os_ptcl3D%set_state(i, 1)
            call proj%os_ptcl3D%set(i, 'eo', real(mod((i+1)/2, 2)))   ! the halves' own even/odd partition, in both fields
            call proj%os_ptcl2D%set(i, 'eo', real(mod((i+1)/2, 2)))
        enddo
        call o_truth%kill
        call fix%kill
        call simple_chdir(gate_root, status)
        ! ---- the training half: the student against the truth map ----
        call enter_stage('train')
        do i = 1, NPTCLS
            if( mod(i,2) == 0 ) call proj%os_ptcl3D%set_state(i, 0)
        enddo
        call proj%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        t0 = tic()
        call cl%set('prg',         'learn_pose_student')
        call cl%set('projfile',    simple_abspath(string(PROJNAME//'.simple')))
        call cl%set('vol1',        truth_abs)
        call cl%set('lp',          LP)
        call cl%set('nspace',      NSPACE)
        call cl%set('mskdiam',     mskdiam)
        call cl%set('pgrp',        'c1')
        call cl%set('nthr',        nthr)
        call cl%set('shift_sigma', SHIFT_SIGMA)
        call cl%set('fname',       MODEL)
        call xlearn%execute(cl)
        call cl%kill
        wall_learn = real(toc(t0))
        model_abs  = simple_abspath(string(MODEL))
        call read_model_record(model_abs)
        call gate%check('model_records_threshold_and_bound', ieee_is_finite(conf_threshold) .and. cone_bound > 0.)
        call gate%report('heldout_fraction_seeded',    heldout_seeded)
        call gate%report('heldout_recovered_fraction', heldout_recovered)
        call gate%report('training_wall_s',            wall_learn)
        call gate%report('model_training_time_s',      time_train)
        call simple_chdir(gate_root, status)
        ! ---- the test half: blank poses, the prediction, the records ----
        call enter_stage('test')
        do i = 1, NPTCLS
            call proj%os_ptcl3D%set_state(i, merge(1, 0, mod(i,2) == 0))
            call proj%os_ptcl3D%set_euler(i, [0.,0.,0.])
            call proj%os_ptcl3D%set_shift(i, [0.,0.])
        enddo
        call proj%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        proj_abs = simple_abspath(string(PROJNAME//'.simple'))
        t0 = tic()
        call cl%set('prg',      'predict_poses')
        call cl%set('projfile', proj_abs)
        call cl%set('infile',   model_abs)
        call cl%set('mskdiam',  mskdiam)
        call cl%set('nthr',     nthr)
        call cl%set('fname',    MODEL)
        call xpredict%execute(cl)
        call cl%kill
        wall_predict = real(toc(t0))
        call pred%read(proj_abs)
        call read_pred_table(string(PRED_TABLE))
        nsel = 0; nseeded = 0; ndense = 0; n_in_cone = 0; n_two = 0
        l_route_ok = .true.; l_cone_ok = .true.; l_dense_ok = .true.; l_fold_ok = .true.
        spacing = huge(spacing)
        do i = 1, NPTCLS
            l_sel(i) = pred%os_ptcl3D%get_state(i) == 1
            if( .not. l_sel(i) ) cycle
            nsel  = nsel + 1
            route(i) = nint(pred%os_ptcl3D%get(i, 'route'))
            if( route(i) /= 1 .and. route(i) /= 2 ) l_route_ok = .false.
            euls   = pred%os_ptcl3D%get_euler(i)
            euls_t = truth%get_euler(i)
            call symop%sym_dists(euls, euls_t, euls_sym, d1, inpld)
            call euler_mirror(euls, euls_m)
            call symop%sym_dists(euls_m, euls_t, euls_sym, d2, inpld)
            if( route(i) == 1 )then
                nseeded = nseeded + 1
                if( .not. (pred%os_ptcl3D%get(i, 'cone') > 0. .and. pred%os_ptcl3D%get(i, 'cone') <= cone_bound + 1.e-3 .and. &
                    &pred%os_ptcl3D%get(i, 'pconf') >= conf_threshold - 1.e-6) ) l_cone_ok = .false.
                spacing = min(spacing, pred%os_ptcl3D%get(i, 'cone'))
                ! the cone's coverage: the truth within the cone about the stored direction or about its mirror
                if( min(d1, d2) <= pred%os_ptcl3D%get(i, 'cone') + 1.e-3 ) n_in_cone = n_in_cone + 1
                if( min(d1, d2) <= TWO_STEPS_DEG ) n_two = n_two + 1
            else
                ndense = ndense + 1
                if( .not. (pred%os_ptcl3D%get(i, 'cone') == 0.) ) l_dense_ok = .false.
            endif
            if( radius(i) > 90. + 1.e-3 ) l_fold_ok = .false.
        enddo
        call gate%check('every_test_row_routed',                 l_route_ok .and. nsel == NPTCLS/2)
        call gate%check('seeded_rows_cone_within_bound_and_confidence', l_cone_ok)
        call gate%check('dense_rows_carry_no_cone',               l_dense_ok)
        call gate%check('folded_cone_radius_at_most_90deg',       l_fold_ok)
        call gate%report('seeded_fraction',                       real(nseeded)/real(max(1,nsel)))
        call gate%report('seeded_truth_inside_cone_fraction',     real(n_in_cone)/real(max(1,nseeded)))
        call gate%report('seeded_mode_within_12deg_of_truth',     real(n_two)/real(max(1,nseeded)))
        call gate%report('prediction_wall_s',                     wall_predict)
        call gate%report('prediction_ms_per_particle',            1000.*wall_predict/real(max(1,nsel)))
        call simple_chdir(gate_root, status)
        ! ---- the neighbourhood search on every test row: the cone where there is one, the cell where not ----
        call enter_stage('search')
        call pred%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        call set_fixed_seed(GATE_SEED)
        t0 = tic()
        call run_search(simple_abspath(string(PROJNAME//'.simple')), 'prob_neigh')
        wall_search = real(toc(t0))
        call searched%read(simple_abspath(string(PROJNAME//'.simple')))
        call simple_chdir(gate_root, status)
        ! ---- the control: the same search without student records (the base path for every row) ----
        call enter_stage('control')
        do i = 1, NPTCLS
            call pred%os_ptcl3D%set(i, 'cone',  0.)
            call pred%os_ptcl3D%set(i, 'pconf', 0.)
            call pred%os_ptcl3D%set(i, 'route', 0.)
        enddo
        call pred%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        call set_fixed_seed(GATE_SEED)
        t0 = tic()
        call run_search(simple_abspath(string(PROJNAME//'.simple')), 'prob_neigh')
        wall_control = real(toc(t0))
        call control%read(simple_abspath(string(PROJNAME//'.simple')))
        call simple_chdir(gate_root, status)
        ! ---- the protocol's fallback: the dense search on every row from the blank poses ----
        call enter_stage('dense')
        call pred%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        call set_fixed_seed(GATE_SEED)
        t0 = tic()
        call run_search(simple_abspath(string(PROJNAME//'.simple')), 'prob')
        wall_dense = real(toc(t0))
        call dense%read(simple_abspath(string(PROJNAME//'.simple')))
        call simple_chdir(gate_root, status)
        ! the reads after the searches: the rows without a cone bit for bit the control's, the seeded rows changed
        ! by their cones, every pose finite, the agreement with the truth per route beside the dense search's on
        ! the same rows
        l_fallback_same = .true.; l_seeded_changed = .false.; l_finite = .true.
        n6 = 0; n12 = 0; nroute = 0; n6d = 0; n12d = 0
        do i = 1, NPTCLS
            if( .not. l_sel(i) ) cycle
            euls   = searched%os_ptcl3D%get_euler(i)
            euls_m = control%os_ptcl3D%get_euler(i)
            l_finite = l_finite .and. all(ieee_is_finite(euls)) .and. all(ieee_is_finite(searched%os_ptcl3D%get_2Dshift(i))) &
                &.and. all(ieee_is_finite(dense%os_ptcl3D%get_euler(i))) .and. all(ieee_is_finite(dense%os_ptcl3D%get_2Dshift(i)))
            if( route(i) == 2 )then
                if( any(euls /= euls_m) .or. any(searched%os_ptcl3D%get_2Dshift(i) /= control%os_ptcl3D%get_2Dshift(i)) ) &
                    &l_fallback_same = .false.
            else
                if( any(euls /= euls_m) .or. any(searched%os_ptcl3D%get_2Dshift(i) /= control%os_ptcl3D%get_2Dshift(i)) ) &
                    &l_seeded_changed = .true.
            endif
            euls_t = truth%get_euler(i)
            j = merge(1, 0, route(i) == 1)
            nroute(j) = nroute(j) + 1
            call symop%sym_dists(euls, euls_t, euls_sym, d1, inpld)
            if( d1 <= AGREE_DEG )     n6(j)  = n6(j)  + 1
            if( d1 <= TWO_STEPS_DEG ) n12(j) = n12(j) + 1
            euls = dense%os_ptcl3D%get_euler(i)
            call symop%sym_dists(euls, euls_t, euls_sym, d1, inpld)
            if( d1 <= AGREE_DEG )     n6d(j)  = n6d(j)  + 1
            if( d1 <= TWO_STEPS_DEG ) n12d(j) = n12d(j) + 1
        enddo
        call gate%check('every_searched_pose_finite',           l_finite)
        call gate%check('rows_without_cone_identical_to_no_student_search', l_fallback_same)
        if( nseeded > 0 )then
            call gate%check('seeded_rows_changed_by_their_cones', l_seeded_changed)
        else
            call gate%report('seeded_rows_changed_by_their_cones', 0.)
        endif
        call gate%report('seeded_within_6deg_of_truth',          real(n6(1))/real(max(1,nroute(1))))
        call gate%report('seeded_within_12deg_of_truth',         real(n12(1))/real(max(1,nroute(1))))
        call gate%report('seeded_rows_by_dense_within_6deg',     real(n6d(1))/real(max(1,nroute(1))))
        call gate%report('unseeded_rows_by_dense_within_6deg',   real(n6d(0))/real(max(1,nroute(0))))
        call gate%report('unseeded_rows_by_dense_within_12deg',  real(n12d(0))/real(max(1,nroute(0))))
        call gate%report('all_rows_by_dense_within_6deg',        real(n6d(0)+n6d(1))/real(max(1,nsel)))
        call gate%report('cone_search_wall_s',                   wall_search)
        call gate%report('control_search_wall_s',                wall_control)
        call gate%report('dense_search_wall_s',                  wall_dense)
        ! ---- the test half reconstructed from the protocol's merge and from the dense search, against the truth ----
        nshell = box/2 - 1
        allocate(corrs(nshell), res(nshell))
        do k = 1, nshell
            res(k) = real(box)*smpd/real(k)
        enddo
        nk12 = max(1, min(nshell, nint(real(box)*smpd/12.)))
        call enter_stage('merged')
        do i = 1, NPTCLS
            if( .not. l_sel(i) ) cycle
            if( route(i) == 1 )then
                call dense%os_ptcl3D%set_euler(i, searched%os_ptcl3D%get_euler(i))
                call dense%os_ptcl3D%set_shift(i, searched%os_ptcl3D%get_2Dshift(i))
            endif
        enddo
        call dense%write(simple_abspath(string(PROJNAME//'.simple'), check_exists=.false.))
        t0 = tic()
        call reconstruct_here(fsc05, fsc0143, corrs(nk12))
        wall_rec = real(toc(t0))
        call gate%report('merged_map_fsc_0143_vs_truth_A', fsc0143)
        call gate%report('merged_map_fsc_05_vs_truth_A',   fsc05)
        call gate%report('merged_map_fsc_at_12A_vs_truth', corrs(nk12))
        call gate%report('reconstruction_wall_s',          wall_rec)
        call simple_chdir(gate_root, status)
        call simple_chdir(string('dense'), status)
        call reconstruct_here(fsc05_d, fsc0143_d, d1)
        call gate%report('dense_map_fsc_0143_vs_truth_A', fsc0143_d)
        call gate%report('dense_map_fsc_05_vs_truth_A',   fsc05_d)
        call gate%report('dense_map_fsc_at_12A_vs_truth', d1)
        call simple_chdir(gate_root, status)
        all_ok = all_ok .and. gate%passed()
        call gate%kill
        call simple_chdir(root, status)
        call vol%kill
        call vol_rec%kill
        call truth%kill
        call symop%kill
        call proj%kill
        call pred%kill
        call searched%kill
        call control%kill
        call dense%kill
        deallocate(cl, proj, pred, searched, control, dense, corrs, res)

    contains

        subroutine enter_stage( dir )
            character(len=*), intent(in) :: dir
            call simple_mkdir(dir)
            call simple_chdir(string(dir), status)
            if( status /= 0 ) THROW_HARD('Could not enter '//dir)
        end subroutine enter_stage

        !> one refine3D iteration against the truth map, in place: refine=prob_neigh (prob_neigh_mode=geom, the
        !! cone or the cell per particle) or refine=prob (the dense search)
        subroutine run_search( projfile, refine )
            class(string),    intent(in) :: projfile
            character(len=*), intent(in) :: refine
            call cl%set('prg',             'refine3D')
            call cl%set('projfile',        projfile)
            call cl%set('mkdir',           'no')
            call cl%set('vol1',            truth_abs)
            call cl%set('refine',          refine)
            if( refine == 'prob_neigh' ) call cl%set('prob_neigh_mode', 'geom')
            call cl%set('objfun',          'euclid')
            call cl%set('pgrp',            'c1')
            call cl%set('mskdiam',         mskdiam)
            call cl%set('lp',              LP)
            call cl%set('center',          'no')
            call cl%set('maxits',          1)
            call cl%set('nthr',            nthr)
            call xrefine3D%execute(cl)
            call cl%kill
        end subroutine run_search

        !> reconstruct3D of the project in the current directory and the map's FSC against the truth: the
        !! resolutions at 0.5 and 0.143 and the correlation at the shell nearest 12 A
        subroutine reconstruct_here( r05, r0143, c12 )
            real, intent(out) :: r05, r0143, c12
            call cl%set('prg',      'reconstruct3D')
            call cl%set('projfile', simple_abspath(string(PROJNAME//'.simple')))
            call cl%set('mkdir',    'no')
            call cl%set('mskdiam',  mskdiam)
            call cl%set('pgrp',     'c1')
            call cl%set('nthr',     nthr)
            call xrec3D%execute(cl)
            call cl%kill
            call vol_rec%new([box,box,box], smpd)
            call vol_rec%read(string('recvol_state01.mrc'))
            call vol_rec%mask3D_soft(mskrad, backgr=0.)
            call vol_rec%fft
            call vol_rec%fsc(vol, corrs)
            call get_resolution(corrs, res, r05, r0143)
            c12 = corrs(nk12)
            call vol_rec%kill
        end subroutine reconstruct_here

        !> the model header's threshold, bound and held-out record
        subroutine read_model_record( fname )
            class(string), intent(in) :: fname
            character(len=512) :: line
            integer :: funit, io_stat, ieq
            conf_threshold = huge(conf_threshold); cone_bound = 0.; heldout_seeded = 0.; heldout_recovered = 0.; time_train = 0.
            call fopen(funit, fname, status='old', action='read', iostat=io_stat)
            if( io_stat /= 0 ) THROW_HARD('cannot open the model '//fname%to_char())
            do
                read(funit,'(A)', iostat=io_stat) line
                if( io_stat /= 0 ) exit
                ieq = index(line, '=')
                if( ieq < 2 ) cycle
                select case(trim(line(1:ieq-1)))
                    case('conf_threshold');    read(line(ieq+1:),*) conf_threshold
                    case('cone_bound_deg');    read(line(ieq+1:),*) cone_bound
                    case('heldout_seeded');    read(line(ieq+1:),*) heldout_seeded
                    case('heldout_recovered'); read(line(ieq+1:),*) heldout_recovered
                    case('time_train_s');      read(line(ieq+1:),*) time_train
                    case default
                end select
            enddo
            call fclose(funit)
        end subroutine read_model_record

        !> the prediction table's cone radius per particle (the folded 99% set's spread before the route's floor)
        subroutine read_pred_table( fname )
            class(string), intent(in) :: fname
            character(len=512) :: line
            real    :: vals(14)
            integer :: funit, io_stat, iptcl
            radius = 0.
            call fopen(funit, fname, status='old', action='read', iostat=io_stat)
            if( io_stat /= 0 ) THROW_HARD('cannot open the prediction table '//fname%to_char())
            do
                read(funit,'(A)', iostat=io_stat) line
                if( io_stat /= 0 ) exit
                if( line(1:1) == '#' ) cycle
                read(line,*, iostat=io_stat) vals
                if( io_stat /= 0 ) cycle
                iptcl = nint(vals(1))
                if( iptcl >= 1 .and. iptcl <= NPTCLS ) radius(iptcl) = vals(8)
            enddo
            call fclose(funit)
        end subroutine read_pred_table

    end subroutine run_pose_student_gate

end submodule simple_commanders_test_highlevel_pose
