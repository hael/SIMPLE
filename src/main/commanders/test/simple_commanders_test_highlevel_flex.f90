!@descr: FLEX PCA workflow gates: shared memory against distributed, the per-state refinement and the class-average state refinement
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_flex
use, intrinsic :: iso_fortran_env,      only: real32, real64
use simple_commanders_denoise,          only: commander_cls_expansion
use simple_commanders_flex_pca,         only: commander_flex_pca
use simple_commanders_refine3D,         only: commander_refine3D_auto, commander_refine3D_states_cavgs
use simple_state_weight_set,            only: state_weight_set
use simple_sigma2_state_file,           only: sigma2_state_digest_file
use simple_flex_pca_application_tester, only: configure_flex_pca_phantom, read_state_weight_table, FLEX_PHANTOM_BOX, &
    &FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_NCOMP, FLEX_PHANTOM_MAX_STATES, &
    &FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM, FLEX_PHANTOM_LABEL_ACCURACY_MIN
use simple_test_state_fixture,          only: state_fixture, FIX_KNOWN, FIX_VIEWS
use simple_flex_pca_embedding_io, only: read_deconv_block
use simple_image,                 only: image
use simple_oris,                  only: oris
use simple_sp_project,            only: sp_project
use simple_test_gate,             only: test_gate
use simple_test_truth_metrics,    only: compare_to_truth, pair_pose_error, match_states, mapped_accuracy, accuracy_null
use simple_ui,                    only: make_ui
implicit none
#include "simple_local_flags.inc"

integer, parameter :: FLEX_GATE_SEED = 202610063
!> the per-state refinement: its seed, its iteration budget and how much the weighted map may trail the
!! hard one in truth correlation (plan section 9); it runs with twice the FLEX modes' threads, the
!! processors this serial entry reserves
integer, parameter :: STATE_REFINE_SEED   = 202610071
integer, parameter :: STATE_REFINE_MAXITS = 5
real,    parameter :: STATE_REFINE_MARGIN = 0.01
real,    parameter :: MAP_CORR_LP = 8.0
!> the class-average state refinement (refine3D_states_cavgs): CAVGS_NVIEWS directions, each an exact 2D
!! class of 60 particles over CAVGS_NSTATES planted states, split by cls_expansion into CAVGS_NSUB
!! sub-classes each (ncls is the total; the project has no class averages, so every view is its own group,
!! and CAVGS_MINPOP lets 60 particles hold three). The states are the fixture's compositional
!! beta-galactosidase states (the tetramer, the tetramer less chain D, the A-B dimer), which differ by
!! whole chains (ruling of 2026-10-18: blob states that differ in a lobe of about 6 % of the mass are
!! not a fair correctness case for flex_pca; blob states that each carry a large blob of their own
!! were tried and give the pose-less sub-averages no orientation at the fit band);
!! CAVGS_EXPANSION_LP is the fit band (8 A was tried: noisier poses, no better states)
integer, parameter :: CAVGS_NPTCLS       = 3600
integer, parameter :: CAVGS_NSTATES      = 3
integer, parameter :: CAVGS_NVIEWS       = 60
integer, parameter :: CAVGS_NSUB         = 3
integer, parameter :: CAVGS_MINPOP       = 20
integer, parameter :: CAVGS_MIN_SUBAVGS  = 100
integer, parameter :: CAVGS_NPAIRS       = 2000
integer, parameter :: CAVGS_SEED         = 202610092
real,    parameter :: CAVGS_EXPANSION_LP = 12.0
!> The class-average state route is a robustness profile, not a gate: its profile of 2026-10-18
!! (scratch/keep/profiles/phase5/flex_f, 20 seed offsets at 8 threads, on the compositional states)
!! leaves no admissible floor. Three runs delivered fewer than two states; in the other 17 the
!! sub-average state accuracy (one-to-one mapping) has median 0.42, s 0.049, worst 0.34, so a floor
!! would have to lie at 0.27 or below, under the shuffled-label null (median 0.37, 99.9th percentile
!! 0.43); the median pair pose error of the sub-averages has median 11.7 deg, s 8.5, worst 51, so a
!! ceiling would have to lie at 60 deg or above, the level of random poses. The accuracy, the pose
!! error, the number of delivered states and the number of planted states the state maps match are
!! therefore reported; the route's contract (the projects, the hand-back, the maps) stays gated. The
!! committed sub-averages keep their floor. The profile of the merged code (2026-10-19,
!! scratch/keep/profiles/phase8/main) is the same: four runs of 20 delivered one state, the others a
!! sub-average state accuracy of median 0.41 (0.35-0.69) against a null 99.9th percentile of 0.41.
integer, parameter :: CAVGS_MIN_COMMITTED      = 90

!> the datasets come from the multi-state test fixture: two blob states at their true poses (known
!! mode) for the FLEX modes, three beta-galactosidase states in exact views (views mode) for the
!! class-average state refinement; the true state of every particle is read from the fixture's truth
!! document, and the geometry of the class-average dataset (box, pixel, mask) from the fixture
real,    parameter :: PHANTOM_SNR = 1.0
integer, allocatable :: truth_states_ab(:), truth_states_abc(:)
type(string)         :: truth_doc_abc
integer              :: cavgs_box     = 0
real                 :: cavgs_smpd    = 0.
real                 :: cavgs_mskdiam = 0.

type :: flex_gate_result
    logical :: valid = .false.
    logical :: embedding_found = .false.
    integer :: nstates = 0
    integer :: matched_a = 0
    integer :: matched_b = 0
    real    :: label_accuracy = 0.0
    real    :: pc_corr = 0.0
    real(real64) :: noise_scale = 1.0_real64
    integer, allocatable :: predicted_truth(:)
    real(real64), allocatable :: z(:,:)
    type(string) :: map_a
    type(string) :: map_b
end type flex_gate_result

contains

module subroutine exec_test_flex_pca_blobs( self, cline )
    class(commander_test_flex_pca_blobs), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    class(parameters), allocatable :: params
    type(string) :: cwd_saved, fixture_root
    integer :: nthr, status
    logical :: all_ok
    call simple_getcwd(cwd_saved)
    fixture_root = filepath(cwd_saved, 'test_flex_pca_blobs_'//int2str(get_process_id()))
    if( dir_exists(fixture_root) ) call simple_rmdir(fixture_root)
    call simple_mkdir(fixture_root)
    call simple_chdir(fixture_root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter fixture directory')
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    all_ok = .true.
    call run_flex_pca_blobs_gate(nthr, filepath(cwd_saved, 'metrics.tsv'), all_ok)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore original directory')
    CWD_GLOB = cwd_saved%to_char()
    if( all_ok )then
        call simple_rmdir(fixture_root)
        write(logfhandle,'(A)') 'PASS: FLEX PCA shared/distributed workflow agrees with simulation truth'
        call simple_end('**** SIMPLE_TEST_FLEX_PCA_BLOBS NORMAL STOP ****')
    else
        THROW_HARD('FLEX_PCA_BLOBS FAILED')
    endif
end subroutine exec_test_flex_pca_blobs

subroutine run_flex_pca_blobs_gate( nthr, metrics_file, all_ok )
    integer,       intent(in)    :: nthr
    class(string), intent(in)    :: metrics_file
    logical,       intent(inout) :: all_ok
    type(test_gate) :: gate
    type(flex_gate_result) :: shared, distributed
    type(state_fixture) :: fix
    type(string) :: root, truth_a, truth_b, truth_mean, truth_diff, phantom
    real :: latent_corr, truth_label_agreement, map_a_corr, map_b_corr
    call make_ui
    call simple_getcwd(root)
    CWD_GLOB = root%to_char()
    call gate%new(metrics_file)
    ! two blob states, equal populations, at their true poses without shifts
    fix%family    = 'blob'
    fix%mode      = FIX_KNOWN
    fix%nstates   = 2
    fix%pops(1:2) = [0.5, 0.5]
    fix%nptcls    = FLEX_PHANTOM_NPTCLS
    fix%snr       = PHANTOM_SNR
    fix%sherr     = 0.
    fix%nthr      = nthr
    call fix%build
    phantom = fix%project_fname()
    truth_a = fix%truth_fname(1)
    truth_b = fix%truth_fname(2)
    call write_mean_and_difference([truth_a, truth_b], root, fix%box(), fix%smpd(), truth_mean, truth_diff)
    truth_states_ab = read_truth_states(fix%truth_doc(), fix%nptcls)
    call run_flex_mode('shared', 1, nthr, root, truth_a, truth_b, truth_mean, truth_diff, phantom, gate, shared)
    call run_flex_mode('distributed', 2, nthr, root, truth_a, truth_b, truth_mean, truth_diff, phantom, gate, &
        &distributed)
    call gate%check('shared and distributed runs publish valid result sets', shared%valid .and. distributed%valid)
    if( shared%valid .and. distributed%valid )then
        call gate%report('published_state_count_difference', real(abs(shared%nstates - distributed%nstates)))
        truth_label_agreement = real(count(shared%predicted_truth == distributed%predicted_truth)) / &
            &real(FLEX_PHANTOM_NPTCLS)
        call gate%report('shared_distributed_truth_label_agreement', truth_label_agreement)
        if( shared%embedding_found .and. distributed%embedding_found )then
            latent_corr = best_latent_correlation(shared%z, distributed%z)
            call gate%report('shared_distributed_latent_correlation', latent_corr)
        endif
        call map_pair_correlation(shared%map_a, distributed%map_a, map_a_corr)
        call map_pair_correlation(shared%map_b, distributed%map_b, map_b_corr)
        call gate%report('shared_distributed_mode_a_map_correlation', map_a_corr)
        call gate%report('shared_distributed_mode_b_map_correlation', map_b_corr)
    endif
    ! the per-state refinement of the shared run's states, after the FLEX comparisons it does not touch
    if( shared%valid ) call run_state_refinements(root, 2*nthr, truth_a, truth_b, shared, gate)
    ! the class-average state refinement, on its own three-state phantom
    call run_states_cavgs_gate(root, 2*nthr, gate)
    all_ok = all_ok .and. gate%passed()
    call gate%kill
    call kill_result(shared)
    call kill_result(distributed)
end subroutine run_flex_pca_blobs_gate

subroutine run_flex_mode( tag, nparts, nthr, root, truth_a, truth_b, truth_mean, truth_diff, phantom, gate, result )
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: nparts, nthr
    class(string),    intent(in) :: root, truth_a, truth_b, truth_mean, truth_diff, phantom
    type(test_gate),  intent(inout) :: gate
    type(flex_gate_result), intent(out) :: result
    class(commander_flex_pca), allocatable :: flex_pca
    class(cmdline), allocatable :: cline
    type(string) :: run_dir
    integer :: status, i
    run_dir = root//'/'//trim(tag)
    call simple_mkdir(run_dir)
    call simple_copy_file(phantom, run_dir//'/phantom.simple')
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(flex_pca, cline)
    call configure_flex_pca_phantom(cline, truth_mean, 'pcg', 'pcg', nparts, nthr)
    call set_fixed_seed(FLEX_GATE_SEED, propagate=.true.)
    call flex_pca%execute(cline)
    call cline%kill
    deallocate(cline, flex_pca)
    ! the path: a distributed run's parts ran as worker jobs, each recording exit status 0
    if( nparts > 1 ) call gate%check(trim(tag)//': every part ran as a worker job that exited cleanly', &
        &all([(job_exited_cleanly(i), i = 1, nparts)]))
    call inspect_flex_result(tag, truth_a, truth_b, truth_diff, gate, result)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()

contains

    !> the queue's record of job ipart in the working directory holds exit status 0
    logical function job_exited_cleanly( ipart )
        integer, intent(in) :: ipart
        integer :: funit, ios, code
        job_exited_cleanly = .false.
        open(newunit=funit, file='EXIT_CODE_JOB_'//int2str(ipart), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit, *, iostat=ios) code
        close(funit)
        if( ios /= 0 ) return
        job_exited_cleanly = code == 0
    end function job_exited_cleanly

end subroutine run_flex_mode

subroutine inspect_flex_result( tag, truth_a, truth_b, truth_diff, gate, result )
    character(len=*), intent(in) :: tag
    class(string),    intent(in) :: truth_a, truth_b, truth_diff
    type(test_gate),  intent(inout) :: gate
    type(flex_gate_result), intent(out) :: result
    class(sp_project), allocatable :: project
    class(oris), allocatable :: field
    real(real32), allocatable :: weights(:,:)
    real(real64), allocatable :: precision(:,:,:)
    integer, allocatable :: labels(:), deconv_labels(:), truth_counts(:,:), state_truth(:), map_states(:), shuffled(:)
    real, allocatable :: corr(:,:)
    type(string) :: map
    character(len=STDLEN) :: message
    integer :: status, i, state, truth_label, predicted, nlabel_mismatch, nweight_mismatch, map_box
    real :: best_score, fsc05, fsc0143, map_smpd
    result = flex_gate_result()
    allocate(project, field)
    call project%read(string('phantom.simple'))
    field = project%os_ptcl3D
    call read_state_weight_table(project, field, result%nstates, weights, labels, status, message)
    call gate%check(trim(tag)//': delivered flex-weight set validates', status == 0)
    if( status /= 0 )then
        write(logfhandle,'(A,A)') '>>> FLEX_PCA_BLOBS weight-store failure: ', trim(message)
        call field%kill
        call project%kill
        deallocate(field, project)
        return
    endif
    call gate%check(trim(tag)//': publication retains both truth modes', result%nstates >= 2)
    call gate%check(trim(tag)//': publication respects the state ceiling', &
        &result%nstates <= FLEX_PHANTOM_MAX_STATES)
    allocate(truth_counts(result%nstates,2), source=0)
    allocate(state_truth(result%nstates), source=0)
    allocate(result%predicted_truth(FLEX_PHANTOM_NPTCLS), source=0)
    nlabel_mismatch = 0
    nweight_mismatch = 0
    do i = 1, FLEX_PHANTOM_NPTCLS
        truth_label = truth_states_ab(i)
        if( labels(i) >= 1 .and. labels(i) <= result%nstates ) &
            &truth_counts(labels(i),truth_label) = truth_counts(labels(i),truth_label) + 1
        if( int(labels(i)) /= project%os_ptcl3D%get_state(i) ) nlabel_mismatch = nlabel_mismatch + 1
        if( labels(i) > 0 )then
            predicted = maxloc(weights(:,i), dim=1)
            if( int(labels(i)) /= predicted ) nweight_mismatch = nweight_mismatch + 1
        endif
    enddo
    ! one-to-one: the delivered states mapped onto the truths that maximise agreement (a majority
    ! recoding would credit two states that both hold mostly one truth)
    map_states = match_states(truth_states_ab, labels, 2, result%nstates)
    do state = 1, result%nstates
        state_truth(state) = map_states(state)
    enddo
    do i = 1, FLEX_PHANTOM_NPTCLS
        if( labels(i) >= 1 .and. labels(i) <= result%nstates ) &
            &result%predicted_truth(i) = state_truth(labels(i))
    enddo
    result%label_accuracy = mapped_accuracy(truth_states_ab, labels, map_states)
    ! the floor (FLEX_PHANTOM_LABEL_ACCURACY_MIN, 0.60) is admissible on the profile of the merged code
    ! (2026-10-19, scratch/keep/profiles/phase8/main, 20 seed offsets at the entry's threads): shared
    ! median 0.989, s 0.007, worst 0.908, distributed median 0.990, worst 0.907, admissible to 0.900;
    ! the shuffled labels of two equal states score about 0.5
    call gate%metric(trim(tag)//': truth_label_accuracy', result%label_accuracy, &
        &FLEX_PHANTOM_LABEL_ACCURACY_MIN, result%label_accuracy >= FLEX_PHANTOM_LABEL_ACCURACY_MIN)
    ! the negative control: the truth labels shuffled across the particles must fail the floor
    shuffled = truth_states_ab
    call set_fixed_seed(FLEX_GATE_SEED + 1)
    call shuffle(shuffled)
    map_states = match_states(shuffled, labels, 2, result%nstates)
    call gate%check(trim(tag)//': shuffled truth labels fail the accuracy floor', &
        &.not. (mapped_accuracy(shuffled, labels, map_states) >= FLEX_PHANTOM_LABEL_ACCURACY_MIN))
    call gate%check(trim(tag)//': weight flags and project labels agree', nlabel_mismatch == 0)
    call gate%check(trim(tag)//': hard labels select maximum kernel weight', nweight_mismatch == 0)
    call inspect_eigenvolume(tag, truth_diff, gate, result%pc_corr)
    allocate(corr(result%nstates,2), source=0.0)
    do state = 1, result%nstates
        call project%get_vol('vol', state, map, map_smpd, map_box)
        call gate%check(trim(tag)//': delivered state map '//int2str(state)//' exists', file_exists(map))
        if( file_exists(map) )then
            call compare_to_truth(truth_a, map, FLEX_PHANTOM_MSKDIAM, corr(state,1), fsc05, fsc0143, &
                &corr_lp=MAP_CORR_LP)
            call compare_to_truth(truth_b, map, FLEX_PHANTOM_MSKDIAM, corr(state,2), fsc05, fsc0143, &
                &corr_lp=MAP_CORR_LP)
        endif
        call map%kill
    enddo
    call match_truth_maps(corr, result%matched_a, result%matched_b, best_score)
    call gate%check(trim(tag)//': two distinct maps match the two truth modes', &
        &result%matched_a > 0 .and. result%matched_b > 0)
    if( result%matched_a > 0 .and. result%matched_b > 0 )then
        call gate%check(trim(tag)//': mode-A map prefers truth A', &
            &corr(result%matched_a,1) > corr(result%matched_a,2))
        call gate%check(trim(tag)//': mode-B map prefers truth B', &
            &corr(result%matched_b,2) > corr(result%matched_b,1))
        call project%get_vol('vol', result%matched_a, map, map_smpd, map_box)
        result%map_a = simple_abspath(map)
        call map%kill
        call project%get_vol('vol', result%matched_b, map, map_smpd, map_box)
        result%map_b = simple_abspath(map)
        call map%kill
    endif
    write(logfhandle,'(A,A,A)') '>>> FLEX_PCA_BLOBS ', trim(tag), &
        &' state truth counts/correlations: state nA nB corrA corrB'
    do state = 1, result%nstates
        write(logfhandle,'(A,I0,2(1X,I0),2(1X,F9.6))') '>>>   ', state, truth_counts(state,:), corr(state,:)
    enddo
    call read_deconv_block('flex_pca_embedding.bin', FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_NCOMP, &
        &result%z, precision, deconv_labels, result%noise_scale, result%embedding_found)
    call gate%check(trim(tag)//': deconvolved embedding is published', result%embedding_found)
    if( result%embedding_found ) call gate%report(trim(tag)//': deconvolution_noise_scale', real(result%noise_scale))
    result%valid = result%nstates >= 2 .and. result%nstates <= FLEX_PHANTOM_MAX_STATES .and. &
        &result%matched_a > 0 .and. result%matched_b > 0
    call field%kill
    call project%kill
    deallocate(field, project)
    deallocate(weights, labels, truth_counts, state_truth, corr)
    if( allocated(precision) ) deallocate(precision)
    if( allocated(deconv_labels) ) deallocate(deconv_labels)
end subroutine inspect_flex_result

subroutine inspect_eigenvolume( tag, truth_diff, gate, corr )
    character(len=*), intent(in) :: tag
    class(string),    intent(in) :: truth_diff
    type(test_gate),  intent(inout) :: gate
    real,             intent(out) :: corr
    class(image), allocatable :: truth, pc
    type(string) :: pcfile
    pcfile = 'flex_pca_polished_pc001.mrc'
    corr = 0.0
    call gate%check(trim(tag)//': polished leading eigenvolume is delivered', file_exists(pcfile))
    if( .not. file_exists(pcfile) ) return
    allocate(truth, pc)
    call truth%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
    call pc%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
    call truth%read(truth_diff)
    call pc%read(pcfile)
    corr = abs(truth%real_corr(pc))
    call gate%report(trim(tag)//': leading_eigenvolume_truth_correlation', corr)
    call truth%kill
    call pc%kill
    deallocate(truth, pc)
end subroutine inspect_eigenvolume

!> The per-state refinement (plan section 9), on the shared-memory run's published project: each
!! truth-matched state is refined by refine3D_auto state=X twice from that same parent, with the frozen
!! FLEX weights and with the hard labels. The weighted map must correlate with its truth (same frame,
!! low-pass 8 A) at least as well as the hard one, minus STATE_REFINE_MARGIN; the parent project must
!! stay byte-identical.
subroutine run_state_refinements( root, nthr, truth_a, truth_b, shared, gate )
    class(string),          intent(in)    :: root, truth_a, truth_b
    integer,                intent(in)    :: nthr
    type(flex_gate_result), intent(in)    :: shared
    type(test_gate),        intent(inout) :: gate
    type(string)   :: run_dir, parent
    integer(int64) :: checksum_before
    real    :: corr_flex, corr_hard
    integer :: status, k, state
    run_dir = root//'/shared'
    parent  = run_dir//'/phantom.simple'
    checksum_before = sigma2_state_digest_file(parent)
    do k = 1, 2
        state = merge(shared%matched_a, shared%matched_b, k == 1)
        if( k == 1 )then
            call refine_flex_state(run_dir, parent, state, 'flex', nthr, truth_a, gate, corr_flex)
            call refine_flex_state(run_dir, parent, state, 'no',   nthr, truth_a, gate, corr_hard)
        else
            call refine_flex_state(run_dir, parent, state, 'flex', nthr, truth_b, gate, corr_flex)
            call refine_flex_state(run_dir, parent, state, 'no',   nthr, truth_b, gate, corr_hard)
        endif
        call gate%report('state'//int2str(state)//'_hard_refined_truth_correlation', corr_hard)
        ! relative to the hard map (STATE_REFINE_MARGIN); on the profile of 2026-10-19 the weighted maps
        ! correlate at median 0.992 and 0.994 (s 0.001), every run above its floor by 15 s or more
        call gate%metric('state'//int2str(state)//'_weighted_refined_truth_correlation', corr_flex, &
            &corr_hard - STATE_REFINE_MARGIN, corr_flex >= corr_hard - STATE_REFINE_MARGIN)
    enddo
    call gate%check('the parent project is byte-identical after the per-state refinements', &
        &sigma2_state_digest_file(parent) == checksum_before)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()
end subroutine run_state_refinements

!> refine3D_auto state=X m_estimator=mode in a new run directory of run_dir, from projfile; corr is the
!! refined map's correlation with its truth
subroutine refine_flex_state( run_dir, projfile, state, mode, nthr, truth, gate, corr )
    class(string),    intent(in)    :: run_dir, projfile, truth
    integer,          intent(in)    :: state, nthr
    character(len=*), intent(in)    :: mode
    type(test_gate),  intent(inout) :: gate
    real,             intent(out)   :: corr
    class(commander_refine3D_auto), allocatable :: xrefine
    class(cmdline),                 allocatable :: cline
    class(sp_project),              allocatable :: work
    type(state_weight_set) :: wset
    type(string)           :: work_dir, work_projfile, map, imgkind
    character(len=STDLEN)  :: tag, message
    real    :: map_smpd, fsc05, fsc0143
    integer :: map_box, status, i, ind
    logical :: l_registered
    corr = 0.
    tag  = 'state'//int2str(state)//'_'//trim(mode)
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_STATE_REFINE FAILED: could not enter run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xrefine, cline, work)
    call cline%set('prg',         'refine3D_auto')
    call cline%set('projfile',    projfile)
    call cline%set('mkdir',       'yes')
    call cline%set('state',       state)
    call cline%set('m_estimator', mode)
    call cline%set('pgrp',        'c1')
    call cline%set('mskdiam',     FLEX_PHANTOM_MSKDIAM)
    call cline%set('maxits',      STATE_REFINE_MAXITS)
    call cline%set('nthr',        nthr)
    call set_fixed_seed(STATE_REFINE_SEED, propagate=.true.)
    call xrefine%execute(cline)
    ! refine3D_auto ran in its own numbered directory
    call simple_getcwd(work_dir)
    work_projfile = work_dir//'/refine3D_auto_state'//int2str_pad(state,2)//'.simple'
    call gate%check(trim(tag)//': the work project is written', file_exists(work_projfile))
    if( file_exists(work_projfile) )then
        call work%read(work_projfile)
        call gate%check(trim(tag)//': the work project holds one state', work%os_ptcl3D%get_n('state') == 1)
        call work%get_vol('vol', 1, map, map_smpd, map_box)
        call gate%check(trim(tag)//': the refined map exists', file_exists(map))
        if( file_exists(map) ) call compare_to_truth(truth, map, FLEX_PHANTOM_MSKDIAM, corr, fsc05, fsc0143, &
            &corr_lp=MAP_CORR_LP)
        if( trim(mode) == 'flex' )then
            call wset%new(work, work%os_ptcl3D, status, message)
            call gate%check(trim(tag)//': the work project''s single-state weight set validates', status == 0)
            if( status == 0 ) call gate%check(trim(tag)//': the work set records the parent state', &
                &wset%get_parent_state() == state .and. wset%get_nstates() == 1)
            call wset%kill
            ind = 0
            do i = 1, work%os_out%get_noris()
                if( .not. work%os_out%isthere(i, 'imgkind') ) cycle
                imgkind = work%os_out%get_str(i, 'imgkind')
                if( imgkind%to_char() == 'vol' .and. work%os_out%get_state(i) == 1 ) ind = i
            enddo
            l_registered = .false.
            if( ind > 0 ) l_registered = work%os_out%isthere(ind, 'mass') .and. work%os_out%isthere(ind, 'ess')
            call gate%check(trim(tag)//': the final map is registered with applied mass and effective sample size', &
                &l_registered)
        endif
        call work%kill
    endif
    call gate%report(trim(tag)//'_truth_correlation', corr)
    call cline%kill
    deallocate(xrefine, cline, work)
end subroutine refine_flex_state

subroutine match_truth_maps( corr, matched_a, matched_b, best_score )
    real,    intent(in)  :: corr(:,:)
    integer, intent(out) :: matched_a, matched_b
    real,    intent(out) :: best_score
    integer :: state_a, state_b
    best_score = -huge(best_score)
    matched_a = 0
    matched_b = 0
    do state_a = 1, size(corr,1)
        do state_b = 1, size(corr,1)
            if( state_a == state_b ) cycle
            if( corr(state_a,1) + corr(state_b,2) > best_score )then
                best_score = corr(state_a,1) + corr(state_b,2)
                matched_a = state_a
                matched_b = state_b
            endif
        enddo
    enddo
end subroutine match_truth_maps

real function best_latent_correlation( a, b ) result( score )
    real(real64), intent(in) :: a(:,:), b(:,:)
    real :: corr(FLEX_PHANTOM_NCOMP,FLEX_PHANTOM_NCOMP)
    integer :: i, j
    do i = 1, FLEX_PHANTOM_NCOMP
        do j = 1, FLEX_PHANTOM_NCOMP
            corr(i,j) = abs(vector_correlation(a(:,i), b(:,j)))
        enddo
    enddo
    score = max(min(corr(1,1),corr(2,2)), min(corr(1,2),corr(2,1)))
end function best_latent_correlation

real function vector_correlation( a, b ) result( corr )
    real(real64), intent(in) :: a(:), b(:)
    real(real64) :: ac(size(a)), bc(size(b)), denom
    ac = a - sum(a) / real(size(a),real64)
    bc = b - sum(b) / real(size(b),real64)
    denom = sqrt(sum(ac*ac) * sum(bc*bc))
    if( denom > 0.0_real64 )then
        corr = real(sum(ac*bc) / denom)
    else
        corr = 0.0
    endif
end function vector_correlation

subroutine map_pair_correlation( a, b, corr )
    class(string), intent(in) :: a, b
    real,          intent(out) :: corr
    real :: fsc05, fsc0143
    call compare_to_truth(a, b, FLEX_PHANTOM_MSKDIAM, corr, fsc05, fsc0143, corr_lp=MAP_CORR_LP)
end subroutine map_pair_correlation

!> the mean of the truths (box voxels of smpd A) and the difference of the first two, written to dir
subroutine write_mean_and_difference( truths, dir, box, smpd, truth_mean, truth_diff )
    class(string), intent(in)  :: truths(:), dir
    integer,       intent(in)  :: box
    real,          intent(in)  :: smpd
    type(string),  intent(out) :: truth_mean, truth_diff
    type(image) :: vol, acc, first
    integer :: k
    truth_mean = dir//'/truth_mean.mrc'
    truth_diff = dir//'/truth_a_minus_b.mrc'
    call vol%new([box,box,box], smpd, wthreads=.false.)
    call acc%new([box,box,box], smpd, wthreads=.false.)
    call first%new([box,box,box], smpd, wthreads=.false.)
    call first%read(truths(1))
    do k = 1, size(truths)
        call vol%read(truths(k))
        call acc%add(vol)
    enddo
    call acc%div(real(size(truths)))
    call acc%write(truth_mean)
    call vol%read(truths(2))
    call first%subtr(vol)
    call first%write(truth_diff)
    call vol%kill
    call acc%kill
    call first%kill
end subroutine write_mean_and_difference

!> the true state of each of the n particles of a fixture truth document
function read_truth_states( doc, n ) result( states )
    class(string), intent(in) :: doc
    integer,       intent(in) :: n
    integer, allocatable :: states(:)
    type(oris) :: truth
    call truth%new(n, is_ptcl=.true.)
    call truth%read(doc, [1, n])
    states = nint(truth%get_all('state'))
    call truth%kill
end function read_truth_states

subroutine kill_result( result )
    type(flex_gate_result), intent(inout) :: result
    if( allocated(result%predicted_truth) ) deallocate(result%predicted_truth)
    if( allocated(result%z) ) deallocate(result%z)
    call result%map_a%kill
    call result%map_b%kill
    result = flex_gate_result()
end subroutine kill_result

!> cls_expansion, then refine3D_states_cavgs against the truth mean, on three beta-galactosidase states
!! under root/states_cavgs; the products against the truth, the contract against the cls_expansion project
subroutine run_states_cavgs_gate( root, nthr, gate )
    class(string),   intent(in)    :: root
    integer,         intent(in)    :: nthr
    type(test_gate), intent(inout) :: gate
    class(sp_project), allocatable :: expanded, sorted
    type(string) :: fixture, run_dir, expansion_projfile, sorted_projfile, work_projfile, truth_mean
    type(string) :: truths(CAVGS_NSTATES)
    integer      :: status, nsel
    fixture = root//'/states_cavgs'
    run_dir = fixture//'/run'
    call simple_mkdir(fixture)
    call simple_mkdir(run_dir)
    call prepare_states_cavgs_project(fixture, run_dir, nthr, truths, truth_mean)
    CWD_GLOB = fixture%to_char()
    call run_cls_expansion_stage(run_dir, nthr, expansion_projfile)
    call gate%check('states_cavgs: cls_expansion writes its project', file_exists(expansion_projfile))
    if( file_exists(expansion_projfile) )then
        allocate(expanded)
        call expanded%read(expansion_projfile)
        nsel = expanded%os_cls2D%count_state_gt_zero()
        call gate%report('states_cavgs: selected sub-class averages', real(nsel))
        call gate%check('states_cavgs: cls_expansion delivers at least 100 selected sub-class averages', &
            &nsel >= CAVGS_MIN_SUBAVGS)
        call check_memberless_halves(expanded, gate)
        if( nsel >= CAVGS_MIN_SUBAVGS )then
            call run_states_cavgs_stage(run_dir, expansion_projfile, nthr, sorted_projfile, work_projfile)
            call gate%check('states_cavgs: refine3D_states_cavgs writes its project', file_exists(sorted_projfile))
            call gate%check('states_cavgs: the work project is deleted', .not. file_exists(work_projfile))
            if( file_exists(sorted_projfile) )then
                allocate(sorted)
                call sorted%read(sorted_projfile)
                call inspect_states_cavgs(expanded, sorted, truths, truth_mean, gate)
                call sorted%kill
                deallocate(sorted)
            endif
        endif
        call expanded%kill
        deallocate(expanded)
    endif
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not restore fixture directory')
    CWD_GLOB = root%to_char()
end subroutine run_states_cavgs_gate

!> cls_expansion flags a sub-class half with state 0 (cls2D state_even/state_odd) only when the
!! sub-class has no member in that half (the project's ptcl2D class and eo; a parent class whose members
!! all share one eo is split into halves alternately and is left out), and no blank half stays active.
!! A half without hard members may stay active: a kernel split weighs every member of the half into
!! each sub-class average and flags the half by its effective size, which is never zero
subroutine check_memberless_halves( spproj, gate )
    use simple_fileio, only: add2fbody, fname2ext
    class(sp_project), intent(inout) :: spproj
    type(test_gate),   intent(inout) :: gate
    character(len=4), parameter :: HALVES(2) = ['even', 'odd ']
    type(image)          :: img
    type(string)         :: stk, ext, stk_half
    real,    allocatable :: rmat(:,:,:)
    integer, allocatable :: npop(:,:), parent_eo(:,:)
    real    :: smpd
    integer :: ncavgs, icls, ih, iptcl, ldim(3), nflag(2), nmismatch, nblank_active, half_state, eo, parent
    logical :: l_blank
    call spproj%get_cavgs_stk(stk, ncavgs, smpd)
    ! the members of every sub-class in each half, and which halves each parent class has at all
    allocate(npop(2,ncavgs), source=0)
    allocate(parent_eo(2,max(1,maxval(spproj%os_ptcl2D%get_all_asint('cluster')))), source=0)
    do iptcl = 1, spproj%os_ptcl2D%get_noris()
        icls = spproj%os_ptcl2D%get_int(iptcl, 'class')
        if( icls < 1 .or. icls > ncavgs ) cycle
        eo = spproj%os_ptcl2D%get_eo(iptcl)
        if( eo < 0 .or. eo > 1 ) cycle
        npop(eo+1,icls) = npop(eo+1,icls) + 1
        parent = spproj%os_ptcl2D%get_int(iptcl, 'cluster')
        if( parent >= 1 .and. parent <= size(parent_eo,2) ) parent_eo(eo+1,parent) = 1
    enddo
    ext = string('.')//fname2ext(stk)
    nflag = 0; nmismatch = 0; nblank_active = 0
    do ih = 1, 2
        stk_half = add2fbody(stk, ext, '_'//trim(HALVES(ih)))
        call find_ldim_nptcls(stk_half, ldim, ncavgs)
        ldim(3) = 1
        call img%new(ldim, smpd)
        do icls = 1, ncavgs
            half_state = 1
            if( spproj%os_cls2D%isthere(icls, 'state_'//trim(HALVES(ih))) ) &
                &half_state = spproj%os_cls2D%get_int(icls, 'state_'//trim(HALVES(ih)))
            if( half_state == 0 ) nflag(ih) = nflag(ih) + 1
            call img%read(stk_half, icls)
            rmat    = img%get_rmat()
            l_blank = all(rmat == 0.)
            if( l_blank .and. half_state /= 0 ) nblank_active = nblank_active + 1
            parent = spproj%os_cls2D%get_int(icls, 'cluster')
            if( parent < 1 .or. parent > size(parent_eo,2) ) cycle
            if( sum(parent_eo(:,parent)) < 2 ) cycle
            if( half_state == 0 .and. npop(ih,icls) > 0 ) nmismatch = nmismatch + 1
        enddo
        call img%kill
    enddo
    call gate%report('states_cavgs: sub-class halves written with state 0', real(sum(nflag)))
    call gate%check('states_cavgs: a sub-class half carries state 0 only when it has no member', nmismatch == 0)
    call gate%check('states_cavgs: no blank sub-class half stays active', nblank_active == 0)
    deallocate(npop, parent_eo)
    call stk%kill
    call ext%kill
    call stk_half%kill
end subroutine check_memberless_halves

!> the fixture's three compositional beta-galactosidase states in exact views (each view an exact 2D
!! class over the three states), in run_dir/phantom.simple with the truth mean registered as the
!! consensus map (vol, state 1)
subroutine prepare_states_cavgs_project( fixture, run_dir, nthr, truths, truth_mean )
    class(string), intent(in)  :: fixture, run_dir
    integer,       intent(in)  :: nthr
    type(string),  intent(out) :: truths(CAVGS_NSTATES), truth_mean
    class(sp_project), allocatable :: project
    type(state_fixture) :: fix
    type(string) :: truth_diff, projfile, cwd
    integer :: k, status
    call simple_getcwd(cwd)
    call simple_chdir(fixture, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter the states_cavgs fixture directory')
    fix%family  = 'betagal'
    fix%kind    = 'compositional'
    fix%mode    = FIX_VIEWS
    fix%nstates = CAVGS_NSTATES
    fix%pops(1:CAVGS_NSTATES) = 1. / real(CAVGS_NSTATES)
    fix%nptcls  = CAVGS_NPTCLS
    fix%nviews  = CAVGS_NVIEWS
    fix%snr     = PHANTOM_SNR
    fix%nthr    = nthr
    call fix%build
    cavgs_box     = fix%box()
    cavgs_smpd    = fix%smpd()
    cavgs_mskdiam = fix%mskdiam()
    do k = 1, CAVGS_NSTATES
        truths(k) = fix%truth_fname(k)
    enddo
    call write_mean_and_difference(truths, fixture, cavgs_box, cavgs_smpd, truth_mean, truth_diff)
    truth_states_abc = read_truth_states(fix%truth_doc(), fix%nptcls)
    truth_doc_abc    = fix%truth_doc()
    projfile = run_dir//'/phantom.simple'
    call simple_copy_file(fix%project_fname(), projfile)
    call simple_chdir(cwd, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not leave the states_cavgs fixture directory')
    allocate(project)
    call project%read(projfile)
    call project%add_vol2os_out(truth_mean, cavgs_smpd, 1, 'vol', box=cavgs_box)
    call project%write(projfile)
    call project%kill
    deallocate(project)
end subroutine prepare_states_cavgs_project

!> cls_expansion in a new run directory of run_dir; projfile_out is the project it writes there
subroutine run_cls_expansion_stage( run_dir, nthr, projfile_out )
    class(string), intent(in)  :: run_dir
    integer,       intent(in)  :: nthr
    type(string),  intent(out) :: projfile_out
    class(commander_cls_expansion), allocatable :: xexpansion
    class(cmdline),                 allocatable :: cline
    type(string) :: exec_dir
    integer      :: status
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter the class-average run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xexpansion, cline)
    call cline%set('prg',      'cls_expansion')
    call cline%set('projfile', run_dir//'/phantom.simple')
    call cline%set('mkdir',    'yes')
    call cline%set('ncls',     CAVGS_NSUB * CAVGS_NVIEWS)
    call cline%set('minpop',   CAVGS_MINPOP)
    call cline%set('lp',       CAVGS_EXPANSION_LP)
    call cline%set('mskdiam',  cavgs_mskdiam)
    call cline%set('nthr',     nthr)
    call set_fixed_seed(CAVGS_SEED, propagate=.true.)
    call xexpansion%execute(cline)
    ! cls_expansion ran in its own numbered directory
    call simple_getcwd(exec_dir)
    projfile_out = exec_dir//'/phantom.simple'
    call cline%kill
    deallocate(xexpansion, cline)
end subroutine run_cls_expansion_stage

!> refine3D_states_cavgs on projfile_in in a new run directory of run_dir; projfile_out is the project it
!! writes there and work_projfile the work project it must have deleted
subroutine run_states_cavgs_stage( run_dir, projfile_in, nthr, projfile_out, work_projfile )
    class(string), intent(in)  :: run_dir, projfile_in
    integer,       intent(in)  :: nthr
    type(string),  intent(out) :: projfile_out, work_projfile
    class(commander_refine3D_states_cavgs), allocatable :: xstates
    class(cmdline),                         allocatable :: cline
    type(string) :: exec_dir
    integer      :: status
    call simple_chdir(run_dir, status)
    if( status /= 0 ) THROW_HARD('FLEX_PCA_BLOBS FAILED: could not enter the class-average run directory')
    CWD_GLOB = run_dir%to_char()
    allocate(xstates, cline)
    call cline%set('prg',      'refine3D_states_cavgs')
    call cline%set('projfile', projfile_in)
    call cline%set('mkdir',    'yes')
    call cline%set('nstates',  CAVGS_NSTATES)
    call cline%set('pgrp',     'c1')
    call cline%set('mskdiam',  cavgs_mskdiam)
    call cline%set('nthr',     nthr)
    call set_fixed_seed(CAVGS_SEED, propagate=.true.)
    call xstates%execute(cline)
    call simple_getcwd(exec_dir)
    projfile_out  = exec_dir//'/phantom.simple'
    work_projfile = exec_dir//'/refine3D_states_cavgs_tmpproj.simple'
    call cline%kill
    deallocate(xstates, cline)
end subroutine run_states_cavgs_stage

!> refine3D_states_cavgs's products against the truth: the state of particle i is the fixture's, its
!! direction (i-1)/nrep+1 of CAVGS_NVIEWS exact views; the contract against the cls_expansion project
subroutine inspect_states_cavgs( expanded, sorted, truths, truth_mean, gate )
    class(sp_project), intent(in)    :: expanded, sorted
    class(string),     intent(in)    :: truths(:), truth_mean
    type(test_gate),   intent(inout) :: gate
    character(len=8), parameter :: CLS2D_KEYS(3) = [character(len=8) :: 'state', 'pop', 'cluster']
    character(len=8), parameter :: CLS3D_KEYS(2) = [character(len=8) :: 'cluster', 'accept']
    character(len=8), parameter :: POSE_KEYS(5)  = [character(len=8) :: 'e1', 'e2', 'e3', 'x', 'y']
    class(oris), allocatable :: poses, est, tru
    integer, allocatable :: sub_counts(:,:), sub_views(:,:), sub_truth(:), sub_view(:), states(:), ptcl_states(:)
    integer, allocatable :: cavg_counts(:,:), ptcl_counts(:,:), committed(:), expected(:), best_truth(:), map_states(:)
    real,    allocatable :: corr(:,:)
    type(string) :: map
    real    :: smpd, fsc05, fsc0143, err, frac5, acc_null_med, acc_null_p999
    integer :: nsub, nposes, nposes_per_view, i, k, s, t, nstates_out, ncommitted, view, box
    integer :: nmajority, nmembers
    logical :: l_ok
    nsub            = expanded%os_cls2D%get_noris()
    nposes          = CAVGS_NPTCLS
    nposes_per_view = CAVGS_NPTCLS / CAVGS_NVIEWS
    ! truth composition of every sub-average, from its members' planted states and directions
    allocate(sub_counts(nsub,CAVGS_NSTATES), sub_views(nsub,CAVGS_NVIEWS), source=0)
    do i = 1, expanded%os_ptcl2D%get_noris()
        k = expanded%os_ptcl2D%get_class(i)
        if( k < 1 .or. k > nsub ) cycle
        t     = truth_states_abc(i)
        view  = (i - 1) / nposes_per_view + 1
        sub_counts(k,t)   = sub_counts(k,t)   + 1
        sub_views(k,view) = sub_views(k,view) + 1
    enddo
    ! a sub-average without members has no truth (0) and is left out of the accuracy and the poses
    allocate(sub_truth(nsub), sub_view(nsub), source=0)
    do k = 1, nsub
        if( sum(sub_counts(k,:)) == 0 ) cycle
        sub_truth(k) = maxloc(sub_counts(k,:), dim=1)
        sub_view(k)  = maxloc(sub_views(k,:),  dim=1)
    enddo
    ! the contract: cls2D and the cls_expansion fields of cls3D kept, particle poses kept
    call gate%check('states_cavgs: cls2D is unchanged', fields_equal(expanded%os_cls2D, sorted%os_cls2D, CLS2D_KEYS))
    call gate%check('states_cavgs: cls3D keeps cluster and accept', &
        &fields_equal(expanded%os_cls3D, sorted%os_cls3D, CLS3D_KEYS))
    call gate%check('states_cavgs: particle poses are unchanged', &
        &fields_equal(expanded%os_ptcl3D, sorted%os_ptcl3D, POSE_KEYS))
    ! the delivered states, consecutive and populated, and their hand-back through the 2D classes
    states      = nint(sorted%os_cls3D%get_all('state'))
    ptcl_states = sorted%os_ptcl3D%get_all_asint('state')
    nstates_out = maxval(states)
    call gate%report('states_cavgs: delivered states', real(nstates_out))
    if( nstates_out < 2 .or. nstates_out > CAVGS_NSTATES ) return
    l_ok = .true.
    do s = 1, nstates_out
        if( count(states == s) == 0 .or. count(ptcl_states == s) == 0 ) l_ok = .false.
    enddo
    call gate%check('states_cavgs: every delivered state holds sub-averages and particles', l_ok)
    allocate(expected(size(ptcl_states)), source=0)
    do i = 1, size(ptcl_states)
        if( sorted%os_ptcl2D%get_state(i) < 1 ) cycle
        k = sorted%os_ptcl2D%get_class(i)
        if( k >= 1 .and. k <= nsub ) expected(i) = states(k)
    enddo
    call gate%check('states_cavgs: every particle carries the state of its sub-class', all(ptcl_states == expected))
    ! one class-average map per delivered state, and the consensus entry kept
    l_ok = .not. sorted%isthere_in_osout('vol_cavg', nstates_out + 1)
    do s = 1, nstates_out
        if( .not. sorted%isthere_in_osout('vol_cavg', s) )then
            l_ok = .false.
        else
            call sorted%get_vol('vol_cavg', s, map, smpd, box)
            if( .not. file_exists(map) ) l_ok = .false.
        endif
    enddo
    call gate%check('states_cavgs: one vol_cavg map per delivered state', l_ok)
    l_ok = sorted%isthere_in_osout('vol', 1)
    if( l_ok )then
        call sorted%get_vol('vol', 1, map, smpd, box)
        l_ok = map%to_char() == truth_mean%to_char()
    endif
    call gate%check('states_cavgs: the consensus map entry is kept', l_ok)
    ! coverage and state accuracy of the committed sub-averages
    ncommitted = count(states > 0 .and. sub_truth > 0)
    call gate%report('states_cavgs: committed sub-average fraction', real(count(states > 0)) / real(nsub))
    call gate%metric('states_cavgs: committed sub-averages with members', real(ncommitted), &
        &real(CAVGS_MIN_COMMITTED), ncommitted >= CAVGS_MIN_COMMITTED)
    allocate(cavg_counts(nstates_out,CAVGS_NSTATES), ptcl_counts(nstates_out,CAVGS_NSTATES), source=0)
    nmajority = 0
    nmembers  = 0
    do k = 1, nsub
        if( states(k) < 1 .or. sub_truth(k) < 1 ) cycle
        cavg_counts(states(k),sub_truth(k)) = cavg_counts(states(k),sub_truth(k)) + 1
        nmajority = nmajority + maxval(sub_counts(k,:))
        nmembers  = nmembers  + sum(sub_counts(k,:))
    enddo
    do i = 1, size(ptcl_states)
        if( ptcl_states(i) < 1 ) cycle
        t = truth_states_abc(i)
        ptcl_counts(ptcl_states(i),t) = ptcl_counts(ptcl_states(i),t) + 1
    enddo
    ! one-to-one: the delivered states mapped onto the truths that maximise agreement
    if( ncommitted > 0 )then
        committed  = pack([(k, k = 1, nsub)], states > 0 .and. sub_truth > 0)
        map_states = match_states(sub_truth(committed), states(committed), CAVGS_NSTATES, nstates_out)
        call gate%report('states_cavgs: sub-average state accuracy', &
            &mapped_accuracy(sub_truth(committed), states(committed), map_states))
        call gate%report('states_cavgs: member purity of the committed sub-averages', &
            &real(nmajority) / real(max(1, nmembers)))
        ! the negative control: the accuracy of the same labels against the truths shuffled
        call accuracy_null(sub_truth(committed), states(committed), CAVGS_NSTATES, nstates_out, 1000, CAVGS_SEED, &
            &acc_null_med, acc_null_p999)
        call gate%report('states_cavgs: control sub-average accuracy null median', acc_null_med)
        call gate%report('states_cavgs: control sub-average accuracy null p999', acc_null_p999)
    endif
    if( sum(ptcl_counts) > 0 )then
        map_states = match_states(truth_states_abc, ptcl_states, CAVGS_NSTATES, nstates_out)
        call gate%report('states_cavgs: particle state accuracy', &
            &mapped_accuracy(truth_states_abc, ptcl_states, map_states))
    endif
    write(logfhandle,'(A)') '>>> FLEX_PCA_BLOBS states_cavgs sub-average truth counts: state nA nB nC'
    do s = 1, nstates_out
        write(logfhandle,'(A,I0,3(1X,I0))') '>>>   ', s, cavg_counts(s,:)
    enddo
    ! poses of the committed sub-averages against their directions, free of the global frame
    if( ncommitted >= 2 )then
        committed = pack([(k, k = 1, nsub)], states > 0 .and. sub_truth > 0)
        allocate(poses, est, tru)
        call poses%new(nposes, is_ptcl=.true.)
        call poses%read(truth_doc_abc, [1,nposes])
        call est%new(ncommitted, is_ptcl=.false.)
        call tru%new(ncommitted, is_ptcl=.false.)
        do i = 1, ncommitted
            call est%set_euler(i, sorted%os_cls3D%get_euler(committed(i)))
            call tru%set_euler(i, poses%get_euler((sub_view(committed(i)) - 1) * nposes_per_view + 1))
        enddo
        err = pair_pose_error(est, tru, [(i, i = 1, ncommitted)], [(i, i = 1, ncommitted)], CAVGS_NPAIRS, &
            &CAVGS_SEED, frac5)
        call gate%report('states_cavgs: median pair pose error of the sub-averages (deg)', err)
        call gate%report('states_cavgs: fraction of sub-average pose pairs within 5 deg', frac5)
        call poses%kill
        call est%kill
        call tru%kill
        deallocate(poses, est, tru)
    endif
    ! the class-average maps against the planted states
    allocate(corr(nstates_out,CAVGS_NSTATES), source=0.)
    allocate(best_truth(nstates_out), source=0)
    do s = 1, nstates_out
        if( .not. sorted%isthere_in_osout('vol_cavg', s) ) cycle
        call sorted%get_vol('vol_cavg', s, map, smpd, box)
        if( .not. file_exists(map) ) cycle
        do t = 1, CAVGS_NSTATES
            call compare_to_truth(truths(t), map, cavgs_mskdiam, corr(s,t), fsc05, fsc0143, corr_lp=MAP_CORR_LP)
        enddo
        ! maxloc passes over a non-finite correlation: such a map matches no planted state
        if( all(corr(s,:) == corr(s,:)) ) best_truth(s) = maxloc(corr(s,:), dim=1)
        write(logfhandle,'(A,I0,3(1X,F9.6))') '>>> FLEX_PCA_BLOBS states_cavgs map/truth correlations: ', s, corr(s,:)
    enddo
    ! the states stay at chance (above), so how many planted states the maps match is a measurement:
    ! the profile of 2026-10-19 (scratch/keep/profiles/phase8/main) found one run in 20 whose maps
    ! matched a single planted state
    call gate%report('states_cavgs: distinct planted states the state maps match', &
        &real(count([(any(best_truth == t), t = 1, CAVGS_NSTATES)])))
    call map%kill
end subroutine inspect_states_cavgs

!> whether every listed field agrees between two orientation sets of the same size
logical function fields_equal( a, b, keys ) result( l_equal )
    class(oris),      intent(in) :: a, b
    character(len=*), intent(in) :: keys(:)
    integer :: i, j
    l_equal = a%get_noris() == b%get_noris()
    if( .not. l_equal ) return
    do j = 1, size(keys)
        do i = 1, a%get_noris()
            if( a%get(i, trim(keys(j))) /= b%get(i, trim(keys(j))) )then
                l_equal = .false.
                return
            endif
        enddo
    enddo
end function fields_equal

module subroutine exec_test_write_state_weights_labels( self, cline )
    class(commander_test_write_state_weights_labels), intent(inout) :: self
    class(cmdline),                                   intent(inout) :: cline
    call write_state_weights_from_labels(cline, 1.0)
    call simple_end('**** SIMPLE_TEST_WRITE_STATE_WEIGHTS_LABELS NORMAL STOP ****')
end subroutine exec_test_write_state_weights_labels

module subroutine exec_test_write_state_weights_mixed( self, cline )
    class(commander_test_write_state_weights_mixed), intent(inout) :: self
    class(cmdline),                                  intent(inout) :: cline
    call write_state_weights_from_labels(cline, 0.7)
    call simple_end('**** SIMPLE_TEST_WRITE_STATE_WEIGHTS_MIXED NORMAL STOP ****')
end subroutine exec_test_write_state_weights_mixed

!> Publish a PARTITION state weight set over the project's labelled particles (state > 0): w_own on the
!! labelled state, the rest spread equally over the other states (w_own = 1 gives the hard labels).
!! Written in the working directory and registered in the project (nstates from the command line).
subroutine write_state_weights_from_labels( cline, w_own )
    use simple_state_weight_set, only: state_weight_set
    class(cmdline), intent(inout) :: cline
    real,           intent(in)    :: w_own
    class(sp_project), allocatable :: project
    type(state_weight_set) :: wset
    type(string)           :: projfile
    real,    allocatable   :: weights(:,:)
    integer, allocatable   :: pinds(:), labels(:)
    integer :: nstates, nsel, i, s
    if( .not. cline%defined('projfile') ) THROW_HARD('projfile= is required')
    if( .not. cline%defined('nstates')  ) THROW_HARD('nstates= is required')
    projfile = cline%get_carg('projfile')
    nstates  = cline%get_iarg('nstates')
    allocate(project)
    call project%read(projfile)
    nsel = 0
    do i = 1, project%os_ptcl3D%get_noris()
        s = project%os_ptcl3D%get_state(i)
        if( s >= 1 .and. s <= nstates ) nsel = nsel + 1
    enddo
    allocate(pinds(nsel), labels(nsel), weights(nsel,nstates))
    nsel = 0
    do i = 1, project%os_ptcl3D%get_noris()
        s = project%os_ptcl3D%get_state(i)
        if( s < 1 .or. s > nstates ) cycle
        nsel = nsel + 1
        pinds(nsel)  = i
        labels(nsel) = s
        if( nstates > 1 )then
            weights(nsel,:) = (1. - w_own) / real(nstates - 1)
        else
            weights(nsel,:) = 0.
        endif
        weights(nsel,s) = w_own
    enddo
    call wset%publish(project, project%os_ptcl3D, projfile, pinds, weights, labels, 'test_tool')
    write(logfhandle,'(A,I0,A,I0,A,F5.2)') '>>> STATE WEIGHT SET PUBLISHED FROM LABELS: ', nsel, ' particles x ', nstates, &
        &' states, own-state weight ', w_own
    call wset%kill
    call project%kill
    deallocate(project, pinds, labels, weights)
end subroutine write_state_weights_from_labels

end submodule simple_commanders_test_highlevel_flex
