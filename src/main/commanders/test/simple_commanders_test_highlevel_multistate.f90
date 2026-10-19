!@descr: the multi-state assignment test of simple_commanders_test_highlevel: refine3D with three states against the known states of a blob phantom
! refine3D nstates=3 (refine=prob, and refine=prob_neigh) on the blob family of the state fixture in
! its search mode (random poses and shifts, populations 50, 30 and 20 %), from the three truths
! low-passed to 20 A and the uniformly random initial states refine3D draws, scored by correlation
! (objfun=cc). The states, the maps and the poses are
! judged against the truth: one mapping of estimated onto true states (match_states), the adjusted
! Rand index, the accuracy against its null and the collapsed predictor, the populations, the map
! mapping and the discrimination of the maps, the poses of the correctly assigned particles. Every
! run also scores the shuffled-label and duplicated-maps controls; suite=identical runs the
! identical-states control (labels and the shuffled-label control only). The starting volumes fix frame and hand, so the maps and poses are
! compared in the truths' frame without docking.
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_multistate
use simple_test_gate,          only: test_gate
use simple_test_state_fixture, only: state_fixture, FIX_SEARCH
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSTATES  = 3
integer, parameter :: NPTCLS   = 1500
integer, parameter :: MAXITS   = 6
integer, parameter :: NSHUFFLE = 1000
real,    parameter :: POPS(NSTATES) = [0.5, 0.3, 0.2]
real,    parameter :: SNR      = 2.0
real,    parameter :: LP_START = 20.  !< A, the low-pass of the starting volumes
real,    parameter :: LP       = 8.   !< A, the matching low-pass of refine3D
real,    parameter :: TRS      = 4.   !< pixels, the shift search range (the truth is within 2.5)
real,    parameter :: DIFF_FRAC   = 0.1 !< the difference mask: where two truths differ by this fraction of the larger maximum,
integer, parameter :: DIFF_DILATE = 1   !< dilated by this many pixels
! The floors (ruling of 2026-10-17: a merge of two states fails; the Phase 3 floors 0.30 and 0.63
! passed one). Under euclid from the noise-model-free starting maps two true states merged in most
! runs; under objfun=cc every state is resolved in all 20 runs of the profile of 2026-10-18
! (keep/profiles/phase5/after, 8 threads), by the floor rule (s = 1.4826 MAD):
! - adjusted Rand index: prob median 0.941, s 0.010, worst 0.928, admissible to 0.912 (floor 0.90);
!   prob_neigh median 0.891, s 0.022, worst 0.849, admissible to 0.824 (floor 0.80); a merge of two
!   states scores 0.45-0.76 and shuffled labels 0.00 (s 0.003);
! - accuracy: prob median 0.978, s 0.004, worst 0.972, to 0.965 (floor 0.95); prob_neigh median
!   0.959, s 0.008, worst 0.943, to 0.934 (floor 0.92); the null's 99.9th percentile is 0.42;
! - smallest per-state recall: prob median 0.969, s 0.008, worst 0.953, to 0.946 (floor 0.93);
!   prob_neigh median 0.949, s 0.010, worst 0.893, to 0.884 (floor 0.85); a merged state's is 0.02;
! - largest population error (ceiling): prob median 0.009, s 0.004, worst 0.013, from 0.021
!   (ceiling 0.03); prob_neigh median 0.015, s 0.007, worst 0.025, from 0.036 (ceiling 0.04);
! - discrimination, smallest margin: prob median 0.345, worst 0.334, to 0.328; prob_neigh median
!   0.355, worst 0.343, to 0.338 (floor 0.30 for both); copies of one map score -0.02;
! - every estimated state populated and the map mapping agreeing with the label mapping: checks.
! The poses and shifts of the correctly assigned particles (ceilings, from the Phase 3 profile under
! euclid, prob_neigh: rotation median 4.75 deg, s 0.37, worst 5.56, admissible from 5.93; shift t1
! median 0.266 px, s 0.038, worst 0.307, from 0.380; under cc now 3.2 deg and 0.26 px); random poses
! err by about 90 deg and zero shifts by about 2 px. The fits are reported; a reported value that is
! not finite (an empty or a duplicated map) still fails its row.
real,    parameter :: MAX_MED_ROT_ERR = 6.0
real,    parameter :: MAX_MED_SH_ERR  = 0.40
real,    parameter :: MAX_POP_ERR_PROB = 0.03, MAX_POP_ERR_NEIGH = 0.04
real,    parameter :: MIN_DISCR_MARGIN = 0.30

contains

!> simple_test_exec test=multistate_assignment [suite=prob|prob_neigh|identical|list] [nthr=<n>]
module subroutine exec_test_multistate_assignment( self, cline )
    use simple_ui, only: make_ui
    class(commander_test_multistate_assignment), intent(inout) :: self
    class(cmdline),                              intent(inout) :: cline
    class(parameters), allocatable :: params
    type(test_gate) :: gate
    type(string)    :: cwd_saved, root, suite
    integer :: nthr, status
    logical :: all_ok
    suite = ''
    if( cline%defined('suite') ) suite = cline%get_carg('suite')
    select case(suite%to_char())
        case('list')
            write(logfhandle,'(a)') 'Available suites for multistate_assignment: prob, prob_neigh, identical'
            return
        case('', 'prob', 'prob_neigh', 'identical')
        case default
            THROW_HARD('no sub-suite '//suite%to_char()//' in multistate_assignment; use suite=list')
    end select
    call make_ui
    call simple_getcwd(cwd_saved)
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))
    root = filepath(cwd_saved, 'test_multistate_assignment_'//int2str(get_process_id()))
    if( dir_exists(root) ) call simple_rmdir(root)
    call simple_mkdir(root)
    select case(suite%to_char())
        case('')
            call run_case(root, 'prob',       'prob',       'distinct', nthr, gate)
            call run_case(root, 'prob_neigh', 'prob_neigh', 'distinct', nthr, gate)
        case('prob', 'prob_neigh')
            call run_case(root, suite%to_char(), suite%to_char(), 'distinct', nthr, gate)
        case('identical')
            call run_case(root, 'identical', 'prob', 'identical', nthr, gate)
    end select
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('could not return to '//cwd_saved%to_char())
    CWD_GLOB = cwd_saved%to_char()
    all_ok = gate%passed()
    call gate%kill
    if( all_ok )then
        call simple_rmdir(root)
        write(logfhandle,'(a)') 'PASS: refine3D recovers the planted states of the blob phantom'
        call simple_end('**** SIMPLE_TEST_MULTISTATE_ASSIGNMENT NORMAL STOP ****')
    else
        THROW_HARD('TEST_MULTISTATE_ASSIGNMENT FAILED')
    endif
end subroutine exec_test_multistate_assignment

!> one refine3D run in root/<tag> and its metrics, named after tag
subroutine run_case( root, tag, refine, kind, nthr, gate )
    use simple_commanders_refine3D, only: commander_refine3D
    use simple_refine3D_fnames,    only: refine3D_state_vol_fname
    use simple_sp_project,         only: sp_project
    use simple_test_utils,         only: set_fixed_seed
    use simple_rnd,                only: shuffle
    use simple_test_truth_metrics, only: match_states, adjusted_rand_index, mapped_accuracy, accuracy_null, &
        &collapsed_accuracy, population_errors, match_maps, discrimination, group_errors, group_operators, compare_to_truth, &
        &finite_min, offdiag, maps_distinct, state_recalls
    class(string),    intent(in)    :: root
    character(len=*), intent(in)    :: tag, refine, kind
    integer,          intent(in)    :: nthr
    type(test_gate),  intent(inout) :: gate
    type(state_fixture) :: fix
    type(commander_refine3D) :: xrefine3D
    class(cmdline),    allocatable :: cl
    class(sp_project), allocatable :: proj
    type(oris)   :: truth
    type(image)  :: est(NSTATES), truths(NSTATES), ordered(NSTATES), mean
    type(string) :: dir, projfile, fname
    integer, allocatable :: tstates(:), estates(:), map(:), mapmaps(:), shuffled(:), mapctrl(:)
    real,    allocatable :: errs(:), ops(:,:,:), re(:,:,:), rt(:,:,:), sherrs(:)
    logical, allocatable :: right(:)
    real    :: ari, acc, med_null, p999, margins(NSTATES), diffcorr(NSTATES,NSTATES), cc(NSTATES,NSTATES), sh(2)
    real    :: corr, fsc05, fsc0143, poperr, ari_shuf, acc_shuf, recall, margin
    real    :: min_ari, min_acc, min_recall, max_poperr
    integer :: k, l, i, e, nmask, status, ldim(3), n
    logical :: l_maps, l_control
    l_control = kind == 'identical'
    ! the floors of the mode (header)
    if( refine == 'prob_neigh' )then
        min_ari    = 0.80
        min_acc    = 0.92
        min_recall = 0.85
        max_poperr = MAX_POP_ERR_NEIGH
    else
        min_ari    = 0.90
        min_acc    = 0.95
        min_recall = 0.93
        max_poperr = MAX_POP_ERR_PROB
    endif
    write(logfhandle,'(a)') '>>> TEST_MULTISTATE_ASSIGNMENT: '//tag
    dir = root//'/'//tag
    call simple_mkdir(dir)
    call simple_chdir(dir, status)
    if( status /= 0 ) THROW_HARD('could not enter '//dir%to_char())
    CWD_GLOB = dir%to_char()
    ! the dataset, from the cache when it is there
    fix%family    = 'blob'
    fix%kind      = kind
    fix%mode      = FIX_SEARCH
    fix%nstates   = NSTATES
    fix%pops(1:NSTATES) = POPS
    fix%nptcls    = NPTCLS
    fix%snr       = SNR
    fix%nthr      = nthr
    call fix%build
    projfile = fix%project_fname()
    ! the project holds every particle in state 1: refine3D draws the initial states itself
    allocate(proj)
    ! the starting volumes: the truths low-passed to 20 A
    allocate(cl)
    do k = 1, NSTATES
        fname = 'startvol_truth'//int2str(k)//'.mrc'
        call fix%write_lowpassed(k, LP_START, fname)
        call cl%set('vol'//int2str(k), fname)
    enddo
    call cl%set('prg',      'refine3D')
    call cl%set('projfile', projfile)
    call cl%set('mkdir',    'no')
    call cl%set('nstates',  NSTATES)
    call cl%set('refine',   refine)
    call cl%set('pgrp',     'c1')
    call cl%set('mskdiam',  fix%mskdiam())
    call cl%set('lp',       LP)
    call cl%set('trs',      TRS)
    call cl%set('maxits',   MAXITS)
    call cl%set('nspace',   2500)
    if( refine == 'prob_neigh' ) call cl%set('nspace_sub', 500)
    ! correlation, as production's state assignment scores it (simple_state_assignment): the
    ! starting maps carry no noise model, and under euclid two true states merged in most runs
    call cl%set('objfun',   'cc')
    call cl%set('nthr',     nthr)
    call xrefine3D%execute(cl)
    call cl%kill
    deallocate(cl)
    ! ---- labels ----
    call proj%read(projfile)
    call truth%new(NPTCLS, is_ptcl=.true.)
    call truth%read(fix%truth_doc(), [1, NPTCLS])
    tstates = nint(truth%get_all('state'))
    estates = nint(proj%os_ptcl3D%get_all('state'))
    map     = match_states(tstates, estates, NSTATES, NSTATES)
    acc     = mapped_accuracy(tstates, estates, map)
    ari     = adjusted_rand_index(tstates, estates)
    call accuracy_null(tstates, estates, NSTATES, NSTATES, NSHUFFLE, 20261018, med_null, p999)
    poperr = maxval(population_errors(tstates, estates, map, NSTATES))
    recall = minval(state_recalls(tstates, estates, map, NSTATES))
    if( l_control )then
        call gate%report(tag//'_ari', ari)
        call gate%report(tag//'_accuracy', acc)
    else
        call gate%metric(tag//'_ari', ari, min_ari, ari >= min_ari)
        call gate%metric(tag//'_accuracy', acc, min_acc, acc >= min_acc)
    endif
    call gate%report(tag//'_accuracy_null_median', med_null)
    call gate%report(tag//'_accuracy_null_p999', p999)
    call gate%report(tag//'_collapsed_accuracy', collapsed_accuracy(tstates, NSTATES))
    ! a merge of two true states: one true state unrecalled, one estimated state empty
    if( l_control )then
        call gate%report(tag//'_largest_population_error', poperr)
        call gate%report(tag//'_min_state_recall', recall)
    else
        call gate%metric(tag//'_largest_population_error', poperr, max_poperr, poperr <= max_poperr, upper=.true.)
        call gate%metric(tag//'_min_state_recall', recall, min_recall, recall >= min_recall)
        call gate%check(tag//'_every_estimated_state_populated', all([(count(estates == k) > 0, k = 1, NSTATES)]))
    endif
    ! the shuffled-label control: truth labels shuffled across particles, the counts kept
    call set_fixed_seed(20261019)
    shuffled = tstates
    call shuffle(shuffled)
    mapctrl = match_states(shuffled, estates, NSTATES, NSTATES)
    ari_shuf = adjusted_rand_index(shuffled, estates)
    acc_shuf = mapped_accuracy(shuffled, estates, mapctrl)
    call gate%report(tag//'_control_shuffled_ari', ari_shuf)
    ! (one draw of the null, so inside it unless the draw is beyond its 99.9th percentile)
    call gate%check(tag//'_control_shuffled_accuracy_within_null', acc_shuf <= p999)
    ! ---- maps, in the truths' frame: the starting volumes, the low-passed truths, fix frame and hand ----
    ldim = [fix%box(), fix%box(), fix%box()]
    l_maps = .true.
    do k = 1, NSTATES
        call truths(k)%new(ldim, fix%smpd())
        call truths(k)%read(fix%truth_fname(k))
    enddo
    do e = 1, NSTATES
        call est(e)%new(ldim, fix%smpd())
        fname = refine3D_state_vol_fname(e)
        if( file_exists(fname) )then
            call est(e)%read(fname)
        else
            l_maps = .false.
        endif
    enddo
    write(logfhandle,'(a)') 'CONFUSION (rows: estimated states; columns: true states)'
    do e = 1, NSTATES
        write(logfhandle,'(8i8)') (count(estates == e .and. tstates == k), k = 1, NSTATES)
    enddo
    ! the map mapping on its own, from the correlations inside the difference mask of the truths
    call discrimination(est, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask, cc)
    write(logfhandle,'(a,i8)') 'MASKED MAP CORRELATIONS (rows: estimated maps; columns: truths); mask voxels: ', nmask
    do e = 1, NSTATES
        write(logfhandle,'(8f8.3)') cc(e,:)
    enddo
    mapmaps = match_maps(cc)
    if( l_control )then
        call gate%report(tag//'_difference_mask_voxels', real(nmask))
    else
        call gate%check(tag//'_every_state_has_a_map', l_maps)
        call gate%check(tag//'_map_mapping_agrees_with_label_mapping', all(mapmaps == map))
    endif
    ! the fit of every mapped map to its truth
    do e = 1, NSTATES
        if( map(e) < 1 ) cycle
        call compare_to_truth(fix%truth_fname(map(e)), refine3D_state_vol_fname(e), fix%mskdiam(), corr, fsc05, fsc0143)
        call gate%report(tag//'_map_corr_state'//int2str(map(e)), corr)
        call gate%report(tag//'_map_fsc0143_state'//int2str(map(e)), fsc0143)
    enddo
    ! the discrimination of the maps in the order of the truths they map onto, and the
    ! duplicated-maps control: every map replaced by the first one, then by their mean
    do k = 1, NSTATES
        e = findloc(map, k, dim=1)
        call ordered(k)%copy(est(max(1, e)))
        if( e < 1 ) call ordered(k)%zero_and_unflag_ft
    enddo
    call discrimination(ordered, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
    if( l_control )then
        ! identical truths leave no difference mask, so no correlation inside it is finite and the
        ! maps cannot count as distinct states whatever refine3D made: this row records that the
        ! discrimination has nothing to discriminate; the duplicated-maps controls of the prob and
        ! prob_neigh runs are the ones that test it
        call gate%check(tag//'_states_not_distinct', .not. maps_distinct(margins, diffcorr))
    else
        margin = finite_min(margins)
        call gate%metric(tag//'_discrimination_min_margin', margin, MIN_DISCR_MARGIN, margin >= MIN_DISCR_MARGIN)
        call gate%report(tag//'_discrimination_min_difference_corr', finite_min(offdiag(diffcorr)))
    endif
    if( .not. l_control )then
        call mean%copy(ordered(1))
        do k = 2, NSTATES
            call mean%add(ordered(k))
        enddo
        call mean%div(real(NSTATES))
        do k = 2, NSTATES
            call ordered(k)%copy(ordered(1))
        enddo
        call discrimination(ordered, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
        write(logfhandle,'(a,2f10.4)') 'COPIES OF ONE MAP: MIN MARGIN, MIN DIFFERENCE CORR ', &
            &finite_min(margins), finite_min(offdiag(diffcorr))
        call gate%check(tag//'_control_copies_not_distinct', .not. maps_distinct(margins, diffcorr))
        do k = 1, NSTATES
            call ordered(k)%copy(mean)
        enddo
        call discrimination(ordered, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
        write(logfhandle,'(a,2f10.4)') 'COPIES OF THE MEAN MAP: MIN MARGIN, MIN DIFFERENCE CORR ', &
            &finite_min(margins), finite_min(offdiag(diffcorr))
        call gate%check(tag//'_control_mean_not_distinct', .not. maps_distinct(margins, diffcorr))
        call mean%kill
    endif
    ! ---- poses and shifts of the particles in their right state (C1 truths, one frame) ----
    allocate(right(NPTCLS))
    do i = 1, NPTCLS
        right(i) = estates(i) >= 1 .and. estates(i) <= NSTATES
        if( right(i) ) right(i) = map(estates(i)) == tstates(i)
    enddo
    n = count(right)
    allocate(re(3,3,n), rt(3,3,n), sherrs(n))
    l = 0
    do i = 1, NPTCLS
        if( .not. right(i) ) cycle
        l = l + 1
        re(:,:,l) = proj%os_ptcl3D%get_mat(i)
        rt(:,:,l) = truth%get_mat(i)
        sh = proj%os_ptcl3D%get_2Dshift(i) - truth%get_2Dshift(i)
        sherrs(l) = sqrt(sum(sh**2))
    enddo
    ops  = group_operators('c1')
    errs = group_errors(re, rt, ops)
    if( n > 0 )then
        if( l_control )then
            call gate%report(tag//'_median_rotation_error_deg', median(errs))
            call gate%report(tag//'_median_shift_error_px', median(sherrs))
        else
            call gate%metric(tag//'_median_rotation_error_deg', median(errs), MAX_MED_ROT_ERR, &
                &median(errs) <= MAX_MED_ROT_ERR, upper=.true.)
            call gate%metric(tag//'_median_shift_error_px', median(sherrs), MAX_MED_SH_ERR, &
                &median(sherrs) <= MAX_MED_SH_ERR, upper=.true.)
        endif
    else
        call gate%check(tag//'_some_particle_in_its_state', .false.)
    endif
    do k = 1, NSTATES
        call est(k)%kill
        call truths(k)%kill
        call ordered(k)%kill
    enddo
    call truth%kill
    call proj%kill
    deallocate(proj)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('could not return to '//root%to_char())
    CWD_GLOB = root%to_char()

end subroutine run_case

end submodule simple_commanders_test_highlevel_multistate
