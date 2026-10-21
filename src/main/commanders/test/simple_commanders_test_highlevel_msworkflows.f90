!@descr: the multi-state workflow tests of simple_commanders_test_highlevel: solve3D and refine3D_states on beta-galactosidase states against their known labels
! Three workflows on the compositional beta-galactosidase states of the state fixture (the tetramer,
! the tetramer less chain D, the A-B dimer; populations 50, 30 and 20 %), each judged against the
! truth with the metrics of the multi-state assignment test:
!   cavgs   solve3D_cavgs nstates=3 on simulated images of the states imported as its class
!           averages (import_cavgs; rulings of 2026-10-17 and 2026-10-18): 600 images of known state,
!           view and shift, populations 50, 30 and 20 %, without CTF, each the average of an even and
!           an odd half that are independent noise realisations; the state and pose of every imported
!           image are judged;
!   shc     solve3D nstates=3 from the particles (random poses, no states), whose stages 1-2 run shc;
!   states  refine3D_states nstates=3 on a single-state solution: the true poses and the consensus
!           map reconstructed from every particle in one state; flex_pca initializes the states and
!           one prob_neigh march refines them;
!   consensus  the route used on real data (ruling of 2026-10-18): single-state solve3D in C1 from
!           the particles, then refine3D_states nstates=3 on its solution, flex_pca initializing the
!           states from it.
! The maps of a de novo run have their own frame and hand, so every estimated map is docked in both
! hands onto every truth: the docked correlations give the map mapping, the docking onto its mapped
! truth the frame in which maps, poses (up to the truth's declared point group) and shifts are
! judged. Every run scores the shuffled-label and duplicated-maps controls. Whether every true state
! is resolved (an estimated state holds most of its particles: state 0 means exclusion) is checked in
! the states and consensus routes; in the de novo multi-state routes (shc, cavgs) it is decided by the
! early stages, which take the best-scoring state outright from random labels, and a single run may
! miss it (ruling of 2026-10-18): there it is a success, which test=profile holds to the fraction of
! 20 seeds below, and the labels are reported. The default (the CTest entry) runs the four routes;
! suite=<route> runs one, suite=<route>_identical a route on identical states, suite=hard_<route>
! on the conformational states (a chain rotated by 10 degrees, low SNR, a small second state), a
! robustness profile.
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_msworkflows
use simple_test_gate,          only: test_gate
use simple_test_state_fixture, only: state_fixture, FIX_SEARCH, FIX_KNOWN
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSTATES  = 3
integer, parameter :: NPTCLS   = 3000
! The class averages (ruling of 2026-10-18): about 600 in all, each with even and odd halves that are
! independent noise realisations of the same state, view and shift, no CTF, small shift errors. Each
! half at SNR 2.5 makes the imported image (their average) SNR 5: 600 images of SNR 5 carry about
! the signal of the particle route's 3 000 particles of SNR 1
integer, parameter :: NCAVGS     = 600  !< the simulated images imported as class averages
real,    parameter :: CAVG_SNR   = 2.5  !< of each half of an imported image
real,    parameter :: CAVG_SHERR = 1.0  !< pixels: the largest shift of an image (both halves share it)
integer, parameter :: NSHUFFLE = 1000
real,    parameter :: POPS(NSTATES) = [0.5, 0.3, 0.2]
real,    parameter :: SNR      = 1.0
! the hard profile: the tetramer and the tetramer with chain D rotated by 10 degrees
integer, parameter :: HARD_NSTATES = 2
real,    parameter :: HARD_POPS(HARD_NSTATES) = [0.8, 0.2]
real,    parameter :: HARD_SNR = 0.5
real,    parameter :: DOCK_HP  = 40.    !< A, the docking band
real,    parameter :: DOCK_LP  = 10.
real,    parameter :: DIFF_FRAC   = 0.1 !< the difference mask: where two truths differ by this fraction of the larger maximum,
integer, parameter :: DIFF_DILATE = 1   !< dilated by this many pixels
! The state labels of the de novo multi-state routes (shc, cavgs) are reported, a robustness profile: a
! single run may merge or lose a state (ruling of 2026-10-18), so whether every true state is resolved
! is a success, held over 20 seeds by test=profile. The labels of the states and consensus routes are
! gated (below).
! Ceilings of the poses and shifts of the shc route from the profile of 2026-10-17
! (scratch/keep/profiles/phase4/current_t{1,8}_shc, 20 seed offsets at 1 and at 8 threads; s is the
! robust spread):
! - median rotation error: t8 median 1.81 deg, s 0.062, worst 2.05, admissible from 2.12; t1 median
!   1.79, s 0.075, worst 2.13, admissible from 2.20; random poses err by about 90 deg;
! - median shift error: t8 median 0.070 px, s 0.008, worst 0.098, admissible from 0.107; t1 median
!   0.077, s 0.012, worst 0.104, admissible from 0.115; zero shifts err by about 2 px.
! Both stay admissible on the profile of the merged code (2026-10-19, msw_t8_shc): rotation median
! 1.83 deg, s 0.065, worst 1.99, admissible from 2.05; shift median 0.075 px, s 0.008, worst 0.082,
! admissible from 0.098.
real,    parameter :: SHC_MAX_ROT_ERR = 2.5
real,    parameter :: SHC_MAX_SH_ERR  = 0.15
! Ceilings of the poses and shifts of the class-average route (cavgs), from its profile of 2026-10-18
! (scratch/keep/profiles/phase5/ms_cavgs, 20 seed offsets at 8 threads, after the rulings of that day):
! median rotation error median 1.79 deg, s 0.35, worst 3.92, admissible from 4.27; median shift error
! median 0.091 px, s 0.013, worst 0.348, admissible from 0.361; random poses err by about 90 deg and
! zero shifts by about 0.6 px (the simulated shifts are uniform within one pixel). Both stay admissible on
! the profile of the merged code (2026-10-19, msw_t8_cavgs): rotation median 1.82 deg, s 0.25, worst
! 2.86, admissible from 3.11; shift median 0.092 px, s 0.012, worst 0.258, admissible from 0.270.
real,    parameter :: CAVGS_MAX_ROT_ERR = 4.5
real,    parameter :: CAVGS_MAX_SH_ERR  = 0.40
! The states route (refine3D_states) resolves every state in all 20 runs of its profile of
! 2026-10-18 (scratch/keep/profiles/phase5/ms_states, 8 threads, after the flex_pca population-floor
! fix), so its labels, maps and poses are gated again (ruling of 2026-10-17), by the floor rule:
! adjusted Rand index median 1.000, s 0.0003, worst 0.9989, admissible to 0.9986 (floor 0.99; shuffled
! labels 0.00); accuracy and smallest per-state recall worst 0.9993-0.9997 (floors 0.99); largest
! population error worst 0.0003 (ceiling 0.01; a lost state errs by 0.2); discrimination margin median
! 0.450, s 0.002, worst 0.4445, to 0.442 (floor 0.40; copies of one map score below 0); median
! rotation error median 1.44 deg, s 0.017, worst 1.455, from 1.494 (ceiling 2.0); median shift error
! median 0.049 px, worst 0.0498, from 0.0507 (ceiling 0.10). Every one stays admissible on the profile of
! the merged code (2026-10-19, msw_t8_states, all 20 runs resolve every state): adjusted Rand index
! worst 0.9977 (s 0.0006), accuracy worst 0.9993, smallest recall worst 0.9987, population error worst
! 0.0007, margin median 0.450, worst 0.446, rotation median 1.44 deg, worst 1.52 (s 0.02), shift worst
! 0.051 px
real,    parameter :: STATES_MIN_ARI      = 0.99, STATES_MIN_ACC = 0.99, STATES_MIN_RECALL = 0.99
real,    parameter :: STATES_MAX_POP_ERR  = 0.01, STATES_MIN_MARGIN = 0.40
real,    parameter :: STATES_MAX_ROT_ERR  = 2.0,  STATES_MAX_SH_ERR = 0.10
! The fractions of 20 seed offsets in which a de novo route must resolve every true state (ruling of
! 2026-10-18): the count k observed in its profile less three binomial standard deviations
! (sqrt(20 p (1-p)), p = k/20), rounded down, over 20. The profiles of 2026-10-19 at 8 threads
! (scratch/keep/profiles/phase8/msw_t8_<route>): shc 13 of 20, so 13 - 3 x 2.13 = 6.6, floor 6 (0.30);
! cavgs 9 of 20, so 9 - 3 x 2.22 = 2.3, floor 2 (0.10). The shuffled-label control resolved every
! state in none of the runs of either. In every run of both profiles the path assertions, the checks
! and the pose ceilings held, so the nightly keeps them as gates
real,    parameter :: SHC_RESOLVED_FLOOR = 0.30, CAVGS_RESOLVED_FLOOR = 0.10
! The consensus route (ruling of 2026-10-18: gate what holds run to run by the floor rule) resolved
! every true state in all 40 runs of its profiles of 2026-10-19 (scratch/keep/profiles/phase8/
! msw_t8_consensus and msw_t1_consensus, 20 seed offsets each at 8 threads and at one thread), so the
! resolution is checked in every run like the states route's. One state is often only partly
! separated from another, so the labels are coarser than the states route's; as median, s, worst and
! the bound of the floor rule over both profiles (8 threads; one thread):
! - adjusted Rand index 0.78, 0.12, 0.67; 0.79, 0.11, 0.69: to 0.427; shuffled labels 0.00 (+3s 0.35):
!   floor 0.40;
! - accuracy 0.927, 0.044, 0.881; 0.932, 0.041, 0.892: to 0.794; the collapsed predictor 0.50 (+3s
!   0.63) and the null's 99.9th percentile 0.43: floor 0.70;
! - largest population error 0.073, 0.042, 0.111; 0.066, 0.040, 0.108: from 0.198; the collapsed
!   predictor's 0.50 (-3s 0.37): ceiling 0.30;
! - discrimination margin 0.439, 0.017, 0.307; 0.436, 0.018, 0.307: to 0.289; copies of one map score
!   about 0 (+3s 0.05): floor 0.20;
! - the smallest per-state recall 0.76, 0.14, 0.62; 0.78, 0.15, 0.64: to 0.342, under what shuffled
!   labels reach, so it is reported (negative floor below);
! - median rotation error 1.63 deg, 0.06, 1.85; 1.72, 0.10, 1.89: from 2.03, random poses err by about
!   90 deg: ceiling 2.5; median shift error 0.073 px, 0.008, 0.086; 0.072, 0.007, 0.086: from 0.096,
!   zero shifts err by about 2 px: ceiling 0.15.
! The maintainer's Mac has not profiled this route yet.
real,    parameter :: CONSENSUS_LABEL_FLOORS(5) = [0.40, 0.70, -1., 0.30, 0.20]
real,    parameter :: CONSENSUS_MAX_ROT_ERR = 2.5, CONSENSUS_MAX_SH_ERR = 0.15
character(len=*), parameter :: SUITES = 'cavgs, shc, states, consensus, cavgs_identical, shc_identical, '//&
    &'states_identical, consensus_identical, hard_cavgs, hard_shc, hard_states, hard_consensus, cavgs_single, '//&
    &'cavgs_as_particles'
!> suite=cavgs_single (by hand): the single-state class-average start, solve3D_cavgs nstates=1 on the
!! class averages of the tetramer alone, its poses and map reported; prob_athres= passed to the test is
!! passed on (prob_athres=90 is the default before the ruling of 2026-10-18, 10 after it).
!! suite=cavgs_as_particles (by hand, the control of the ruling of 2026-10-18): the same simulated
!! halves imported as 2*NCAVGS particles and solved by particle solve3D nstates=3
real :: cavgs_prob_athres = 0.

contains

!> simple_test_exec test=multistate_workflows [suite=<route>|<route>_identical|hard_<route>|cavgs_single|
!! cavgs_as_particles|list]
!! [nthr=<n>] [prob_athres=<deg>]
module subroutine exec_test_multistate_workflows( self, cline )
    use simple_ui, only: make_ui
    class(commander_test_multistate_workflows), intent(inout) :: self
    class(cmdline),                             intent(inout) :: cline
    class(parameters), allocatable :: params
    type(test_gate) :: gate
    type(string)    :: cwd_saved, root, suite_arg
    character(len=:), allocatable :: suite
    integer :: nthr, status
    logical :: all_ok
    suite_arg = ''
    if( cline%defined('suite') ) suite_arg = cline%get_carg('suite')
    suite = suite_arg%to_char()
    select case(suite)
        case('list')
            write(logfhandle,'(a)') 'Available suites for multistate_workflows: '//SUITES
            return
        case('', 'cavgs', 'shc', 'states', 'consensus', 'cavgs_identical', 'shc_identical', 'states_identical', &
            &'consensus_identical', 'hard_cavgs', 'hard_shc', 'hard_states', 'hard_consensus', 'cavgs_single', &
            &'cavgs_as_particles')
        case default
            THROW_HARD('no sub-suite '//suite//' in multistate_workflows; use suite=list')
    end select
    call make_ui
    call simple_getcwd(cwd_saved)
    if( cline%defined('prob_athres') ) cavgs_prob_athres = cline%get_rarg('prob_athres')
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))
    root = filepath(cwd_saved, 'test_multistate_workflows_'//int2str(get_process_id()))
    if( dir_exists(root) ) call simple_rmdir(root)
    call simple_mkdir(root)
    select case(suite)
        case('')
            call run_case(root, 'cavgs',  'cavgs',  'compositional', nthr, gate)
            call run_case(root, 'shc',    'shc',    'compositional', nthr, gate)
            call run_case(root, 'states', 'states', 'compositional', nthr, gate)
            call run_case(root, 'consensus', 'consensus', 'compositional', nthr, gate)
        case('cavgs', 'shc', 'states', 'consensus')
            call run_case(root, suite, suite, 'compositional', nthr, gate)
        case('cavgs_identical', 'shc_identical', 'states_identical', 'consensus_identical')
            call run_case(root, suite, suite(1:index(suite, '_') - 1), 'identical', nthr, gate)
        case('hard_cavgs', 'hard_shc', 'hard_states', 'hard_consensus')
            call run_case(root, suite, suite(6:), 'conformational', nthr, gate)
        case('cavgs_single')
            call run_case(root, suite, 'cavgs', 'identical', nthr, gate)
        case('cavgs_as_particles')
            call run_case(root, suite, 'cavgs', 'compositional', nthr, gate)
    end select
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('could not return to '//cwd_saved%to_char())
    CWD_GLOB = cwd_saved%to_char()
    all_ok = gate%passed()
    call gate%kill
    if( all_ok )then
        call simple_rmdir(root)
        write(logfhandle,'(a)') 'PASS: the multi-state workflows meet their floors against the planted states'
        call simple_end('**** SIMPLE_TEST_MULTISTATE_WORKFLOWS NORMAL STOP ****')
    else
        THROW_HARD('TEST_MULTISTATE_WORKFLOWS FAILED')
    endif
end subroutine exec_test_multistate_workflows

!> one workflow (route cavgs, shc, states or consensus) in root/<tag> on the betagal states of kind, and its
!! metrics, named after tag
subroutine run_case( root, tag, route, kind, nthr, gate )
    use simple_bench_io, only: bench_report_fname
    class(string),    intent(in)    :: root
    character(len=*), intent(in)    :: tag, route, kind
    integer,          intent(in)    :: nthr
    type(test_gate),  intent(inout) :: gate
    type(state_fixture) :: fix
    type(string) :: dir, projfile, bench
    integer :: status, nstates_run
    logical :: l_control
    write(logfhandle,'(a)') '>>> TEST_MULTISTATE_WORKFLOWS: '//tag
    l_control = kind == 'identical'
    dir = root//'/'//tag
    call simple_mkdir(dir)
    call simple_chdir(dir, status)
    if( status /= 0 ) THROW_HARD('could not enter '//dir%to_char())
    CWD_GLOB = dir%to_char()
    ! the dataset, from the cache when it is there
    fix%family = 'betagal'
    fix%kind   = kind
    fix%nptcls = NPTCLS
    fix%nthr   = nthr
    if( kind == 'conformational' )then
        fix%nstates = HARD_NSTATES
        fix%pops(1:HARD_NSTATES) = HARD_POPS
        fix%snr = HARD_SNR
    else if( tag == 'cavgs_single' )then
        ! one state: the tetramer
        fix%nstates = 1
        fix%pops(1) = 1.
        fix%snr = SNR
    else
        fix%nstates = NSTATES
        fix%pops(1:NSTATES) = POPS
        fix%snr = SNR
    endif
    ! the states the programs are asked for: refine3D_states takes a ceiling of at least three
    nstates_run = fix%nstates
    if( route == 'states' .or. route == 'consensus' ) nstates_run = max(3, fix%nstates)
    select case(route)
        case('cavgs')
            ! the truth maps of the fixture; the images are simulated here
            fix%mode = FIX_SEARCH
            call fix%build
            projfile = string(trim(dir%to_char())//'/cavgs_imported.simple')
            call run_cavgs_route(fix, projfile, nstates_run, nthr, tag, l_control, gate, tag == 'cavgs_as_particles')
        case('shc')
            fix%mode = FIX_SEARCH
            call fix%build
            projfile = fix%project_fname()
            call run_solve3D(fix, projfile, nstates_run, nthr, .false.)
            ! stages 1-2 shc over 500 and 1 000 directions, 3-4 prob over 1 000, 5 prob at the cap
            call check_path(tag, ['shc       ', 'shc       ', 'prob      ', 'prob      '], [500, 1000, 1000, 2500], &
                &nstates_run, gate)
            call evaluate(fix, projfile, tag, l_control, gate, gated=kind == 'compositional', de_novo=.true., &
                &resolved_floor=SHC_RESOLVED_FLOOR)
        case('states')
            fix%mode = FIX_KNOWN
            call fix%build
            projfile = fix%project_fname()
            call run_states_route(fix, projfile, nstates_run, nthr)
            ! the neighbourhood stage: prob_neigh over the multi-state cap of 2 500 directions, with
            ! the states flex_pca delivered (at least two)
            call check_path(tag, ['prob_neigh'], [2500], -2, gate)
            call evaluate(fix, projfile, tag, l_control, gate, gated=kind == 'compositional', &
                &max_rot_err=STATES_MAX_ROT_ERR, max_sh_err=STATES_MAX_SH_ERR, labels_gated=kind == 'compositional')
        case('consensus')
            fix%mode = FIX_SEARCH
            call fix%build
            projfile = fix%project_fname()
            ! the single-state solution in the point group the states share (C1 for these)
            call run_solve3D(fix, projfile, 1, nthr, .false.)
            ! the strategy's record of solve3D's iterations is kept aside, so that the path checked is
            ! refine3D_states' own
            bench = bench_report_fname('REFINE3D_STRATEGY')
            if( file_exists(bench) ) call simple_rename(bench%to_char(), bench%to_char()//'.solve3D')
            call run_refine3D_states(fix, projfile, nstates_run, nthr)
            ! the neighbourhood stage of refine3D_states, with the states flex_pca delivered (at least two)
            call check_path(tag, ['prob_neigh'], [2500], -2, gate)
            call evaluate(fix, projfile, tag, l_control, gate, gated=kind == 'compositional', &
                &max_rot_err=CONSENSUS_MAX_ROT_ERR, max_sh_err=CONSENSUS_MAX_SH_ERR, &
                &labels_gated=kind == 'compositional', label_floors=CONSENSUS_LABEL_FLOORS)
    end select
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('could not return to '//root%to_char())
    CWD_GLOB = root%to_char()
end subroutine run_case

!> NCAVGS simulated images of the states (simulate_particles: random views, the fixture's
!! populations, no CTF, shifts of at most CAVG_SHERR pixels), each the average of an even and an odd
!! half simulated at the same views, states and shifts with independent noise, imported as class averages (import_cavgs) into a new project, then
!! solve3D_cavgs nstates=nstates_run; the state and pose solve3D_cavgs gives every imported image
!! (the class-average segment cls3D) against its truth (cavgs_truth.txt). With l_as_particles, the
!! control: the even and odd halves imported as particles (import_particles) and solved by particle
!! solve3D, each half judged against its image's truth
subroutine run_cavgs_route( fix, projfile, nstates_run, nthr, tag, l_control, gate, l_as_particles )
    use simple_commanders_sim,       only: commander_simulate_particles
    use simple_commanders_project_cls, only: commander_import_cavgs
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_commanders_solve3D,   only: commander_solve3D_cavgs
    use simple_sp_project,           only: sp_project
    use simple_test_utils,           only: set_fixed_seed
    class(state_fixture), intent(in)    :: fix
    class(string),        intent(in)    :: projfile
    integer,              intent(in)    :: nstates_run, nthr
    character(len=*),     intent(in)    :: tag
    logical,              intent(in)    :: l_control
    type(test_gate),      intent(inout) :: gate
    logical,              intent(in)    :: l_as_particles
    character(len=*), parameter :: TRUTH_DOC = 'cavgs_truth.txt', HALVES_DOC = 'halves_truth.txt'
    class(commander_simulate_particles), allocatable :: xsim
    class(commander_import_cavgs),       allocatable :: ximport
    class(commander_solve3D_cavgs),      allocatable :: xcavgs
    class(cmdline),    allocatable :: cl
    class(sp_project), allocatable :: proj
    type(cmdline) :: cline_env
    type(image)   :: even, odd
    character(len=256) :: pops
    character(len=16)  :: buf
    integer :: k, i
    ! the even halves draw the views, the states and the shifts, the odd halves reuse them
    call set_fixed_seed(fix%seed + 1, propagate=.true.)
    allocate(xsim, cl)
    call simulate_half('cavgs_even.mrc', '')
    call simulate_half('cavgs_odd.mrc', TRUTH_DOC)
    deallocate(xsim)
    ! the imported image: the average of its halves
    call even%new([fix%box(), fix%box(), 1], fix%smpd())
    call odd%new([fix%box(), fix%box(), 1], fix%smpd())
    do i = 1, NCAVGS
        call even%read(string('cavgs_even.mrc'), i)
        call odd%read(string('cavgs_odd.mrc'), i)
        call even%add(odd)
        call even%div(2.)
        call even%write(string('cavgs.mrc'), i)
    enddo
    call even%kill
    call odd%kill
    ! a new project, then the images as its class averages
    allocate(proj)
    call proj%update_projinfo(projfile)
    call proj%update_compenv(cline_env)
    call proj%write(projfile)
    call proj%kill
    deallocate(proj)
    if( l_as_particles )then
        call solve_halves_as_particles
        deallocate(cl)
        return
    endif
    allocate(ximport)
    call cl%set('prg',      'import_cavgs')
    call cl%set('projfile', projfile)
    call cl%set('stk',      'cavgs.mrc')
    call cl%set('smpd',     fix%smpd())
    call cl%set('mkdir',    'no')
    call ximport%execute(cl)
    call cl%kill
    deallocate(ximport)
    allocate(xcavgs)
    call cl%set('prg',           'solve3D_cavgs')
    call cl%set('projfile',      projfile)
    call cl%set('mkdir',         'no')
    call cl%set('nstates',       nstates_run)
    call cl%set('pgrp',          'c1')
    call cl%set('mskdiam',       fix%mskdiam())
    ! imported images have no class FRCs to plan the resolution ladder from
    call cl%set('lpstart_ini3D', 20.)
    call cl%set('lpstop_ini3D',  8.)
    call cl%set('rank_cavgs',    'no')
    if( cavgs_prob_athres > 0. ) call cl%set('prob_athres', cavgs_prob_athres)
    call cl%set('nthr',          nthr)
    call xcavgs%execute(cl)
    call cl%kill
    deallocate(xcavgs, cl)
    ! stages 1-2 shc and 3 prob over 500 directions, 4 prob over 1 000, 5 prob and 6-7 prob_neigh
    ! at the cap (multi-state stages 1-2 run shc as particle solve3D does: ruling of 2026-10-18;
    ! prob_neigh from stage 6 as in single state: 2026-10-21); the single-state start is a
    ! measurement, its path not asserted
    if( nstates_run > 1 ) call check_path(tag, ['shc       ', 'prob      ', 'prob      ', 'prob      ', &
        &'prob_neigh'], [500, 500, 1000, 2500, 2500], nstates_run, gate)
    call evaluate(fix, projfile, tag, l_control, gate, gated=.not. l_control .and. nstates_run > 1, de_novo=.true., &
        &images_truth=TRUTH_DOC, nimages=NCAVGS, max_rot_err=CAVGS_MAX_ROT_ERR, max_sh_err=CAVGS_MAX_SH_ERR, &
        &resolved_floor=CAVGS_RESOLVED_FLOOR)

  contains

    !> the control: the even halves then the odd halves as 2*NCAVGS particles without CTF, each with
    !! its image's truth, through particle solve3D
    subroutine solve_halves_as_particles
        class(commander_import_particles), allocatable :: xptcls
        type(oris)  :: tin, tout
        type(image) :: img
        call img%new([fix%box(), fix%box(), 1], fix%smpd())
        do i = 1, NCAVGS
            call img%read(string('cavgs_even.mrc'), i)
            call img%write(string('halves.mrc'), i)
        enddo
        do i = 1, NCAVGS
            call img%read(string('cavgs_odd.mrc'), i)
            call img%write(string('halves.mrc'), NCAVGS + i)
        enddo
        call img%kill
        call tin%new(NCAVGS, is_ptcl=.true.)
        call tin%read(string(TRUTH_DOC), [1, NCAVGS])
        call tout%new(2 * NCAVGS, is_ptcl=.true.)
        do i = 1, NCAVGS
            call tout%transfer_ori(i, tin, i)
            call tout%transfer_ori(NCAVGS + i, tin, i)
        enddo
        call tout%write(string(HALVES_DOC), [1, 2 * NCAVGS])
        call tin%kill
        call tout%kill
        allocate(xptcls)
        call cl%set('prg',      'import_particles')
        call cl%set('projfile', projfile)
        call cl%set('stk',      'halves.mrc')
        call cl%set('smpd',     fix%smpd())
        call cl%set('ctf',      'no')
        call cl%set('mkdir',    'no')
        call xptcls%execute(cl)
        call cl%kill
        deallocate(xptcls)
        call run_solve3D(fix, projfile, nstates_run, nthr, .false.)
        call evaluate(fix, projfile, tag, l_control, gate, images_truth=HALVES_DOC, nimages=2 * NCAVGS, &
            &images_are_particles=.true.)
    end subroutine solve_halves_as_particles

    !> one half: NCAVGS images; with oritab the views and states of the even half
    subroutine simulate_half( outstk, oritab )
        character(len=*), intent(in) :: outstk, oritab
        do k = 1, fix%nstates
            call cl%set('vol'//int2str(k), fix%truth_fname(k))
        enddo
        if( fix%nstates > 1 ) call cl%set('nstates', fix%nstates)
        if( len_trim(oritab) > 0 )then
            call cl%set('oritab', oritab)
        else if( fix%nstates > 1 )then
            pops = ''
            do k = 1, fix%nstates
                write(buf,'(f8.5)') fix%pops(k)
                pops = trim(pops)//trim(adjustl(buf))
                if( k < fix%nstates ) pops = trim(pops)//','
            enddo
            call cl%set('state_pops', trim(pops))
        endif
        call cl%set('prg',     'simulate_particles')
        call cl%set('mkdir',   'no')
        call cl%set('smpd',    fix%smpd())
        call cl%set('mskdiam', fix%mskdiam())
        call cl%set('nptcls',  NCAVGS)
        call cl%set('pgrp',    'c1')
        call cl%set('ctf',     'no')
        call cl%set('snr',     CAVG_SNR)
        ! the even half draws the shifts; the odd half takes them from the oritab
        call cl%set('sherr',   merge(CAVG_SHERR, 0., len_trim(oritab) == 0))
        call cl%set('bfac',    0.)
        call cl%set('outstk',  outstk)
        call cl%set('outfile', trim(merge('cavgs_truth_even.txt', 'cavgs_odd_oris.txt  ', len_trim(oritab) == 0)))
        call cl%set('nthr',    nthr)
        call xsim%execute(cl)
        call cl%kill
        if( len_trim(oritab) == 0 ) call simple_rename('cavgs_truth_even.txt', TRUTH_DOC)
    end subroutine simulate_half

end subroutine run_cavgs_route

!> solve3D nstates=nstates_run in the working directory, from the particles or (cavg_ini_ext) from
!! the class-average solution in the project
subroutine run_solve3D( fix, projfile, nstates_run, nthr, l_cavg_ini_ext )
    use simple_commanders_solve3D, only: commander_solve3D
    class(state_fixture), intent(in) :: fix
    class(string),        intent(in) :: projfile
    integer,              intent(in) :: nstates_run, nthr
    logical,              intent(in) :: l_cavg_ini_ext
    class(commander_solve3D), allocatable :: xsolve
    class(cmdline),           allocatable :: cl
    allocate(xsolve, cl)
    call cl%set('prg',         'solve3D')
    call cl%set('projfile',    projfile)
    call cl%set('mkdir',       'no')
    call cl%set('nstates',     nstates_run)
    call cl%set('pgrp',        'c1')
    call cl%set('mskdiam',     fix%mskdiam())
    call cl%set('print_states', 'yes')
    call cl%set('nthr',        nthr)
    if( l_cavg_ini_ext ) call cl%set('cavg_ini_ext', 'yes')
    call xsolve%execute(cl)
    call cl%kill
    deallocate(xsolve, cl)
end subroutine run_solve3D

!> the single-state solution (the consensus map of every particle at its true pose, registered as
!! the project's state 1 map: reconstruct3D with mkdir=no registers nothing), then refine3D_states.
!! On identical states flex_pca may find fewer than two states and stop the run: run the
!! states_identical control as a suite of its own
subroutine run_states_route( fix, projfile, nstates_run, nthr )
    use simple_commanders_rec,      only: commander_rec3D
    use simple_refine3D_fnames,     only: refine3D_state_vol_fname, refine3D_fsc_fname
    use simple_sp_project,          only: sp_project
    use simple_syslib,              only: simple_abspath
    class(state_fixture), intent(in) :: fix
    class(string),        intent(in) :: projfile
    integer,              intent(in) :: nstates_run, nthr
    class(commander_rec3D), allocatable :: xrec
    class(cmdline),         allocatable :: cl
    class(sp_project), allocatable :: proj
    type(string) :: vol, fsc
    allocate(xrec, cl)
    call cl%set('prg',      'reconstruct3D')
    call cl%set('projfile', projfile)
    call cl%set('mkdir',    'no')
    call cl%set('oritype',  'ptcl3D')
    call cl%set('nstates',  1)
    call cl%set('pgrp',     'c1')
    call cl%set('mskdiam',  fix%mskdiam())
    call cl%set('nthr',     nthr)
    call xrec%execute(cl)
    call cl%kill
    deallocate(xrec)
    vol = refine3D_state_vol_fname(1)
    fsc = refine3D_fsc_fname(1)
    if( .not. file_exists(vol) .or. .not. file_exists(fsc) ) THROW_HARD('the consensus reconstruction is missing')
    allocate(proj)
    call proj%read(projfile)
    call proj%add_vol2os_out(vol, fix%smpd(), 1, 'vol', box=fix%box())
    call proj%add_fsc2os_out(simple_abspath(fsc), 1, fix%box())
    call proj%write(projfile)
    call proj%kill
    deallocate(proj, cl)
    call run_refine3D_states(fix, projfile, nstates_run, nthr)
end subroutine run_states_route

!> refine3D_states nstates=nstates_run in C1 on the single-state solution of the project; its
!! particles carry no state labels, so flex_pca initializes the states
subroutine run_refine3D_states( fix, projfile, nstates_run, nthr )
    use simple_commanders_refine3D, only: commander_refine3D_states
    class(state_fixture), intent(in) :: fix
    class(string),        intent(in) :: projfile
    integer,              intent(in) :: nstates_run, nthr
    class(commander_refine3D_states), allocatable :: xstates
    class(cmdline),                   allocatable :: cl
    allocate(xstates, cl)
    call cl%set('prg',      'refine3D_states')
    call cl%set('projfile', projfile)
    call cl%set('mkdir',    'no')
    call cl%set('nstates',  nstates_run)
    call cl%set('pgrp',     'c1')
    call cl%set('mskdiam',  fix%mskdiam())
    call cl%set('nthr',     nthr)
    call xstates%execute(cl)
    call cl%kill
    deallocate(xstates, cl)
end subroutine run_refine3D_states

!> The path the refinements took, from the refine3D strategy's record of every iteration (refine
!! mode, directions, states): the blocks of consecutive iterations with one mode and one number of
!! directions must contain the expected ones in order; every iteration runs nstates_exp states
!! (-n: at least n) and no more than the multi-state cap of 2 500 directions
subroutine check_path( tag, modes, nspaces, nstates_exp, gate )
    use simple_bench_io, only: bench_report_fname
    character(len=*), intent(in)    :: tag
    character(len=*), intent(in)    :: modes(:)
    integer,          intent(in)    :: nspaces(:), nstates_exp
    type(test_gate),  intent(inout) :: gate
    character(len=32), allocatable :: bmodes(:)
    integer,           allocatable :: bnspaces(:)
    character(len=512) :: line
    character(len=32)  :: mode
    type(string) :: fname
    integer :: funit, ios, nspace, nst, nrec, j, k
    logical :: l_states, l_cap, l_order
    fname = bench_report_fname('REFINE3D_STRATEGY')
    if( .not. file_exists(fname) )then
        call gate%check(tag//'_path_recorded', .false.)
        return
    endif
    allocate(bmodes(0), bnspaces(0))
    mode     = ''
    nspace   = 0
    nrec     = 0
    l_states = .true.
    l_cap    = .true.
    open(newunit=funit, file=fname%to_char(), status='old', action='read')
    do
        read(funit, '(a)', iostat=ios) line
        if( ios /= 0 ) exit
        if( index(line, 'refine3D refine mode') == 1 ) mode = adjustl(line(index(line, ':') + 1:))
        if( index(line, 'refine3D nspace ') == 1 ) read(line(index(line, ':') + 1:), *) nspace
        if( index(line, 'refine3D nstates') == 1 )then
            ! the last field of a record's context: the record is complete
            read(line(index(line, ':') + 1:), *) nst
            nrec = nrec + 1
            if( nstates_exp > 0 )then
                l_states = l_states .and. nst == nstates_exp
            else
                l_states = l_states .and. nst >= -nstates_exp
            endif
            l_cap = l_cap .and. nspace <= 2500
            j = size(bmodes)
            if( j == 0 )then
                bmodes   = [mode]
                bnspaces = [nspace]
            else if( bmodes(j) /= mode .or. bnspaces(j) /= nspace )then
                bmodes   = [bmodes, mode]
                bnspaces = [bnspaces, nspace]
            endif
        endif
    enddo
    close(funit)
    write(logfhandle,'(a)') 'PATH '//tag//' (refine mode, directions):'
    do j = 1, size(bmodes)
        write(logfhandle,'(4x,a,i6)') trim(bmodes(j)), bnspaces(j)
    enddo
    ! the expected blocks as an ordered subsequence
    k = 1
    do j = 1, size(bmodes)
        if( k > size(modes) ) exit
        if( trim(bmodes(j)) == trim(modes(k)) .and. bnspaces(j) == nspaces(k) ) k = k + 1
    enddo
    l_order = k > size(modes)
    call gate%check(tag//'_path_recorded', nrec > 0)
    call gate%check(tag//'_path_modes_and_directions', l_order)
    call gate%check(tag//'_path_states', l_states)
    call gate%check(tag//'_path_direction_cap', l_cap)
end subroutine check_path

!> labels, maps (docked in both hands onto every truth), discrimination, poses and shifts of the
!! finished run in the working directory against the truth; the controls
subroutine evaluate( fix, projfile, tag, l_control, gate, gated, images_truth, nimages, images_are_particles, &
    &max_rot_err, max_sh_err, labels_gated, de_novo, resolved_floor, label_floors )
    use simple_sp_project,         only: sp_project
    use simple_rnd,                only: shuffle
    use simple_test_utils,         only: set_fixed_seed
    use simple_test_truth_metrics, only: match_states, adjusted_rand_index, mapped_accuracy, accuracy_null, &
        &collapsed_accuracy, population_errors, match_maps, discrimination, group_errors, group_operators, &
        &compare_to_truth, dock_both_hands, carry_pose, finite_min, offdiag, maps_distinct, state_recalls
    use simple_oris,               only: oris
    class(state_fixture), intent(in)    :: fix
    class(string),        intent(in)    :: projfile
    character(len=*),     intent(in)    :: tag
    logical,              intent(in)    :: l_control
    type(test_gate),      intent(inout) :: gate
    logical, optional,    intent(in)    :: gated !< the poses and shifts held to their ceilings
    !> imported images judged instead of particles: their truth document and their number (the
    !! states and poses are those of the class-average segment cls3D)
    character(len=*), optional, intent(in) :: images_truth
    integer,          optional, intent(in) :: nimages
    logical,          optional, intent(in) :: images_are_particles !< the judged images are the particle segment
    real,             optional, intent(in) :: max_rot_err, max_sh_err !< the ceilings when gated (the shc route's by default)
    logical,          optional, intent(in) :: labels_gated !< the labels and the discrimination held to the states route's floors
    !> a de novo route: whether every true state is resolved is a success, held by test=profile to the
    !! fraction resolved_floor of its runs, not a check of the run
    logical,          optional, intent(in) :: de_novo
    real,             optional, intent(in) :: resolved_floor
    !> the floors of gated labels: adjusted Rand index, accuracy, smallest recall (negative: reported),
    !! the ceiling of the population error and the floor of the discrimination margin; the states
    !! route's by default
    real,             optional, intent(in) :: label_floors(5)
    class(sp_project), allocatable, target :: proj
    class(oris), pointer :: os => null()
    type(string) :: truth_doc
    type(oris)   :: truth
    type(image), allocatable :: docked(:), truths(:)
    type(image)  :: mean
    type(string) :: fname
    integer, allocatable :: tstates(:), estates(:), map(:), mapmaps(:), shuffled(:), mapctrl(:)
    real,    allocatable :: cc(:,:), eul(:,:,:), sh3(:,:,:), margins(:), diffcorr(:,:), ops(:,:,:)
    real,    allocatable :: roterr(:), sherr(:), rt(:,:,:), rc(:,:,:)
    logical, allocatable :: mirrored(:,:), has_map(:), l_mapped(:)
    real    :: ari, acc, med_null, p999, poperr, ccd, ccm, corr, fsc05, fsc0143, shc(2), r(3,3), d(2)
    integer :: ntrue, nest, k, e, i, n, nmask, ldim(3), nobs
    logical :: l_gated, l_labels, l_de_novo
    real    :: rot_ceil, sh_ceil, recall, margin, lf(5)
    l_gated = .false.
    if( present(gated) ) l_gated = gated
    l_de_novo = .false.
    if( present(de_novo) ) l_de_novo = de_novo .and. .not. l_control
    lf = [STATES_MIN_ARI, STATES_MIN_ACC, STATES_MIN_RECALL, STATES_MAX_POP_ERR, STATES_MIN_MARGIN]
    if( present(label_floors) ) lf = label_floors
    if( l_de_novo .and. .not. present(resolved_floor) ) THROW_HARD('a de novo route needs the floor of its successes')
    l_labels = .false.
    if( present(labels_gated) ) l_labels = labels_gated
    rot_ceil = SHC_MAX_ROT_ERR
    sh_ceil  = SHC_MAX_SH_ERR
    if( present(max_rot_err) ) rot_ceil = max_rot_err
    if( present(max_sh_err)  ) sh_ceil  = max_sh_err
    nobs      = fix%nptcls
    truth_doc = fix%truth_doc()
    if( present(images_truth) .and. present(nimages) )then
        nobs      = nimages
        truth_doc = images_truth
    endif
    ntrue = fix%nstates
    ldim  = [fix%box(), fix%box(), fix%box()]
    ! ---- labels ----
    allocate(proj)
    call proj%read(projfile)
    os => proj%os_ptcl3D
    if( present(images_truth) ) os => proj%os_cls3D
    if( present(images_are_particles) )then
        if( images_are_particles ) os => proj%os_ptcl3D
    endif
    call truth%new(nobs, is_ptcl=.true.)
    call truth%read(truth_doc, [1, nobs])
    tstates = nint(truth%get_all('state'))
    estates = nint(os%get_all('state'))
    if( size(estates) /= nobs ) THROW_HARD('the project does not hold one record per judged image')
    nest    = max(1, maxval(estates))
    write(logfhandle,'(a,i3,a,i3)') 'TRUE STATES ', ntrue, '; ESTIMATED STATES ', nest
    map = match_states(tstates, estates, ntrue, nest)
    acc = mapped_accuracy(tstates, estates, map)
    ari = adjusted_rand_index(tstates, estates)
    call accuracy_null(tstates, estates, ntrue, nest, NSHUFFLE, 20261018, med_null, p999)
    poperr = maxval(population_errors(tstates, estates, map, ntrue))
    call gate%report(tag//'_estimated_states', real(nest))
    ! particles a workflow leaves in state 0 (inactive) count as wrong in the accuracy; the adjusted
    ! Rand index sees them as one more cluster
    call gate%report(tag//'_unassigned_fraction', real(count(estates < 1)) / real(nobs))
    if( l_labels )then
        call gate%metric(tag//'_ari', ari, lf(1), ari >= lf(1))
    else
        call gate%report(tag//'_ari', ari)
    endif
    ! a merge of two true states leaves one unrecalled and an estimated state empty; state 0 means
    ! exclusion, so a true state is resolved when an estimated state holds most of its particles
    recall = minval(state_recalls(tstates, estates, map, ntrue))
    if( l_labels .and. lf(3) >= 0. )then
        call gate%metric(tag//'_min_state_recall', recall, lf(3), recall >= lf(3))
    else
        call gate%report(tag//'_min_state_recall', recall)
    endif
    if( l_de_novo )then
        call gate%success(tag//'_every_true_state_resolved', all(state_recalls(tstates, estates, map, ntrue) > 0.5), &
            &resolved_floor)
    else if( .not. l_control )then
        call gate%check(tag//'_every_true_state_resolved', all(state_recalls(tstates, estates, map, ntrue) > 0.5))
        call gate%check(tag//'_no_estimated_state_empty', all([(count(estates == e) > 0, e = 1, nest)]))
    endif
    if( l_labels )then
        call gate%metric(tag//'_accuracy', acc, lf(2), acc >= lf(2))
    else
        call gate%report(tag//'_accuracy', acc)
    endif
    call gate%report(tag//'_accuracy_null_median', med_null)
    call gate%report(tag//'_accuracy_null_p999', p999)
    call gate%report(tag//'_collapsed_accuracy', collapsed_accuracy(tstates, ntrue))
    if( l_labels )then
        call gate%metric(tag//'_largest_population_error', poperr, lf(4), poperr <= lf(4), &
            &upper=.true.)
    else
        call gate%report(tag//'_largest_population_error', poperr)
    endif
    ! the shuffled-label control: truth labels shuffled across particles, the counts kept
    call set_fixed_seed(20261019)
    shuffled = tstates
    call shuffle(shuffled)
    mapctrl = match_states(shuffled, estates, ntrue, nest)
    call gate%report(tag//'_control_shuffled_ari', adjusted_rand_index(shuffled, estates))
    ! the success count of the shuffled labels, which the floor of the true labels must clear
    if( l_de_novo ) call gate%success(tag//'_control_shuffled_every_true_state_resolved', &
        &all(state_recalls(shuffled, estates, mapctrl, ntrue) > 0.5), 0.)
    call gate%check(tag//'_control_shuffled_accuracy_within_null', mapped_accuracy(shuffled, estates, mapctrl) <= p999)
    write(logfhandle,'(a)') 'CONFUSION (rows: estimated states; columns: true states)'
    do e = 1, nest
        write(logfhandle,'(8i8)') (count(estates == e .and. tstates == k), k = 1, ntrue)
    enddo
    ! ---- maps: every estimated map docked in both hands onto every truth ----
    allocate(cc(nest,ntrue), eul(3,nest,ntrue), sh3(3,nest,ntrue), mirrored(nest,ntrue), has_map(nest))
    cc = 0.
    eul = 0.
    sh3 = 0.
    mirrored = .false.
    do e = 1, nest
        fname = 'rec_final_state'//int2str_pad(e, 2)//'.mrc'
        has_map(e) = file_exists(fname)
        if( .not. has_map(e) ) cycle
        do k = 1, ntrue
            call dock_both_hands(fix%truth_fname(k), fname, fix%mskdiam(), DOCK_HP, DOCK_LP, &
                &'dock_e'//int2str(e)//'_t'//int2str(k), string(docked_fname(e, k)), ccd, ccm, &
                &eul(:,e,k), sh3(:,e,k), mirrored(e,k))
            cc(e,k) = max(ccd, ccm)
        enddo
    enddo
    write(logfhandle,'(a)') 'DOCKED MAP CORRELATIONS (rows: estimated maps; columns: truths)'
    do e = 1, nest
        write(logfhandle,'(8f8.3)') cc(e,:)
    enddo
    mapmaps = match_maps(cc)
    ! every populated state has a map (an emptied state keeps its number and loses its map)
    call gate%check(tag//'_every_populated_state_has_a_map', all(has_map .or. [(count(estates == e) == 0, e = 1, nest)]))
    call gate%report(tag//'_empty_states', real(count([(count(estates == e) == 0, e = 1, nest)])))
    call gate%report(tag//'_map_mapping_agrees_with_label_mapping', merge(1., 0., all(mapmaps == map)))
    ! the fit of every mapped map to its truth, in the truth's frame
    do e = 1, nest
        k = map(e)
        if( k < 1 .or. .not. has_map(e) ) cycle
        call compare_to_truth(fix%truth_fname(k), string(docked_fname(e, k)), fix%mskdiam(), corr, fsc05, fsc0143)
        call gate%report(tag//'_map_corr_state'//int2str(k), corr)
        call gate%report(tag//'_map_fsc0143_state'//int2str(k), fsc0143)
    enddo
    ! the discrimination of the docked maps in the order of the truths they map onto, and the
    ! duplicated-maps control: every map replaced by the first one, then by their mean
    allocate(docked(ntrue), truths(ntrue), margins(ntrue), diffcorr(ntrue,ntrue))
    allocate(l_mapped(ntrue), source=.false.)
    do k = 1, ntrue
        call truths(k)%new(ldim, fix%smpd())
        call truths(k)%read(fix%truth_fname(k))
        call docked(k)%new(ldim, fix%smpd())
        e = findloc(map, k, dim=1)
        if( e > 0 ) l_mapped(k) = has_map(e)
        if( l_mapped(k) ) call docked(k)%read(string(docked_fname(e, k)))
    enddo
    if( ntrue > 1 )then
        call discrimination(docked, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
        if( l_control )then
            ! identical truths leave no difference mask: nothing to discriminate
            call gate%check(tag//'_states_not_distinct', .not. maps_distinct(margins, diffcorr))
        else if( all(l_mapped) )then
            margin = finite_min(margins)
            if( l_labels )then
                call gate%metric(tag//'_discrimination_min_margin', margin, lf(5), margin >= lf(5))
            else
                call gate%report(tag//'_discrimination_min_margin', margin)
            endif
            call gate%report(tag//'_discrimination_min_difference_corr', finite_min(offdiag(diffcorr)))
        else
            ! a truth without a map discriminates nothing: the lowest values, not a non-finite
            ! correlation that would fail the row of a reported (ungated) metric
            call gate%report(tag//'_discrimination_min_margin', -2.)
            call gate%report(tag//'_discrimination_min_difference_corr', -1.)
        endif
        if( .not. l_control )then
            call mean%copy(docked(1))
            do k = 2, ntrue
                call mean%add(docked(k))
            enddo
            call mean%div(real(ntrue))
            do k = 2, ntrue
                call docked(k)%copy(docked(1))
            enddo
            call discrimination(docked, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
            call gate%check(tag//'_control_copies_not_distinct', .not. maps_distinct(margins, diffcorr))
            do k = 1, ntrue
                call docked(k)%copy(mean)
            enddo
            call discrimination(docked, truths, DIFF_FRAC, DIFF_DILATE, margins, diffcorr, nmask)
            call gate%check(tag//'_control_mean_not_distinct', .not. maps_distinct(margins, diffcorr))
            call mean%kill
        endif
    endif
    ! ---- poses and shifts of the particles in their right state, carried into the truth's frame ----
    allocate(roterr(0), sherr(0))
    do k = 1, ntrue
        e = findloc(map, k, dim=1)
        if( e < 1 ) cycle
        if( .not. has_map(e) ) cycle
        n = count(estates == e .and. tstates == k)
        if( n == 0 ) cycle
        allocate(rc(3,3,n), rt(3,3,n))
        n = 0
        do i = 1, nobs
            if( estates(i) /= e .or. tstates(i) /= k ) cycle
            n = n + 1
            r = os%get_mat(i)
            call carry_pose(r, os%get_2Dshift(i), mirrored(e,k), fix%box(), euler2m(eul(:,e,k)), &
                &sh3(:,e,k), rc(:,:,n), shc)
            rt(:,:,n) = truth%get_mat(i)
            d = shc - truth%get_2Dshift(i)
            sherr = [sherr, sqrt(sum(d**2))]
        enddo
        ops    = group_operators(trim(fix%pgrps(k)))
        roterr = [roterr, group_errors(rc, rt, ops)]
        deallocate(rc, rt)
    enddo
    if( size(roterr) > 0 )then
        if( l_gated )then
            call gate%metric(tag//'_median_rotation_error_deg', median(roterr), rot_ceil, &
                &median(roterr) <= rot_ceil, upper=.true.)
            call gate%metric(tag//'_median_shift_error_px', median(sherr), sh_ceil, &
                &median(sherr) <= sh_ceil, upper=.true.)
        else
            call gate%report(tag//'_median_rotation_error_deg', median(roterr))
            call gate%report(tag//'_median_shift_error_px', median(sherr))
        endif
    else
        call gate%check(tag//'_some_particle_in_its_state', .false.)
    endif
    do k = 1, ntrue
        call docked(k)%kill
        call truths(k)%kill
    enddo
    call truth%kill
    call fname%kill
    nullify(os)
    call proj%kill
    deallocate(proj)

  contains

    function docked_fname( ie, ik ) result( f )
        integer, intent(in) :: ie, ik
        character(len=:), allocatable :: f
        f = 'docked_e'//int2str(ie)//'_t'//int2str(ik)//'.mrc'
    end function docked_fname

end subroutine evaluate

end submodule simple_commanders_test_highlevel_msworkflows
