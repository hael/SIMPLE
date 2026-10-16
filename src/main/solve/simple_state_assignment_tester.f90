!@descr: unit tests of the multistate class-average assignment rules shared by stream p07 and solve3D_pyramid (simple_state_assignment)
! The rules on fixed numbers (the check of the reference volumes, the renumbering, the sub-class budget,
! the mask diameter, the rejection by fit, the ambiguity flag, the first set's cut at whole stacks), the
! class averages grouped by state on fixture images, a state's own solve3D_cavgs (its classes selected,
! its poses mapped onto its particles) on fixture projects, and the jobs' command lines. The steps on projects
! are tested through stream p07's stage (simple_stream_stage_solve3D_tester), the docking in its
! library suite.
module simple_state_assignment_tester
use simple_string,           only: string
use simple_fileio,           only: file_exists
use simple_imghead,          only: find_ldim_nptcls
use simple_cmdline,          only: cmdline
use simple_image,            only: image
use simple_sp_project,       only: sp_project
use simple_state_assignment, only: NCLS2D_FIRST, LPSTART_CAVGS3D, LPSTOP_CAVGS3D, LP_ASSIGN, NSPACE_ASSIGN,&
                                  &refs_usable, reject_by_fit, flag_ambiguous, compact_states, expansion_budget,&
                                  &count_active_classes, fit_mskdiam, cut_whole_stacks, write_state_cavgs,&
                                  &set_solve2D_route_cline, set_clsexp_cline, set_cavgs3D_cline, set_assign_cline,&
                                  &set_refine_cont_cline, set_state_cavgs3D_cline, select_state_classes, map_state_poses
use simple_test_utils
implicit none
private
public :: run_all_state_assignment_tests

integer, parameter :: IMG_BOX  = 16
real,    parameter :: IMG_SMPD = 2.0

contains

    subroutine run_all_state_assignment_tests()
        write(*,'(A)') '**** running all state assignment tests ****'
        call test_check_rules()
        call test_fit_rules()
        call test_cut_whole_stacks()
        call test_state_cavgs()
        call test_state_runs()
        call test_job_clines()
    end subroutine run_all_state_assignment_tests

    !> the check of the reference volumes, and the renumbering, sub-class budget and mask diameter
    !! the route uses
    subroutine test_check_rules()
        integer, allocatable :: states(:), newind(:)
        integer :: nstates_held
        write(*,'(A)') 'test_check_rules'
        ! the check: even/odd angular distance, state agreement, distinct states, two states of 5%
        call assert_true(refs_usable(12., 85., [40, 30, 30], 0.5),   'halves that agree, three distinct states')
        call assert_true(refs_usable(50., 70., [50, 50], 0.97),      'at the thresholds')
        call assert_false(refs_usable(55., 85., [40, 30, 30], 0.5),  'halves too far apart')
        call assert_false(refs_usable(12., 60., [40, 30, 30], 0.5),  'states that disagree between the halves')
        call assert_false(refs_usable(-1., 85., [40, 30, 30], 0.5),  'an angular distance not measured')
        call assert_false(refs_usable(12., -1., [40, 30, 30], 0.5),  'an agreement not measured')
        call assert_false(refs_usable(12., 85., [40, 30, 30], 0.99), 'two near-duplicate states')
        call assert_false(refs_usable(12., 85., [96, 4], 0.5),       'a state under 5%: one state counts')
        call assert_true(refs_usable(12., 85., [95, 5], 0.5),        'two states of 5%')
        call assert_false(refs_usable(12., 85., [100], -1.),         'one state')
        ! the states renumbered from 1
        states = [0, 3, 1, 3, 0]
        call compact_states(states, newind, nstates_held)
        call assert_int(2, nstates_held, 'two states are held')
        call assert_true(all(states == [0, 2, 1, 2, 0]), 'renumbered from 1 in order')
        call assert_true(all(newind == [1, 0, 2]),        'state 2, held by none, gets no number')
        states = [0, 0]
        call compact_states(states, newind, nstates_held)
        call assert_int(0, nstates_held, 'none when no state is held')
        ! cls_expansion's total of sub-classes, 500 and never fewer than the classes
        call assert_int(500, expansion_budget(150), '150 classes: 500 sub-classes in all')
        call assert_int(500, expansion_budget(0),   'no class: 500 (cls_expansion then finds nothing to split)')
        call assert_int(600, expansion_budget(600), '600 classes: one each at least')
        ! the classes cls_expansion splits: distinct labels of the selected particles
        call assert_int(2, count_active_classes([1, 1, 0, 1, 1], [3, 7, 9, 3, 0]), 'classes 3 and 7; a deselected or unlabelled row is not counted')
        call assert_int(0, count_active_classes([0, 0], [1, 2]),                  'none without a selected particle')
        ! the mask diameter that fits the class averages
        call assert_real(100., fit_mskdiam(100.,  64, 2.0), 1.e-4, 'a mask diameter that fits is kept')
        call assert_real(116., fit_mskdiam(400.,  64, 2.0), 1.e-4, 'one too large for the box gets the box default')
        call assert_real(116., fit_mskdiam(0.,    64, 2.0), 1.e-4, 'none gets the box default')
        call assert_real(116., fit_mskdiam(-7.8,  64, 2.0), 1.e-4, 'a negative one too')
    end subroutine test_check_rules

    !> the rejection by fit (below the median less 3 MADs, the worst first, at most 30% of the
    !! candidates, never a class that is no candidate) and the ambiguity flag
    subroutine test_fit_rules()
        logical, allocatable :: l_rej(:), l_cand(:)
        real    :: med, mad, thres
        integer :: nmax
        real,    parameter   :: CORRS_TWO(10)  = [0.50, 0.52, 0.48, 0.51, 0.49, 0.50, 0.53, 0.47, 0.10, 0.20]
        real,    parameter   :: CORRS_FOUR(10) = [0.50, 0.51, 0.49, 0.50, 0.52, 0.48, 0.10, 0.12, 0.11, 0.13]
        write(*,'(A)') 'test_fit_rules'
        allocate(l_cand(10), source=.true.)
        ! two outliers, under the cap; the statistics for the log
        call reject_by_fit(CORRS_TWO, l_cand, l_rej, med, mad, thres, nmax)
        call assert_int(2, count(l_rej),           'two outliers are rejected')
        call assert_true(l_rej(9) .and. l_rej(10), 'the two that fit worst')
        call assert_real(0.495, med,   1.e-5, 'the candidates'' median')
        call assert_real(0.02,  mad,   1.e-5, 'their median absolute deviation')
        call assert_real(0.435, thres, 1.e-5, 'the threshold, 3 MADs below the median')
        call assert_int(3, nmax,                   'the cap, 30% of 10')
        ! four outliers: the cap of 3 of 10, the worst first
        call reject_by_fit(CORRS_FOUR, l_cand, l_rej)
        call assert_int(3, count(l_rej),                         'at most 30% of the candidates')
        call assert_true(l_rej(7) .and. l_rej(9) .and. l_rej(8), 'the worst three')
        call assert_false(l_rej(10),                             'the fourth is kept')
        ! a class that is no candidate is never rejected, and the statistics are the candidates'
        l_cand(7) = .false.
        call reject_by_fit(CORRS_FOUR, l_cand, l_rej)
        call assert_false(l_rej(7),     'a class that is no candidate is not rejected')
        call assert_int(2, count(l_rej), 'two of the nine candidates, the cap')
        call assert_true(l_rej(9) .and. l_rej(8), 'the worst two candidates')
        ! none when every class fits
        call reject_by_fit([0.50, 0.51, 0.49, 0.50], [.true., .true., .true., .true.], l_rej)
        call assert_int(0, count(l_rej), 'none when every class fits')
        ! no candidate: nothing, and zero statistics
        call reject_by_fit([0.5, 0.1], [.false., .false.], l_rej, med, mad, thres, nmax)
        call assert_int(0, count(l_rej), 'none without candidates')
        call assert_int(0, nmax,         'and no cap')
        ! the margin flags the kept classes below 0.02
        call assert_true(all(flag_ambiguous([0.01, 0.05, 0.0, 0.019], [.true., .true., .false., .true.]) .eqv.&
            &[.true., .false., .false., .true.]), 'small margins of kept classes are ambiguous')
    end subroutine test_fit_rules

    !> the first set: whole stacks in row order while they fit; the stack that crosses the count
    !! ends the cut; a first stack alone above the count gives its first particles
    subroutine test_cut_whole_stacks()
        integer, parameter :: STKIND(9) = [1, 1, 1, 2, 2, 3, 3, 3, 3] ! stacks of 3, 2 and 4 rows
        logical, parameter :: ALL_SEL(9) = .true.
        logical, parameter :: SOME_SEL(9) = [.true., .false., .true., .true., .true., .false., .true., .true., .true.]
        logical :: l_in(9)
        write(*,'(A)') 'test_cut_whole_stacks'
        l_in = cut_whole_stacks(ALL_SEL, STKIND, 5)
        call assert_true(all(l_in .eqv. [.true., .true., .true., .true., .true., .false., .false., .false., .false.]),&
            &'two whole stacks fill the count')
        l_in = cut_whole_stacks(ALL_SEL, STKIND, 6)
        call assert_int(5, count(l_in), 'the stack that crosses the count is left out')
        call assert_false(any(l_in(6:9)), 'and every later one')
        l_in = cut_whole_stacks(ALL_SEL, STKIND, 2)
        call assert_true(all(l_in .eqv. [.true., .true., .false., .false., .false., .false., .false., .false., .false.]),&
            &'a first stack above the count: its first particles')
        ! selected per stack: 2, 2 and 3
        l_in = cut_whole_stacks(SOME_SEL, STKIND, 4)
        call assert_int(4, count(l_in), 'selected particles are counted, not rows')
        call assert_true(l_in(1) .and. l_in(3) .and. l_in(4) .and. l_in(5), 'the selected ones of stacks 1 and 2')
        call assert_false(l_in(2), 'a deselected row stays out')
        l_in = cut_whole_stacks(SOME_SEL, STKIND, 3)
        call assert_int(2, count(l_in), 'stack 2 crosses 3: stack 1''s two only')
        l_in = cut_whole_stacks(SOME_SEL, STKIND, 100)
        call assert_true(all(l_in .eqv. SOME_SEL), 'a count above the selection: every selected particle')
        l_in = cut_whole_stacks(ALL_SEL, STKIND, 0)
        call assert_false(any(l_in), 'none for a count of 0')
        l_in = cut_whole_stacks([.false., .false., .false., .true., .true., .false., .false., .false., .false.], STKIND, 1)
        call assert_true(l_in(4) .and. count(l_in) == 1, 'a first stack without selected particles is passed over')
    end subroutine test_cut_whole_stacks

    !> class averages grouped by state: a stack and a JPEG per state, one of the rejected and one of
    !! the ambiguous ones (the ambiguous also in their state's); a state without class averages
    !! writes nothing
    subroutine test_state_cavgs()
        type(image),  allocatable :: cavgs(:)
        real,         allocatable :: rmat(:,:,:)
        type(string)       :: cwd_saved, root
        integer, parameter :: STATES(5) = [1, 2, 1, 0, 2]
        logical, parameter :: L_REJ(5)  = [.false., .false., .false., .true., .false.]
        logical, parameter :: L_AMB(5)  = [.false., .true., .false., .false., .false.]
        integer :: nfail0, i, x, y, ldim(3), n
        write(*,'(A)') 'test_state_cavgs'
        nfail0 = tests_failed
        call enter_fixture('state_assignment_cavgs', cwd_saved, root)
        allocate(cavgs(size(STATES)), rmat(IMG_BOX,IMG_BOX,1))
        do i = 1,size(STATES)
            do y = 1,IMG_BOX
                do x = 1,IMG_BOX
                    rmat(x,y,1) = real(i) * real(x + 2 * y)
                enddo
            enddo
            call cavgs(i)%new([IMG_BOX, IMG_BOX, 1], IMG_SMPD)
            call cavgs(i)%set_rmat(rmat, .false.)
        enddo
        call write_state_cavgs(cavgs, STATES, 3, string('.'), 'TEST', l_rej=L_REJ, l_amb=L_AMB)
        call find_ldim_nptcls(string('cavgs_state01.mrc'), ldim, n)
        call assert_int(2, n, 'state 1''s two class averages')
        call find_ldim_nptcls(string('cavgs_state02.mrc'), ldim, n)
        call assert_int(2, n, 'state 2''s two, the ambiguous one included')
        call assert_false(file_exists(string('cavgs_state03.mrc')), 'a state without class averages writes nothing')
        call find_ldim_nptcls(string('cavgs_rejected.mrc'), ldim, n)
        call assert_int(1, n, 'the rejected one')
        call find_ldim_nptcls(string('cavgs_ambiguous.mrc'), ldim, n)
        call assert_int(1, n, 'the ambiguous one')
        call assert_true(file_exists(string('cavgs_state01.jpg')) .and. file_exists(string('cavgs_state02.jpg')),&
            &'a JPEG per state with class averages')
        call assert_true(file_exists(string('cavgs_rejected.jpg')) .and. file_exists(string('cavgs_ambiguous.jpg')),&
            &'and of the rejected and the ambiguous')
        do i = 1,size(cavgs)
            call cavgs(i)%kill
        enddo
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_state_cavgs

    !> a state's own solve3D_cavgs: the classes it holds selected, the others not; its class poses
    !! mapped onto its particles only, a row of another state or none keeping its pose; a result of
    !! other rows not taken
    subroutine test_state_runs()
        class(sp_project), allocatable :: spproj, res
        integer, parameter :: ROW_STATES(4) = [1, 2, 1, 0], ROW_CLASSES(4) = [1, 2, 1, 3]
        integer :: i, nsel, nposes
        write(*,'(A)') 'test_state_runs'
        allocate(spproj, res)
        ! the assignment's project: three classes in states 1, 2 and none (rejected)
        call res%os_cls2D%new(3, is_ptcl=.false.)
        call res%os_cls3D%new(3, is_ptcl=.false.)
        do i = 1,3
            call res%os_cls2D%set_state(i, 1)
        enddo
        call res%os_cls3D%set_state(1, 1)
        call res%os_cls3D%set_state(2, 2)
        call res%os_cls3D%set_state(3, 0)
        nsel = select_state_classes(res, 1)
        call assert_int(1, nsel, 'state 1 holds one class')
        call assert_int(1, res%os_cls2D%get_state(1), 'its class is selected')
        call assert_int(0, res%os_cls2D%get_state(2), 'another state''s class is not')
        call assert_int(0, res%os_cls2D%get_state(3), 'nor a rejected one')
        ! the state's solve3D_cavgs result: its class with a pose, the others deselected
        call res%os_ptcl2D%new(4, is_ptcl=.true.)
        call res%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1,4
            call res%os_ptcl2D%set_class(i, ROW_CLASSES(i))
            call res%os_ptcl2D%set(i, 'corr', 0.5)
            call res%os_ptcl2D%set_state(i, 1)
        enddo
        call res%os_cls3D%set_euler(1, [10., 45., 0.])
        call res%os_cls3D%set_euler(2, [20., 90., 0.])
        call res%os_cls3D%set_state(2, 0)
        ! the project's rows: rows 1 and 3 in state 1, row 2 in state 2, row 4 deselected
        call spproj%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1,4
            call spproj%os_ptcl3D%set_state(i, ROW_STATES(i))
            call spproj%os_ptcl3D%set_euler(i, [1., 2., 3.])
        enddo
        nposes = map_state_poses(spproj, res, 1)
        call assert_int(2, nposes, 'the state''s two particles take its poses')
        call assert_real(45., spproj%os_ptcl3D%get(1, 'e2'), 1.e-3, 'a particle takes its class''s pose')
        call assert_real(10., spproj%os_ptcl3D%get(3, 'e1'), 1.e-3, 'every particle of the state')
        call assert_real(2.,  spproj%os_ptcl3D%get(2, 'e2'), 1.e-3, 'a particle of another state keeps its pose')
        call assert_real(2.,  spproj%os_ptcl3D%get(4, 'e2'), 1.e-3, 'a deselected particle too')
        call assert_int(2, spproj%os_ptcl3D%get_state(2), 'and the states stay')
        ! a result of other rows
        call res%os_ptcl2D%new(3, is_ptcl=.true.)
        call assert_int(-1, map_state_poses(spproj, res, 1), 'a result of other rows is not taken')
        call spproj%kill
        call res%kill
        deallocate(spproj, res)
    end subroutine test_state_runs

    !> the jobs' command lines: each starts afresh in the job's folder and sets the keys of the
    !! method; a mask diameter not positive is left to the job
    subroutine test_job_clines()
        class(cmdline), allocatable :: cline
        type(string) :: vols(2)
        write(*,'(A)') 'test_job_clines'
        allocate(cline)
        call cline%set('stray', 'yes')
        call set_solve2D_route_cline(cline, 'route.simple', 120.)
        call assert_false(cline%defined('stray'), 'a command line starts afresh')
        call assert_true(cline%get_carg('prg') == 'solve2D', 'the route''s solve2D')
        call assert_true(cline%get_carg('mkdir') == 'no' .and. cline%get_carg('projfile') == 'route.simple',&
            &'in the job''s folder, on its project')
        call assert_int(NCLS2D_FIRST, cline%get_iarg('ncls'), 'into NCLS2D_FIRST classes')
        call assert_real(120., cline%get_rarg('mskdiam'), 1.e-4, 'with the mask diameter')
        call assert_true(cline%get_carg('cache') == 'yes', 'cache=yes')
        call set_solve2D_route_cline(cline, 'route.simple', 0.)
        call assert_false(cline%defined('mskdiam'), 'no mask diameter: the job''s default')
        call set_clsexp_cline(cline, 'exp.simple', 500, 110.)
        call assert_true(cline%get_carg('prg') == 'cls_expansion', 'cls_expansion')
        call assert_int(500, cline%get_iarg('ncls'), 'the sub-classes in all')
        call set_cavgs3D_cline(cline, 'c3D.simple', 3, 0.)
        call assert_true(cline%get_carg('prg') == 'solve3D_cavgs', 'solve3D_cavgs')
        call assert_int(3, cline%get_iarg('nstates'), 'the states')
        call assert_true(cline%get_carg('pgrp') == 'c1' .and. cline%get_carg('prune') == 'no', 'c1, no pruning')
        call assert_real(LPSTART_CAVGS3D, cline%get_rarg('lpstart'), 1.e-4, 'the fixed ramp''s start')
        call assert_real(LPSTOP_CAVGS3D,  cline%get_rarg('lpstop'),  1.e-4, 'and stop')
        call assert_true(cline%get_carg('force_lp_range') == 'yes' .and. cline%get_carg('prob_athres_taper') == 'yes',&
            &'forced, with the angular taper')
        vols(1) = 'state_refs/vol_state01.mrc'
        vols(2) = 'state_refs/vol_state02.mrc'
        call set_assign_cline(cline, 'match.simple', vols, 90.)
        call assert_true(cline%get_carg('prg') == 'assign_cavgs_states', 'assign_cavgs_states')
        call assert_int(2, cline%get_iarg('nstates'), 'a state per reference')
        call assert_true(cline%get_carg('vol2') == 'state_refs/vol_state02.mrc', 'each reference as a volume')
        call assert_real(LP_ASSIGN, cline%get_rarg('lp'), 1.e-4, 'at the matching''s limit')
        call assert_int(NSPACE_ASSIGN, cline%get_iarg('nspace'), 'over its directions')
        call set_state_cavgs3D_cline(cline, 'state.simple', 120.)
        call assert_true(cline%get_carg('prg') == 'solve3D_cavgs', 'a state''s own solve3D_cavgs')
        call assert_int(1, cline%get_iarg('nstates'), 'one state')
        call assert_false(cline%defined('lpstart') .or. cline%defined('lpstop') .or. cline%defined('force_lp_range'),&
            &'its limits from the class FRCs, no fixed ramp')
        call set_refine_cont_cline(cline, 'ref.simple', 2, 2500, 90.)
        call assert_true(cline%get_carg('prg') == 'refine3D' .and. cline%get_carg('refine') == 'cont', 'refine3D refine=cont')
        call assert_int(1, cline%get_iarg('maxits'), 'one iteration')
        call assert_int(2, cline%get_iarg('nstates'), 'the states held')
        call assert_int(2500, cline%get_iarg('nspace'), 'the seeds'' grid')
        call assert_true(cline%get_carg('objfun') == 'cc' .and. cline%get_carg('pgrp') == 'c1', 'objfun=cc, c1')
        call cline%kill
        deallocate(cline)
    end subroutine test_job_clines

end module simple_state_assignment_tester
