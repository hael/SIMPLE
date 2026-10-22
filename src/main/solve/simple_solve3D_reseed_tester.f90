!@descr: unit tests of solve3D's stage rules: the axis-search stage, stage skipping, symmetry randomization, and solve3D_cavgs' state reseeding
! In-memory labels and scores of the even (1..ncavgs) and odd (ncavgs+1..2*ncavgs) halves of the
! class averages; no volume or job. The symmetry randomization test writes one small project file.
module simple_solve3D_reseed_tester
use simple_solve3D_utils, only: reseed_state_labels, solve3D_symsrch_stage, solve3D_stage_skipping, &
    &symmetry_randomize_poses
use simple_test_utils
implicit none
private
public :: run_all_solve3D_reseed_tests

contains

    subroutine run_all_solve3D_reseed_tests()
        write(*,'(A)') '**** running all solve3D_cavgs stage rule tests ****'
        call test_symsrch_stage()
        call test_stage_skipping()
        call test_symmetry_randomize()
        call test_no_weak_state()
        call test_one_empty_state()
        call test_all_in_one_state()
        call test_weak_by_fraction()
        call test_too_few_to_reseed()
    end subroutine run_all_solve3D_reseed_tests

    !> the symmetry axis search: after stage 4 in a run that searches an axis, stage 3 otherwise
    subroutine test_symsrch_stage()
        write(*,'(A)') 'test_symsrch_stage'
        call assert_int(3, solve3D_symsrch_stage(.false.), 'no axis search: stage 3')
        call assert_int(4, solve3D_symsrch_stage(.true.),  'axis search after stage 4')
    end subroutine test_symsrch_stage

    !> stage omission and the mid-stage resolution exit: single-state runs only
    subroutine test_stage_skipping()
        write(*,'(A)') 'test_stage_skipping'
        call assert_true(solve3D_stage_skipping(1),        'single state: stages may be skipped')
        call assert_false(solve3D_stage_skipping(2),       'two states: no stage skipped')
        call assert_false(solve3D_stage_skipping(4),       'four states: no stage skipped')
    end subroutine test_stage_skipping

    !> the randomization acts on the poses the stage wrote into the project file, not on the in-memory
    !! copy of the previous stage: every pose in the file ends as itself or its c2 mate
    subroutine test_symmetry_randomize()
        use simple_sp_project, only: sp_project
        use simple_sym,        only: sym
        use simple_ori,        only: ori
        use simple_string,     only: string
        integer,          parameter :: NPTCLS = 12
        character(len=*), parameter :: PROJFILE = 'symrand.simple'
        class(sp_project), allocatable :: proj_file, proj_stale, proj_back
        type(sym)    :: symobj
        type(ori)    :: o_ref, o_mate, o_back, o_stale
        type(string) :: cwd, root
        integer      :: i, nfail0, nkept, nmoved, nstale
        real         :: mat_back(3,3), mat_stale(3,3,NPTCLS)
        write(*,'(A)') 'test_symmetry_randomize'
        nfail0 = tests_failed
        call enter_fixture('symrand', cwd, root)
        allocate(proj_file, proj_stale, proj_back)
        call proj_file%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            call proj_file%os_ptcl3D%set_euler(i, [real(mod(37*i,360)), 20.+6.*real(i), real(mod(53*i,360))])
        enddo
        call proj_file%update_projinfo(string(PROJFILE))
        call proj_file%write(string(PROJFILE))
        ! the in-memory project is the previous stage's: other poses than the file's
        call proj_stale%read(string(PROJFILE))
        do i = 1, NPTCLS
            call proj_stale%os_ptcl3D%set_euler(i, [real(mod(91*i+7,360)), 100.-5.*real(i), real(mod(29*i+3,360))])
        enddo
        do i = 1, NPTCLS
            call proj_stale%os_ptcl3D%get_ori(i, o_stale)
            mat_stale(:,:,i) = o_stale%get_mat()
        enddo
        call symobj%new('c2')
        call symmetry_randomize_poses(symobj, proj_stale, string(PROJFILE))
        call proj_back%read(string(PROJFILE))
        nkept = 0; nmoved = 0; nstale = 0
        do i = 1, NPTCLS
            call proj_file%os_ptcl3D%get_ori(i, o_ref)
            call symobj%apply(o_ref, 2, o_mate)
            call proj_back%os_ptcl3D%get_ori(i, o_back)
            mat_back = o_back%get_mat()
            if( maxval(abs(mat_back - o_ref%get_mat())) < 1.e-4 )  nkept  = nkept  + 1
            if( maxval(abs(mat_back - o_mate%get_mat())) < 1.e-4 ) nmoved = nmoved + 1
            if( maxval(abs(mat_back - mat_stale(:,:,i))) < 1.e-4 ) nstale = nstale + 1
        enddo
        call assert_int(NPTCLS, nkept + nmoved, 'each pose of the file ends as itself or its mate')
        call assert_int(0, nstale, 'the in-memory poses of the previous stage do not reach the file')
        call proj_file%kill; call proj_stale%kill; call proj_back%kill
        deallocate(proj_file, proj_stale, proj_back)
        call symobj%kill
        call leave_fixture(cwd, root, nfail0)
    end subroutine test_symmetry_randomize

    !> populated states, and a single state, are left as they are
    subroutine test_no_weak_state()
        integer, parameter :: NCAVGS = 12
        integer :: states(2*NCAVGS), nmoved
        real    :: scores(2*NCAVGS)
        integer :: i
        write(*,'(A)') 'test_no_weak_state'
        do i = 1, NCAVGS
            states(i)        = mod(i - 1, 3) + 1
            states(NCAVGS+i) = states(i)
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(0, nmoved, 'three populated states: nothing moves')
        states = 1
        call reseed_state_labels(states, scores, NCAVGS, 1, nmoved)
        call assert_int(0, nmoved, 'one state: nothing to reseed')
    end subroutine test_no_weak_state

    !> an empty state takes the donor's worst-fitting classes, both halves together, up to half
    !! the donor; deselected classes and the other state stay as they are
    subroutine test_one_empty_state()
        integer, parameter :: NCAVGS = 12
        real,    parameter :: SCORES1(6) = [0.9, 0.2, 0.8, 0.1, 0.7, 0.3] ! classes 4, 2 and 6 fit worst
        integer :: states(2*NCAVGS), expected(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_one_empty_state'
        ! classes 1-6 in state 1, 7-11 in state 2, 12 deselected, state 3 empty
        states(1:6)   = 1
        states(7:11)  = 2
        states(12)    = 0
        states(NCAVGS+1:2*NCAVGS) = states(1:NCAVGS)
        scores        = 0.5
        scores(1:6)   = SCORES1
        scores(NCAVGS+1:NCAVGS+6) = SCORES1
        expected = states
        do i = 1, NCAVGS
            if( any(i == [2, 4, 6]) )then
                expected(i)        = 3
                expected(NCAVGS+i) = 3
            endif
        enddo
        ! 22 selected halves: an equal share is 7, half the donor's 12 is 6
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(6, nmoved, 'half the donor moves')
        call assert_true(all(states == expected), 'the three worst-fitting classes of state 1, both halves')
        call assert_true(states(12) == 0 .and. states(2*NCAVGS) == 0, 'the deselected class stays deselected')
    end subroutine test_one_empty_state

    !> every class in one state: the two empty states each take an equal share, the worst-fitting first
    subroutine test_all_in_one_state()
        integer, parameter :: NCAVGS = 12
        integer :: states(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_all_in_one_state'
        states = 1
        do i = 1, NCAVGS
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(16, nmoved, 'two equal shares of 8 move')
        call assert_true(all(states(1:4) == 2) .and. all(states(NCAVGS+1:NCAVGS+4) == 2),&
            &'state 2 takes the four worst-fitting classes')
        call assert_true(all(states(5:8) == 3) .and. all(states(NCAVGS+5:NCAVGS+8) == 3),&
            &'state 3 the next four')
        call assert_true(all(states(9:12) == 1) .and. all(states(NCAVGS+9:NCAVGS+12) == 1),&
            &'state 1 keeps the best-fitting four')
    end subroutine test_all_in_one_state

    !> a state above the search's floor but under 2% of the selected halves is reseeded
    subroutine test_weak_by_fraction()
        integer, parameter :: NCAVGS = 200
        integer :: states(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_weak_by_fraction'
        states(1:100)   = 1
        states(101:197) = 2
        states(198:200) = 3 ! 6 halves of 400
        states(NCAVGS+1:2*NCAVGS) = states(1:NCAVGS)
        do i = 1, NCAVGS
            scores(i)        = 0.5 + real(i) / 1000.
            scores(NCAVGS+i) = scores(i)
        enddo
        ! an equal share is 133, half the donor's 200 is 100
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(100, nmoved, 'half the donor moves')
        call assert_int(106, count(states == 3), 'state 3 holds its own and the moved halves')
        call assert_true(all(states(1:50) == 3) .and. all(states(51:100) == 1),&
            &'the worst-fitting 50 classes of state 1 move')
    end subroutine test_weak_by_fraction

    !> too few classes to leave a reseeded state above the search's floor: nothing moves, so no
    !! state is left empty while another is relabelled
    subroutine test_too_few_to_reseed()
        integer, parameter :: NCAVGS = 4
        integer :: states(2*NCAVGS), states_in(2*NCAVGS), nmoved, i
        real    :: scores(2*NCAVGS)
        write(*,'(A)') 'test_too_few_to_reseed'
        states = 1
        do i = 1, NCAVGS
            scores(i)        = real(i) / 100.
            scores(NCAVGS+i) = scores(i)
        enddo
        states_in = states
        call reseed_state_labels(states, scores, NCAVGS, 3, nmoved)
        call assert_int(0, nmoved, 'nothing moves')
        call assert_true(all(states == states_in), 'the labels are unchanged')
    end subroutine test_too_few_to_reseed

end module simple_solve3D_reseed_tester
