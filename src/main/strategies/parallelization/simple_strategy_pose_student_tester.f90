!@descr: unit tests of the pose student strategy's helpers: the held-out split leaves both cohorts non-empty, the resampled design comparison of zero groups, the confidence threshold at the end of a tie, the model's physical geometry against the project's, a part's range of the selected particles
submodule (simple_strategy_pose_student) simple_strategy_pose_student_tester
use simple_test_utils
implicit none

contains

    module subroutine run_all_strategy_pose_student_tests()
        write(*,'(A)') '**** running all pose student strategy tests ****'
        call test_holdout_split()
        call test_boot_quantile()
        call test_conf_threshold()
        call test_geometry()
        call test_part_range()
    end subroutine run_all_strategy_pose_student_tests

    !> The held-out cohort is never empty: two to nine stacks and selections without a tenth stack fall back to every
    !! tenth selected particle; ten or more selected stacks hold out every tenth of them in order of stkind, whatever
    !! their numbers; fewer than ten particles on one stack, or on a few, are refused
    subroutine test_holdout_split()
        integer, parameter :: N = 120
        integer :: stk(N), group_of(N), i, nstk, nheld_stk
        logical :: held(N), l_by_stack, l_ok
        character(len=:), allocatable :: msg
        write(*,'(A)') 'test_holdout_split'
        do nstk = 2,9
            do i = 1,N
                stk(i) = 1 + mod(i-1, nstk)
            enddo
            call student_holdout_split(N, stk, held, group_of, l_by_stack, msg)
            l_ok = len(msg) == 0 .and. .not. l_by_stack .and. count(held) == N / 10 .and. count(.not. held) > 0
            call assert_true(l_ok, 'stacks 1 to '//int2str(nstk)//': every tenth selected particle is held out')
        enddo
        ! fifteen selected stacks, none a multiple of ten (11..25 without 20): every tenth by rank, one stack
        do i = 1,N
            stk(i) = 11 + mod(i-1, 15)
            if( stk(i) == 20 ) stk(i) = 26
        enddo
        call student_holdout_split(N, stk, held, group_of, l_by_stack, msg)
        nheld_stk = 0
        do i = 1,N
            if( held(i) .and. stk(i) /= 21 ) nheld_stk = nheld_stk + 1
        enddo
        call assert_true(len(msg) == 0 .and. l_by_stack, 'fifteen selected stacks are held out by stack')
        call assert_true(count(held) > 0 .and. nheld_stk == 0, 'the tenth selected stack (stkind 21), although no stkind is a multiple of ten')
        call assert_true(all(group_of == stk), 'the groups are the stacks')
        ! a single stack with nine selected particles, and nine on several stacks: refused
        stk = 1
        call student_holdout_split(9, stk(1:9), held(1:9), group_of(1:9), l_by_stack, msg)
        call assert_true(len(msg) > 0, 'nine particles on one stack are too few to split')
        do i = 1,9
            stk(i) = i
        enddo
        call student_holdout_split(9, stk(1:9), held(1:9), group_of(1:9), l_by_stack, msg)
        call assert_true(len(msg) > 0, 'nine particles on nine stacks are too few to split')
        stk = 1
        call student_holdout_split(10, stk(1:10), held(1:10), group_of(1:10), l_by_stack, msg)
        call assert_true(len(msg) == 0 .and. count(held(1:10)) == 1, 'ten particles on one stack hold out one')
    end subroutine test_holdout_split

    !> The resampled design comparison: zero groups give no gain (zero, not zero divided by zero), equal gains their
    !! value, a mixed set a quantile within its range
    subroutine test_boot_quantile()
        real :: gain(4), q
        write(*,'(A)') 'test_boot_quantile'
        gain = 0.
        q = student_boot_quantile(gain, 0, 20261126_int64)
        call assert_real(0., q, 0., 'zero groups give no gain')
        gain = 0.25
        q = student_boot_quantile(gain, 4, 20261126_int64)
        call assert_real(0.25, q, 1.e-6, 'equal gains give their value')
        gain = [-1., 0., 0.5, 1.]
        q = student_boot_quantile(gain, 4, 20261126_int64)
        call assert_true(q >= -1. .and. q <= 1., 'a quantile of resampled means lies within the gains')
    end subroutine test_boot_quantile

    !> The threshold ends a group of equal confidences: eighteen recovered and two missed rows at confidence one give
    !! none (the tie recovers 90%); the tie below a recovered row ends the prefix before it; a recovered tie is taken
    !! whole; the rows at or above the threshold are the prefix; ineligible rows are passed over
    subroutine test_conf_threshold()
        integer, parameter :: N = 21
        real    :: conf(N), thr
        logical :: elig(N), recov(N)
        integer :: nsel, nrec
        write(*,'(A)') 'test_conf_threshold'
        conf  = 1.
        elig  = .true.
        recov = .true.
        recov(19:20) = .false.
        call student_conf_threshold(20, conf(1:20), elig(1:20), recov(1:20), 0.95, thr, nsel, nrec)
        call assert_int(0, nsel, 'a tie of eighteen recovered and two missed rows qualifies no prefix')
        call assert_true(thr > 1., 'and the threshold seeds nobody')
        ! one recovered row above the tie
        conf(1:20)  = 0.98
        conf(21)    = 0.999
        recov(21)   = .true.
        call student_conf_threshold(N, conf, elig, recov, 0.95, thr, nsel, nrec)
        call assert_int(1, nsel, 'the prefix ends before the tie it cannot take whole')
        call assert_real(0.999, thr, 0., 'at the higher confidence')
        call assert_int(nsel, count(elig .and. conf >= thr), 'the rows at or above the threshold are the prefix')
        ! the tie fully recovered: it is taken whole
        recov(19:20) = .true.
        call student_conf_threshold(N, conf, elig, recov, 0.95, thr, nsel, nrec)
        call assert_int(N, nsel, 'a recovered tie is taken whole')
        call assert_int(nsel, count(elig .and. conf >= thr), 'and the threshold selects it whole')
        ! an ineligible row above all is passed over
        elig(21) = .false.
        call student_conf_threshold(N, conf, elig, recov, 0.95, thr, nsel, nrec)
        call assert_int(20, nsel, 'an ineligible row is not seeded')
        call assert_int(nsel, count(elig .and. conf >= thr), 'the threshold selects the eligible prefix')
    end subroutine test_conf_threshold

    !> The model's physical geometry against the project's with the same box and half ring: another sampling, crop or
    !! mask is refused; values within the header's decimals pass
    subroutine test_geometry()
        character(len=:), allocatable :: msg
        write(*,'(A)') 'test_geometry'
        msg = student_geometry_mismatch(160, 256, 256, 1.275, 180., .true., 160, 256, 256, 1.275, 180., .true.)
        call assert_true(len(msg) == 0, 'the same geometry passes')
        msg = student_geometry_mismatch(160, 256, 256, 1.27500, 180.000, .true., 160, 256, 256, 1.2750004, 180.0004, .true.)
        call assert_true(len(msg) == 0, 'a difference within the header''s decimals passes')
        msg = student_geometry_mismatch(160, 256, 256, 4., 200., .true., 160, 256, 256, 2., 200., .true.)
        call assert_true(len(msg) > 0, 'another sampling with the same box and half ring is refused')
        msg = student_geometry_mismatch(160, 256, 256, 1.275, 200., .true., 160, 256, 256, 1.275, 100., .true.)
        call assert_true(len(msg) > 0, 'another mask with the same box and half ring is refused')
        msg = student_geometry_mismatch(160, 256, 256, 1.275, 180., .true., 160, 256, 200, 1.275, 180., .true.)
        call assert_true(len(msg) > 0, 'another crop with the same box and half ring is refused')
        msg = student_geometry_mismatch(160, 256, 256, 1.275, 180., .true., 160, 256, 256, 1.275, 180., .false.)
        call assert_true(len(msg) > 0, 'another CTF treatment is refused')
        msg = student_geometry_mismatch(160, 256, 256, 1.275, 180., .true., 150, 256, 256, 1.275, 180., .true.)
        call assert_true(len(msg) > 0, 'another half ring is refused')
    end subroutine test_geometry

    !> A part's range of the selected particles: uneven ranges over a field cover the selection disjointly, a range of
    !! unselected rows holds none (and is no error), a range outside the field or a part outside nparts is refused
    subroutine test_part_range()
        integer, parameter :: NP = 50, NPARTS = 4, RANGES(2,NPARTS) = reshape([1,7, 8,40, 41,41, 42,50], [2,NPARTS])
        integer, allocatable :: sel(:), pinds(:)
        integer :: covered(NP), i, ipart, nsel
        character(len=:), allocatable :: msg
        logical :: l_ok
        write(*,'(A)') 'test_part_range'
        ! the selection: every even particle (41 is odd: part 3 holds only an unselected row)
        nsel = NP / 2
        allocate(sel(nsel))
        do i = 1,nsel
            sel(i) = 2 * i
        enddo
        covered = 0
        l_ok    = .true.
        do ipart = 1,NPARTS
            call student_part_pinds(sel, ipart, NPARTS, RANGES(1,ipart), RANGES(2,ipart), NP, pinds, msg)
            if( len(msg) > 0 ) l_ok = .false.
            do i = 1,size(pinds)
                covered(pinds(i)) = covered(pinds(i)) + 1
                if( pinds(i) < RANGES(1,ipart) .or. pinds(i) > RANGES(2,ipart) ) l_ok = .false.
            enddo
            if( ipart == 3 ) call assert_int(0, size(pinds), 'a part of only unselected rows holds no particle')
        enddo
        call assert_true(l_ok, 'every part keeps to its range without error')
        call assert_true(all(covered(sel) == 1), 'the uneven parts cover every selected particle once')
        call assert_int(nsel, sum(covered), 'and nothing else')
        call student_part_pinds(sel, 1, NPARTS, 0, 7, NP, pinds, msg)
        call assert_true(len(msg) > 0, 'a range starting before the first particle is refused')
        call student_part_pinds(sel, 1, NPARTS, 42, NP + 1, NP, pinds, msg)
        call assert_true(len(msg) > 0, 'a range past the last particle is refused')
        call student_part_pinds(sel, 1, NPARTS, 9, 8, NP, pinds, msg)
        call assert_true(len(msg) > 0, 'an inverted range is refused')
        call student_part_pinds(sel, NPARTS + 1, NPARTS, 1, 7, NP, pinds, msg)
        call assert_true(len(msg) > 0, 'a part beyond nparts is refused')
        deallocate(sel)
    end subroutine test_part_range

end submodule simple_strategy_pose_student_tester
