!@descr: regression of the coarse candidates retained by the polar volume symmetry-axis search
submodule (simple_volpft_symsrch) simple_volpft_symsrch_tester
use simple_test_utils, only: assert_int, assert_real, assert_true
implicit none

contains

    module subroutine run_all_volpft_symsrch_tests()
        write(*,'(A)') '**** running all symmetry-axis candidate tests ****'
        call test_distinct_scored_candidates()
    end subroutine run_all_volpft_symsrch_tests

    subroutine test_distinct_scored_candidates()
        integer, parameter :: NAXES = 24, NSCORED = 20
        ! Binary fractions give exact scores; the descending indices are ordered by hand.
        real, parameter :: SCORES(NSCORED) = real([1,17,3,15,5,13,20,7,14,4,12,10,19,2,11,9,18,6,16,8]) / 32.
        integer, parameter :: EXPECTED_TOP(10) = [7,13,17,2,19,4,9,6,11,15]
        class(oris), allocatable :: espace, candidates
        integer, allocatable     :: ranked(:)
        integer                  :: i
        real                     :: euler(3), input_axes(3,NAXES)
        write(*,'(A)') 'test_distinct_scored_candidates'
        allocate(espace, candidates)
        call espace%new(NAXES, is_ptcl=.false.)
        do i = 1, NAXES
            call espace%set_euler(i, real([7*i,10+i,3*i]))
            ! set_euler canonicalizes through a rotation matrix; retain the stored input exactly.
            input_axes(:,i) = espace%get_euler(i)
        enddo
        call fill_symaxis_candidates(espace, [1,NSCORED], SCORES, candidates)
        do i = 1, NSCORED
            euler = candidates%get_euler(i)
            call assert_true(all(euler == input_axes(:,i)), 'each scored candidate keeps its own axis')
            call assert_real(SCORES(i), candidates%get(i,'corr'), 0., 'each candidate keeps its own exact score')
        enddo
        do i = NSCORED+1, NAXES
            call assert_real(-1., candidates%get(i,'corr'), 0., 'unsearched axes remain excluded from refinement')
        enddo
        ranked = candidates%order()
        call assert_int(size(EXPECTED_TOP), NBEST, 'the symmetry search refines ten coarse candidates')
        do i = 1, size(EXPECTED_TOP)
            call assert_int(EXPECTED_TOP(i), ranked(i), 'refinement starts are the ten highest-scoring axes')
        enddo
        call candidates%kill
        call espace%kill
        deallocate(espace, candidates, ranked)
    end subroutine test_distinct_scored_candidates

end submodule simple_volpft_symsrch_tester
