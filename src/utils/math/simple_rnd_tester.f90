!@descr: unit tests for the random draws of simple_rnd (shuffle, partial_shuffle, multinomal)
! Full and partial Fisher-Yates shuffles keep every value and give a duplicate-free
! prefix; the multinomial draw reproduces its probabilities within four binomial
! standard deviations (1000 draws: 0.013 at p = 0.8, 0.009 at p = 0.1) and never draws
! an entry of probability zero. Every test starts from a fixed seed.
module simple_rnd_tester
use simple_rnd,        only: shuffle, partial_shuffle, multinomal
use simple_math,       only: hpsort
use simple_test_utils
implicit none
private
public :: run_all_rnd_tests

integer, parameter :: N = 100, NPARTIAL = 17, NDRAWS = 1000

contains

    subroutine run_all_rnd_tests()
        write(*,'(A)') '**** running all random draw tests ****'
        call test_shuffle()
        call test_partial_shuffle()
        call test_multinomal()
    end subroutine run_all_rnd_tests

    !> a full shuffle is a permutation of its input, and it does permute
    subroutine test_shuffle()
        integer :: iarr(N), iarr_copy(N), expected(N), i
        real    :: rarr(N), rarr_copy(N), rexpected(N)
        write(*,'(A)') 'test_shuffle'
        call set_fixed_seed(20260929)
        expected  = [(i,i=1,N)]
        rexpected = real(expected)
        iarr = expected
        call shuffle(iarr)
        call assert_true(any(iarr /= expected), 'integer full shuffle changes the order')
        iarr_copy = iarr
        call hpsort(iarr_copy)
        call assert_true(all(iarr_copy == expected), 'integer full shuffle preserves the input permutation')
        rarr = rexpected
        call shuffle(rarr)
        call assert_true(any(rarr /= rexpected), 'real full shuffle changes the order')
        rarr_copy = rarr
        call hpsort(rarr_copy)
        call assert_true(all(rarr_copy == rexpected), 'real full shuffle preserves the input permutation')
    end subroutine test_shuffle

    !> a partial shuffle selects an ordered, duplicate-free prefix and keeps every value
    subroutine test_partial_shuffle()
        integer :: iarr(N), iarr_copy(N), expected(N), i
        real    :: rarr(N), rarr_copy(N), rexpected(N)
        write(*,'(A)') 'test_partial_shuffle'
        call set_fixed_seed(20260930)
        expected  = [(i,i=1,N)]
        rexpected = real(expected)
        iarr = expected
        call partial_shuffle(iarr, 0)
        call assert_true(all(iarr == expected), 'zero-length integer partial shuffle is a no-op')
        call partial_shuffle(iarr, NPARTIAL)
        iarr_copy = iarr
        call hpsort(iarr_copy)
        call assert_true(all(iarr_copy == expected), 'integer partial shuffle preserves the input permutation')
        iarr_copy(:NPARTIAL) = iarr(:NPARTIAL)
        call hpsort(iarr_copy(:NPARTIAL))
        call assert_true(all(iarr_copy(2:NPARTIAL) /= iarr_copy(1:NPARTIAL-1)), &
            &'integer partial shuffle prefix contains no duplicates')
        rarr = rexpected
        call partial_shuffle(rarr, NPARTIAL)
        rarr_copy = rarr
        call hpsort(rarr_copy)
        call assert_true(all(rarr_copy == rexpected), 'real partial shuffle preserves the input permutation')
        ! selecting the complete prefix is a full Fisher-Yates permutation
        iarr = expected
        call partial_shuffle(iarr, N)
        iarr_copy = iarr
        call hpsort(iarr_copy)
        call assert_true(all(iarr_copy == expected), 'complete partial shuffle preserves the input permutation')
    end subroutine test_partial_shuffle

    !> the draw frequencies follow the probabilities; a zero probability is never drawn
    subroutine test_multinomal()
        real    :: pvec(10), pzero(4), freq
        integer :: i, cnt, which
        logical :: l_zero_drawn
        write(*,'(A)') 'test_multinomal'
        call set_fixed_seed(20260926)
        pvec(1)  = 0.8
        pvec(2:) = 0.2/9.
        cnt = 0
        do i = 1, NDRAWS
            if( multinomal(pvec) == 1 ) cnt = cnt + 1
        end do
        freq = real(cnt)/real(NDRAWS)
        call assert_real(0.8, freq, 0.05, 'the dominant entry (p = 0.8) is drawn at its probability')
        pvec = 0.1
        cnt  = 0
        do i = 1, NDRAWS
            if( multinomal(pvec) == 1 ) cnt = cnt + 1
        end do
        freq = real(cnt)/real(NDRAWS)
        call assert_real(0.1, freq, 0.04, 'one of ten equal entries (p = 0.1) is drawn at its probability')
        pzero = [0., 0.5, 0., 0.5]
        l_zero_drawn = .false.
        do i = 1, NDRAWS
            which = multinomal(pzero)
            if( which == 1 .or. which == 3 ) l_zero_drawn = .true.
        end do
        call assert_false(l_zero_drawn, 'an entry of probability zero is never drawn')
    end subroutine test_multinomal

end module simple_rnd_tester
