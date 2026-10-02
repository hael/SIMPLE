!@descr: unit tests for simple_rnd integer draws, shuffles and multinomial sampling
! Full and partial Fisher-Yates shuffles keep every value and give a duplicate-free
! prefix; the multinomial draw reproduces its probabilities within four binomial
! standard deviations (1000 draws: 0.013 at p = 0.8, 0.009 at p = 0.1) and never draws
! an entry of probability zero. Every test starts from a fixed seed.
module simple_rnd_tester
use simple_defs,       only: dp
use simple_rnd,        only: irnd_uni, shuffle, partial_shuffle, multinomal
use simple_math,       only: hpsort
use simple_test_utils
implicit none
private
public :: run_all_rnd_tests

integer, parameter :: N = 100, NPARTIAL = 17, NDRAWS = 1000

contains

    subroutine run_all_rnd_tests()
        write(*,'(A)') '**** running all random draw tests ****'
        call test_irnd_uni()
        call test_shuffle()
        call test_partial_shuffle()
        call test_multinomal()
    end subroutine run_all_rnd_tests

    !> Guard against single-precision rounding that makes valid indices unreachable for large NP.
    !! Check bounds around 2**23 and 2**24, and huge(0), by comparing each sampled index with
    !! the probability interval containing the same underlying double-precision random draw.
    subroutine test_irnd_uni()
        integer, parameter :: BOUNDS(7) = [2**23-1, 2**23, 2**23+1, 2**24-1, 2**24, 2**24+1, huge(0)]
        integer, parameter :: NBOUND_DRAWS = 256
        integer, allocatable :: seed(:)
        integer :: nseed, ibound, idraw, np, which
        real(dp) :: harvest, lower, upper
        logical :: in_range, correct_bins
        character(len=128) :: msg
        write(*,'(A)') 'test_irnd_uni'
        call set_fixed_seed(20261002)
        call assert_int(1, irnd_uni(1), 'a singleton integer range returns one')
        call random_seed(size=nseed)
        allocate(seed(nseed))
        do ibound = 1,size(BOUNDS)
            np = BOUNDS(ibound)
            in_range = .true.
            correct_bins = .true.
            do idraw = 1,NBOUND_DRAWS
                call random_seed(get=seed)
                which = irnd_uni(np)
                ! Replay the same draw: index k must correspond to a value in [(k-1)/NP, k/NP).
                ! This detects rounding errors without waiting to randomly hit a rare endpoint.
                call random_seed(put=seed)
                call random_number(harvest)
                if( which < 1 .or. which > np )then
                    in_range = .false.
                    correct_bins = .false.
                    cycle
                endif
                lower = real(which-1,dp) / real(np,dp)
                upper = real(which,dp)   / real(np,dp)
                if( harvest < lower .or. harvest >= upper ) correct_bins = .false.
            enddo
            write(msg,'(A,I0)') 'integer draws stay within [1, NP] for NP=', np
            call assert_true(in_range, trim(msg))
            write(msg,'(A,I0)') 'integer draws match equal-width probability bins for NP=', np
            call assert_true(correct_bins, trim(msg))
        enddo
        deallocate(seed)
    end subroutine test_irnd_uni

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
