!@descr: unit tests for simple_rnd integer draws, shuffles, multinomial sampling and the keyed generator
! Full and partial Fisher-Yates shuffles keep every value and give a duplicate-free
! prefix; the multinomial draw reproduces its probabilities within four binomial
! standard deviations (1000 draws: 0.013 at p = 0.8, 0.009 at p = 0.1) and never draws
! an entry of probability zero; the Poisson deviates of both branches (means 4 and 250) have
! their mean and variance within four standard errors. Every test starts from a fixed seed.
! The keyed generator of the probabilistic alignment (rnd_keyed, prob_align_batch_pipeline.md section 4):
! ten pinned vectors (key, 64-bit word, unit variate) from an independent Python reference, the
! unit-interval mapping's endpoints, independence of the order of draws, the sensitivity of the word to
! every field of the key, uniformity and no serial correlation over consecutive keys, the keyed integer
! draw's bins, and the keyed partial Fisher-Yates permutation's distribution (every subset and every order
! equally likely over many keys). The vectors below were produced on the Dell (GCC 15.2.1, 2026-11-04) and
! agree between the Debug and the Release build; they are recorded in the plan's Progress for the Mac check.
module simple_rnd_tester
use simple_defs,       only: dp
use simple_rnd,        only: irnd_uni, shuffle, partial_shuffle, multinomal, seed_rnd, seed_rnd_fixed, has_global_seed, &
    &poisson_rnd, ran3, rnd_key, rnd_keyed, rnd_keyed_bits, rnd_keyed_int, rnd_keyed_partial_shuffle, rnd_key_lanes, &
    &rnd_bits_to_unit, rnd_scope, RND_KEY_NLANES, RND_FAM_INPL, RND_FAM_DISK, RND_FAM_CAND, RND_FAM_ASSIGN, &
    &RND_PROG_REFINE3D, RND_PROG_REFINE2D, RND_PASS_ALIGN, RND_KEYED_VERSION, rnd_run_seed, ensure_run_seed, &
    &rnd_scope_key, rnd_draw_key, get_rnd_state, put_rnd_state
use simple_math,       only: hpsort
use simple_string_utils, only: int2str
use simple_error,        only: simple_exception
use simple_test_utils
use, intrinsic :: iso_fortran_env, only: int64
implicit none
private
public :: run_all_rnd_tests
#include "simple_local_flags.inc"

integer, parameter :: N = 100, NPARTIAL = 17, NDRAWS = 1000
! the pinned vectors of the keyed generator: seed, (version, scope, iteration, particle, family, reference,
! ordinal), the 64-bit word and the unit variate (scratch/keep/rnd/vectors_python.txt of the run)
integer, parameter :: NVEC = 10
integer(int64), parameter :: VEC_SEED(NVEC) = [0_int64, 20261103_int64, 20261103_int64, 20261103_int64, 20261103_int64, &
    &20261103_int64, -1_int64, 4611686018427400249_int64, 123456789012345_int64, 20261103_int64]
integer, parameter :: VEC_FIELDS(7,NVEC) = reshape([ &
    &1,  0, 0,      0, 0,     0,          0, &
    &1, 17, 1,      1, 2,     1,          0, &
    &1, 17, 1,      2, 2,     1,          0, &
    &1, 17, 1,      1, 2,     2,          0, &
    &1, 17, 2,      1, 2,     1,          0, &
    &1, 17, 1, 750000, 5,     0,          0, &
    &1, 33, 3,  12345, 3,     0,          7, &
    &1, 18, 1,      1, 4, 20000, 2147483647, &
    &1, 17, 25, 35645, 1,     0,      35645, &
    &1, 17, 1,      1, 6,     0,          1], [7,NVEC])
integer(int64), parameter :: VEC_WORD(NVEC) = [-4131892048949355934_int64, -3402735607585271042_int64, &
    &4682640351440764835_int64, -2904353384885242822_int64, 6790283504847523914_int64, -4401584271307322743_int64, &
    &-2048942294402738840_int64, 3588558501286480495_int64, -3713091948123483349_int64, -5700306898616924037_int64]
real(dp), parameter :: VEC_UNIT(NVEC) = [0.7760096832026763_dp, 0.8155373330931133_dp, 0.25384644210001805_dp, &
    &0.8425546875220895_dp, 0.36810200638740853_dp, 0.7613896385335287_dp, 0.8889266156555558_dp, &
    &0.19453614615930637_dp, 0.7987128821603042_dp, 0.6909857438342711_dp]
integer, parameter :: RND_SCOPE_EXPECTED = 17 !< rnd_scope(RND_PROG_REFINE3D, RND_PASS_ALIGN) = 16 + 1

contains

    subroutine run_all_rnd_tests()
        write(*,'(A)') '**** running all random draw tests ****'
        call test_global_seed()
        call test_irnd_uni()
        call test_shuffle()
        call test_partial_shuffle()
        call test_multinomal()
        call test_poisson()
        call test_keyed_vectors()
        call test_keyed_mapping()
        call test_keyed_order()
        call test_keyed_sensitivity()
        call test_keyed_int()
        call test_keyed_partial_shuffle()
        call test_keyed_scope()
        call test_rnd_state_restore()
    end subroutine run_all_rnd_tests

    !> a routine that reseeds (a parameters object of its own) and restores the saved state leaves its caller's
    !! stream as if it had not run: the next draw and the reseed count are the caller's
    subroutine test_rnd_state_restore()
        integer, allocatable :: seed(:)
        integer :: ncalls
        real    :: r_ref, r_after
        write(*,'(A)') 'test_rnd_state_restore'
        call seed_rnd_fixed(31337)
        r_ref = ran3()
        call seed_rnd_fixed(31337)
        call get_rnd_state(seed, ncalls)
        call seed_rnd_fixed(4242)
        r_after = ran3()
        call put_rnd_state(seed, ncalls)
        r_after = ran3()
        call assert_true(r_after == r_ref, 'put_rnd_state restores the caller''s next draw')
    end subroutine test_rnd_state_restore

    !> Library reseeding and intrinsic draws must replay from the propagated root seed.
    subroutine test_global_seed()
        real :: first(16), replay(16), next(16), changed(16)
        character(len=64) :: seed_value
        call seed_rnd_fixed(1729, propagate=.true.)
        call assert_true(has_global_seed(), 'propagated seed enables reproducible library mode')
        call get_environment_variable('SIMPLE_SEED', seed_value)
        call assert_char('1729', trim(seed_value), 'library seed is inherited by child processes')
        call seed_rnd_fixed(99)
        call get_environment_variable('SIMPLE_SEED', seed_value)
        call assert_char('1729', trim(seed_value), 'local seeding preserves the global root seed')
        call seed_rnd_fixed(100, propagate=.false.)
        call get_environment_variable('SIMPLE_SEED', seed_value)
        call assert_char('1729', trim(seed_value), 'disabled propagation preserves the global root seed')
        call seed_rnd()
        call random_number(first)
        call seed_rnd()
        call random_number(next)
        call seed_rnd_fixed(1729, propagate=.true.)
        call seed_rnd()
        call random_number(replay)
        call assert_true(all(first == replay), 'global seed replays intrinsic random draws')
        call seed_rnd()
        call random_number(replay)
        call assert_true(all(next == replay), 'global seed replays successive library reseeds')
        call assert_true(any(first /= next), 'successive library reseeds use distinct streams')
        call seed_rnd_fixed(1730, propagate=.true.)
        call seed_rnd()
        call random_number(changed)
        call assert_true(any(first /= changed), 'changing the global seed changes the stream')
        call fix_test_seed()
    end subroutine test_global_seed

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

    !> Poisson deviates: zero at mean zero, never negative, and the sample mean and variance of NPOIS
    !! draws within four standard errors of the mean (sqrt(mu/n) for the mean, about mu sqrt(2/n) for
    !! the variance) on both branches, multiplication below 30 and the transformed rejection above
    subroutine test_poisson()
        integer, parameter :: NPOIS = 20000
        real,    parameter :: MEANS(2) = [4.0, 250.0]
        real(dp) :: s1, s2, smean, svar
        real     :: tol
        integer  :: i, j, k, nneg
        call fix_test_seed()
        call assert_int(0, poisson_rnd(0.0), 'a Poisson deviate of mean zero is zero')
        do j = 1, size(MEANS)
            s1 = 0.0_dp; s2 = 0.0_dp; nneg = 0
            do i = 1, NPOIS
                k = poisson_rnd(MEANS(j))
                if( k < 0 ) nneg = nneg + 1
                s1 = s1 + real(k,dp)
                s2 = s2 + real(k,dp)**2
            enddo
            smean = s1 / real(NPOIS,dp)
            svar  = s2 / real(NPOIS,dp) - smean**2
            call assert_int(0, nneg, 'Poisson deviates are never negative')
            call assert_true(abs(smean - MEANS(j)) < 4.0 * sqrt(MEANS(j) / real(NPOIS)), &
                &'Poisson deviates have the requested mean')
            tol = 4.0 * sqrt(2.0 / real(NPOIS)) + 4.0 / sqrt(MEANS(j) * real(NPOIS))
            call assert_true(abs(svar / MEANS(j) - 1.0_dp) < tol, 'Poisson deviates have a variance equal to their mean')
        enddo
    end subroutine test_poisson

    ! ---- the keyed generator (prob_align_batch_pipeline.md, section 4; level 1 of its section 8) ----

    !> The ten pinned vectors: the key, its 64-bit word and its unit variate, from the independent Python
    !! reference (scratch/keep/rnd/rnd_keyed_ref.py of the run prob_pipeline, phase 0, 2026-11-04). The words
    !! are printed so that a build's values can be recorded; Debug and Release of GCC 15.2 on the Dell agree.
    subroutine test_keyed_vectors()
        type(rnd_key)  :: key
        integer(int64) :: word
        real(dp)       :: u
        integer        :: i
        character(len=256) :: line
        write(*,'(A)') 'test_keyed_vectors'
        write(*,'(A)') '  seed version scope iteration particle family reference ordinal | word | unit variate'
        do i = 1, NVEC
            key  = vec_key(i)
            word = rnd_keyed_bits(key)
            u    = rnd_keyed(key)
            write(line,'(i0,7(1x,i0),a,i0,a,es25.17)') key%seed, key%version, key%scope, key%iteration, key%particle,&
                &key%family, key%reference, key%ordinal, ' | ', word, ' | ', u
            write(*,'(2x,a)') trim(line)
            call assert_true(word == VEC_WORD(i), 'keyed vector '//int2str(i)//' gives the pinned 64-bit word')
            call assert_true(u == VEC_UNIT(i),    'keyed vector '//int2str(i)//' gives the pinned unit variate')
            call assert_true(u == rnd_bits_to_unit(word), 'the unit variate is the word''s top 53 bits over 2^53')
            call assert_true(u >= 0._dp .and. u < 1._dp, 'the unit variate lies in [0,1)')
        enddo
    end subroutine test_keyed_vectors

    !> the ith pinned key
    function vec_key( i ) result( key )
        integer, intent(in) :: i
        type(rnd_key) :: key
        key%seed      = VEC_SEED(i)
        key%version   = VEC_FIELDS(1,i)
        key%scope     = VEC_FIELDS(2,i)
        key%iteration = VEC_FIELDS(3,i)
        key%particle  = VEC_FIELDS(4,i)
        key%family    = VEC_FIELDS(5,i)
        key%reference = VEC_FIELDS(6,i)
        key%ordinal   = VEC_FIELDS(7,i)
    end function vec_key

    !> The unit-interval mapping's endpoints: the all-zero word gives 0, the all-one word 1 - 2^-53 (below 1),
    !! the low eleven bits are dropped (2^11 is the smallest word above zero, 2^11 - 1 still zero), and the
    !! lanes are the key's fields as unsigned 32-bit words in the pinned order (a negative seed gives its
    !! two's-complement words)
    subroutine test_keyed_mapping()
        type(rnd_key)  :: key
        integer(int64) :: lanes(RND_KEY_NLANES)
        write(*,'(A)') 'test_keyed_mapping'
        call assert_true(rnd_bits_to_unit(0_int64) == 0._dp, 'the zero word maps to 0')
        call assert_true(rnd_bits_to_unit(-1_int64) == 1._dp - 2._dp**(-53), 'the all-one word maps to 1 - 2^-53')
        call assert_true(rnd_bits_to_unit(-1_int64) < 1._dp, 'the largest variate is below 1')
        call assert_true(rnd_bits_to_unit(2048_int64) == 2._dp**(-53), 'the word 2^11 maps to the smallest positive variate')
        call assert_true(rnd_bits_to_unit(2047_int64) == 0._dp, 'the low eleven bits of the word are dropped')
        call assert_true(rnd_bits_to_unit(huge(0_int64)) == 0.5_dp - 2._dp**(-53), 'the largest signed word maps below one half')
        key%seed = -1_int64; key%version = 7; key%scope = 33; key%iteration = 3; key%particle = 12345
        key%family = RND_FAM_DISK; key%reference = 2; key%ordinal = 7
        lanes = rnd_key_lanes(key)
        call assert_true(lanes(1) == 4294967295_int64 .and. lanes(2) == 4294967295_int64, 'a negative seed gives two all-one lanes')
        call assert_true(lanes(3) == 7_int64 .and. lanes(4) == 33_int64 .and. lanes(5) == 3_int64 .and. &
            &lanes(6) == 12345_int64 .and. lanes(7) == 3_int64 .and. lanes(8) == 2_int64 .and. lanes(9) == 7_int64, &
            &'the lanes follow the key''s fields in the pinned order')
        key%seed = 4611686018427400249_int64
        lanes = rnd_key_lanes(key)
        call assert_true(lanes(1) == 12345_int64 .and. lanes(2) == 1073741824_int64, 'the seed splits into its low and high words')
        call assert_int(RND_SCOPE_EXPECTED, rnd_scope(RND_PROG_REFINE3D, RND_PASS_ALIGN), 'the scope packs the program and the pass kind')
    end subroutine test_keyed_mapping

    !> Independence of order: the variates of a key sequence are the same drawn in order, in reverse, and
    !! interleaved with intrinsic draws and reseeds; the keyed draws leave the intrinsic stream untouched
    subroutine test_keyed_order()
        integer, parameter :: NKEYS = 200
        type(rnd_key) :: key
        real(dp) :: inorder(NKEYS), reversed(NKEYS), interleaved(NKEYS)
        real     :: before(8), after(8), replay(8)
        integer  :: i
        write(*,'(A)') 'test_keyed_order'
        key%seed = 20261104_int64; key%scope = rnd_scope(RND_PROG_REFINE2D, RND_PASS_ALIGN); key%iteration = 2
        key%family = RND_FAM_ASSIGN
        do i = 1, NKEYS
            key%particle = i
            inorder(i) = rnd_keyed(key)
        enddo
        do i = NKEYS, 1, -1
            key%particle = i
            reversed(i) = rnd_keyed(key)
        enddo
        call assert_true(all(inorder == reversed), 'a key sequence gives the same variates drawn in reverse')
        call set_fixed_seed(20261104)
        do i = 1, NKEYS
            key%particle = i
            if( mod(i, 7) == 0 ) call set_fixed_seed(i)
            before(1) = ran3()
            interleaved(i) = rnd_keyed(key)
        enddo
        call assert_true(all(inorder == interleaved), 'intrinsic draws and reseeds between keyed draws change nothing')
        call set_fixed_seed(20261105)
        call random_number(before)
        call set_fixed_seed(20261105)
        do i = 1, NKEYS
            key%particle = i
            interleaved(i) = rnd_keyed(key)
        enddo
        call random_number(after)
        call assert_true(all(before == after), 'keyed draws leave the intrinsic stream untouched')
        call set_fixed_seed(20261105)
        call random_number(replay)
        call assert_true(all(replay == before), 'the intrinsic stream replays from its seed')
    end subroutine test_keyed_order

    !> Every field of the key moves the word (version, scope, iteration, particle, family, reference, ordinal and
    !! both seed words); consecutive particle keys give distinct variates; and over a long run of consecutive
    !! keys the variates are uniform (64 bins, chi-square of 63 degrees of freedom below 110, which a uniform
    !! source exceeds with probability 2e-4) with no serial correlation (|lag-1| below 0.02 where its standard
    !! error is 0.003)
    subroutine test_keyed_sensitivity()
        integer, parameter :: NBASE = 20, NRUN = 100000, NBINS = 64
        type(rnd_key)  :: key, kmod
        integer(int64) :: word, words(NRUN)
        integer        :: i, ifield, nsame, keys31(NRUN), inds(NRUN), counts(NBINS), ndup
        real(dp)       :: u(NRUN), mean, num, den, chi2, expected
        write(*,'(A)') 'test_keyed_sensitivity'
        nsame = 0
        key%seed = 20261103_int64; key%scope = rnd_scope(RND_PROG_REFINE3D, RND_PASS_ALIGN); key%iteration = 1
        key%family = RND_FAM_INPL
        do i = 1, NBASE
            key%particle  = 10 * i
            key%reference = 3 * i
            word = rnd_keyed_bits(key)
            do ifield = 1, 9
                kmod = key
                select case(ifield)
                    case(1); kmod%seed      = kmod%seed + 1_int64
                    case(2); kmod%seed      = kmod%seed + 4294967296_int64
                    case(3); kmod%version   = kmod%version + 1
                    case(4); kmod%scope     = kmod%scope + 1
                    case(5); kmod%iteration = kmod%iteration + 1
                    case(6); kmod%particle  = kmod%particle + 1
                    case(7); kmod%family    = kmod%family + 1
                    case(8); kmod%reference = kmod%reference + 1
                    case(9); kmod%ordinal   = kmod%ordinal + 1
                end select
                if( rnd_keyed_bits(kmod) == word ) nsame = nsame + 1
            enddo
        enddo
        call assert_int(0, nsame, 'changing any field of the key changes the word')
        key%reference = 7
        do i = 1, NRUN
            key%particle = i
            words(i)  = rnd_keyed_bits(key)
            u(i)      = rnd_bits_to_unit(words(i))
            keys31(i) = int(shiftr(words(i), 33))
            inds(i)   = i
        enddo
        ! duplicates: sort by the top 31 bits, compare full words within runs of equal keys
        call hpsort(keys31, inds)
        ndup = 0
        do i = 2, NRUN
            if( keys31(i) /= keys31(i-1) ) cycle
            if( words(inds(i)) == words(inds(i-1)) ) ndup = ndup + 1
        enddo
        call assert_int(0, ndup, 'consecutive particle keys give distinct words')
        counts = 0
        do i = 1, NRUN
            counts(min(NBINS, 1 + int(u(i) * NBINS))) = counts(min(NBINS, 1 + int(u(i) * NBINS))) + 1
        enddo
        expected = real(NRUN, dp) / real(NBINS, dp)
        chi2     = sum((real(counts, dp) - expected)**2) / expected
        call assert_true(chi2 < 110._dp, 'the variates of consecutive keys are uniform over 64 bins')
        mean = sum(u) / real(NRUN, dp)
        num  = sum((u(1:NRUN-1) - mean) * (u(2:NRUN) - mean))
        den  = sum((u - mean)**2)
        call assert_true(abs(num / den) < 0.02_dp, 'consecutive keys'' variates have no serial correlation')
        call assert_true(abs(mean - 0.5_dp) < 0.005_dp, 'the variates of consecutive keys have mean one half')
    end subroutine test_keyed_sensitivity

    !> The keyed integer draw stays in range and follows the variate's equal-width bin, around 2^24 and at huge
    subroutine test_keyed_int()
        integer, parameter :: BOUNDS(5) = [2**23, 2**24-1, 2**24+1, 2**30, huge(0)]
        type(rnd_key) :: key
        integer  :: ibound, i, np, which
        real(dp) :: u
        logical  :: in_range, correct_bins
        write(*,'(A)') 'test_keyed_int'
        key%seed = 1_int64; key%family = RND_FAM_CAND
        call assert_int(1, rnd_keyed_int(key, 1), 'a singleton keyed integer range returns one')
        do ibound = 1, size(BOUNDS)
            np = BOUNDS(ibound)
            in_range = .true.; correct_bins = .true.
            do i = 1, 64
                key%ordinal = i
                which = rnd_keyed_int(key, np)
                u     = rnd_keyed(key)
                if( which < 1 .or. which > np ) in_range = .false.
                if( u < real(which-1,dp) / real(np,dp) .or. u >= real(which,dp) / real(np,dp) ) correct_bins = .false.
            enddo
            call assert_true(in_range,     'keyed integer draws stay within [1, NP] for NP='//int2str(np))
            call assert_true(correct_bins, 'keyed integer draws match equal-width probability bins for NP='//int2str(np))
        enddo
    end subroutine test_keyed_int

    !> The keyed partial permutation: a duplicate-free prefix of a permutation of the input; the same key
    !! gives the same permutation whatever was drawn before, another key another; every subset of size three of
    !! six is equally likely over many keys (20 subsets, chi-square of 19 degrees of freedom below 60) and so is
    !! every order of a full permutation of four (24 orders, chi-square of 23 degrees of freedom below 70; a
    !! uniform source exceeds either bound with probability below 1e-5)
    subroutine test_keyed_partial_shuffle()
        integer, parameter :: NSUB = 20000, NPERM = 24000
        type(rnd_key) :: key
        integer :: arr(N), copy(N), expected(N), six(6), four(4), i, k, code, counts(20), perms(24), nchanged
        real(dp) :: chi2
        write(*,'(A)') 'test_keyed_partial_shuffle'
        expected = [(i, i=1,N)]
        key%seed = 20261103_int64; key%scope = rnd_scope(RND_PROG_REFINE2D, RND_PASS_ALIGN); key%iteration = 4
        key%family = RND_FAM_CAND; key%particle = 11
        arr = expected
        call rnd_keyed_partial_shuffle(key, arr, 0)
        call assert_true(all(arr == expected), 'a zero-length keyed partial shuffle is a no-op')
        call rnd_keyed_partial_shuffle(key, arr, NPARTIAL)
        copy = arr
        call hpsort(copy)
        call assert_true(all(copy == expected), 'the keyed partial shuffle preserves the input permutation')
        copy(:NPARTIAL) = arr(:NPARTIAL)
        call hpsort(copy(:NPARTIAL))
        call assert_true(all(copy(2:NPARTIAL) /= copy(1:NPARTIAL-1)), 'the keyed partial shuffle prefix contains no duplicates')
        call assert_true(any(arr(:NPARTIAL) /= expected(:NPARTIAL)), 'the keyed partial shuffle changes the order')
        copy = expected
        call set_fixed_seed(5)
        do i = 1, 3
            k = int(10. * ran3())
        enddo
        call rnd_keyed_partial_shuffle(key, copy, NPARTIAL)
        call assert_true(all(copy == arr), 'the same key gives the same partial shuffle whatever was drawn before')
        nchanged = 0
        do i = 1, 50
            key%particle = 11 + i
            copy = expected
            call rnd_keyed_partial_shuffle(key, copy, NPARTIAL)
            if( any(copy(:NPARTIAL) /= arr(:NPARTIAL)) ) nchanged = nchanged + 1
        enddo
        call assert_int(50, nchanged, 'another key gives another partial shuffle')
        ! subsets of size three of six, coded by their sorted members
        counts = 0
        do i = 1, NSUB
            key%particle = i
            six = [1,2,3,4,5,6]
            call rnd_keyed_partial_shuffle(key, six, 3)
            call hpsort(six(1:3))
            code = subset_code(six(1), six(2), six(3))
            counts(code) = counts(code) + 1
        enddo
        chi2 = sum((real(counts, dp) - real(NSUB,dp)/20._dp)**2) / (real(NSUB,dp)/20._dp)
        call assert_true(all(counts > 0), 'every subset of size three of six is drawn')
        call assert_true(chi2 < 60._dp, 'the subsets of size three of six are equally likely over many keys')
        ! the orders of a full permutation of four
        perms = 0
        do i = 1, NPERM
            key%particle = i
            four = [1,2,3,4]
            call rnd_keyed_partial_shuffle(key, four, 4)
            code = perm_code(four)
            perms(code) = perms(code) + 1
        enddo
        chi2 = sum((real(perms, dp) - real(NPERM,dp)/24._dp)**2) / (real(NPERM,dp)/24._dp)
        call assert_true(all(perms > 0), 'every order of a permutation of four is drawn')
        call assert_true(chi2 < 70._dp, 'the orders of a full keyed permutation of four are equally likely')
        ! the prefix's ordinals: swapping the key's own ordinal changes nothing (the swap index overrides it)
        key%particle = 11; key%ordinal = 99
        copy = expected
        call rnd_keyed_partial_shuffle(key, copy, NPARTIAL)
        call assert_true(all(copy == arr), 'the key''s ordinal is overridden by the swap index')

    contains

        !> 1..20 for a sorted triple of 1..6
        integer function subset_code( a, b, c )
            integer, intent(in) :: a, b, c
            integer :: ia, ib, ic
            subset_code = 0
            do ia = 1, 6
                do ib = ia+1, 6
                    do ic = ib+1, 6
                        subset_code = subset_code + 1
                        if( ia == a .and. ib == b .and. ic == c ) return
                    enddo
                enddo
            enddo
            THROW_HARD('not a subset of size three of six; subset_code')
        end function subset_code

        !> 1..24 for a permutation of 1..4 (its Lehmer code)
        integer function perm_code( p )
            integer, intent(in) :: p(4)
            integer :: ia, ib, lesser
            perm_code = 0
            do ia = 1, 3
                lesser = 0
                do ib = ia+1, 4
                    if( p(ib) < p(ia) ) lesser = lesser + 1
                enddo
                perm_code = perm_code * (5 - ia) + lesser
            enddo
            perm_code = perm_code + 1
        end function perm_code

    end subroutine test_keyed_partial_shuffle

    !> The run seed and the scope: ensure_run_seed keeps a set SIMPLE_SEED and rnd_run_seed reads it back; the
    !! scope key packs the seed, the version, the program and pass kind and the iteration with the draw's
    !! fields zero; a draw key within it carries the particle, family, reference and ordinal and nothing else
    subroutine test_keyed_scope()
        type(rnd_key)  :: scope, key
        integer(int64) :: seed
        write(*,'(A)') 'test_keyed_scope'
        call seed_rnd_fixed(20261104, propagate=.true.)
        seed = ensure_run_seed()
        call assert_true(seed == 20261104_int64, 'ensure_run_seed keeps a set SIMPLE_SEED')
        call assert_true(rnd_run_seed() == 20261104_int64, 'rnd_run_seed reads SIMPLE_SEED')
        scope = rnd_scope_key(seed, RND_PROG_REFINE2D, RND_PASS_ALIGN, 7)
        call assert_true(scope%seed == seed .and. scope%version == RND_KEYED_VERSION .and. &
            &scope%scope == rnd_scope(RND_PROG_REFINE2D, RND_PASS_ALIGN) .and. scope%iteration == 7 .and. &
            &scope%particle == 0 .and. scope%family == 0 .and. scope%reference == 0 .and. scope%ordinal == 0, &
            &'the scope key carries the seed, version, program, pass and iteration only')
        key = rnd_draw_key(scope, 35645, RND_FAM_INPL, reference=2500, ordinal=3)
        call assert_true(key%seed == seed .and. key%version == scope%version .and. key%scope == scope%scope .and. &
            &key%iteration == 7 .and. key%particle == 35645 .and. key%family == RND_FAM_INPL .and. &
            &key%reference == 2500 .and. key%ordinal == 3, 'a draw key completes the scope with the draw''s fields')
        key = rnd_draw_key(scope, 1, RND_FAM_ASSIGN)
        call assert_true(key%reference == 0 .and. key%ordinal == 0, 'a draw key without reference and ordinal has them zero')
        call assert_true(rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_ASSIGN)) /= rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_INPL)), &
            &'the family separates the draws of one particle')
        call assert_true(rnd_keyed(rnd_draw_key(scope, 1, RND_FAM_ASSIGN)) /= &
            &rnd_keyed(rnd_draw_key(rnd_scope_key(seed, RND_PROG_REFINE3D, RND_PASS_ALIGN, 7), 1, RND_FAM_ASSIGN)), &
            &'the program separates the draws of one particle')
        call fix_test_seed()
    end subroutine test_keyed_scope

end module simple_rnd_tester
