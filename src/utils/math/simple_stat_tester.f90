!@descr: unit test routines for the statistics utilities (simple_stat)
module simple_stat_tester
use simple_test_utils ! assertions etc.
use simple_defs       ! TINY
use simple_type_defs  ! weighting criteria enumerators
use simple_stat,      only: corrs2weights, conv2rank_weights, rank_sum_weights, rank_centroid_weights,&
                           &rank_exponent_weights, rank_inverse_weights, median, median_nocopy, calc_stats
implicit none
private
public :: run_all_stat_tests

real, parameter :: EPS = 1.0e-5

contains

    subroutine run_all_stat_tests()
        write(*,'(A)') '**** running all stat tests ****'
        call test_rank_weight_kernels()
        call test_conv2rank_weights()
        call test_corrs2weights_corrw()
        call test_corrs2weights_other_criteria()
        call test_median()
        call test_calc_stats()
    end subroutine run_all_stat_tests

    !---------------- median ----------------

    ! two arrays on which the selection routine's two-element final partition matters (they
    ! gave 14.5 and 23 with the `ir-1` typo in selec, 2026-09-23); the true medians by sorting
    subroutine test_median()
        real, parameter :: EVEN(10) = [4., 15., 3., 36., 9., 19., 27., 10., 35., 8.]
        real, parameter :: ODD(11)  = [25., 23., 20., 28., 6., 4., 31., 13., 24., 35., 29.]
        real :: work(11)
        write(*,'(A)') 'test_median'
        call assert_real(12.5, median(EVEN), 0., 'median of an even count is the mean of the two middle values')
        call assert_real(24.0, median(ODD),  0., 'median of an odd count is the middle value')
        call assert_real(2.0,  median([2.]), 0., 'median of one value')
        call assert_real(2.5,  median([3., 2.]), 0., 'median of two values')
        call assert_real(7.0,  median([7., 7., 7., 7.]), 0., 'median of equal values')
        work(1:10) = EVEN
        call assert_real(12.5, median_nocopy(work(1:10)), 0., 'median_nocopy gives the same value')
        work = ODD
        call assert_real(24.0, median_nocopy(work), 0., 'median_nocopy on the odd array')
    end subroutine test_median

    ! calc_stats on the even array above and on its masked half
    subroutine test_calc_stats()
        real, parameter :: EVEN(10) = [4., 15., 3., 36., 9., 19., 27., 10., 35., 8.]
        type(stats_struct) :: st
        logical :: mask(10)
        write(*,'(A)') 'test_calc_stats'
        call calc_stats(EVEN, st)
        call assert_real(16.6,     st%avg,  1.e-4, 'calc_stats: mean')
        call assert_real(12.5,     st%med,  0.,    'calc_stats: median')
        call assert_real(12.24926, st%sdev, 1.e-3, 'calc_stats: standard deviation (n-1 normalisation)')
        call assert_real(3.,       st%minv, 0.,    'calc_stats: minimum')
        call assert_real(36.,      st%maxv, 0.,    'calc_stats: maximum')
        mask = [.true., .false., .true., .false., .true., .false., .true., .false., .true., .false.]
        call calc_stats(EVEN, st, mask)
        call assert_real(15.6,     st%avg,  1.e-4, 'calc_stats with mask: mean of the odd entries')
        call assert_real(9.,       st%med,  0.,    'calc_stats with mask: median of the odd entries')
        call assert_real(3.,       st%minv, 0.,    'calc_stats with mask: minimum')
        call assert_real(35.,      st%maxv, 0.,    'calc_stats with mask: maximum')
    end subroutine test_calc_stats

    !---------------- the four rank-weight kernels ----------------

    ! each kernel maps rank 1..n (1 = best) to weights that are positive,
    ! non-increasing in rank and sum to one
    subroutine test_rank_weight_kernels()
        integer, parameter :: N = 200
        real :: w(N)
        write(*,'(A)') 'test_rank_weight_kernels'
        call rank_sum_weights(N, w)
        call check_kernel(w, 'rank_sum_weights')
        call assert_real(2.0/real(N+1), w(1), EPS, 'rank_sum_weights first weight 2/(n+1)')
        call rank_centroid_weights(N, w)
        call check_kernel(w, 'rank_centroid_weights')
        call assert_real(1.0/real(N)**2, w(N), 1.0e-7, 'rank_centroid_weights last weight 1/n^2')
        call rank_exponent_weights(N, 10.0, w)
        call check_kernel(w, 'rank_exponent_weights')
        call rank_exponent_weights(N, 1.0, w)
        call assert_real(2.0/real(N+1), w(1), EPS, 'rank_exponent_weights p=1 equals rank_sum_weights')
        call rank_inverse_weights(N, w)
        call check_kernel(w, 'rank_inverse_weights')
        call assert_real(0.5, w(2)/w(1), EPS, 'rank_inverse_weights second/first = 1/2')

        contains

            subroutine check_kernel( w, name )
                real,             intent(in) :: w(:)
                character(len=*), intent(in) :: name
                integer :: i
                logical :: nonincreasing
                nonincreasing = .true.
                do i = 2,size(w)
                    if( w(i) > w(i-1) + EPS ) nonincreasing = .false.
                end do
                call assert_true(all(w > 0.0),           name//': all weights positive')
                call assert_true(nonincreasing,          name//': weights non-increasing in rank')
                call assert_real(1.0, sum(w), 1.0e-4,    name//': weights sum to one')
            end subroutine check_kernel

    end subroutine test_rank_weight_kernels

    !---------------- conv2rank_weights ----------------

    subroutine test_conv2rank_weights()
        real    :: w(6), w_in(6)
        integer(kind=kind(ENUM_WCRIT)) :: crit_list(4)
        integer :: ic
        character(len=16) :: crit_name(4)
        write(*,'(A)') 'test_conv2rank_weights'
        ! one non-zero entry takes all the weight
        w = [0.0, 0.0, 0.3, 0.0, 0.0, 0.0]
        call conv2rank_weights(6, w, RANK_SUM_CRIT)
        call assert_real(1.0, w(3), EPS,        'conv2rank_weights: single non-zero gets weight 1')
        call assert_real(0.0, sum(w) - w(3), EPS, 'conv2rank_weights: single non-zero, others 0')
        ! nothing to rank
        w = 0.0
        call conv2rank_weights(6, w, RANK_SUM_CRIT)
        call assert_true(all(w == 0.0),         'conv2rank_weights: all-zero input stays zero')
        ! general case: zeros stay zero, the rest follow the rank of the input values
        w_in = [0.2, 0.0, 0.9, 0.5, 0.0, 0.7]
        crit_list = [RANK_SUM_CRIT, RANK_CEN_CRIT, RANK_EXP_CRIT, RANK_INV_CRIT]
        crit_name = ['RANK_SUM_CRIT   ', 'RANK_CEN_CRIT   ', 'RANK_EXP_CRIT   ', 'RANK_INV_CRIT   ']
        do ic = 1,4
            w = w_in
            call conv2rank_weights(6, w, crit_list(ic), p=2.0)
            call assert_real(0.0, w(2), EPS,    trim(crit_name(ic))//': zero input stays zero (2)')
            call assert_real(0.0, w(5), EPS,    trim(crit_name(ic))//': zero input stays zero (5)')
            call assert_real(1.0, sum(w), 1.0e-4, trim(crit_name(ic))//': non-zero weights sum to one')
            call assert_true(w(3) > w(6) .and. w(6) > w(4) .and. w(4) > w(1),&
                &trim(crit_name(ic))//': weights follow the input ranking')
        end do
        ! the largest weight of the four-entry ranking is the kernel's first weight
        w = w_in
        call conv2rank_weights(6, w, RANK_INV_CRIT)
        call assert_real(1.0/(1.0 + 0.5 + 1.0/3.0 + 0.25), w(3), EPS,&
            &'RANK_INV_CRIT: best entry gets the first inverse-rank weight')
    end subroutine test_conv2rank_weights

    !---------------- corrs2weights, CORRW_CRIT ----------------

    subroutine test_corrs2weights_corrw()
        real, parameter   :: corrs(12) = [-1.0, 0.0, 0.005, 0.1, 0.2, 0.3, 0.4, 0.5, 0.51, 0.52, 0.53, 0.6]
        real, allocatable :: w(:)
        integer :: i
        logical :: nondecreasing
        write(*,'(A)') 'test_corrs2weights_corrw'
        w = corrs2weights(corrs, CORRW_CRIT)
        call assert_int(12, size(w),               'corrs2weights returns one weight per correlation')
        call assert_real(1.0, sum(w), 1.0e-4,      'CORRW_CRIT: weights sum to one')
        call assert_real(0.0, w(1), EPS,           'CORRW_CRIT: negative correlation gets weight 0')
        call assert_real(0.0, w(2), EPS,           'CORRW_CRIT: zero correlation gets weight 0')
        call assert_true(all(w(3:) > 0.0),         'CORRW_CRIT: positive correlations get positive weight')
        nondecreasing = .true.
        do i = 4,12
            if( w(i) < w(i-1) - EPS ) nondecreasing = .false.
        end do
        call assert_true(nondecreasing,            'CORRW_CRIT: weights are monotone in the correlation')
        call assert_true(w(12) > 2.0 * w(3),       'CORRW_CRIT: the sigmoid spread separates best from worst')
        ! without the sigmoid normalisation the weights are exp(corr) up to a constant
        w = corrs2weights([0.2, 0.7], CORRW_CRIT, norm_sigm=.false.)
        call assert_real(exp(0.5), w(2)/w(1), 1.0e-4, 'CORRW_CRIT, norm_sigm=.false.: weight ratio is exp(dcorr)')
        ! a tight spread (max/min below the threshold) skips the sigmoid: same ratio rule
        w = corrs2weights([0.50, 0.55, 0.60], CORRW_CRIT)
        call assert_real(exp(0.1), w(3)/w(1), 1.0e-4, 'CORRW_CRIT: no sigmoid below the max/min threshold')
        ! nothing positive: all weights zero
        w = corrs2weights([-0.5, -0.1, 0.0], CORRW_CRIT)
        call assert_true(all(w == 0.0),            'CORRW_CRIT: no positive correlation gives all-zero weights')
    end subroutine test_corrs2weights_corrw

    !---------------- corrs2weights, the other criteria ----------------

    subroutine test_corrs2weights_other_criteria()
        real, parameter   :: corrs(5) = [0.2, -0.3, 0.9, 0.5, 0.0]
        real, allocatable :: w(:)
        write(*,'(A)') 'test_corrs2weights_other_criteria'
        w = corrs2weights(corrs, UNIFORM_CRIT)
        call assert_true(all(abs(w - 0.2) < EPS),  'UNIFORM_CRIT: every weight is 1/n')
        w = corrs2weights(corrs, RANK_SUM_CRIT)
        call assert_real(1.0, sum(w), 1.0e-4,      'RANK_SUM_CRIT: weights sum to one')
        call assert_real(0.0, w(2), EPS,           'RANK_SUM_CRIT: negative correlation gets weight 0')
        call assert_real(0.0, w(5), EPS,           'RANK_SUM_CRIT: zero correlation gets weight 0')
        call assert_true(w(3) > w(4) .and. w(4) > w(1), 'RANK_SUM_CRIT: weights follow the correlation ranking')
        call assert_real(0.5, w(3), EPS,           'RANK_SUM_CRIT: best of three gets 2/(3+1)')
        w = corrs2weights(corrs, RANK_EXP_CRIT, p=2.0)
        call assert_real(1.0, sum(w), 1.0e-4,      'RANK_EXP_CRIT: weights sum to one')
        call assert_real(9.0/14.0, w(3), EPS,      'RANK_EXP_CRIT p=2: best of three gets 9/14')
        w = corrs2weights(corrs, CORRW_ZSCORE_CRIT)
        call assert_true(all(w >= 0.0),            'CORRW_ZSCORE_CRIT: weights are non-negative')
        call assert_real(0.0, w(2), EPS,           'CORRW_ZSCORE_CRIT: negative correlation gets weight 0')
        call assert_real(0.0, w(5), EPS,           'CORRW_ZSCORE_CRIT: zero correlation gets weight 0')
        call assert_true(w(3) > w(4) .and. w(4) > w(1), 'CORRW_ZSCORE_CRIT: weights follow the correlation ranking')
    end subroutine test_corrs2weights_other_criteria

end module simple_stat_tester
