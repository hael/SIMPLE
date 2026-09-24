!@descr: unit tests for online moments estimation (simple_online_var): closed-form means and sample variances
! Replaces the in-module test_online_var, which compared with simple_stat%moment and printed its
! verdict without recording it, so it could not fail.
module simple_online_var_tester
use simple_test_utils
use simple_defs
use simple_online_var, only: online_var
use simple_stat,       only: moment
use simple_rnd,        only: gasdev
implicit none
private
public :: run_all_online_var_tests

contains

    subroutine run_all_online_var_tests()
        write(*,'(A)') '**** running all online_var tests ****'
        call test_closed_forms()
        call test_few_samples()
        call test_large_offset()
        call test_against_moment()
    end subroutine run_all_online_var_tests

    ! 2,4,4,4,5,5,7,9: mean 5, squared deviations 32, sample variance 32/7
    subroutine test_closed_forms()
        real, parameter  :: X(8) = [2., 4., 4., 4., 5., 5., 7., 9.]
        type(online_var) :: ov
        integer :: i
        write(*,'(A)') 'test_closed_forms'
        do i = 1,size(X)
            call ov%add(X(i))
        end do
        call assert_real(5.,      ov%get_mean(), 1.e-6, 'mean of 2,4,4,4,5,5,7,9 is 5')
        call assert_real(32./7.,  ov%get_var(),  1.e-6, 'sample variance of 2,4,4,4,5,5,7,9 is 32/7 (n-1)')
    end subroutine test_closed_forms

    subroutine test_few_samples()
        type(online_var) :: ov, empty
        write(*,'(A)') 'test_few_samples'
        call assert_real(0., empty%get_mean(), 0., 'no observations: mean 0')
        call assert_real(0., empty%get_var(),  0., 'no observations: variance 0')
        call ov%add(3.5)
        call assert_real(3.5, ov%get_mean(), 0., 'one observation: the mean is the observation')
        call assert_real(0.,  ov%get_var(),  0., 'one observation: variance 0')
    end subroutine test_few_samples

    ! 10001..10004 are exact in single precision; mean 10002.5, sample variance 5/3. A single-precision
    ! sum of squares loses this completely (10002.5**2 needs 27 bits)
    subroutine test_large_offset()
        type(online_var) :: ov
        integer :: i
        write(*,'(A)') 'test_large_offset'
        do i = 1,4
            call ov%add(10000. + real(i))
        end do
        call assert_real(10002.5, ov%get_mean(), 0.,   'mean of 10001..10004')
        call assert_real(5./3.,   ov%get_var(),  1.e-5, 'sample variance of 10001..10004 is 5/3 despite the offset')
    end subroutine test_large_offset

    ! 10000 seeded Gaussian draws: the online estimate agrees with the two-pass moment()
    subroutine test_against_moment()
        integer, parameter :: N = 10000
        type(online_var) :: ov
        real    :: samples(N), ave, sdev, var
        integer :: i
        logical :: err
        write(*,'(A)') 'test_against_moment'
        do i = 1,N
            samples(i) = gasdev(5., 2.)
            call ov%add(samples(i))
        end do
        call moment(samples, ave, sdev, var, err)
        call assert_false(err, 'moment() succeeds on the draws')
        call assert_real(ave, ov%get_mean(), 1.e-4, 'online mean equals the two-pass mean')
        call assert_real(var, ov%get_var(),  1.e-3 * var, 'online sample variance equals the two-pass variance')
    end subroutine test_against_moment

end module simple_online_var_tester
