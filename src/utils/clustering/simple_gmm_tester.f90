!@descr: unit tests for the tied-covariance Gaussian mixture (simple_gmm): recovery, order, floors, separation
! Two planar Gaussians with a shared diagonal covariance (1400 points at (0,0), 600 at (6,0), standard
! deviations 1 and 0.5) are fitted from initial means at (5,0) and (1,0): means, mixing proportions and covariance
! fall within a few sampling errors, the larger component is labelled 1 and the hard labels match the
! truth; the respawn rule leaves well-separated components alone, so the fit with it is the same; floored
! responsibilities are exact zeros or at least the floor and sum to one; a mixing-proportion floor of 0.2
! lifts a 10 % component; the Mahalanobis separation of the means is about 6 and the ICL is not below
! the BIC.
module simple_gmm_tester
use simple_core_module_api, only: dp
use simple_test_utils
use simple_rnd, only: gasdev
use simple_gmm, only: gmm
implicit none
private
public :: run_all_gmm_tests

contains

    subroutine run_all_gmm_tests()
        write(*,'(A)') '**** running all gmm tests ****'
        call test_gmm_recovery()
        call test_gmm_pi_floor()
    end subroutine run_all_gmm_tests

    !> n1 points around c1 and n2 around c2 with standard deviations sd(1:2) per axis, in that order
    subroutine two_blobs( X, n1, c1, n2, c2, sd )
        real(dp), allocatable, intent(out) :: X(:,:)
        integer,               intent(in)  :: n1, n2
        real(dp),              intent(in)  :: c1(2), c2(2), sd(2)
        integer :: i
        allocate(X(n1+n2,2))
        do i = 1, n1 + n2
            X(i,1) = real(gasdev(),dp)*sd(1)
            X(i,2) = real(gasdev(),dp)*sd(2)
            if( i <= n1 )then
                X(i,:) = X(i,:) + c1
            else
                X(i,:) = X(i,:) + c2
            endif
        end do
    end subroutine two_blobs

    subroutine test_gmm_recovery()
        real(dp), allocatable :: X(:,:), resp(:,:)
        real(dp) :: means(2,2), cov(2,2), pival(2), sep(2,2), rowsum
        integer  :: labels(2000), labels2(2000), i
        logical  :: ok, l_floor_ok
        type(gmm) :: gm
        write(*,'(A)') 'test_gmm_recovery'
        call set_fixed_seed(20261008)
        call two_blobs(X, 1400, [0._dp,0._dp], 600, [6._dp,0._dp], [1._dp,0.5_dp])
        means(:,1) = [5._dp, 0._dp]
        means(:,2) = [1._dp, 0._dp]
        call gm%new(X, 2, means, reg=1.e-6_dp, tol=1.e-8_dp, maxits=200, resp_floor=1.e-3_dp, respawn=.false.)
        call gm%fit(ok)
        call assert_true(ok, 'gmm fits two separated Gaussians')
        if( .not. ok )then
            call gm%kill
            deallocate(X)
            return
        endif
        call gm%get_means(means)
        call gm%get_cov(cov)
        call gm%get_pi(pival)
        call gm%get_labels(labels)
        call assert_true(maxval(abs(means(:,1) - [0._dp,0._dp])) < 0.15_dp, 'gmm component 1 is the larger Gaussian, mean recovered')
        call assert_true(maxval(abs(means(:,2) - [6._dp,0._dp])) < 0.2_dp,  'gmm component 2 mean recovered')
        call assert_true(abs(pival(1) - 0.7_dp) < 0.04_dp, 'gmm mixing proportions recovered')
        call assert_true(abs(cov(1,1) - 1._dp) < 0.15_dp .and. abs(cov(2,2) - 0.25_dp) < 0.05_dp .and. &
            &abs(cov(1,2)) < 0.1_dp, 'gmm tied covariance recovered')
        call assert_true(count(labels(1:1400) == 1) + count(labels(1401:2000) == 2) >= 1980, &
            &'gmm hard labels match the truth')
        allocate(resp(2000,2))
        call gm%get_resp(resp)
        l_floor_ok = .true.
        do i = 1, 2000
            rowsum = sum(resp(i,:))
            if( abs(rowsum - 1._dp) > 1.e-12_dp ) l_floor_ok = .false.
            if( any(resp(i,:) > 0._dp .and. resp(i,:) < 1.e-3_dp) ) l_floor_ok = .false.
        end do
        call assert_true(l_floor_ok, 'gmm floored responsibilities are zero or at least the floor and sum to one')
        call gm%get_pairsep(sep, ok)
        call assert_true(ok .and. abs(sep(1,2) - 6._dp) < 0.5_dp .and. sep(1,2) == sep(2,1), &
            &'gmm separation of the means is the Mahalanobis distance')
        call assert_true(gm%get_icl() >= gm%get_bic(), 'gmm ICL is not below the BIC')
        ! the respawn rule leaves separated components alone: the same fit
        means(:,1) = [5._dp, 0._dp]
        means(:,2) = [1._dp, 0._dp]
        call gm%new(X, 2, means, reg=1.e-6_dp, tol=1.e-8_dp, maxits=200, respawn=.true.)
        call gm%fit(ok)
        labels2 = 0
        if( ok ) call gm%get_labels(labels2)
        call assert_true(ok .and. all(labels2 == labels), 'gmm respawn leaves well-separated components alone')
        call gm%kill
        deallocate(X, resp)
    end subroutine test_gmm_recovery

    !> 1800 points at (0,0) and 200 at (8,0): the minor proportion is about 0.1, and a floor of 0.2 lifts it to
    !! about 0.2/1.1 after renormalization
    subroutine test_gmm_pi_floor()
        real(dp), allocatable :: X(:,:)
        real(dp) :: means(2,2), pival(2)
        logical  :: ok
        type(gmm) :: gm
        write(*,'(A)') 'test_gmm_pi_floor'
        call set_fixed_seed(20261009)
        call two_blobs(X, 1800, [0._dp,0._dp], 200, [8._dp,0._dp], [0.5_dp,0.5_dp])
        means(:,1) = [0._dp, 0._dp]
        means(:,2) = [8._dp, 0._dp]
        call gm%new(X, 2, means, reg=1.e-6_dp, tol=1.e-8_dp, maxits=200)
        call gm%fit(ok)
        call gm%get_pi(pival)
        call assert_true(ok .and. pival(2) < 0.12_dp, 'gmm without a floor keeps the minor proportion near 0.1')
        call gm%new(X, 2, means, reg=1.e-6_dp, tol=1.e-8_dp, maxits=200, pi_floor=0.2_dp)
        call gm%fit(ok)
        call gm%get_pi(pival)
        call assert_true(ok .and. pival(2) > 0.16_dp, 'gmm mixing-proportion floor lifts the minor component')
        call gm%kill
        deallocate(X)
    end subroutine test_gmm_pi_floor

end module simple_gmm_tester
