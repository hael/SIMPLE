!@descr: unit tests for extreme deconvolution (simple_xd_gmm): noise-free mixture recovery, K selection, posterior
! 3000 planar points v from 0.6 N((0,0), 0.25 I) + 0.4 N((5,0), 0.25 I), each observed as x = r v + e with
! r in {1, 0.8} and e ~ N(0, s^2 I), s in {0.3, 0.6}, all known. Fitted with two components the means,
! proportions and the noise-free covariances (0.25 I, not the inflated observed spread) fall within a few
! sampling errors, the heavier component first. Held-out K selection prefers two components clearly
! over one. The posterior responsibilities sum to one and the posterior covariances are smaller on
! average than the noise of the underlying points.
module simple_xd_gmm_tester
use simple_core_module_api, only: dp
use simple_test_utils
use simple_rnd,    only: gasdev
use simple_xd_gmm, only: xd_gmm, xd_select_k
implicit none
private
public :: run_all_xd_gmm_tests

integer, parameter :: NPTS = 3000

contains

    subroutine run_all_xd_gmm_tests()
        write(*,'(A)') '**** running all xd_gmm tests ****'
        call test_xd_gmm_recovery()
        call test_xd_gmm_select_k()
    end subroutine run_all_xd_gmm_tests

    !> the fixture: X(NPTS,2), R(2,2,NPTS), Nz(2,2,NPTS), the first 1800 points from component 1
    subroutine fixture( X, R, Nz, vnoise )
        real(dp), intent(out) :: X(NPTS,2), R(2,2,NPTS), Nz(2,2,NPTS)
        real(dp), intent(out) :: vnoise   !< mean trace of the noise covariance of v, (sd/rs)^2 per axis
        real(dp) :: v(2), rs, sd
        integer  :: i
        call set_fixed_seed(20261010)
        vnoise = 0._dp
        do i = 1, NPTS
            v = 0.5_dp*[real(gasdev(),dp), real(gasdev(),dp)]
            if( i > 1800 ) v(1) = v(1) + 5._dp
            rs = merge(1._dp, 0.8_dp, mod(i,2) == 0)
            sd = merge(0.3_dp, 0.6_dp, mod(i,3) == 0)
            R(:,:,i)  = 0._dp
            R(1,1,i)  = rs
            R(2,2,i)  = rs
            Nz(:,:,i) = 0._dp
            Nz(1,1,i) = sd*sd
            Nz(2,2,i) = sd*sd
            X(i,:)    = rs*v + sd*[real(gasdev(),dp), real(gasdev(),dp)]
            vnoise    = vnoise + 2._dp*(sd/rs)**2
        end do
        vnoise = vnoise/real(NPTS,dp)
    end subroutine fixture

    subroutine test_xd_gmm_recovery()
        real(dp), allocatable :: X(:,:), R(:,:,:), Nz(:,:,:), xhat(:,:), xcov(:,:,:), resp(:,:)
        real(dp) :: mu(2,2), Sig(2,2,2), pik(2), vnoise, trpost
        integer  :: i
        type(xd_gmm) :: xd
        write(*,'(A)') 'test_xd_gmm_recovery'
        allocate(X(NPTS,2), R(2,2,NPTS), Nz(2,2,NPTS))
        call fixture(X, R, Nz, vnoise)
        call xd%new(2, 2, 150, 1.e-6_dp, 1.e-4_dp)
        call xd%init(X, Nz)
        call xd%fit(X, R, Nz)
        call xd%get_means(mu)
        call xd%get_covs(Sig)
        call xd%get_pi(pik)
        call assert_true(abs(pik(1) - 0.6_dp) < 0.04_dp, 'xd_gmm recovers the proportions, the heavier component first')
        call assert_true(maxval(abs(mu(:,1) - [0._dp,0._dp])) < 0.15_dp, 'xd_gmm recovers the mean of component 1')
        call assert_true(maxval(abs(mu(:,2) - [5._dp,0._dp])) < 0.15_dp, 'xd_gmm recovers the mean of component 2')
        call assert_true(abs(Sig(1,1,1) - 0.25_dp) < 0.08_dp .and. abs(Sig(2,2,1) - 0.25_dp) < 0.08_dp .and. &
            &abs(Sig(1,1,2) - 0.25_dp) < 0.08_dp .and. abs(Sig(2,2,2) - 0.25_dp) < 0.08_dp, &
            &'xd_gmm recovers the noise-free covariances')
        call assert_true(abs(Sig(1,2,1)) < 0.06_dp .and. abs(Sig(1,2,2)) < 0.06_dp, 'xd_gmm covariances stay diagonal')
        allocate(xhat(NPTS,2), xcov(2,2,NPTS), resp(NPTS,2))
        call xd%posterior(X, R, Nz, xhat, xcov, resp)
        call assert_true(maxval(abs(sum(resp, dim=2) - 1._dp)) < 1.e-10_dp, 'xd_gmm responsibilities sum to one')
        call assert_true(count(maxloc(resp(1:1800,:), dim=2) == 1) + count(maxloc(resp(1801:NPTS,:), dim=2) == 2) &
            &>= NPTS - 30, 'xd_gmm responsibilities assign the points to their components')
        trpost = 0._dp
        do i = 1, NPTS
            trpost = trpost + xcov(1,1,i) + xcov(2,2,i)
        end do
        trpost = trpost/real(NPTS,dp)
        call assert_true(trpost > 0._dp .and. trpost < vnoise, 'xd_gmm posterior covariances shrink the noise')
        call xd%kill
        deallocate(X, R, Nz, xhat, xcov, resp)
    end subroutine test_xd_gmm_recovery

    subroutine test_xd_gmm_select_k()
        real(dp), allocatable :: X(:,:), R(:,:,:), Nz(:,:,:)
        real(dp) :: scores(4), vnoise
        integer  :: kbest
        write(*,'(A)') 'test_xd_gmm_select_k'
        allocate(X(NPTS,2), R(2,2,NPTS), Nz(2,2,NPTS))
        call fixture(X, R, Nz, vnoise)
        call xd_select_k(X, R, Nz, 4, NPTS, 150, 1.e-6_dp, 1.e-4_dp, kbest, scores)
        call assert_true(kbest >= 2, 'xd_select_k chooses more than one component for a two-component mixture')
        call assert_true(scores(2) > scores(1) + 0.1_dp, 'xd_select_k held-out likelihood rises clearly from K = 1 to 2')
        deallocate(X, R, Nz)
    end subroutine test_xd_gmm_select_k

end module simple_xd_gmm_tester
