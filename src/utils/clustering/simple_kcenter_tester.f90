!@descr: unit tests for farthest-point k-center clustering (simple_kcenter): coverage, centers, labels, edges
! A dense ring of 50 points at the origin, a far singleton at (100,0) and three points near (0,100): the
! first center is the singleton (farthest from the mean), the second the farthest of the three, the third
! a ring point, so every region gets a center; labels follow population (ring, triple, singleton). A
! repeated call gives the same result; k = 1 has the farthest point from the mean as its center.
module simple_kcenter_tester
use simple_core_module_api, only: dp
use simple_test_utils
use simple_kcenter, only: kcenter
implicit none
private
public :: run_all_kcenter_tests

contains

    subroutine run_all_kcenter_tests()
        write(*,'(A)') '**** running all kcenter tests ****'
        call test_kcenter_coverage()
        call test_kcenter_edges()
    end subroutine run_all_kcenter_tests

    !> 50 ring points at the origin (radius 0.5), the singleton (100,0) as point 51 and the triple
    !! (0,100), (1,100), (0,101) as points 52-54
    subroutine fixture( X )
        real(dp), intent(out) :: X(54,2)
        real(dp), parameter :: TWOPI = 6.283185307179586_dp
        integer :: i
        do i = 1, 50
            X(i,1) = 0.5_dp*cos(TWOPI*real(i,dp)/50._dp)
            X(i,2) = 0.5_dp*sin(TWOPI*real(i,dp)/50._dp)
        end do
        X(51,:) = [100._dp,   0._dp]
        X(52,:) = [  0._dp, 100._dp]
        X(53,:) = [  1._dp, 100._dp]
        X(54,:) = [  0._dp, 101._dp]
    end subroutine fixture

    subroutine test_kcenter_coverage()
        real(dp) :: X(54,2)
        integer  :: centers(3), labels(54), centers2(3), labels2(54)
        type(kcenter) :: kc
        write(*,'(A)') 'test_kcenter_coverage'
        call fixture(X)
        call kc%new(X, 3)
        call kc%cluster(centers, labels)
        call assert_true(all(labels(1:50) == 1), 'kcenter gives the dense ring one cluster, labelled 1')
        call assert_true(all(labels(52:54) == 2), 'kcenter gives the triple one cluster, labelled 2')
        call assert_int(3, labels(51), 'kcenter gives the far singleton its own cluster')
        call assert_int(51, centers(3), 'kcenter takes the point farthest from the mean as a center')
        call assert_int(54, centers(2), 'kcenter takes the farthest point of the triple as its center')
        call assert_true(centers(1) >= 1 .and. centers(1) <= 50, 'kcenter places the third center on the ring')
        call kc%cluster(centers2, labels2)
        call assert_true(all(labels2 == labels) .and. all(centers2 == centers), 'kcenter is deterministic')
        call kc%kill
    end subroutine test_kcenter_coverage

    subroutine test_kcenter_edges()
        real(dp) :: X(54,2)
        integer  :: centers(1), labels(54)
        type(kcenter) :: kc
        write(*,'(A)') 'test_kcenter_edges'
        call fixture(X)
        call kc%new(X, 1)
        call kc%cluster(centers, labels)
        call assert_true(all(labels == 1), 'kcenter with k = 1 gives one cluster')
        call assert_int(51, centers(1), 'kcenter with k = 1 centers on the point farthest from the mean')
        call kc%kill
    end subroutine test_kcenter_edges

end module simple_kcenter_tester
