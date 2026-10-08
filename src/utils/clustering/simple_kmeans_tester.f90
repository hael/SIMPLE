!@descr: unit tests for k-means clustering (simple_kmeans): recovery, population order, weighted metric, edges
! Three well-separated planar groups of 30, 20 and 10 points are recovered exactly, labelled by decreasing
! population, with their means as centroids; four corner points split along the long side under unit
! weights and along the short side when that dimension is weighted 1000-fold (seeding and assignment
! worked by hand); a repeated call gives the same labels; k = 1 gives the mean and k = n the identity.
module simple_kmeans_tester
use simple_core_module_api, only: dp
use simple_test_utils
use simple_kmeans, only: kmeans
implicit none
private
public :: run_all_kmeans_tests

contains

    subroutine run_all_kmeans_tests()
        write(*,'(A)') '**** running all kmeans tests ****'
        call test_kmeans_separated_groups()
        call test_kmeans_weighted_metric()
        call test_kmeans_edges()
    end subroutine run_all_kmeans_tests

    !> npts points on a small deterministic ring around (cx,cy)
    subroutine ring( X, i0, npts, cx, cy )
        real(dp), intent(inout) :: X(:,:)
        integer,  intent(in)    :: i0, npts
        real(dp), intent(in)    :: cx, cy
        real(dp), parameter :: TWOPI = 6.283185307179586_dp
        integer :: i
        do i = 1, npts
            X(i0+i-1,1) = cx + 0.5_dp*cos(TWOPI*real(i,dp)/real(npts,dp))
            X(i0+i-1,2) = cy + 0.5_dp*sin(TWOPI*real(i,dp)/real(npts,dp))
        end do
    end subroutine ring

    !> groups of 10 at (0,0), 30 at (20,0) and 20 at (0,20), listed in that order: labels follow population
    !! (the 30 first), every group is one cluster and the centroids are the group means; a second call on
    !! the same construction gives the same labels
    subroutine test_kmeans_separated_groups()
        real(dp) :: X(60,2), cen(2,3), gmean(2)
        integer  :: labels(60), labels2(60)
        type(kmeans) :: km
        write(*,'(A)') 'test_kmeans_separated_groups'
        call ring(X,  1, 10,  0._dp,  0._dp)
        call ring(X, 11, 30, 20._dp,  0._dp)
        call ring(X, 41, 20,  0._dp, 20._dp)
        call km%new(X, 3)
        call km%cluster(labels, cen)
        call assert_true(all(labels(11:40) == 1), 'kmeans labels the 30-point group 1')
        call assert_true(all(labels(41:60) == 2), 'kmeans labels the 20-point group 2')
        call assert_true(all(labels(1:10)  == 3), 'kmeans labels the 10-point group 3')
        gmean = [sum(X(11:40,1)), sum(X(11:40,2))] / 30._dp
        call assert_true(maxval(abs(cen(:,1) - gmean)) < 1.e-12_dp, 'kmeans centroid 1 is the mean of its group')
        gmean = [sum(X(1:10,1)), sum(X(1:10,2))] / 10._dp
        call assert_true(maxval(abs(cen(:,3) - gmean)) < 1.e-12_dp, 'kmeans centroid 3 is the mean of its group')
        call assert_int(0, km%get_nchanged(), 'kmeans converges on separated groups')
        call km%cluster(labels2, cen)
        call assert_true(all(labels2 == labels), 'kmeans is deterministic')
        call km%kill
    end subroutine test_kmeans_separated_groups

    !> corners (0,0), (0,10), (1,0), (1,10), k = 2. Unit weights: every corner is equally near the mean, so
    !! the first seeds; the farthest from it is (1,10); the split is along y, {1,3} and {2,4}. Weights
    !! (1000,1): the farthest from (0,0) is again (1,10), but x now dominates and the split is {1,2}, {3,4}
    subroutine test_kmeans_weighted_metric()
        real(dp) :: X(4,2), cen(2,2)
        integer  :: labels(4)
        type(kmeans) :: km
        write(*,'(A)') 'test_kmeans_weighted_metric'
        X(:,1) = [0._dp, 0._dp, 1._dp,  1._dp]
        X(:,2) = [0._dp, 10._dp, 0._dp, 10._dp]
        call km%new(X, 2)
        call km%cluster(labels, cen)
        call assert_true(all(labels == [1,2,1,2]), 'kmeans with unit weights splits along the long side')
        call km%new(X, 2, wdim=[1000._dp, 1._dp])
        call km%cluster(labels, cen)
        call assert_true(all(labels == [1,1,2,2]), 'kmeans with a weighted dimension splits along it')
        call assert_true(maxval(abs(cen(:,1) - [0._dp, 5._dp])) < 1.e-12_dp, 'kmeans weighted centroids stay in the input frame')
        call km%kill
    end subroutine test_kmeans_weighted_metric

    !> k = 1: one cluster, centroid the mean; k = n on distinct points: every point its own cluster,
    !! labelled by index (all populations 1, ties by member index)
    subroutine test_kmeans_edges()
        real(dp) :: X(5,1), cen1(1,1), cen5(1,5)
        integer  :: labels(5)
        type(kmeans) :: km
        write(*,'(A)') 'test_kmeans_edges'
        X(:,1) = [0._dp, 1._dp, 3._dp, 7._dp, 8._dp]
        call km%new(X, 1)
        call km%cluster(labels, cen1)
        call assert_true(all(labels == 1), 'kmeans with k = 1 gives one cluster')
        call assert_true(abs(cen1(1,1) - 3.8_dp) < 1.e-12_dp, 'kmeans with k = 1 has the mean as centroid')
        call km%new(X, 5)
        call km%cluster(labels, cen5)
        call assert_true(all(labels == [1,2,3,4,5]), 'kmeans with k = n is the identity')
        call assert_true(maxval(abs(cen5(1,:) - X(:,1))) < 1.e-12_dp, 'kmeans with k = n has the points as centroids')
        call km%kill
    end subroutine test_kmeans_edges

end module simple_kmeans_tester
