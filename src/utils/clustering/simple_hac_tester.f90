!@descr: unit tests for hierarchical agglomerative clustering (simple_hac): average and complete linkage, stop rules, mask
! Points on a line, distance |x_i - x_j|; every expectation follows from the merge order by hand.
! Average linkage: separated groups are recovered with labels by decreasing population and their medoids;
! a dense group of twelve items stays one cluster beside four far singletons; a five-point case separates
! the population-weighted update (UPGMA) from the unweighted one (WPGMA). Complete linkage: a four-point
! case where the maximum cross distance sends the second merge elsewhere than the average, with the merge
! history. Threshold stop under both linkages, the mask, and nclust = n and nclust = 1 as the edges, run
! on one object that is reconstructed in between.
module simple_hac_tester
use simple_core_module_api, only: dp
use simple_test_utils
use simple_hac, only: hac
implicit none
private
public :: run_all_hac_tests

contains

    subroutine run_all_hac_tests()
        write(*,'(A)') '**** running all hac tests ****'
        call test_hac_separated_groups()
        call test_hac_dense_group_stays_whole()
        call test_hac_population_weighted_update()
        call test_hac_complete_linkage()
        call test_hac_threshold_stop()
        call test_hac_mask()
        call test_hac_edges()
    end subroutine run_all_hac_tests

    !> clusters the points x on a line into nclust groups
    subroutine run_line( x, linkage, nclust, i_medoids, labels )
        real,                 intent(in)    :: x(:)
        character(len=*),     intent(in)    :: linkage
        integer,              intent(in)    :: nclust
        integer, allocatable, intent(inout) :: i_medoids(:), labels(:)
        type(hac) :: hc
        call hc%new(size(x), line_dmat(x), linkage, nclust=nclust)
        call hc%cluster(i_medoids, labels)
        call hc%kill
    end subroutine run_line

    !> |x_i - x_j| for points on a line
    function line_dmat( x ) result( dmat )
        real, intent(in)      :: x(:)
        real(dp), allocatable :: dmat(:,:)
        integer :: i, j
        allocate(dmat(size(x),size(x)))
        do j = 1, size(x)
            do i = 1, size(x)
                dmat(i,j) = real(abs(x(i) - x(j)),dp)
            enddo
        enddo
    end function line_dmat

    !> groups {0,1,2}, {50,51}, {100,101,102,103}: six within-group merges (distances <= 3) precede
    !! any merge across groups (>= 48); labels by decreasing population, the medoid is the member with
    !! the smallest summed distance, the first one on a tie (101 before 102 in the largest group)
    subroutine test_hac_separated_groups()
        real, parameter      :: X(9) = [0., 1., 2., 50., 51., 100., 101., 102., 103.]
        integer, allocatable :: i_medoids(:), labels(:)
        write(*,'(A)') 'test_hac_separated_groups'
        call run_line(X, 'average', 3, i_medoids, labels)
        call assert_int(3, size(i_medoids), 'hac forms the requested number of clusters')
        call assert_true(all(labels == [2,2,2,3,3,1,1,1,1]), 'hac recovers the groups, labelled by decreasing population')
        call assert_true(all(i_medoids == [7,2,4]), 'hac medoids are the minimal-sum members, first on a tie')
    end subroutine test_hac_separated_groups

    !> twelve items at 0..11 and four at 1000..4000: the eleven merges that take sixteen items to five
    !! clusters all fall inside the dense group (mean distances <= 11 against >= 989), so it is one
    !! cluster; its medoid is the first median item (x = 5)
    subroutine test_hac_dense_group_stays_whole()
        real    :: x(16)
        integer, allocatable :: i_medoids(:), labels(:)
        integer :: i
        write(*,'(A)') 'test_hac_dense_group_stays_whole'
        x(1:12)  = [(real(i), i=0,11)]
        x(13:16) = [1000., 2000., 3000., 4000.]
        call run_line(x, 'average', 5, i_medoids, labels)
        call assert_true(all(labels(1:12) == 1), 'hac keeps the dense group in one cluster')
        call assert_true(all(labels(13:16) == [2,3,4,5]), 'hac leaves the far items as singletons, by index')
        call assert_int(6, i_medoids(1), 'hac medoid of the dense group is its first median item')
    end subroutine test_hac_dense_group_stays_whole

    !> x = 0, 0, 10, 4.6, 18. The merges (1,2) at 0 and ({1,2},4) at 4.6 give T = {1,2,4}. Population
    !! weighting gives d(T,3) = (2*10 + 5.4)/3 = 8.47 > d(3,5) = 8, so 3 joins 5; the unweighted mean
    !! (10 + 5.4)/2 = 7.7 would have joined 3 to T instead
    subroutine test_hac_population_weighted_update()
        real, parameter      :: X(5) = [0., 0., 10., 4.6, 18.]
        integer, allocatable :: i_medoids(:), labels(:)
        write(*,'(A)') 'test_hac_population_weighted_update'
        call run_line(X, 'average', 2, i_medoids, labels)
        call assert_true(all(labels == [1,1,2,1,2]), 'hac weights the average-linkage update by cluster population')
    end subroutine test_hac_population_weighted_update

    !> x = 0, 1, 3, 5.8. Both linkages first join (1,2) at 1. Then average gives d({1,2},3) = 2.5 < d(3,4)
    !! = 2.8, so 3 joins {1,2}; complete gives d({1,2},3) = max(3,2) = 3 > 2.8, so 3 joins 4. The complete
    !! history is (1,2) at 1 and (3,4) at 2.8; the tie of populations 2 and 2 labels by representative
    subroutine test_hac_complete_linkage()
        real, parameter       :: X(4) = [0., 1., 3., 5.8]
        integer,  allocatable :: i_medoids(:), labels(:), pairs(:,:)
        real(dp), allocatable :: heights(:)
        type(hac) :: hc
        write(*,'(A)') 'test_hac_complete_linkage'
        call run_line(X, 'average', 2, i_medoids, labels)
        call assert_true(all(labels == [1,1,1,2]), 'hac average linkage joins the middle point to the pair')
        call assert_true(all(i_medoids == [2,4]), 'hac average-linkage medoids')
        call hc%new(size(X), line_dmat(X), 'complete', nclust=2)
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == [1,1,2,2]), 'hac complete linkage joins the middle point to the far point')
        call hc%get_history(pairs, heights)
        call assert_int(2, size(heights), 'hac history holds one entry per merge')
        if( size(heights) == 2 )then
            call assert_true(all(pairs(:,1) == [1,2]) .and. all(pairs(:,2) == [3,4]), 'hac history names the joined clusters')
            call assert_true(abs(heights(1) - 1._dp) < 1.e-5_dp .and. abs(heights(2) - 2.8_dp) < 1.e-5_dp, &
                &'hac history holds the complete-linkage distances')
        endif
        call hc%kill
    end subroutine test_hac_complete_linkage

    !> x = 0, 1, 3, 7, 8 with threshold 2.7. Both linkages join (1,2) and (4,5) at 1; average then joins 3 to
    !! {1,2} at 2.5 and stops at 6.17 (two clusters), complete would join it at 3 > 2.7 and stops (three)
    subroutine test_hac_threshold_stop()
        real, parameter      :: X(5) = [0., 1., 3., 7., 8.]
        integer, allocatable :: i_medoids(:), labels(:)
        type(hac) :: hc
        write(*,'(A)') 'test_hac_threshold_stop'
        call hc%new(size(X), line_dmat(X), 'average', thres=2.7_dp)
        call hc%cluster(i_medoids, labels)
        call assert_int(2, hc%get_nclust(), 'hac average linkage at threshold 2.7 forms two clusters')
        call assert_true(all(labels == [1,1,1,2,2]), 'hac average-linkage threshold labels')
        call hc%new(size(X), line_dmat(X), 'complete', thres=2.7_dp)
        call hc%cluster(i_medoids, labels)
        call assert_int(3, hc%get_nclust(), 'hac complete linkage at threshold 2.7 forms three clusters')
        call assert_true(all(labels == [1,1,3,2,2]), 'hac complete-linkage threshold labels, population ties by representative')
        call assert_int(3, size(i_medoids), 'hac returns a medoid per cluster under the threshold stop')
        call hc%kill
    end subroutine test_hac_threshold_stop

    !> x = 0, 0.5, 1, 10 at threshold 2. Unmasked: (1,2) at 0.5 and 3 at 0.75 join; masking out entry 2
    !! leaves it a singleton and (1,3) join at 1
    subroutine test_hac_mask()
        real, parameter      :: X(4) = [0., 0.5, 1., 10.]
        integer, allocatable :: i_medoids(:), labels(:)
        type(hac) :: hc
        write(*,'(A)') 'test_hac_mask'
        call hc%new(size(X), line_dmat(X), 'average', thres=2._dp)
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == [1,1,1,2]), 'hac without a mask joins the three close points')
        call hc%new(size(X), line_dmat(X), 'average', thres=2._dp, mask=[.true., .false., .true., .true.])
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == [1,2,1,3]), 'hac leaves a masked-out entry a singleton')
        call assert_int(3, hc%get_nclust(), 'hac counts the masked singleton as a cluster')
        call hc%kill
    end subroutine test_hac_mask

    !> nclust = n: every item its own cluster, labelled by index; nclust = 1: one cluster whose medoid
    !! is the global medoid (x = 1 of 0, 1, 5). One object serves both: new on a live object replaces
    !! it, and cluster may be called twice on the same construction
    subroutine test_hac_edges()
        real, parameter      :: X(3) = [0., 1., 5.]
        integer, allocatable :: i_medoids(:), labels(:)
        type(hac)            :: hc
        write(*,'(A)') 'test_hac_edges'
        call hc%new(3, line_dmat(X), 'average', nclust=3)
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == [1,2,3]) .and. all(i_medoids == [1,2,3]), 'hac with nclust = n is the identity')
        call hc%new(3, line_dmat(X), 'average', nclust=1)
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == 1), 'hac with nclust = 1 gives one cluster')
        call assert_int(2, i_medoids(1), 'hac with nclust = 1 has the global medoid')
        call hc%cluster(i_medoids, labels)
        call assert_true(all(labels == 1) .and. i_medoids(1) == 2, 'hac reclustering one construction repeats the result')
        call hc%kill
    end subroutine test_hac_edges

end module simple_hac_tester
