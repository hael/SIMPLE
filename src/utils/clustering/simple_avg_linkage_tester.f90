!@descr: unit tests for average-linkage hierarchical clustering (simple_avg_linkage): groups, medoids, UPGMA update and edges
! Points on a line, distance |x_i - x_j|; every expectation follows from the UPGMA merge order by hand.
! Separated groups are recovered with labels by decreasing population and their medoids; a dense group of
! twelve items stays one cluster beside four far singletons; a five-point case separates the
! population-weighted update (UPGMA) from the unweighted one (WPGMA); nclust = n and nclust = 1 are the
! edges, run on one object that is reconstructed in between.
module simple_avg_linkage_tester
use simple_test_utils
use simple_avg_linkage, only: avg_linkage
implicit none
private
public :: run_all_avg_linkage_tests

contains

    subroutine run_all_avg_linkage_tests()
        write(*,'(A)') '**** running all avg_linkage tests ****'
        call test_avglink_separated_groups()
        call test_avglink_dense_group_stays_whole()
        call test_avglink_population_weighted_update()
        call test_avglink_edges()
    end subroutine run_all_avg_linkage_tests

    !> clusters the points x on a line into nclust groups
    subroutine run_line( x, nclust, i_medoids, labels )
        real,                 intent(in)    :: x(:)
        integer,              intent(in)    :: nclust
        integer, allocatable, intent(inout) :: i_medoids(:), labels(:)
        type(avg_linkage) :: avglink
        call avglink%new(size(x), line_dmat(x), nclust)
        call avglink%cluster(i_medoids, labels)
        call avglink%kill
    end subroutine run_line

    !> |x_i - x_j| for points on a line
    function line_dmat( x ) result( dmat )
        real, intent(in)  :: x(:)
        real, allocatable :: dmat(:,:)
        integer :: i, j
        allocate(dmat(size(x),size(x)))
        do j = 1, size(x)
            do i = 1, size(x)
                dmat(i,j) = abs(x(i) - x(j))
            enddo
        enddo
    end function line_dmat

    !> groups {0,1,2}, {50,51}, {100,101,102,103}: six within-group merges (distances <= 3) precede
    !! any merge across groups (>= 48); labels by decreasing population, the medoid is the member with
    !! the smallest summed distance, the first one on a tie (101 before 102 in the largest group)
    subroutine test_avglink_separated_groups()
        real, parameter      :: X(9) = [0., 1., 2., 50., 51., 100., 101., 102., 103.]
        integer, allocatable :: i_medoids(:), labels(:)
        write(*,'(A)') 'test_avglink_separated_groups'
        call run_line(X, 3, i_medoids, labels)
        call assert_int(3, size(i_medoids), 'avglink forms the requested number of clusters')
        call assert_true(all(labels == [2,2,2,3,3,1,1,1,1]), 'avglink recovers the groups, labelled by decreasing population')
        call assert_true(all(i_medoids == [7,2,4]), 'avglink medoids are the minimal-sum members, first on a tie')
    end subroutine test_avglink_separated_groups

    !> twelve items at 0..11 and four at 1000..4000: the eleven merges that take sixteen items to five
    !! clusters all fall inside the dense group (mean distances <= 11 against >= 989), so it is one
    !! cluster; its medoid is the first median item (x = 5)
    subroutine test_avglink_dense_group_stays_whole()
        real    :: x(16)
        integer, allocatable :: i_medoids(:), labels(:)
        integer :: i
        write(*,'(A)') 'test_avglink_dense_group_stays_whole'
        x(1:12)  = [(real(i), i=0,11)]
        x(13:16) = [1000., 2000., 3000., 4000.]
        call run_line(x, 5, i_medoids, labels)
        call assert_true(all(labels(1:12) == 1), 'avglink keeps the dense group in one cluster')
        call assert_true(all(labels(13:16) == [2,3,4,5]), 'avglink leaves the far items as singletons, by index')
        call assert_int(6, i_medoids(1), 'avglink medoid of the dense group is its first median item')
    end subroutine test_avglink_dense_group_stays_whole

    !> x = 0, 0, 10, 4.6, 18. The merges (1,2) at 0 and ({1,2},4) at 4.6 give T = {1,2,4}. Population
    !! weighting gives d(T,3) = (2*10 + 5.4)/3 = 8.47 > d(3,5) = 8, so 3 joins 5; the unweighted mean
    !! (10 + 5.4)/2 = 7.7 would have joined 3 to T instead
    subroutine test_avglink_population_weighted_update()
        real, parameter      :: X(5) = [0., 0., 10., 4.6, 18.]
        integer, allocatable :: i_medoids(:), labels(:)
        write(*,'(A)') 'test_avglink_population_weighted_update'
        call run_line(X, 2, i_medoids, labels)
        call assert_true(all(labels == [1,1,2,1,2]), 'avglink weights the linkage update by cluster population')
    end subroutine test_avglink_population_weighted_update

    !> nclust = n: every item its own cluster, labelled by index; nclust = 1: one cluster whose medoid
    !! is the global medoid (x = 1 of 0, 1, 5). One object serves both: new on a live object replaces
    !! it, and cluster may be called twice on the same construction
    subroutine test_avglink_edges()
        real, parameter      :: X(3) = [0., 1., 5.]
        integer, allocatable :: i_medoids(:), labels(:)
        type(avg_linkage)    :: avglink
        write(*,'(A)') 'test_avglink_edges'
        call avglink%new(3, line_dmat(X), 3)
        call avglink%cluster(i_medoids, labels)
        call assert_true(all(labels == [1,2,3]) .and. all(i_medoids == [1,2,3]), 'avglink with nclust = n is the identity')
        call avglink%new(3, line_dmat(X), 1)
        call avglink%cluster(i_medoids, labels)
        call assert_true(all(labels == 1), 'avglink with nclust = 1 gives one cluster')
        call assert_int(2, i_medoids(1), 'avglink with nclust = 1 has the global medoid')
        call avglink%cluster(i_medoids, labels)
        call assert_true(all(labels == 1) .and. i_medoids(1) == 2, 'avglink reclustering one construction repeats the result')
        call avglink%kill
    end subroutine test_avglink_edges

end module simple_avg_linkage_tester
