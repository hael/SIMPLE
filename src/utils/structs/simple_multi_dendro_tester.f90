!@descr: unit test routines for simple_multi_dendro class
module simple_multi_dendro_tester
use simple_multi_dendro
use simple_test_utils
implicit none

public :: run_all_multi_dendro_tests
private

contains

    subroutine run_all_multi_dendro_tests()
        write(*,'(A)') '**** running multi_dendro tests ****'
        call test_initial_state()
        call test_new_counts_and_getters()
        call test_build_multi_dendro_singletons()
        call test_build_multi_dendro_pairs_medoids()
        call test_build_multi_dendro_two_clusters_medoids()
        call test_get_left_right_idxs_not_found()
        call test_kill_and_reuse()
        ! call report_summary()
        write(*,'(A)') '**** finished multi_dendro tests ****'
    end subroutine run_all_multi_dendro_tests

    !===========================
    ! Helpers
    !===========================

    subroutine make_distmat(n, dist)
        integer, intent(in) :: n
        real,    intent(out) :: dist(n,n)
        integer :: i,j
        dist = 0.0
        do i=1,n
            dist(i,i) = 0.0
        end do
        do i=1,n
            do j=i+1,n
                dist(i,j) = real(abs(i-j))
                dist(j,i) = dist(i,j)
            end do
        end do
    end subroutine make_distmat

    subroutine make_distmat_two_triples(dist)
        ! Designed so medoids are deterministic:
        ! cluster {1,2,3} -> medoid 2
        ! cluster {4,5,6} -> medoid 5
        real, intent(out) :: dist(6,6)
        integer :: i,j
        dist = 0.0
        do i=1,6
            dist(i,i) = 0.0
        end do

        ! cluster 1
        dist(1,2)=1.0; dist(2,1)=1.0
        dist(2,3)=1.0; dist(3,2)=1.0
        dist(1,3)=2.0; dist(3,1)=2.0

        ! cluster 2
        dist(4,5)=1.0; dist(5,4)=1.0
        dist(5,6)=1.0; dist(6,5)=1.0
        dist(4,6)=2.0; dist(6,4)=2.0

        ! cross cluster large
        do i=1,3
            do j=4,6
                dist(i,j)=10.0
                dist(j,i)=10.0
            end do
        end do
    end subroutine make_distmat_two_triples

    subroutine make_distmat_two_pairs(dist)
        ! Two clusters of size 2:
        !  {1,2} -> tie => medoid should be 1 (your binary_tree medoid logic picks first on tie)
        !  {3,4} -> tie => medoid should be 3
        real, intent(out) :: dist(4,4)
        integer :: i,j
        dist = 0.0
        do i=1,4
            dist(i,i)=0.0
        end do
        dist(1,2)=1.0; dist(2,1)=1.0
        dist(3,4)=1.0; dist(4,3)=1.0
        do i=1,2
            do j=3,4
                dist(i,j)=10.0
                dist(j,i)=10.0
            end do
        end do
    end subroutine make_distmat_two_pairs

    subroutine assert_pair(expected_l, expected_r, l, r, msg)
        integer, intent(in) :: expected_l, expected_r, l, r
        character(*), intent(in) :: msg
        call assert_int(expected_l, l, trim(msg)//' left')
        call assert_int(expected_r, r, trim(msg)//' right')
    end subroutine assert_pair

    !===========================
    ! Tests
    !===========================

    subroutine test_initial_state()
        type(multi_dendro) :: md
        write(*,'(A)') 'test_initial_state'
        call assert_int(0, md%get_n_trees(), 'initial get_n_trees=0')
        call assert_int(0, md%get_n_refs(),  'initial get_n_refs=0')
        call assert_int(0, md%get_medoid(1), 'initial get_medoid(1)=0')
        call assert_int(0, md%get_cls_pop(1),'initial get_cls_pop(1)=0')
    end subroutine test_initial_state

    subroutine test_new_counts_and_getters()
        type(multi_dendro) :: md
        real :: dist(5,5)
        integer :: labels(5)

        write(*,'(A)') 'test_new_counts_and_getters'

        call make_distmat(5, dist)
        labels = [1,2,1,3,2]  ! {1,3}, {2,5}, {4}

        call md%new(dist, labels)

        call assert_int(3, md%get_n_trees(), 'n_trees from labels')
        call assert_int(5, md%get_n_refs(),  'n_refs from labels')

        call assert_int(2, md%get_cls_pop(1), 'cls_pop(1)=2')
        call assert_int(2, md%get_cls_pop(2), 'cls_pop(2)=2')
        call assert_int(1, md%get_cls_pop(3), 'cls_pop(3)=1')

        call assert_int(0, md%get_cls_pop(0), 'cls_pop(0)=0')
        call assert_int(0, md%get_cls_pop(4), 'cls_pop(4)=0')
        call assert_int(0, md%get_medoid(0),  'medoid(0)=0')
        call assert_int(0, md%get_medoid(4),  'medoid(4)=0')

        call md%kill()
    end subroutine test_new_counts_and_getters

    subroutine test_build_multi_dendro_singletons()
        type(multi_dendro) :: md
        real :: dist(4,4)
        integer :: labels(4)
        integer :: l, r, i

        write(*,'(A)') 'test_build_multi_dendro_singletons'

        call make_distmat(4, dist)
        labels = [1,2,3,4]

        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)

        call assert_int(4, md%get_n_trees(), 'n_trees=4')
        call assert_int(4, md%get_n_refs(),  'n_refs=4')

        do i=1,4
            call assert_int(1, md%get_cls_pop(i), 'singleton cls_pop=1')
            call assert_int(i, md%get_medoid(i),  'singleton medoid equals ref')
            call md%get_left_right_idxs(i, l, r)
            call assert_pair(0,0, l, r, 'singleton left/right are 0,0')
        end do

        call md%kill()
    end subroutine test_build_multi_dendro_singletons

    subroutine test_build_multi_dendro_pairs_medoids()
        type(multi_dendro) :: md
        real :: dist(4,4)
        integer :: labels(4)

        write(*,'(A)') 'test_build_multi_dendro_pairs_medoids'

        call make_distmat_two_pairs(dist)
        labels = [1,1,2,2]   ! {1,2} and {3,4}

        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)

        call assert_int(2, md%get_n_trees(), 'pairs: n_trees=2')
        call assert_int(4, md%get_n_refs(),  'pairs: n_refs=4')
        call assert_int(2, md%get_cls_pop(1),'pairs: cls_pop(1)=2')
        call assert_int(2, md%get_cls_pop(2),'pairs: cls_pop(2)=2')

        ! With tie distances in each pair, binary_tree medoid selection picks first
        call assert_int(1, md%get_medoid(1), 'pair cluster1 medoid=1 (tie -> first)')
        call assert_int(3, md%get_medoid(2), 'pair cluster2 medoid=3 (tie -> first)')

        call md%kill()
    end subroutine test_build_multi_dendro_pairs_medoids

    subroutine test_build_multi_dendro_two_clusters_medoids()
        type(multi_dendro) :: md
        real :: dist(6,6)
        integer :: labels(6)
        integer :: l, r

        write(*,'(A)') 'test_build_multi_dendro_two_clusters_medoids'

        call make_distmat_two_triples(dist)
        labels = [1,1,1, 2,2,2]

        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)

        call assert_int(2, md%get_n_trees(), 'n_trees=2')
        call assert_int(6, md%get_n_refs(),  'n_refs=6')
        call assert_int(3, md%get_cls_pop(1),'cls_pop(1)=3')
        call assert_int(3, md%get_cls_pop(2),'cls_pop(2)=3')

        call assert_int(2, md%get_medoid(1), 'cluster1 medoid=2')
        call assert_int(5, md%get_medoid(2), 'cluster2 medoid=5')

        ! leaf-hit behavior (binary_tree finds leaves first)
        call md%get_left_right_idxs(1, l, r)
        call assert_pair(0,0, l, r, 'leaf ref=1 left/right 0,0')
        call md%get_left_right_idxs(6, l, r)
        call assert_pair(0,0, l, r, 'leaf ref=6 left/right 0,0')

        call md%kill()
    end subroutine test_build_multi_dendro_two_clusters_medoids

    subroutine test_get_left_right_idxs_not_found()
        type(multi_dendro) :: md
        real :: dist(5,5)
        integer :: labels(5)
        integer :: l, r

        write(*,'(A)') 'test_get_left_right_idxs_not_found'

        call make_distmat(5, dist)
        labels = [1,2,1,3,2]
        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)

        call md%get_left_right_idxs(999, l, r)
        call assert_pair(0,0, l, r, 'unknown ref returns 0,0')

        call md%kill()
    end subroutine test_get_left_right_idxs_not_found

    subroutine test_kill_and_reuse()
        type(multi_dendro) :: md
        real :: dist(6,6)
        integer :: labels(6)

        write(*,'(A)') 'test_kill_and_reuse'

        call make_distmat_two_triples(dist)
        labels = [1,1,1,2,2,2]
        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)
        call md%kill()

        call assert_int(0, md%get_n_trees(), 'after kill n_trees=0')
        call assert_int(0, md%get_n_refs(),  'after kill n_refs=0')
        call assert_int(0, md%get_medoid(1), 'after kill get_medoid(1)=0')
        call assert_int(0, md%get_cls_pop(1),'after kill get_cls_pop(1)=0')

        labels = [1,2,2,2,2,2]  ! singleton + larger cluster
        call md%new(dist, labels)
        call md%build_multi_dendro(linkage=1)

        call assert_int(2, md%get_n_trees(), 'reuse n_trees=2')
        call assert_int(6, md%get_n_refs(),  'reuse n_refs=6')
        call assert_int(1, md%get_cls_pop(1),'reuse cls_pop(1)=1')
        call assert_int(5, md%get_cls_pop(2),'reuse cls_pop(2)=5')
        call assert_int(1, md%get_medoid(1), 'reuse singleton medoid=1')

        call md%kill()
    end subroutine test_kill_and_reuse

end module simple_multi_dendro_tester
