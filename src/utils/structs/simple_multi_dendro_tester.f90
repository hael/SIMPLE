!@descr: unit test routines for simple_multi_dendro class (refactored: per-tree build from externally provided sub_distmat)
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
        call test_get_tree_refs()
        call test_build_tree_singletons()
        call test_build_tree_pairs_medoids()
        call test_build_tree_two_clusters_medoids()
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
        !  {1,2} -> tie => medoid should be 1 (binary_tree medoid logic picks first on tie)
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

    subroutine assert_int_in_set(x, setv, msg)
        integer, intent(in) :: x
        integer, intent(in) :: setv(:)
        character(*), intent(in) :: msg
        call assert_true(any(setv == x), msg)
    end subroutine assert_int_in_set

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
        integer :: labels(5)
        write(*,'(A)') 'test_new_counts_and_getters'
        labels = [1,2,1,3,2]  ! {1,3}, {2,5}, {4}
        call md%new(labels)
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

    subroutine test_get_tree_refs()
        type(multi_dendro) :: md
        integer :: labels(6)
        integer, allocatable :: refs(:)
        write(*,'(A)') 'test_get_tree_refs'
        labels = [1,1,1, 2,2,2]
        call md%new(labels)
        refs = md%get_tree_refs(1)
        call assert_int(3, size(refs), 'tree1 refs size=3')
        call assert_int_in_set(1, refs, 'tree1 contains ref 1')
        call assert_int_in_set(2, refs, 'tree1 contains ref 2')
        call assert_int_in_set(3, refs, 'tree1 contains ref 3')
        deallocate(refs)
        refs = md%get_tree_refs(2)
        call assert_int(3, size(refs), 'tree2 refs size=3')
        call assert_int_in_set(4, refs, 'tree2 contains ref 4')
        call assert_int_in_set(5, refs, 'tree2 contains ref 5')
        call assert_int_in_set(6, refs, 'tree2 contains ref 6')
        deallocate(refs)
        refs = md%get_tree_refs(0)
        call assert_int(0, size(refs), 'get_tree_refs(0) returns empty')
        deallocate(refs)
        refs = md%get_tree_refs(3)
        call assert_int(0, size(refs), 'get_tree_refs(out-of-range) returns empty')
        deallocate(refs)
        call md%kill()
    end subroutine test_get_tree_refs

    subroutine test_build_tree_singletons()
        type(multi_dendro) :: md
        integer :: labels(4)
        integer :: itree
        integer, allocatable :: refs(:)
        real, allocatable :: submat(:,:)
        integer :: l, r
        write(*,'(A)') 'test_build_tree_singletons'
        labels = [1,2,3,4]
        call md%new(labels)
        call assert_int(4, md%get_n_trees(), 'n_trees=4')
        call assert_int(4, md%get_n_refs(),  'n_refs=4')
        do itree = 1, 4
            refs = md%get_tree_refs(itree)
            call assert_int(1, size(refs), 'singleton refs size=1')
            allocate(submat(1,1))
            submat = 0.0
            call md%build_tree_from_subdist(itree, refs, submat, linkage=1)
            call assert_int(1, md%get_cls_pop(itree), 'singleton cls_pop=1')
            call assert_int(refs(1), md%get_medoid(itree), 'singleton medoid equals only ref')
            call md%get_left_right_idxs(refs(1), l, r)
            call assert_pair(0,0, l, r, 'singleton left/right are 0,0')
            deallocate(submat)
            deallocate(refs)
        end do
        call md%kill()
    end subroutine test_build_tree_singletons

    subroutine test_build_tree_pairs_medoids()
        type(multi_dendro) :: md
        integer :: labels(4)
        integer, allocatable :: refs(:)
        real :: dist(4,4)
        real, allocatable :: submat(:,:)
        write(*,'(A)') 'test_build_tree_pairs_medoids'
        call make_distmat_two_pairs(dist)
        labels = [1,1,2,2]   ! {1,2} and {3,4}
        call md%new(labels)
        ! Tree 1
        refs = md%get_tree_refs(1)           ! expected [1,2]
        call assert_int(2, size(refs), 'tree1 refs size=2')
        allocate(submat(2,2))
        submat(1,1)=0.0; submat(2,2)=0.0
        submat(1,2)=dist(refs(1), refs(2))
        submat(2,1)=submat(1,2)
        call md%build_tree_from_subdist(1, refs, submat, linkage=1)
        ! tie -> first
        call assert_int(refs(1), md%get_medoid(1), 'pair tree1 medoid is first (tie)')
        deallocate(submat); deallocate(refs)
        ! Tree 2
        refs = md%get_tree_refs(2)           ! expected [3,4]
        call assert_int(2, size(refs), 'tree2 refs size=2')
        allocate(submat(2,2))
        submat(1,1)=0.0; submat(2,2)=0.0
        submat(1,2)=dist(refs(1), refs(2))
        submat(2,1)=submat(1,2)
        call md%build_tree_from_subdist(2, refs, submat, linkage=1)
        call assert_int(refs(1), md%get_medoid(2), 'pair tree2 medoid is first (tie)')
        deallocate(submat); deallocate(refs)
        call assert_int(2, md%get_n_trees(), 'pairs: n_trees=2')
        call assert_int(4, md%get_n_refs(),  'pairs: n_refs=4')
        call assert_int(2, md%get_cls_pop(1),'pairs: cls_pop(1)=2')
        call assert_int(2, md%get_cls_pop(2),'pairs: cls_pop(2)=2')
        call md%kill()
    end subroutine test_build_tree_pairs_medoids

    subroutine test_build_tree_two_clusters_medoids()
        type(multi_dendro) :: md
        integer :: labels(6)
        integer, allocatable :: refs(:)
        real :: dist(6,6)
        real, allocatable :: submat(:,:)
        integer :: i, j
        integer :: l, r
        write(*,'(A)') 'test_build_tree_two_clusters_medoids'
        call make_distmat_two_triples(dist)
        labels = [1,1,1, 2,2,2]
        call md%new(labels)
        ! Build tree 1 with external submatrix
        refs = md%get_tree_refs(1)   ! [1,2,3]
        call assert_int(3, size(refs), 'tree1 refs size=3')
        allocate(submat(3,3))
        submat = 0.0
        do i=1,3
            do j=i+1,3
                submat(i,j) = dist(refs(i), refs(j))
                submat(j,i) = submat(i,j)
            end do
        end do
        call md%build_tree_from_subdist(1, refs, submat, linkage=1)
        call assert_int(2, md%get_medoid(1), 'cluster1 medoid=2')
        deallocate(submat); deallocate(refs)
        ! Build tree 2 with external submatrix
        refs = md%get_tree_refs(2)   ! [4,5,6]
        call assert_int(3, size(refs), 'tree2 refs size=3')
        allocate(submat(3,3))
        submat = 0.0
        do i=1,3
            do j=i+1,3
                submat(i,j) = dist(refs(i), refs(j))
                submat(j,i) = submat(i,j)
            end do
        end do
        call md%build_tree_from_subdist(2, refs, submat, linkage=1)
        call assert_int(5, md%get_medoid(2), 'cluster2 medoid=5')
        deallocate(submat); deallocate(refs)
        call assert_int(2, md%get_n_trees(), 'n_trees=2')
        call assert_int(6, md%get_n_refs(),  'n_refs=6')
        call assert_int(3, md%get_cls_pop(1),'cls_pop(1)=3')
        call assert_int(3, md%get_cls_pop(2),'cls_pop(2)=3')
        ! leaf-hit behavior (binary_tree finds leaves first)
        call md%get_left_right_idxs(1, l, r)
        call assert_pair(0,0, l, r, 'leaf ref=1 left/right 0,0')
        call md%get_left_right_idxs(6, l, r)
        call assert_pair(0,0, l, r, 'leaf ref=6 left/right 0,0')
        call md%kill()
    end subroutine test_build_tree_two_clusters_medoids

    subroutine test_get_left_right_idxs_not_found()
        type(multi_dendro) :: md
        integer :: labels(5)
        integer, allocatable :: refs(:)
        real, allocatable :: submat(:,:)
        integer :: itree, n
        integer :: l, r, i, j
        write(*,'(A)') 'test_get_left_right_idxs_not_found'
        ! labels: 3 trees: {1,3}, {2,5}, {4}
        labels = [1,2,1,3,2]
        call md%new(labels)
        ! Build all trees externally using abs(i-j) distances (simple deterministic)
        do itree = 1, md%get_n_trees()
            refs = md%get_tree_refs(itree)
            n = size(refs)
            if (n == 0) then
                deallocate(refs)
                cycle
            end if
            allocate(submat(n,n))
            submat = 0.0
            if (n > 1) then
                do i = 1, n
                    do j = i+1, n
                        submat(i,j) = real(abs(refs(i) - refs(j)))
                        submat(j,i) = submat(i,j)
                    end do
                end do
            end if
            call md%build_tree_from_subdist(itree, refs, submat, linkage=1)
            deallocate(submat)
            deallocate(refs)
        end do
        call md%get_left_right_idxs(999, l, r)
        call assert_pair(0,0, l, r, 'unknown ref returns 0,0')
        call md%kill()
    end subroutine test_get_left_right_idxs_not_found

    subroutine test_kill_and_reuse()
        type(multi_dendro) :: md
        integer :: labels(6)
        integer, allocatable :: refs(:)
        real :: dist(6,6)
        real, allocatable :: submat(:,:)
        integer :: itree, n, i, j
        write(*,'(A)') 'test_kill_and_reuse'
        call make_distmat_two_triples(dist)
        labels = [1,1,1,2,2,2]
        call md%new(labels)
        ! Build trees externally
        do itree = 1, md%get_n_trees()
            refs = md%get_tree_refs(itree)
            n = size(refs)
            allocate(submat(n,n))
            submat = 0.0
            do i=1,n
                do j=i+1,n
                    submat(i,j) = dist(refs(i), refs(j))
                    submat(j,i) = submat(i,j)
                end do
            end do
            call md%build_tree_from_subdist(itree, refs, submat, linkage=1)
            deallocate(submat)
            deallocate(refs)
        end do
        call md%kill()
        call assert_int(0, md%get_n_trees(), 'after kill n_trees=0')
        call assert_int(0, md%get_n_refs(),  'after kill n_refs=0')
        call assert_int(0, md%get_medoid(1), 'after kill get_medoid(1)=0')
        call assert_int(0, md%get_cls_pop(1),'after kill get_cls_pop(1)=0')
        ! Reuse: singleton + larger cluster (still use dist from make_distmat_two_triples)
        labels = [1,2,2,2,2,2]
        call md%new(labels)
        do itree = 1, md%get_n_trees()
            refs = md%get_tree_refs(itree)
            n = size(refs)
            allocate(submat(n,n))
            submat = 0.0
            do i=1,n
                do j=i+1,n
                    submat(i,j) = dist(refs(i), refs(j))
                    submat(j,i) = submat(i,j)
                end do
            end do
            call md%build_tree_from_subdist(itree, refs, submat, linkage=1)
            deallocate(submat)
            deallocate(refs)
        end do
        call assert_int(2, md%get_n_trees(), 'reuse n_trees=2')
        call assert_int(6, md%get_n_refs(),  'reuse n_refs=6')
        call assert_int(1, md%get_cls_pop(1),'reuse cls_pop(1)=1')
        call assert_int(5, md%get_cls_pop(2),'reuse cls_pop(2)=5')
        call assert_int(1, md%get_medoid(1), 'reuse singleton medoid=1')
        call md%kill()
    end subroutine test_kill_and_reuse

end module simple_multi_dendro_tester
