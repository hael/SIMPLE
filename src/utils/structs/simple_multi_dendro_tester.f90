!@descr: unit test routines for simple_multi_dendro class (refactored)
module simple_multi_dendro_tester
use simple_multi_dendro
use simple_test_utils
use simple_binary_tree, only: bt_node
implicit none

public :: run_all_multi_dendro_tests
private

contains

    subroutine run_all_multi_dendro_tests()
        write(*,'(A)') '**** running multi_dendro tests ****'
        call test_initial_state()
        call test_new_counts_and_getters()
        call test_get_tree_refs()
        call test_build_tree_singletons_and_node_api()
        call test_build_tree_pairs_medoids_and_node_api()
        call test_build_tree_two_clusters_medoids_and_node_api()
        call test_local_to_global_wrapper()
        call test_kill_and_reuse()
        write(*,'(A)') '**** finished multi_dendro tests ****'
    end subroutine run_all_multi_dendro_tests

    !===========================
    ! Helpers
    !===========================

    subroutine make_distmat_two_triples(dist)
        ! Deterministic medoids:
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
        ! {1,2} -> tie => medoid is first (1)
        ! {3,4} -> tie => medoid is first (3)
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

    subroutine build_external_tree(md, itree, dist_full)
        type(multi_dendro), intent(inout) :: md
        integer,            intent(in)    :: itree
        real,               intent(in)    :: dist_full(:,:)
        integer, allocatable :: refs(:)
        real,    allocatable :: submat(:,:)
        integer :: n, i, j

        refs = md%get_tree_refs(itree)
        n = size(refs)
        if (n == 0) then
            deallocate(refs)
            return
        end if

        allocate(submat(n,n))
        submat = 0.0
        if (n > 1) then
            do i=1,n
                do j=i+1,n
                    submat(i,j) = dist_full(refs(i), refs(j))
                    submat(j,i) = submat(i,j)
                end do
            end do
        end if

        call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)

        deallocate(submat)
        deallocate(refs)
    end subroutine build_external_tree

    !===========================
    ! Tests
    !===========================

    subroutine test_initial_state()
        type(multi_dendro) :: md
        write(*,'(A)') 'test_initial_state'
        call assert_int(0, md%get_n_trees(), 'initial get_n_trees=0')
        call assert_int(0, md%get_n_refs(),  'initial get_n_refs=0')
        call assert_int(0, md%get_medoid(1), 'initial get_medoid(1)=0')
        call assert_int(0, md%get_tree_pop(1),'initial get_tree_pop(1)=0')
        call assert_int(0, md%get_root_node_idx(1), 'initial root_node_idx=0')
    end subroutine test_initial_state

    subroutine test_new_counts_and_getters()
        type(multi_dendro) :: md
        integer :: labels(5)
        write(*,'(A)') 'test_new_counts_and_getters'
        labels = [1,2,1,3,2]  ! {1,3}, {2,5}, {4}
        call md%new(labels)
        call assert_int(3, md%get_n_trees(), 'n_trees from labels')
        call assert_int(5, md%get_n_refs(),  'n_refs from labels')
        call assert_int(2, md%get_tree_pop(1), 'tree_pop(1)=2')
        call assert_int(2, md%get_tree_pop(2), 'tree_pop(2)=2')
        call assert_int(1, md%get_tree_pop(3), 'tree_pop(3)=1')
        call assert_int(0, md%get_tree_pop(0), 'tree_pop(0)=0')
        call assert_int(0, md%get_tree_pop(4), 'tree_pop(4)=0')
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

    subroutine test_build_tree_singletons_and_node_api()
        type(multi_dendro) :: md
        integer :: labels(4)
        integer :: itree
        integer, allocatable :: refs(:)
        real, allocatable :: submat(:,:)
        integer :: root, lidx, ridx
        type(bt_node) :: node
        write(*,'(A)') 'test_build_tree_singletons_and_node_api'

        labels = [1,2,3,4]
        call md%new(labels)
        call assert_int(4, md%get_n_trees(), 'n_trees=4')
        call assert_int(4, md%get_n_refs(),  'n_refs=4')

        do itree = 1, 4
            refs = md%get_tree_refs(itree)
            call assert_int(1, size(refs), 'singleton refs size=1')

            allocate(submat(1,1))
            submat = 0.0
            call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)

            call assert_int(1, md%get_tree_pop(itree), 'singleton tree_pop=1')
            call assert_int(refs(1), md%get_medoid(itree), 'singleton medoid equals only ref')

            root = md%get_root_node_idx(itree)
            call assert_int(1, root, 'singleton root node idx=1')
            call assert_true(md%is_leaf(itree, root), 'singleton root is leaf')
            call assert_int(1, md%get_subset_size(itree, root), 'singleton subset size=1')

            node = md%get_node_by_idx(itree, root)
            call assert_int(refs(1), node%ref_idx, 'singleton node ref is only ref')

            call md%get_children_idx(itree, root, lidx, ridx)
            call assert_pair(0,0, lidx, ridx, 'singleton children idx are 0,0')

            deallocate(submat)
            deallocate(refs)
        end do

        call md%kill()
    end subroutine test_build_tree_singletons_and_node_api

    subroutine test_build_tree_pairs_medoids_and_node_api()
        type(multi_dendro) :: md
        integer :: labels(4)
        integer, allocatable :: refs(:)
        real :: dist(4,4)
        real, allocatable :: submat(:,:)
        integer :: itree, root, lidx, ridx
        integer :: lref, rref
        type(bt_node) :: node
        write(*,'(A)') 'test_build_tree_pairs_medoids_and_node_api'

        call make_distmat_two_pairs(dist)
        labels = [1,1,2,2]   ! {1,2} and {3,4}
        call md%new(labels)

        do itree = 1, 2
            refs = md%get_tree_refs(itree)
            call assert_int(2, size(refs), 'pair refs size=2')

            allocate(submat(2,2))
            submat = 0.0
            submat(1,2) = dist(refs(1), refs(2))
            submat(2,1) = submat(1,2)

            call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)

            ! tie -> first
            call assert_int(refs(1), md%get_medoid(itree), 'pair medoid is first (tie)')

            root = md%get_root_node_idx(itree)
            call assert_int(3, root, 'pair tree root idx = 3 (2*n-1)')
            call assert_true(.not. md%is_leaf(itree, root), 'pair root is not leaf')
            call assert_int(2, md%get_subset_size(itree, root), 'pair root subset size=2')

            node = md%get_node_by_idx(itree, root)
            call assert_int(refs(1), node%ref_idx, 'pair root node medoid ref matches')

            call md%get_children_idx(itree, root, lidx, ridx)
            call assert_pair(1,2, lidx, ridx, 'pair root children idx are leaves 1,2')

            call md%get_children_ref(itree, root, lref, rref)
            call assert_int(refs(1), lref, 'pair left child ref is refs(1)')
            call assert_int(refs(2), rref, 'pair right child ref is refs(2)')

            deallocate(submat)
            deallocate(refs)
        end do

        call md%kill()
    end subroutine test_build_tree_pairs_medoids_and_node_api

    subroutine test_build_tree_two_clusters_medoids_and_node_api()
        type(multi_dendro) :: md
        integer :: labels(6)
        real :: dist(6,6)
        integer :: itree, root, lidx, ridx
        integer :: med
        type(bt_node) :: node
        write(*,'(A)') 'test_build_tree_two_clusters_medoids_and_node_api'

        call make_distmat_two_triples(dist)
        labels = [1,1,1, 2,2,2]
        call md%new(labels)

        ! build both trees
        do itree = 1, md%get_n_trees()
            call build_external_tree(md, itree, dist)
        end do

        call assert_int(2, md%get_n_trees(), 'n_trees=2')
        call assert_int(6, md%get_n_refs(),  'n_refs=6')
        call assert_int(3, md%get_tree_pop(1),'tree_pop(1)=3')
        call assert_int(3, md%get_tree_pop(2),'tree_pop(2)=3')

        ! medoids deterministic
        call assert_int(2, md%get_medoid(1), 'cluster1 medoid=2')
        call assert_int(5, md%get_medoid(2), 'cluster2 medoid=5')

        ! Node API checks on tree 1
        itree = 1
        root = md%get_root_node_idx(itree)
        call assert_int(5, root, 'tree1 root idx = 2*n-1 = 5')
        call assert_true(.not. md%is_leaf(itree, root), 'tree1 root is not leaf')
        call assert_int(3, md%get_subset_size(itree, root), 'tree1 root subset size=3')
        node = md%get_node_by_idx(itree, root)
        call assert_int(2, node%ref_idx, 'tree1 root node medoid=2')

        call md%get_children_idx(itree, root, lidx, ridx)
        call assert_true(lidx /= 0 .and. ridx /= 0, 'tree1 root has two children')

        ! Node API checks on tree 2
        itree = 2
        root = md%get_root_node_idx(itree)
        call assert_int(5, root, 'tree2 root idx = 5')
        call assert_int(3, md%get_subset_size(itree, root), 'tree2 root subset size=3')
        node = md%get_node_by_idx(itree, root)
        call assert_int(5, node%ref_idx, 'tree2 root node medoid=5')

        call md%kill()
    end subroutine test_build_tree_two_clusters_medoids_and_node_api

    subroutine test_local_to_global_wrapper()
        type(multi_dendro) :: md
        integer :: labels(4)
        integer, allocatable :: refs(:)
        real :: dist(4,4)
        real, allocatable :: submat(:,:)
        integer :: g
        write(*,'(A)') 'test_local_to_global_wrapper'

        call make_distmat_two_pairs(dist)
        labels = [1,1,2,2]
        call md%new(labels)

        ! build tree 1 so mapping exists in binary_tree
        refs = md%get_tree_refs(1)   ! expected [1,2]
        allocate(submat(2,2))
        submat = 0.0
        submat(1,2) = dist(refs(1), refs(2))
        submat(2,1) = submat(1,2)
        call md%build_tree_from_subdistmat(1, refs, submat, linkage=1)
        deallocate(submat)

        ! In this test refs=[1,2] so local->global is identity
        g = md%local_to_global_ref(1, 1)
        call assert_int(1, g, 'local_to_global_ref(tree1,1)=1')
        g = md%local_to_global_ref(1, 2)
        call assert_int(2, g, 'local_to_global_ref(tree1,2)=2')

        ! out of range should return 0
        g = md%local_to_global_ref(1, 0)
        call assert_int(0, g, 'local_to_global_ref(tree1,0)=0')
        g = md%local_to_global_ref(1, 999)
        call assert_int(0, g, 'local_to_global_ref(tree1,999)=0')

        deallocate(refs)
        call md%kill()
    end subroutine test_local_to_global_wrapper

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

        do itree = 1, md%get_n_trees()
            call build_external_tree(md, itree, dist)
        end do

        call md%kill()
        call assert_int(0, md%get_n_trees(), 'after kill n_trees=0')
        call assert_int(0, md%get_n_refs(),  'after kill n_refs=0')
        call assert_int(0, md%get_medoid(1), 'after kill get_medoid(1)=0')
        call assert_int(0, md%get_tree_pop(1),'after kill get_tree_pop(1)=0')
        call assert_int(0, md%get_root_node_idx(1), 'after kill root_node_idx=0')

        ! Reuse: singleton + larger cluster
        labels = [1,2,2,2,2,2]
        call md%new(labels)

        do itree = 1, md%get_n_trees()
            refs = md%get_tree_refs(itree)
            n = size(refs)
            allocate(submat(n,n))
            submat = 0.0
            if (n > 1) then
                do i=1,n
                    do j=i+1,n
                        submat(i,j) = dist(refs(i), refs(j))
                        submat(j,i) = submat(i,j)
                    end do
                end do
            end if
            call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)
            deallocate(submat)
            deallocate(refs)
        end do

        call assert_int(2, md%get_n_trees(), 'reuse n_trees=2')
        call assert_int(6, md%get_n_refs(),  'reuse n_refs=6')
        call assert_int(1, md%get_tree_pop(1),'reuse tree_pop(1)=1')
        call assert_int(5, md%get_tree_pop(2),'reuse tree_pop(2)=5')
        call assert_int(1, md%get_medoid(1), 'reuse singleton medoid=1')

        call md%kill()
    end subroutine test_kill_and_reuse

end module simple_multi_dendro_tester