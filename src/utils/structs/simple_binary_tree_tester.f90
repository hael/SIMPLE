!@descr: unit test routines for refactored simple_binary_tree class (no subsets; node-only getters)
module simple_binary_tree_tester
use simple_binary_tree
use simple_test_utils
implicit none

public :: run_all_tree_tests
private

contains

    subroutine run_all_tree_tests()
        write(*,'(A)') '**** running binary_tree tests ****'
        call test_initial_state()
        call test_build_from_hclust_basic_invariants()
        call test_get_node_copy_semantics()
        call test_kill_and_rebuild()
        write(*,'(A)') '**** finished binary_tree tests ****'
    end subroutine run_all_tree_tests

    !===========================
    ! Helpers
    !===========================

    subroutine build_test_tree(t)
        type(binary_tree), intent(inout) :: t
        integer, allocatable :: merge_mat(:,:)
        integer :: refs(4)
        real, allocatable :: dist(:,:)

        ! Leaves: 1..4, internal: 5..7
        ! merges: (1,2)->5, (3,4)->6, (5,6)->7
        allocate(merge_mat(2,3))
        merge_mat(1,:) = [1, 3, 5]
        merge_mat(2,:) = [2, 4, 6]

        refs = [1,2,3,4]

        allocate(dist(4,4))
        dist = 0.0
        dist(1,2)=1.0;  dist(2,1)=1.0
        dist(3,4)=1.0;  dist(4,3)=1.0
        dist(1,3)=10.0; dist(3,1)=10.0
        dist(1,4)=10.0; dist(4,1)=10.0
        dist(2,3)=9.0;  dist(3,2)=9.0
        dist(2,4)=9.0;  dist(4,2)=9.0

        call t%build_from_hclust(merge_mat, refs, dist)

        deallocate(merge_mat)
        deallocate(dist)
    end subroutine build_test_tree

    subroutine assert_idx_valid(idx, n, msg)
        integer, intent(in) :: idx, n
        character(*), intent(in) :: msg
        call assert_true(idx >= 1 .and. idx <= n, trim(msg))
    end subroutine assert_idx_valid

    !===========================
    ! Tests
    !===========================

    subroutine test_initial_state()
        type(binary_tree) :: t
        write(*,'(A)') 'test_initial_state'
        call assert_int(0, t%get_n_nodes(), 'new tree has 0 nodes')
        call assert_int(0, t%get_nrefs(),   'new tree nrefs=0')
        call assert_int(0, t%get_height(),  'new tree height=0')

        ! root node getter should return zero-initialized node
        block
            type(bt_node) :: r
            r = t%get_root_node()
            call assert_int(0, r%node_idx, 'new tree root node_idx=0')
            call assert_int(0, r%ref_idx,  'new tree root ref_idx=0')
            call assert_int(0, r%left_idx, 'new tree root left_idx=0')
            call assert_int(0, r%right_idx,'new tree root right_idx=0')
            call assert_int(0, r%parent_idx,'new tree root parent_idx=0')
        end block
    end subroutine test_initial_state

    subroutine test_build_from_hclust_basic_invariants()
        type(binary_tree) :: t
        type(bt_node) :: root, n1, n2, n3, n4, n5, n6
        integer :: nrefs, n_total

        write(*,'(A)') 'test_build_from_hclust_basic_invariants'
        call build_test_tree(t)

        nrefs = t%get_nrefs()
        call assert_int(4, nrefs, 'nrefs stored')
        n_total = 2*nrefs - 1
        call assert_int(n_total, t%get_n_nodes(), 'n_nodes = 2*nrefs-1')

        ! For a full 4-leaf binary merge tree, height must be 3 (nodes along root->leaf path)
        call assert_int(3, t%get_height(), 'height=3 for 4-leaf full merge')

        root = t%get_root_node()
        call assert_int(n_total, root%node_idx, 'root node_idx is last slot')
        call assert_int(0, root%parent_idx, 'root parent_idx=0')
        call assert_true(.not. t%is_leaf(root%node_idx), 'root is not leaf')
        call assert_idx_valid(root%left_idx,  n_total, 'root left_idx valid')
        call assert_idx_valid(root%right_idx, n_total, 'root right_idx valid')
        call assert_true(root%left_idx /= 0 .and. root%right_idx /= 0, 'root has two children')

        ! Leaves: nodes 1..4 are leaves with refs 1..4 (GLOBAL ids)
        n1 = t%get_node(1)
        call assert_true(t%is_leaf(1), 'node1 is leaf')
        call assert_int(1, n1%ref_idx, 'leaf1 ref_idx=1')
        call assert_int(0, n1%left_idx, 'leaf1 left=0')
        call assert_int(0, n1%right_idx,'leaf1 right=0')
        call assert_int(5, n1%parent_idx,'leaf1 parent=5')

        n2 = t%get_node(2)
        call assert_true(t%is_leaf(2), 'node2 is leaf')
        call assert_int(2, n2%ref_idx, 'leaf2 ref_idx=2')
        call assert_int(5, n2%parent_idx,'leaf2 parent=5')

        n3 = t%get_node(3)
        call assert_true(t%is_leaf(3), 'node3 is leaf')
        call assert_int(3, n3%ref_idx, 'leaf3 ref_idx=3')
        call assert_int(6, n3%parent_idx,'leaf3 parent=6')

        n4 = t%get_node(4)
        call assert_true(t%is_leaf(4), 'node4 is leaf')
        call assert_int(4, n4%ref_idx, 'leaf4 ref_idx=4')
        call assert_int(6, n4%parent_idx,'leaf4 parent=6')

        ! Internal node 5 merges (1,2)
        n5 = t%get_node(5)
        call assert_true(.not. t%is_leaf(5), 'node5 is internal')
        call assert_int(1, n5%left_idx,  'node5 left=1')
        call assert_int(2, n5%right_idx, 'node5 right=2')
        call assert_int(7, n5%parent_idx,'node5 parent=7')
        call assert_int(1, n5%level,     'node5 level=1')
        ! Deterministic medoid for {1,2} tie -> local 1 -> global 1
        call assert_int(1, n5%ref_idx,   'node5 ref_idx (medoid)=1')

        ! Internal node 6 merges (3,4)
        n6 = t%get_node(6)
        call assert_true(.not. t%is_leaf(6), 'node6 is internal')
        call assert_int(3, n6%left_idx,  'node6 left=3')
        call assert_int(4, n6%right_idx, 'node6 right=4')
        call assert_int(7, n6%parent_idx,'node6 parent=7')
        call assert_int(2, n6%level,     'node6 level=2')
        ! Deterministic medoid for {3,4} tie -> local 3 -> global 3
        call assert_int(3, n6%ref_idx,   'node6 ref_idx (medoid)=3')

        ! Root node 7 merges (5,6)
        call assert_int(5, root%left_idx,  'root left=5')
        call assert_int(6, root%right_idx, 'root right=6')
        call assert_int(3, root%level,     'root level=3')
        ! Deterministic medoid for {1,2,3,4} is global 2
        call assert_int(2, root%ref_idx,   'root ref_idx (medoid)=2')

        call t%kill()
    end subroutine test_build_from_hclust_basic_invariants

    subroutine test_get_node_copy_semantics()
        type(binary_tree) :: t
        type(bt_node) :: a, b
        write(*,'(A)') 'test_get_node_copy_semantics'
        call build_test_tree(t)

        a = t%get_node(5)
        call assert_int(5, a%node_idx, 'get_node returns correct node_idx')
        call assert_int(1, a%ref_idx,  'get_node returns correct ref_idx')

        ! mutate local copy; tree must be unchanged
        a%ref_idx = 999
        b = t%get_node(5)
        call assert_int(1, b%ref_idx, 'tree unaffected by local mutation')

        call t%kill()
    end subroutine test_get_node_copy_semantics

    subroutine test_kill_and_rebuild()
        type(binary_tree) :: t
        type(bt_node) :: root
        integer :: n_total

        write(*,'(A)') 'test_kill_and_rebuild'
        call build_test_tree(t)

        call assert_true(t%get_n_nodes() > 0, 'tree built has nodes')
        call t%kill()
        call assert_int(0, t%get_n_nodes(), 'kill deallocates nodes')
        call assert_int(0, t%get_nrefs(),   'kill resets nrefs')
        call assert_int(0, t%get_height(),  'kill => height 0')

        call build_test_tree(t)
        call assert_int(4, t%get_nrefs(), 'rebuild nrefs ok')
        n_total = 2*t%get_nrefs() - 1
        call assert_int(n_total, t%get_n_nodes(), 'rebuild n_nodes ok')

        root = t%get_root_node()
        call assert_int(n_total, root%node_idx, 'rebuild root in last slot')
        call assert_int(0, root%parent_idx,     'rebuild root parent_idx=0')
        call assert_true(.not. t%is_leaf(root%node_idx), 'rebuild root is internal')
        call assert_int(3, t%get_height(), 'rebuild height=3')
        call assert_int(2, root%ref_idx,   'rebuild root medoid=2')

        call t%kill()
    end subroutine test_kill_and_rebuild

end module simple_binary_tree_tester