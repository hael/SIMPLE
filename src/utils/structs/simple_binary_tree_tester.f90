!@descr: unit test routines for simple_binary_tree class
module simple_binary_tree_tester
use simple_binary_tree
use simple_test_utils
implicit none

public :: run_all_tree_tests
private

!-----------------------------------
! Visitor capture utilities
!-----------------------------------
integer, allocatable, save :: visited_nodes(:)
integer, save :: nvisited = 0

contains

    subroutine run_all_tree_tests()
        write(*,'(A)') '**** running binary_tree tests ****'
        call test_initial_state()
        call test_empty_traversals_noop()
        call test_build_from_hclust_basic_structure()
        call test_getters_and_copy_semantics()
        call test_new_helpers_and_mapping()
        call test_traversal_orders()
        call test_kill_and_rebuild()
        write(*,'(A)') '**** finished binary_tree tests ****'
    end subroutine run_all_tree_tests

    !===========================
    ! Helpers
    !===========================

    subroutine reset_visit_buffer(nmax)
        integer, intent(in) :: nmax
        if (allocated(visited_nodes)) deallocate(visited_nodes)
        allocate(visited_nodes(nmax))
        visited_nodes = 0
        nvisited = 0
    end subroutine reset_visit_buffer

    subroutine record_visit(node)
        type(bt_node), intent(in) :: node
        nvisited = nvisited + 1
        if (.not. allocated(visited_nodes)) stop "visit buffer not allocated"
        if (nvisited > size(visited_nodes)) stop "visit buffer overflow"
        visited_nodes(nvisited) = node%node_idx
    end subroutine record_visit

    subroutine assert_visit_sequence(expected, msg)
        integer, intent(in) :: expected(:)
        character(*), intent(in) :: msg
        integer :: i
        call assert_int(size(expected), nvisited, trim(msg)//' (count)')
        do i = 1, size(expected)
            call assert_int(expected(i), visited_nodes(i), trim(msg)//' (seq)')
        end do
    end subroutine assert_visit_sequence

    subroutine build_test_tree(t)
        type(binary_tree), intent(inout) :: t
        integer, allocatable :: merge_mat(:,:)
        integer :: refs(4)
        real, allocatable :: dist(:,:)
        integer :: i

        ! Leaves: 1..4, internal: 5..7
        ! merges:
        !   (1,2)->5, (3,4)->6, (5,6)->7
        allocate(merge_mat(2,3))
        merge_mat(1,:) = [1, 3, 5]
        merge_mat(2,:) = [2, 4, 6]

        refs = [1,2,3,4]

        allocate(dist(4,4))
        dist = 0.0

        ! pair (1,2): tie => choose local 1 => global 1
        ! pair (3,4): tie => choose local 3 => global 3
        ! root (1,2,3,4): best is local 2 => global 2
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

    !===========================
    ! Tests
    !===========================

    subroutine test_initial_state()
        type(binary_tree) :: t
        write(*,'(A)') 'test_initial_state'
        call assert_int(0, t%n_nodes(), 'new tree has 0 nodes')
        call assert_int(0, t%get_root_idx(), 'new tree root_idx=0')
        call assert_int(0, t%get_nref(), 'new tree nref=0')
        call assert_int(0, t%get_medoid(), 'new tree medoid=0')
    end subroutine test_initial_state

    subroutine test_empty_traversals_noop()
        type(binary_tree) :: t
        write(*,'(A)') 'test_empty_traversals_noop'
        call reset_visit_buffer(10)
        call t%traverse_preorder(record_visit)
        call assert_int(0, nvisited, 'preorder on empty visits 0')
        call t%traverse_inorder(record_visit)
        call assert_int(0, nvisited, 'inorder on empty visits 0')
        call t%traverse_postorder(record_visit)
        call assert_int(0, nvisited, 'postorder on empty visits 0')
        call t%traverse_levelorder(record_visit)
        call assert_int(0, nvisited, 'levelorder on empty visits 0')
    end subroutine test_empty_traversals_noop

    subroutine test_build_from_hclust_basic_structure()
        type(binary_tree) :: t
        type(bt_node) :: n
        integer :: lidx, ridx
        integer :: lref, rref

        write(*,'(A)') 'test_build_from_hclust_basic_structure'
        call build_test_tree(t)

        call assert_int(7, t%n_nodes(), 'n_nodes = 2*nref-1')
        call assert_int(7, t%get_root_idx(), 'root is last node')
        call assert_int(4, t%get_nref(), 'nref stored in tree')

        ! Leaf nodes 1..4
        n = t%get_node(1)
        call assert_int(1, n%node_idx, 'leaf1 node_idx')
        call assert_int(1, n%ref_idx,  'leaf1 ref (global)')
        call assert_int(0, n%left_idx, 'leaf1 left=0')
        call assert_int(0, n%right_idx,'leaf1 right=0')
        call assert_int(5, n%parent_idx,'leaf1 parent=5')
        call assert_int(1, size(n%subset), 'leaf1 subset size')
        call assert_int(1, n%subset(1), 'leaf1 subset local=1')

        n = t%get_node(2)
        call assert_int(5, n%parent_idx,'leaf2 parent=5')
        call assert_int(2, n%subset(1), 'leaf2 subset local=2')

        n = t%get_node(3)
        call assert_int(6, n%parent_idx,'leaf3 parent=6')
        call assert_int(3, n%subset(1), 'leaf3 subset local=3')

        n = t%get_node(4)
        call assert_int(6, n%parent_idx,'leaf4 parent=6')
        call assert_int(4, n%subset(1), 'leaf4 subset local=4')

        ! Internal node 5 merges (1,2)
        n = t%get_node(5)
        call assert_int(1, n%left_idx,  'node5 left')
        call assert_int(2, n%right_idx, 'node5 right')
        call assert_int(7, n%parent_idx, 'node5 parent=7')
        call assert_int(1, n%level, 'node5 level=1')
        call assert_int(2, size(n%subset), 'node5 subset size=2')
        call assert_int(1, n%subset(1), 'node5 subset local(1)=1')
        call assert_int(2, n%subset(2), 'node5 subset local(2)=2')
        call assert_int(1, n%ref_idx, 'node5 medoid global ref=1')

        ! Internal node 6 merges (3,4)
        n = t%get_node(6)
        call assert_int(3, n%left_idx,  'node6 left')
        call assert_int(4, n%right_idx, 'node6 right')
        call assert_int(7, n%parent_idx,'node6 parent=7')
        call assert_int(2, n%level, 'node6 level=2')
        call assert_int(2, size(n%subset), 'node6 subset size=2')
        call assert_int(3, n%subset(1), 'node6 subset local(1)=3')
        call assert_int(4, n%subset(2), 'node6 subset local(2)=4')
        call assert_int(3, n%ref_idx, 'node6 medoid global ref=3')

        ! Root node 7 merges (5,6)
        n = t%get_node(7)
        call assert_int(5, n%left_idx,  'root left=5')
        call assert_int(6, n%right_idx, 'root right=6')
        call assert_int(0, n%parent_idx,'root parent=0')
        call assert_int(3, n%level, 'root level=3')
        call assert_int(4, size(n%subset), 'root subset size=4')
        call assert_int(1, n%subset(1), 'root subset local(1)=1')
        call assert_int(2, n%subset(2), 'root subset local(2)=2')
        call assert_int(3, n%subset(3), 'root subset local(3)=3')
        call assert_int(4, n%subset(4), 'root subset local(4)=4')
        call assert_int(2, n%ref_idx, 'root medoid global ref=2')
        call assert_int(2, t%get_medoid(), 'tree get_medoid() matches root medoid')

        ! Sanity-check new child helpers on the root
        call t%get_children_idx(t%get_root_idx(), lidx, ridx)
        call assert_int(5, lidx, 'get_children_idx(root) left=5')
        call assert_int(6, ridx, 'get_children_idx(root) right=6')

        call t%get_children_ref(t%get_root_idx(), lref, rref)
        call assert_int(1, lref, 'get_children_ref(root) left medoid=1')
        call assert_int(3, rref, 'get_children_ref(root) right medoid=3')

        call t%kill()
    end subroutine test_build_from_hclust_basic_structure

    subroutine test_getters_and_copy_semantics()
        type(binary_tree) :: t
        type(bt_node) :: a, b
        write(*,'(A)') 'test_getters_and_copy_semantics'
        call build_test_tree(t)

        call assert_int(7, t%get_root_idx(), 'get_root_idx ok')
        a = t%get_node(1)
        call assert_int(1, a%ref_idx, 'get_node returns correct ref')

        ! Ensure get_node returns a copy (mutating local should not affect tree)
        a%ref_idx = 999
        b = t%get_node(1)
        call assert_int(1, b%ref_idx, 'tree unaffected by local mutation')

        call t%kill()
    end subroutine test_getters_and_copy_semantics

    subroutine test_new_helpers_and_mapping()
        type(binary_tree) :: t
        integer :: root, lidx, ridx
        integer :: g

        write(*,'(A)') 'test_new_helpers_and_mapping'
        call build_test_tree(t)

        call assert_int(4, t%get_nref(), 'get_nref=4')

        ! local_to_global_ref is identity in this test because refs=[1,2,3,4]
        g = t%local_to_global_ref(1); call assert_int(1, g, 'local_to_global_ref(1)=1')
        g = t%local_to_global_ref(4); call assert_int(4, g, 'local_to_global_ref(4)=4')
        g = t%local_to_global_ref(0); call assert_int(0, g, 'local_to_global_ref(0)=0 out-of-range')
        g = t%local_to_global_ref(999); call assert_int(0, g, 'local_to_global_ref(999)=0 out-of-range')

        root = t%get_root_idx()
        call assert_int(0, merge(1,0,t%is_leaf(root)), 'root is not leaf')  ! avoid assert_logical dependency
        call assert_int(4, t%get_subset_size(root), 'root subset size=4')
        call assert_int(2, t%get_node_ref(root), 'get_node_ref(root)=2')

        call t%get_children_idx(root, lidx, ridx)
        call assert_int(5, lidx, 'children_idx(root) left=5')
        call assert_int(6, ridx, 'children_idx(root) right=6')

        call assert_int(2, t%get_subset_size(5), 'subset size(node5)=2')
        call assert_int(1, t%get_node_ref(5), 'node_ref(node5)=1')

        call assert_int(1, merge(1,0,t%is_leaf(1)), 'leaf1 is leaf')
        call assert_int(1, t%get_subset_size(1), 'leaf1 subset size=1')

        call t%kill()
    end subroutine test_new_helpers_and_mapping

    subroutine test_traversal_orders()
        type(binary_tree) :: t
        integer, allocatable :: seq(:)
        write(*,'(A)') 'test_traversal_orders'
        call build_test_tree(t)

        ! Preorder: root, left, right
        call reset_visit_buffer(t%n_nodes())
        call t%traverse_preorder(record_visit)
        seq = [7,5,1,2,6,3,4]
        call assert_visit_sequence(seq, 'preorder node_idx order')

        ! Inorder: left, root, right
        call reset_visit_buffer(t%n_nodes())
        call t%traverse_inorder(record_visit)
        seq = [1,5,2,7,3,6,4]
        call assert_visit_sequence(seq, 'inorder node_idx order')

        ! Postorder: left, right, root
        call reset_visit_buffer(t%n_nodes())
        call t%traverse_postorder(record_visit)
        seq = [1,2,5,3,4,6,7]
        call assert_visit_sequence(seq, 'postorder node_idx order')

        ! Levelorder: BFS
        call reset_visit_buffer(t%n_nodes())
        call t%traverse_levelorder(record_visit)
        seq = [7,5,6,1,2,3,4]
        call assert_visit_sequence(seq, 'levelorder node_idx order')

        call t%kill()
    end subroutine test_traversal_orders

    subroutine test_kill_and_rebuild()
        type(binary_tree) :: t
        write(*,'(A)') 'test_kill_and_rebuild'
        call build_test_tree(t)

        call assert_int(7, t%n_nodes(), 'built tree has nodes')
        call assert_int(4, t%get_nref(), 'built tree has nref')

        call t%kill()
        call assert_int(0, t%n_nodes(), 'kill deallocates nodes')
        call assert_int(0, t%get_root_idx(), 'kill resets root')
        call assert_int(0, t%get_nref(), 'kill resets nref')

        ! Rebuild after kill should work
        call build_test_tree(t)
        call assert_int(7, t%n_nodes(), 'rebuild after kill works')
        call assert_int(7, t%get_root_idx(), 'rebuild root ok')
        call assert_int(4, t%get_nref(), 'rebuild nref ok')

        call t%kill()
    end subroutine test_kill_and_rebuild

end module simple_binary_tree_tester