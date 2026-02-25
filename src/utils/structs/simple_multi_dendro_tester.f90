!@descr: unit test routines for simple_multi_dendro class (matches current minimal API)
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
      call test_build_tree_singletons()
      call test_build_tree_pairs_basic_structure()
      call test_build_tree_two_triples_basic_structure()
      call test_kill_and_reuse()
      write(*,'(A)') '**** finished multi_dendro tests ****'
   end subroutine run_all_multi_dendro_tests

   !===========================
   ! Helpers
   !===========================

   subroutine make_distmat_two_triples(dist)
      ! Two clusters: {1,2,3} and {4,5,6}. Cross distances large.
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

   subroutine assert_int_in_set(x, setv, msg)
      integer, intent(in) :: x
      integer, intent(in) :: setv(:)
      character(*), intent(in) :: msg
      call assert_true(any(setv == x), msg)
   end subroutine assert_int_in_set

   subroutine build_tree_from_full_dist(md, itree, dist_full)
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
         do i = 1, n
            do j = i+1, n
               submat(i,j) = dist_full(refs(i), refs(j))
               submat(j,i) = submat(i,j)
            end do
         end do
      end if

      call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)

      deallocate(submat)
      deallocate(refs)
   end subroutine build_tree_from_full_dist

   !===========================
   ! Tests
   !===========================

   subroutine test_initial_state()
      type(multi_dendro) :: md
      type(bt_node) :: root
      write(*,'(A)') 'test_initial_state'
      call assert_int(0, md%get_n_trees(), 'initial get_n_trees=0')
      call assert_int(0, md%get_n_refs(),  'initial get_n_refs=0')
      call assert_int(0, md%get_tree_pop(1),'initial get_tree_pop(1)=0')
      call assert_int(0, md%get_tree_height(1), 'initial get_tree_height(1)=0 (out-of-range safe)')
      root = md%get_root_node(1)
      call assert_int(0, root%node_idx, 'initial get_root_node(1) returns zeroed node')
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
      type(bt_node) :: root
      write(*,'(A)') 'test_build_tree_singletons'

      labels = [1,2,3,4]
      call md%new(labels)

      do itree = 1, md%get_n_trees()
         refs = md%get_tree_refs(itree)
         call assert_int(1, size(refs), 'singleton refs size=1')

         allocate(submat(1,1))
         submat = 0.0
         call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)
         deallocate(submat)

         call assert_int(1, md%get_tree_pop(itree), 'singleton tree_pop=1')
         call assert_int(1, md%get_tree_height(itree), 'singleton height=1')

         root = md%get_root_node(itree)
         call assert_int(1, root%node_idx, 'singleton root node_idx=1')
         call assert_true(md%is_leaf(itree, root%node_idx), 'singleton root is leaf')
         call assert_int(refs(1), root%ref_idx, 'singleton root ref == only ref')
         call assert_int(0, root%left_idx,  'singleton root left=0')
         call assert_int(0, root%right_idx, 'singleton root right=0')
         call assert_int(0, root%parent_idx,'singleton root parent=0')

         deallocate(refs)
      end do

      call md%kill()
   end subroutine test_build_tree_singletons

   subroutine test_build_tree_pairs_basic_structure()
      type(multi_dendro) :: md
      integer :: labels(4)
      real :: dist(4,4)
      integer :: itree
      integer, allocatable :: refs(:)
      real, allocatable :: submat(:,:)
      type(bt_node) :: root, lch, rch
      integer :: n

      write(*,'(A)') 'test_build_tree_pairs_basic_structure'
      call make_distmat_two_pairs(dist)
      labels = [1,1,2,2]   ! {1,2} and {3,4}
      call md%new(labels)

      do itree = 1, md%get_n_trees()
         refs = md%get_tree_refs(itree)
         n = size(refs)
         call assert_int(2, n, 'pair refs size=2')

         allocate(submat(n,n))
         submat = 0.0
         submat(1,2) = dist(refs(1), refs(2))
         submat(2,1) = submat(1,2)

         call md%build_tree_from_subdistmat(itree, refs, submat, linkage=1)
         deallocate(submat)

         call assert_int(2, md%get_tree_pop(itree), 'pair tree_pop=2')
         call assert_int(2, md%get_tree_height(itree), 'pair height=2')

         root = md%get_root_node(itree)
         call assert_int(3, root%node_idx, 'pair root node_idx = 2*n-1 = 3')
         call assert_true(.not. md%is_leaf(itree, root%node_idx), 'pair root is not leaf')
         call assert_int(0, root%parent_idx, 'pair root parent=0')
         call assert_true(root%left_idx  /= 0, 'pair root left_idx nonzero')
         call assert_true(root%right_idx /= 0, 'pair root right_idx nonzero')

         lch = md%get_node(itree, root%left_idx)
         rch = md%get_node(itree, root%right_idx)
         call assert_true(md%is_leaf(itree, lch%node_idx), 'pair left child is leaf')
         call assert_true(md%is_leaf(itree, rch%node_idx), 'pair right child is leaf')

         ! leaves carry the physical/global ref ids from refs(:)
         call assert_true(lch%ref_idx == refs(1) .or. lch%ref_idx == refs(2), 'pair left leaf ref in refs')
         call assert_true(rch%ref_idx == refs(1) .or. rch%ref_idx == refs(2), 'pair right leaf ref in refs')
         call assert_true(lch%ref_idx /= rch%ref_idx, 'pair leaves have distinct refs')

         deallocate(refs)
      end do

      call md%kill()
   end subroutine test_build_tree_pairs_basic_structure

   subroutine test_build_tree_two_triples_basic_structure()
      type(multi_dendro) :: md
      integer :: labels(6)
      real :: dist(6,6)
      integer :: itree
      type(bt_node) :: root
      integer :: n

      write(*,'(A)') 'test_build_tree_two_triples_basic_structure'
      call make_distmat_two_triples(dist)
      labels = [1,1,1, 2,2,2]
      call md%new(labels)

      do itree = 1, md%get_n_trees()
         call build_tree_from_full_dist(md, itree, dist)
         n = md%get_tree_pop(itree)
         call assert_int(3, n, 'triple tree_pop=3')
         call assert_int(3, md%get_tree_height(itree), 'triple height=3 (3 leaves)')
         root = md%get_root_node(itree)
         call assert_int(5, root%node_idx, 'triple root node_idx=2*n-1=5')
         call assert_true(.not. md%is_leaf(itree, root%node_idx), 'triple root is not leaf')
         call assert_int(0, root%parent_idx, 'triple root parent=0')
         call assert_true(root%left_idx  /= 0, 'triple root left_idx nonzero')
         call assert_true(root%right_idx /= 0, 'triple root right_idx nonzero')
      end do

      call md%kill()
   end subroutine test_build_tree_two_triples_basic_structure

   subroutine test_kill_and_reuse()
      type(multi_dendro) :: md
      integer :: labels(6)
      real :: dist(6,6)
      integer :: itree

      write(*,'(A)') 'test_kill_and_reuse'
      call make_distmat_two_triples(dist)
      labels = [1,1,1,2,2,2]
      call md%new(labels)

      do itree = 1, md%get_n_trees()
         call build_tree_from_full_dist(md, itree, dist)
      end do

      call md%kill()
      call assert_int(0, md%get_n_trees(), 'after kill n_trees=0')
      call assert_int(0, md%get_n_refs(),  'after kill n_refs=0')

      ! Reuse with new labels
      labels = [1,2,2,2,2,2]
      call md%new(labels)
      call assert_int(2, md%get_n_trees(), 'reuse n_trees=2')
      call assert_int(6, md%get_n_refs(),  'reuse n_refs=6')
      call assert_int(1, md%get_tree_pop(1),'reuse tree_pop(1)=1')
      call assert_int(5, md%get_tree_pop(2),'reuse tree_pop(2)=5')

      call md%kill()
   end subroutine test_kill_and_reuse

end module simple_multi_dendro_tester
