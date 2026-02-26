!@descr: simple multi-dendrogram structure to hold multiple hierarchical clusterings
module simple_multi_dendro
use simple_core_module_api
use simple_srch_sort_loc, only: mask2inds
use simple_hclust,        only: hclust
use simple_binary_tree,   only: binary_tree, bt_node
implicit none

public :: multi_dendro
private
#include "simple_local_flags.inc"

type :: multi_dendro
   private
   type(binary_tree), allocatable :: trees(:)
   integer,           allocatable :: tree_pops(:)
   logical,           allocatable :: tree_map(:,:)
   integer :: n_trees = 0
   integer :: n_refs  = 0
contains
   procedure :: new
   procedure :: get_n_trees
   procedure :: get_n_nodes
   procedure :: get_n_refs
   procedure :: get_tree_pop
   procedure :: get_tree_refs
   procedure :: get_tree_height
   procedure :: build_tree_from_subdistmat
   procedure :: get_root_node           ! returns bt_node
   procedure :: get_node                ! returns bt_node
   procedure :: is_leaf
   procedure :: kill
end type multi_dendro

contains

   subroutine new(self, labels)
      class(multi_dendro), intent(inout) :: self
      integer,             intent(in)    :: labels(:)
      integer :: i, j
      call self%kill()
      if (any(labels == 0)) THROW_HARD('0 labels not allowed')
      self%n_trees = maxval(labels)
      self%n_refs  = size(labels)
      allocate(self%tree_pops(self%n_trees), source=0)
      do i = 1, self%n_trees
         self%tree_pops(i) = count(labels == i)
      end do
      allocate(self%tree_map(self%n_trees, self%n_refs))
      do i = 1, self%n_trees
         do j = 1, self%n_refs
            self%tree_map(i,j) = (labels(j) == i)
         end do
      end do
      allocate(self%trees(self%n_trees))
   end subroutine new

   pure integer function get_n_trees(self) result(n)
      class(multi_dendro), intent(in) :: self
      n = self%n_trees
   end function get_n_trees

   pure integer function get_n_nodes(self, itree) result(n)
      class(multi_dendro), intent(in) :: self
      integer,             intent(in) :: itree
      n = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%trees)) return
      n = self%trees(itree)%get_n_nodes()
   end function get_n_nodes

   pure integer function get_n_refs(self) result(n)
      class(multi_dendro), intent(in) :: self
      n = self%n_refs
   end function get_n_refs

   pure integer function get_tree_pop(self, itree) result(pop)
      class(multi_dendro), intent(in) :: self
      integer, intent(in) :: itree
      pop = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%tree_pops)) return
      pop = self%tree_pops(itree)
   end function get_tree_pop

   !--- Return refs for a given tree (GLOBAL ref IDs).
   function get_tree_refs(self, itree) result(refs)
      class(multi_dendro), intent(in) :: self
      integer,             intent(in) :: itree
      integer, allocatable :: refs(:)
      integer :: n
      if (itree < 1 .or. itree > self%n_trees) then
         allocate(refs(0))
         return
      end if
      n = count(self%tree_map(itree,:))
      allocate(refs(n))
      refs = mask2inds(self%tree_map(itree,:))
   end function get_tree_refs

   pure integer function get_tree_height(self, itree) result(h)
      class(multi_dendro), intent(in) :: self
      integer,             intent(in) :: itree
      h = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%trees)) return
      h = self%trees(itree)%get_height()
   end function get_tree_height

   !--- Build a single tree from an externally provided sub-distance matrix.
   subroutine build_tree_from_subdistmat(self, itree, refs, sub_distmat, linkage)
      class(multi_dendro), intent(inout) :: self
      integer,             intent(in)    :: itree
      integer,             intent(in)    :: refs(:)            ! GLOBAL ref ids (tree members)
      real,                intent(in)    :: sub_distmat(:,:)   ! (n,n) LOCAL to refs ordering
      integer,             intent(in)    :: linkage
      integer, allocatable :: merge_mat(:,:)
      real,    allocatable :: height(:)
      type(hclust) :: hc
      integer :: n
      if (itree < 1 .or. itree > self%n_trees) then
         THROW_HARD('build_tree_from_subdistmat: itree out of range')
      end if
      n = size(refs)
      if (n == 0) THROW_HARD('build_tree_from_subdistmat: empty refs')
      if (size(sub_distmat,1) /= n .or. size(sub_distmat,2) /= n) then
         THROW_HARD('build_tree_from_subdistmat: sub_distmat must be (n,n) with n=size(refs)')
      end if
      call self%trees(itree)%kill()
      if (n == 1) then
         ! For n=1, merge_mat must be (2, n-1) = (2,0)
         allocate(merge_mat(2,0))
         call self%trees(itree)%build_from_hclust(merge_mat, refs, sub_distmat)
         deallocate(merge_mat)
         return
      end if
      allocate(merge_mat(2, n-1), height(n-1))
      call hc%new(n, sub_distmat, linkage)
      call hc%cluster(merge_mat, height)
      call hc%kill()
      call self%trees(itree)%build_from_hclust(merge_mat, refs, sub_distmat)
      deallocate(merge_mat, height)
   end subroutine build_tree_from_subdistmat

   pure function get_root_node(self, itree) result(root_node)
      class(multi_dendro), intent(in) :: self
      integer,            intent(in) :: itree
      type(bt_node) :: root_node
      root_node%left_idx   = 0
      root_node%right_idx  = 0
      root_node%parent_idx = 0
      root_node%level      = 0
      root_node%ref_idx    = 0
      root_node%node_idx   = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%trees)) return
      root_node = self%trees(itree)%get_root_node()
   end function get_root_node

   pure function get_node(self, itree, inode) result(node)
      class(multi_dendro), intent(in) :: self
      integer,             intent(in) :: itree, inode
      type(bt_node) :: node
      node%left_idx   = 0
      node%right_idx  = 0
      node%parent_idx = 0
      node%level      = 0
      node%ref_idx    = 0
      node%node_idx   = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%trees)) return
      node = self%trees(itree)%get_node(inode)
   end function get_node

   pure logical function is_leaf(self, itree, inode) result(leaf)
      class(multi_dendro), intent(in) :: self
      integer,            intent(in) :: itree, inode
      leaf = .false.
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%trees)) return
      leaf = self%trees(itree)%is_leaf(inode)
   end function is_leaf

   subroutine kill(self)
      class(multi_dendro), intent(inout) :: self
      integer :: i
      if (allocated(self%trees)) then
         do i = 1, size(self%trees)
            call self%trees(i)%kill()
         end do
         deallocate(self%trees)
      end if
      if (allocated(self%tree_pops))  deallocate(self%tree_pops)
      if (allocated(self%tree_map))   deallocate(self%tree_map)
      self%n_trees = 0
      self%n_refs  = 0
   end subroutine kill

end module simple_multi_dendro