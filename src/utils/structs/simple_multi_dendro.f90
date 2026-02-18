!@descr: simple multi-dendrogram structure to hold multiple hierarchical clusterings of different reference chunks (from coarse clustering step)
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
   integer,           allocatable :: cls_pops(:)
   logical,           allocatable :: cls_map(:,:)
   integer,           allocatable :: medoids(:)
   integer :: n_trees = 0
   integer :: n_refs  = 0
   logical :: exists  = .false.
contains
   procedure :: new
   procedure :: get_n_trees
   procedure :: get_n_refs
   procedure :: get_medoid
   procedure :: get_cls_pop
   procedure :: get_tree_refs
   procedure :: build_tree_from_subdist
   procedure :: get_left_right_idxs
   procedure :: get_tree_indx
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
      allocate(self%cls_pops(self%n_trees), self%medoids(self%n_trees), source=0)
      do i = 1, self%n_trees
         self%cls_pops(i) = count(labels == i)
      end do
      allocate(self%cls_map(self%n_trees,self%n_refs))
      do i = 1, self%n_trees
         do j = 1, self%n_refs
            self%cls_map(i,j) = (labels(j) == i)
         end do
      end do
      allocate(self%trees(self%n_trees))
      self%exists = .true.
   end subroutine new

   pure integer function get_n_trees(self) result(n)
      class(multi_dendro), intent(in) :: self
      n = self%n_trees
   end function

   pure integer function get_n_refs(self) result(n)
      class(multi_dendro), intent(in) :: self
      n = self%n_refs
   end function

   pure integer function get_medoid(self, itree) result(m)
      class(multi_dendro), intent(in) :: self
      integer, intent(in) :: itree
      m = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%medoids)) return
      m = self%medoids(itree)
   end function

   pure integer function get_cls_pop(self, itree) result(pop)
      class(multi_dendro), intent(in) :: self
      integer, intent(in) :: itree
      pop = 0
      if (itree < 1 .or. itree > self%n_trees) return
      if (.not. allocated(self%cls_pops)) return
      pop = self%cls_pops(itree)
   end function get_cls_pop

   !--- Return refs for a given tree (global ref IDs).
   function get_tree_refs(self, itree) result(refs)
        class(multi_dendro), intent(in) :: self
        integer,             intent(in) :: itree
        integer, allocatable :: refs(:)
        integer :: n
        if (itree < 1 .or. itree > self%n_trees) then
            allocate(refs(0))
            return
        end if
        n = count(self%cls_map(itree,:))
        allocate(refs(n))
        refs = mask2inds(self%cls_map(itree,:))
   end function get_tree_refs

   !--- Build a single tree from an externally provided sub-distance matrix.
   !    refs must be the same list returned by get_tree_refs(itree, ...)
   !    sub_distmat must be n x n, symmetric (upper triangle accepted).
   subroutine build_tree_from_subdist(self, itree, refs, sub_distmat, linkage)
      class(multi_dendro), intent(inout) :: self
      integer,            intent(in)    :: itree
      integer,            intent(in)    :: refs(:)            ! global ref ids
      real,               intent(in)    :: sub_distmat(:,:)   ! n x n symmetric
      integer,            intent(in)    :: linkage
      ! local vars
      integer :: n, i, j
      integer, allocatable :: merge_mat(:,:)
      real,    allocatable :: height(:)
      type(hclust)  :: hc
      type(bt_node) :: node
      ! Basic validations
      if (itree < 1 .or. itree > self%n_trees) then
         THROW_HARD('build_tree_from_subdist: itree out of range')
      end if
      n = size(refs)
      if (n == 0) then
         THROW_HARD('build_tree_from_subdist: empty refs')
      end if
      if (size(sub_distmat,1) /= n .or. size(sub_distmat,2) /= n) then
         THROW_HARD('build_tree_from_subdist: sub_distmat must be (n,n) with n=size(refs)')
      end if
      ! If the tree already exists, clear it first
      if (allocated(self%trees)) then
         call self%trees(itree)%kill()
      end if
      ! For n==1: create degenerate single-node tree
      if (n == 1) then
         ! Create a one-node tree (build_from_hclust supports this as before)
         allocate(merge_mat(2,1))  ! dummy but not used by build logic for n=1
         merge_mat = 1   ! value doesn't matter
         call self%trees(itree)%build_from_hclust(merge_mat, refs, sub_distmat)
         self%medoids(itree) = refs(1)
         deallocate(merge_mat)
         return
      end if
      ! Run hclust on the provided local submatrix
      allocate(merge_mat(2, n-1), height(n-1))
      call hc%new(n, sub_distmat, linkage)
      call hc%cluster(merge_mat, height)
      call hc%kill()
      ! Build binary_tree from merge matrix using local submatrix (local builder uses sub_distmat)
      call self%trees(itree)%build_from_hclust(merge_mat, refs, sub_distmat)
      ! set medoid for this tree
      node = self%trees(itree)%get_node(self%trees(itree)%get_root_idx())
      self%medoids(itree) = node%ref_idx
      deallocate(merge_mat, height)
   end subroutine build_tree_from_subdist

   pure subroutine get_left_right_idxs(self, ref_idx, left_ref_idx, right_ref_idx)
      class(multi_dendro), intent(in)  :: self
      integer,             intent(in)  :: ref_idx
      integer,             intent(out) :: left_ref_idx, right_ref_idx
      integer :: itree
      left_ref_idx  = 0
      right_ref_idx = 0
      itree = self%get_tree_indx(ref_idx)
      if (itree == 0) return
      call self%trees(itree)%get_left_right_ref(ref_idx, left_ref_idx, right_ref_idx)
   end subroutine get_left_right_idxs

   pure integer function get_tree_indx(self, ref_idx) result(tree_idx)
      class(multi_dendro), intent(in) :: self
      integer,             intent(in) :: ref_idx
      integer :: itree
      tree_idx = 0
      if (self%n_trees == 0) return ! should never happen, but just in case
      if (self%n_trees == 1) then
         if (self%trees(1)%find_node_by_ref(ref_idx) /= 0) tree_idx = 1
         return
      end if
      do itree = 1, self%n_trees
         if (self%trees(itree)%find_node_by_ref(ref_idx) /= 0) then
            tree_idx = itree
            return
         end if
      end do
   end function get_tree_indx

   subroutine kill(self)
      class(multi_dendro), intent(inout) :: self
      integer :: i
      if (.not. self%exists) return
      if (allocated(self%trees)) then
         do i = 1, size(self%trees)
            call self%trees(i)%kill()
         end do
         deallocate(self%trees)
      end if
      if (allocated(self%cls_pops))  deallocate(self%cls_pops)
      if (allocated(self%cls_map))   deallocate(self%cls_map)
      if (allocated(self%medoids))   deallocate(self%medoids)
      self%n_trees = 0
      self%n_refs  = 0
      self%exists  = .false.
   end subroutine kill

end module simple_multi_dendro
