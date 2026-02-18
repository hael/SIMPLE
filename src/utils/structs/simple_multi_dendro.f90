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
   real,              allocatable :: dist_mat(:,:)
   integer,           allocatable :: cls_pops(:)
   logical,           allocatable :: cls_map(:,:)
   integer,           allocatable :: medoids(:)
   integer :: n_trees = 0
   integer :: n_refs  = 0
   logical :: exists  = .false.
contains
   procedure          :: new
   procedure          :: get_n_trees
   procedure          :: get_n_refs
   procedure          :: get_medoid
   procedure          :: get_cls_pop
   procedure          :: build_multi_dendro
   procedure          :: get_left_right_idxs
   procedure, private :: get_tree_indx
   procedure          :: kill
end type multi_dendro

contains

   subroutine new(self, dist_mat, labels)
      class(multi_dendro), intent(inout) :: self
      real,                intent(in)    :: dist_mat(:,:)
      integer,             intent(in)    :: labels(:)
      integer :: i, j
      call self%kill()
      if (any(labels == 0)) THROW_HARD('0 labels not allowed')
      self%n_trees = maxval(labels)
      self%n_refs  = size(labels)
      allocate(self%dist_mat(self%n_refs,self%n_refs), source=dist_mat)
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

   subroutine build_multi_dendro(self, linkage)
      class(multi_dendro), intent(inout) :: self
      integer,             intent(in)    :: linkage
      real,    allocatable :: sub_distmat(:,:), height(:)
      integer, allocatable :: refs(:), merge_mat(:,:)
      integer :: itree, i, j, nref_sub, root_idx
      type(hclust)  :: hc
      type(bt_node) :: node
      do itree = 1, self%n_trees
         nref_sub = count(self%cls_map(itree,:))
         refs     = mask2inds(self%cls_map(itree,:))
         if (nref_sub == 1) then
            ! degenerate tree: single node
            call self%trees(itree)%build_from_hclust(reshape([1,1],[2,1]), refs, self%dist_mat)
            self%medoids(itree) = refs(1)
            deallocate(refs)
            cycle
         end if
         allocate(sub_distmat(nref_sub, nref_sub))
         sub_distmat = 1.0
         do i = 1, nref_sub
            do j = i + 1, nref_sub
               sub_distmat(i,j) = self%dist_mat(refs(i), refs(j))
               sub_distmat(j,i) = sub_distmat(i,j)
            end do
         end do
         allocate(merge_mat(2, nref_sub-1), height(nref_sub-1))
         call hc%new(nref_sub, sub_distmat, linkage)
         call hc%cluster(merge_mat, height)
         call hc%kill()
         call self%trees(itree)%build_from_hclust(merge_mat, refs, self%dist_mat)
         root_idx = self%trees(itree)%get_root_idx()
         node     = self%trees(itree)%get_node(root_idx)
         self%medoids(itree) = node%ref_idx
         deallocate(sub_distmat, refs, merge_mat, height)
      end do
   end subroutine build_multi_dendro

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
      if (allocated(self%dist_mat))  deallocate(self%dist_mat)
      if (allocated(self%cls_pops))  deallocate(self%cls_pops)
      if (allocated(self%cls_map))   deallocate(self%cls_map)
      if (allocated(self%medoids))   deallocate(self%medoids)
      self%n_trees = 0
      self%n_refs  = 0
      self%exists  = .false.
   end subroutine kill

end module simple_multi_dendro
