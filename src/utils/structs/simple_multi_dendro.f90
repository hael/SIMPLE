!@descr: simple multi-dendrogram structure to hold multiple hierarchical clusterings
module simple_multi_dendro
use simple_core_module_api
use simple_srch_sort_loc, only: mask2inds
use simple_hclust,        only: hclust
use simple_binary_tree,   only: binary_tree, bt_node, serialize_tree, deserialize_tree
implicit none

public :: multi_dendro, serialize_multi_dendro, deserialize_multi_dendro
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

   ! -------------------------------------------------------------------
   ! Serialize / Deserialize routines for multi_dendro
   ! -------------------------------------------------------------------
   subroutine serialize_multi_dendro(self, mat, trees_meta, offsets, lengths, map_int, tree_pops_out)
      class(multi_dendro), intent(in)  :: self
      integer,           allocatable, intent(out) :: mat(:,:)        ! (total_nodes,6)
      integer,           allocatable, intent(out) :: trees_meta(:,:) ! (n_trees,3): root,nrefs,exists_flag
      integer,           allocatable, intent(out) :: offsets(:)     ! 1-based start idx in mat, 0 if length=0
      integer,           allocatable, intent(out) :: lengths(:)     ! number rows per tree (can be 0)
      integer,           allocatable, intent(out) :: map_int(:,:)   ! (n_trees, n_refs) 0/1
      integer,           allocatable, intent(out) :: tree_pops_out(:)
      integer :: i, n_trees_loc, n_refs_loc, total_nodes, pos, n_nodes, tmeta(3), tmpmeta(3)
      integer, allocatable :: tmat(:,:), tmpmat(:,:)
      ! initialize outputs
      if (.not. allocated(self%trees) .and. self%n_trees > 0) then
         ! tree array not allocated but n_trees set: treat as empty trees
      end if
      n_trees_loc = self%n_trees
      n_refs_loc  = self%n_refs
      if (n_trees_loc < 0 .or. n_refs_loc < 0) then
         THROW_HARD("serialize_multi_dendro: invalid n_trees/n_refs")
      end if
      allocate(trees_meta(max(0,n_trees_loc),3))
      allocate(offsets(max(0,n_trees_loc)))
      allocate(lengths(max(0,n_trees_loc)))
      allocate(map_int(max(0,n_trees_loc), max(0,n_refs_loc)))
      allocate(tree_pops_out(max(0,n_trees_loc)))
      if (n_trees_loc == 0) then
         ! empty multi-dendro
         allocate(mat(0,6))
         return
      end if
      ! First pass: gather per-tree sizes (do not serialize full trees yet)
      total_nodes = 0
      do i = 1, n_trees_loc
         if (allocated(self%trees)) then
               n_nodes = self%trees(i)%get_n_nodes()
         else
               n_nodes = 0
         end if
         lengths(i) = n_nodes
         total_nodes = total_nodes + n_nodes
         tree_pops_out(i) = 0
         if (allocated(self%tree_pops)) tree_pops_out(i) = self%tree_pops(i)
         ! tree_map
         if (allocated(self%tree_map)) then
               do pos = 1, n_refs_loc
                  map_int(i,pos) = merge(1,0, self%tree_map(i,pos))
               end do
         else
               if (n_refs_loc > 0) map_int(i,1:n_refs_loc) = 0
         end if
      end do
      ! allocate big mat and fill by serializing each tree
      if (total_nodes == 0) then
         allocate(mat(0,6))
      else
         allocate(mat(total_nodes,6))
      end if
      pos = 1
      do i = 1, n_trees_loc
         offsets(i) = 0
         if (lengths(i) > 0) then
              
               call serialize_tree(self%trees(i), tmat, tmeta)
               ! sanity: tmat rows must equal lengths(i)
               if (size(tmat,1) /= lengths(i)) then
                  ! fallback: adjust length to actual serialized rows
                  lengths(i) = size(tmat,1)
               end if
               offsets(i) = pos
               if (lengths(i) > 0) then
                  mat(pos:pos+lengths(i)-1, :) = tmat(:, :)
                  pos = pos + lengths(i)
               end if
               trees_meta(i,1) = tmeta(1)
               trees_meta(i,2) = tmeta(2)
               trees_meta(i,3) = tmeta(3)
               if (allocated(tmat)) deallocate(tmat)
         else
               ! empty tree: metadata still set (from object's saved values if any)
               if (allocated(self%trees)) then
                  call serialize_tree(self%trees(i), tmpmat, tmpmeta)
                  trees_meta(i,1) = tmpmeta(1)
                  trees_meta(i,2) = tmpmeta(2)
                  trees_meta(i,3) = tmpmeta(3)
                  if (allocated(tmpmat)) deallocate(tmpmat)
               else
                  trees_meta(i,1:3) = 0
               end if
               offsets(i) = 0
         end if
      end do
   end subroutine serialize_multi_dendro

   subroutine deserialize_multi_dendro(self, mat, trees_meta, offsets, lengths, map_int, tree_pops_in)
      class(multi_dendro), intent(inout) :: self
      integer,           intent(in) :: mat(:,:)         ! (total_nodes,6) or (0,6)
      integer,           intent(in) :: trees_meta(:,:)  ! (n_trees,3)
      integer,           intent(in) :: offsets(:)       ! start idx per tree (1-based) or 0
      integer,           intent(in) :: lengths(:)       ! rows per tree
      integer,           intent(in) :: map_int(:,:)     ! (n_trees, n_refs) 0/1
      integer,           intent(in) :: tree_pops_in(:)
      integer :: i, n_trees_in, n_refs_in, start, len, total_nodes, tmpmeta(3)
      integer, allocatable :: submat(:,:)
      ! validate dims
      n_trees_in = size(map_int,1)
      n_refs_in  = size(map_int,2)
      if (size(trees_meta,1) /= n_trees_in .or. size(offsets) /= n_trees_in .or. size(lengths) /= n_trees_in .or. size(tree_pops_in) /= n_trees_in) then
         THROW_HARD("deserialize_multi_dendro: inconsistent input dimensions")
      end if
      ! clear existing
      call self%kill()
      ! set counts and allocate
      self%n_trees = n_trees_in
      self%n_refs  = n_refs_in
      if (n_trees_in > 0) then
         allocate(self%trees(n_trees_in))
         allocate(self%tree_map(n_trees_in, n_refs_in))
         allocate(self%tree_pops(n_trees_in))
      end if
      ! set tree_map and tree_pops
      do i = 1, n_trees_in
         self%tree_pops(i) = tree_pops_in(i)
         do start = 1, n_refs_in
               self%tree_map(i,start) = merge(.true., .false., map_int(i,start) /= 0)
         end do
      end do
      total_nodes = size(mat,1)
      ! reconstruct each tree from its slice
      do i = 1, n_trees_in
         start = offsets(i)
         len   = lengths(i)
         if (len < 0) then
               THROW_HARD("deserialize_multi_dendro: negative length")
         end if
         if (len == 0) then
               ! leave tree empty (already killed by kill())
               cycle
         end if
         if (start < 1 .or. start+len-1 > total_nodes) then
               THROW_HARD("deserialize_multi_dendro: offsets/lengths out of range")
         end if
         ! extract submatrix (Fortran creates a view)
         
         allocate(submat(len, size(mat,2)))
         submat(:, :) = mat(start:start+len-1, :)
         tmpmeta(:) = trees_meta(i, 1:3)
         call deserialize_tree(self%trees(i), submat, tmpmeta)
         deallocate(submat)
      end do
   end subroutine deserialize_multi_dendro

end module simple_multi_dendro