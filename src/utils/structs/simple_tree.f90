!@descr: binary tree data structure
module simple_tree  
use simple_core_module_api
use simple_srch_sort_loc, only: hpsort, locate, mask2inds
use simple_hclust,        only: hclust
implicit none

public :: multi_dendro
private
#include "simple_local_flags.inc"

type  :: s2_node 
   type(s2_node), pointer     :: left => null(), right => null(), parent => null()
   integer,       allocatable :: subset(:) ! array of refs in node
   integer     :: level      = 0           ! where we are in the dendrogram, 0 is bottom
   logical     :: visit      = .false.     ! have visited this one
   logical     :: is_pop     = .false.     ! is populated 
   integer     :: ref_idx    = 0, node_idx = 0
end type s2_node

type node_storage
   type(s2_node), allocatable :: nodes(:)
   integer :: root_idx = 0
end type node_storage

type  :: multi_dendro
   private
   type(node_storage), allocatable :: node_store(:)
   real,               allocatable :: dist_mat(:,:) ! full_distmat
   integer,            allocatable :: cls_pops(:) 
   logical,            allocatable :: cls_map(:,:)  ! AP cluster mapping
   integer,            allocatable :: medoids(:)   
   integer :: n_trees = 0                           ! number of AP clusters
   integer :: n_refs  = 0                           ! number of AP clustered references
   integer :: imedoid                               ! current medoid being tracked
   logical :: exists  = .false.
contains
   ! constructor
   procedure          :: new
   ! tree builder
   procedure          :: build_multi_dendro
   ! setters 
   procedure          :: set_imedoid
   ! left/right search
   procedure          :: get_left_right_idxs
   ! private accesors
   procedure, private :: get_cls_indx
   ! destructor
   procedure          :: kill
end type multi_dendro

contains

   subroutine new( self, dist_mat, labels )
      class(multi_dendro), intent(inout) :: self
      real,                intent(in)    :: dist_mat(:,:)
      integer,             intent(in)    :: labels(:)
      integer :: i, j
      call self%kill
      if( any(labels == 0 )) THROW_HARD('0 labels not allowed')
      self%n_trees = maxval(labels)
      self%n_refs  = size(labels)
      allocate(self%dist_mat(self%n_refs,self%n_refs), source=dist_mat)
      allocate(self%cls_pops(self%n_trees), source=0)
      do i = 1, self%n_trees
         self%cls_pops(i) = count(labels == i)
      end do
      allocate(self%cls_map(self%n_trees,self%n_refs))
      do i = 1, self%n_trees
         do j = 1,self%n_refs
            self%cls_map(i,j) = labels(j) == i
         end do
      end do
      allocate(self%node_store(self%n_trees))
      self%exists = .true.
   end subroutine new

   subroutine build_multi_dendro(self, linkage)
      class(multi_dendro), intent(inout) :: self
      integer,             intent(in)    :: linkage
      real,    allocatable :: sub_distmat(:,:), height(:)
      integer, allocatable :: refs(:), merge_mat(:,:)
      integer  :: icls, i, j, ncls_ap, nref_sub
      type(hclust) :: hc
      do icls = 1,self%n_trees
         nref_sub = count(self%cls_map(icls,:))
         refs     = mask2inds(self%cls_map(icls,:))
         ! if cluster size is 1, done 
         if(nref_sub == 1) then 
            allocate(self%node_store(icls)%nodes(1))
            self%node_store(icls)%nodes(1)%ref_idx = refs(1)
            allocate(self%node_store(icls)%nodes(1)%subset(1))
            self%node_store(icls)%nodes(1)%subset = [refs(1)]
            self%node_store(icls)%root_idx = 1
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
         ! Aggl. Clustering
         allocate(merge_mat(2, nref_sub-1))
         allocate(height(nref_sub - 1))
         call hc%new(nref_sub, sub_distmat, linkage)
         call hc%cluster(merge_mat, height)
         call hc%kill()
         ! Tracking Merges
         call gen_tree4ap_cluster(merge_mat, height, refs, self%node_store(icls)%nodes, self%node_store(icls)%root_idx)
         deallocate(sub_distmat, refs, merge_mat, height)
      end do

      contains
      
         subroutine gen_tree4ap_cluster(merge_mat, height, refs, nodes, root_idx)
            integer,                            intent(in)    :: merge_mat(:,:)
            real,                               intent(in)    :: height(:)
            integer,                            intent(in)    :: refs(:)
            type(s2_node), allocatable, target, intent(inout) :: nodes(:)
            integer,                            intent(out)   :: root_idx
            integer, allocatable :: tmp(:)
            integer :: p, k, s, l, r, m, best, n_nodes, n 
            real    :: best_sum, sum
            n_nodes = 2*nref_sub - 1
            if( allocated(nodes) ) deallocate(nodes)
            allocate(nodes(n_nodes))
            ! allocating nodes 
            do k = 1, n_nodes
               nullify(nodes(k)%left, nodes(k)%right, nodes(k)%parent)
               nodes(k)%level      = 0
               nodes(k)%visit      = .false.
               nodes(k)%is_pop     = .false.
               nodes(k)%ref_idx    = 0
               nodes(k)%node_idx        = k 
               if (allocated(nodes(k)%subset)) deallocate(nodes(k)%subset)
               ! can set all leaves
               if (k <= nref_sub) then
                  nodes(k)%ref_idx  = refs(k)
                  nodes(k)%is_pop   = .true.
                  nodes(k)%level    = 0 
                  allocate(nodes(k)%subset(1))
                  nodes(k)%subset = [refs(k)]
               end if
            end do
            ! set internal nodes / root. 
            m = size(refs)
            do s = 1, nref_sub - 1
               l = merge_mat(1, s)
               r = merge_mat(2, s)
               p = m + s
               ! assigning pointers 
               nodes(p)%left   => nodes(l)
               nodes(p)%right  => nodes(r)
               nodes(l)%parent => nodes(p)
               nodes(r)%parent => nodes(p)
               nodes(p)%level  =  s
               nodes(p)%is_pop =  .true.
               ! parent subset is union of children 
               allocate(tmp(size(nodes(l)%subset) + size(nodes(r)%subset)))
               tmp = [nodes(l)%subset, nodes(r)%subset]
               call hpsort(tmp)
               call move_alloc(tmp, nodes(p)%subset)
               ! Calculate a new medoid in merged set 
               best = -1
               best_sum = huge(1.0)
               do i = 1, size(nodes(p)%subset)
                  sum = 0.0
                  do j = 1, size(nodes(p)%subset)
                     sum = sum + self%dist_mat(nodes(p)%subset(i), nodes(p)%subset(j))
                  end do
                  if (sum < best_sum) then
                     best_sum = sum 
                     best = nodes(p)%subset(i)
                  end if
               end do
               nodes(p)%ref_idx = best
            end do
            root_idx = n_nodes
         end subroutine gen_tree4ap_cluster

   end subroutine build_multi_dendro

   subroutine set_imedoid(self, ref_idx)
      class(multi_dendro), intent(inout)  :: self
      integer,             intent(in)  :: ref_idx
      integer  :: imedoid 
      self%imedoid = self%get_cls_indx(ref_idx)
   end subroutine set_imedoid
   
   ! getter to return left and right indices
   function get_left_right_idxs(self, node_idx) result(idxs)
      class(multi_dendro), intent(inout) :: self 
      integer,             intent(inout) :: node_idx(2)
      type(s2_node), pointer  :: rootp, refp
      integer :: idxs(2), imedoid, ref_idx
      imedoid = self%imedoid
      if(node_idx(1) == 2*self%cls_pops(imedoid) - 1) then 
         rootp = self%node_store(imedoid)%nodes(self%node_store(imedoid)%root_idx)
      end if  
      node_idx(1) = search_tree4ref(self, rootp, ref_idx)
      refp = self%node_store(imedoid)%nodes(node_idx(1))
      idxs = 0
      if ( associated(refp%left)  ) idxs(1) = refp%left%ref_idx ; node_idx(1)  = refp%left%node_idx
      if ( associated(refp%right) ) idxs(2) = refp%right%ref_idx ; node_idx(2) = refp%right%node_idx
   end function get_left_right_idxs

   ! find the cluster to which an aribitrary reference belongs
   pure function get_cls_indx(self, ref_idx) result(tree_idx)
      class(multi_dendro), intent(in) :: self 
      integer,             intent(in) :: ref_idx
      integer :: tree_idx , nrefs_in_cls, ncls, icls, iref
      ncls = self%n_trees
      if(ncls == 1) then 
         tree_idx = 1
         return 
      end if 
      do icls = 1, ncls
         associate(refinds_in_cls => self%node_store(icls)%nodes(self%node_store(icls)%root_idx)%subset)
            nrefs_in_cls = size(refinds_in_cls)
            do iref = 1, nrefs_in_cls
               if( refinds_in_cls(iref) == ref_idx )then
                  tree_idx = icls
                  return
               endif
            end do
         end associate
      end do
   end function get_cls_indx

   ! private helper to search for node with inputted reference index
   recursive function search_tree4ref(self, root, ref_idx) result(node_index)
      type(multi_dendro), intent(inout)  :: self 
      type(s2_node), pointer, intent(in) :: root
      integer,                intent(in) :: ref_idx
      type(s2_node), pointer :: final
      type(s2_node), pointer :: cur
      integer,   allocatable :: tmp(:)
      integer :: j, n, node_index
      logical :: in_left
      ! root should be immutable 
      final => null()
      cur   => root
      do while (associated(cur))
         ! return if ref_idx found   
         if (cur%ref_idx == ref_idx) then
            final => cur
            node_index = final%node_idx
            return
         end if
         in_left = .false.
         ! heap sort + binary search to check if ref_idx is in the subset
         if (associated(cur%left)) then
            n = size(cur%left%subset)
            if (n > 0) then
               allocate(tmp(n), source=cur%left%subset)
               j = locate(tmp, n, ref_idx)
               if (j >= 1 .and. j <= n) then
                  if (tmp(j) == ref_idx .or. tmp(j+1) == ref_idx) in_left = .true.
               end if
               deallocate(tmp)
            end if
         end if
         if (in_left) then
            cur => cur%left
         else
            cur => cur%right
         end if
      end do
      final => null()
      node_index = final%node_idx
   end function search_tree4ref
   
   subroutine kill( self )
      class(multi_dendro), intent(inout) :: self
      integer  :: i, j
      if( self%exists )then
         if( allocated(self%node_store) )then
            do i = 1, self%n_trees
               if( allocated(self%node_store(i)%nodes) )then
                  do j = 1, size(self%node_store(i)%nodes)
                     nullify(self%node_store(i)%nodes(j)%left,self%node_store(i)%nodes(j)%right,self%node_store(i)%nodes(j)%parent)
                     if( allocated(self%node_store(i)%nodes(j)%subset) ) deallocate(self%node_store(i)%nodes(j)%subset)
                     self%node_store(i)%nodes(j)%level      = 0
                     self%node_store(i)%nodes(j)%visit      = .false.
                     self%node_store(i)%nodes(j)%is_pop     = .false.
                     self%node_store(i)%nodes(j)%ref_idx    = 0
                  end do
                  deallocate(self%node_store(i)%nodes)
               endif
               self%node_store(i)%root_idx = 0
            end do
            deallocate(self%node_store)
         endif
         if( allocated(self%dist_mat) ) deallocate(self%dist_mat)
         if( allocated(self%cls_pops) ) deallocate(self%cls_pops)
         if( allocated(self%cls_map)  ) deallocate(self%cls_map)
         if( allocated(self%medoids)  ) deallocate(self%medoids)
      endif
   end subroutine kill

end module simple_tree
