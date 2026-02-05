!@descr: binary tree data structure
module simple_tree  
use simple_core_module_api
use simple_srch_sort_loc, only: hpsort, locate, mask2inds
use simple_hclust,        only: hclust
implicit none

public :: multi_dendro
public :: test_simple_tree
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
   logical :: exists  = .false.
contains
   ! constructor
   procedure          :: new
   ! tree builder
   procedure          :: build_multi_dendro
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
      allocate(self%medoids(self%n_trees))
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
         self%medoids(icls) = self%node_store(icls)%nodes(self%node_store(icls)%root_idx)%ref_idx
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

   ! getter to return left and right indices
   subroutine get_left_right_idxs(self, node_idx, ref_idx)
      class(multi_dendro), intent(in)    :: self
      integer,             intent(inout) :: node_idx(2), ref_idx(2)
      integer :: imedoid, n_nodes, cur
      integer :: root_idx
      type(s2_node), pointer :: lp, rp
      imedoid  = get_cls_indx(self, ref_idx(1))
      n_nodes  = 2*self%cls_pops(imedoid) - 1
      root_idx = self%node_store(imedoid)%root_idx
      if (node_idx(1) == n_nodes) then
         cur = search_tree4ref_idx(self%node_store(imedoid)%nodes, ref_idx(1))
         if (cur == 0) THROW_HARD('search_tree4ref_idx failed in get_left_right_idxs')
         node_idx(1) = cur
      else
         cur = node_idx(1)
         if (cur < 1 .or. cur > n_nodes) THROW_HARD('node_idx out of range in get_left_right_idxs')
      end if
      ! Left child
      lp => self%node_store(imedoid)%nodes(cur)%left
      if (associated(lp)) then
         ref_idx(1)  = lp%ref_idx
         node_idx(1) = lp%node_idx
      else
         ref_idx(1)  = self%node_store(imedoid)%nodes(cur)%ref_idx
         node_idx(1) = self%node_store(imedoid)%nodes(cur)%node_idx
      end if
      ! Right child
      rp => self%node_store(imedoid)%nodes(cur)%right
      if (associated(rp)) then
         ref_idx(2)  = rp%ref_idx
         node_idx(2) = rp%node_idx
      else
         ref_idx(2)  = self%node_store(imedoid)%nodes(cur)%ref_idx
         node_idx(2) = self%node_store(imedoid)%nodes(cur)%node_idx
      end if
   end subroutine get_left_right_idxs

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

   ! private helper to search for node idx with inputted reference index
   pure function search_tree4ref_idx(nodes, ref_idx) result(found_idx)
      type(s2_node), intent(in) :: nodes(:)
      integer,       intent(in) :: ref_idx
      integer :: found_idx
      integer :: cur_idx, j, n
      logical :: in_left
      found_idx = 0
      cur_idx   = size(nodes)
      do
         if (cur_idx < 1 .or. cur_idx > size(nodes)) return
         if (nodes(cur_idx)%ref_idx == ref_idx) then
            found_idx = cur_idx
            return
         end if
         in_left = .false.
         if (associated(nodes(cur_idx)%left)) then
            n = size(nodes(cur_idx)%left%subset)
            if (n == 1) then
               in_left = (nodes(cur_idx)%left%subset(1) == ref_idx)
            else if (n > 1) then
               j = locate(nodes(cur_idx)%left%subset, n, ref_idx)
               if (j >= 1 .and. j <= n-1) then
                  in_left = (nodes(cur_idx)%left%subset(j) == ref_idx) .or. &
                            (nodes(cur_idx)%left%subset(j+1) == ref_idx)
               end if
            end if
         end if
         if (in_left) then
            if (.not. associated(nodes(cur_idx)%left)) return
            cur_idx = nodes(cur_idx)%left%node_idx
         else
            if (.not. associated(nodes(cur_idx)%right)) return
            cur_idx = nodes(cur_idx)%right%node_idx
         end if
      end do
   end function search_tree4ref_idx
   
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

   subroutine test_simple_tree()
      type(multi_dendro) :: md
      integer, allocatable :: labels(:)
      real,    allocatable :: dist(:,:)
      integer :: n, i, j, icls
      integer :: ref_in, n_nodes
      integer :: node_idxs(2), ref_idxs(2)
      integer :: found_idx
      integer :: imed

      n = 6
      allocate(dist(n,n), labels(n))

      ! Two clusters: {1,2,3} and {4,5,6}
      labels = [1,1,1, 2,2,2]

      dist = 0.0
      do i = 1, n
         do j = 1, n
            if (i == j) then
               dist(i,j) = 0.0
            else if (labels(i) == labels(j)) then
               dist(i,j) = 1.0 + 0.1*abs(real(i-j))
            else
               dist(i,j) = 10.0 + 0.1*abs(real(i-j))
            end if
         end do
      end do

      call md%new(dist, labels)
      call md%build_multi_dendro(1)

      do icls = 1, md%n_trees
         if (2*md%cls_pops(icls) - 1 /= size(md%node_store(icls)%nodes)) then
            print *, 'TEST FAILED: TREE ASSEMBLED INCORRECT'
            return
         end if
      end do
      print *, 'TREE ASSEMBLED CORRECTLY'

      ! pick a reference to test
      ref_in = 6

      ! verify get_cls_indx works (private helper)
      imed = get_cls_indx(md, ref_in)
      if (imed == labels(ref_in)) then
         print *, 'GET_CLS_INDX IS CORRECT'
      else
         print *, 'GET_CLS_INDX IS INCORRECT: imed=', imed, 'labels(ref_in)=', labels(ref_in)
         return
      end if

      ! Show the root node (access by index, no copying)
      print *, 'root(ref_idx,node_idx)=', &
               md%node_store(imed)%nodes(md%node_store(imed)%root_idx)%ref_idx, &
               md%node_store(imed)%nodes(md%node_store(imed)%root_idx)%node_idx
      print *, 'ROOT NODE SET (INDEX ACCESS)'

      ! Check search_tree4ref_idx (index-based)
      found_idx = search_tree4ref_idx( md%node_store(imed)%nodes, ref_in )
      if (found_idx > 0 .and. found_idx < 2*md%cls_pops(imed) - 1 .and. &
          md%node_store(imed)%nodes(found_idx)%ref_idx == ref_in) then
         print *, 'SEARCH_TREE4REF_IDX IS CORRECT'
      else
         print *, 'SEARCH_TREE4REF_IDX IS INCORRECT. found_idx=', found_idx
         return
      end if

      ! ------------------------------------------------------------
      ! Single-call sanity for get_left_right_idxs
      ! We must set ref_idxs(1) so get_left_right_idxs can infer the cluster
      ! starting at sentinel root: node_idxs(1) = n_nodes, ref_idxs(1) = medoid
      ! ------------------------------------------------------------
      n_nodes = 2*md%cls_pops(imed) - 1
      node_idxs = [n_nodes, 0]
      ref_idxs  = [ md%medoids(imed), 0 ]   ! IMPORTANT: seed ref_idx(1) with medoid

      call get_left_right_idxs(md, node_idxs, ref_idxs)
      print *, '>>> GETTING LEFT/RIGHT (single call)'
      print *, 'node_idxs', node_idxs, 'ref_idxs', ref_idxs

      ! ============================================================
      ! Robust traversal test: repeatedly call get_left_right_idxs to traverse
      ! entire tree from root, like in production usage. The routine
      ! infers the tree from ref_idx(1) on each call.
      ! ============================================================
      block
         integer :: top, cur_node, cur_ref
         integer, allocatable :: node_stack(:), ref_stack(:)
         logical, allocatable :: seen(:)
         integer :: node_idxs2(2), ref_idxs2(2)
         integer :: expansions, max_expansions
         logical :: failed

         n_nodes = 2*md%cls_pops(imed) - 1

         allocate(seen(n_nodes), source=.false.)
         allocate(node_stack(n_nodes), ref_stack(n_nodes))

         top = 1
         node_stack(top) = n_nodes
         ref_stack(top)  = md%medoids(imed)

         expansions     = 0
         max_expansions = n_nodes
         failed         = .false.

         do while (top > 0)
            cur_node = node_stack(top)
            cur_ref  = ref_stack(top)
            top = top - 1

            node_idxs2 = [cur_node, 0]
            ref_idxs2  = [cur_ref,  0]

            ! call the version that infers the imedoid from ref_idxs2(1)
            call get_left_right_idxs(md, node_idxs2, ref_idxs2)
            expansions = expansions + 1

            if (expansions > max_expansions) then
               print *, 'TEST FAILED: get_left_right loop exceeded max expansions (possible cycle)'
               failed = .true.
               exit
            end if

            if (node_idxs2(1) < 1 .or. node_idxs2(1) > n_nodes) then
               print *, 'TEST FAILED: left child node index out of range:', node_idxs2(1)
               failed = .true.
               exit
            end if
            if (node_idxs2(2) < 1 .or. node_idxs2(2) > n_nodes) then
               print *, 'TEST FAILED: right child node index out of range:', node_idxs2(2)
               failed = .true.
               exit
            end if

            ! Mark the node we are expanding (counts root too)
            if (cur_node >= 1 .and. cur_node <= n_nodes) seen(cur_node) = .true.

            ! Mark returned children
            seen(node_idxs2(1)) = .true.
            seen(node_idxs2(2)) = .true.

            ! Optional stronger invariant: children subsets partition parent subset
            ! (uncomment if you want extra checking)
            ! if (.not. (node_idxs2(1) == node_idxs2(2) .and. ref_idxs2(1) == ref_idxs2(2))) then
            !    if ( size(md%node_store(imed)%nodes(node_idxs2(1))%subset) + &
            !         size(md%node_store(imed)%nodes(node_idxs2(2))%subset) /= &
            !         size(md%node_store(imed)%nodes(cur_node)%subset) ) then
            !       print *, 'TEST FAILED: child subsets do not sum to parent subset size at node ', cur_node
            !       failed = .true.
            !       exit
            !    end if
            ! end if

            ! Leaf convention: if no children, routine returns same node for left/right
            if (.not. (node_idxs2(1) == node_idxs2(2) .and. ref_idxs2(1) == ref_idxs2(2))) then
               if (top + 2 > n_nodes) then
                  print *, 'TEST FAILED: traversal stack overflow (unexpected for a tree)'
                  failed = .true.
                  exit
               end if

               top = top + 1
               node_stack(top) = node_idxs2(1)
               ref_stack(top)  = ref_idxs2(1)

               top = top + 1
               node_stack(top) = node_idxs2(2)
               ref_stack(top)  = ref_idxs2(2)
            end if
         end do

         if (.not. failed) then
            if (count(seen) /= n_nodes) then
               print *, 'TEST FAILED: traversal did not visit all nodes. Seen=', count(seen), 'Expected=', n_nodes
            else
               print *, 'GET_LEFT_RIGHT LOOP TEST PASSED: visited all nodes via repeated calls'
            end if
         end if

         deallocate(seen, node_stack, ref_stack)
      end block

      call md%kill()
      print *, 'TREE KILLED'
      deallocate(dist, labels)
   end subroutine test_simple_tree

end module simple_tree
