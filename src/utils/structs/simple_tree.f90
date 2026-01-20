module simple_tree  
use simple_srch_sort_loc
use simple_hclust
use simple_linked_list
implicit none 
#include "simple_local_flags.inc"

type  :: s2_node 
   type(s2_node), pointer   :: left => null(), right => null(), parent => null()
   real        :: F_ptcl2ref  
   integer     :: level
   logical     :: visit = .false.
   logical     :: is_pop   = .false.
   integer     :: ref_idx = 0
   integer, allocatable :: subset(:) ! array of refs remaining 
end type s2_node
type :: node_storage
   type(s2_node), pointer :: nodes(:) => null()
   integer :: root_idx = 0
end type node_storage
type  :: multi_dendro 
   type(s2_node), allocatable :: root_array(:)
   type(node_storage), allocatable :: node_store(:)
   real, allocatable          :: dist_mat(:,:) ! full_distmat
   integer, allocatable       :: cls_pops(:) 
   integer, allocatable       :: subsets(:,:) 
   integer, allocatable       :: medoids(:)
   integer, allocatable       :: heights(:)
   integer  :: n_trees ! number of AP clusters
contains 
   ! Setters / Getters
   procedure   :: get_heights
   procedure   :: set_cls_pops  
   procedure   :: set_subsets
   procedure   :: set_distmat
   ! Constructor
   procedure   :: build_multi_dendro
end type multi_dendro  
contains 

   pure function get_heights(self) result(k)
      class(multi_dendro ), intent(in) :: self
      integer  :: k(self%n_trees)
      k = self%heights
   end function get_heights

   pure subroutine set_cls_pops(self, labels)
      class(multi_dendro), intent(inout) :: self
      integer, intent(in)               :: labels(:)
      integer :: i
      if (allocated(self%cls_pops)) deallocate(self%cls_pops)
      self%n_trees = maxval(labels)
      allocate(self%cls_pops(self%n_trees), source=0)
      do i = 1, self%n_trees
         self%cls_pops(i) = count(labels == i)
      end do
   end subroutine set_cls_pops

   pure subroutine set_subsets(self, labels)
      class(multi_dendro ), intent(inout)  :: self 
      integer, intent(inout)           :: labels(:) 
      integer, allocatable          :: tmp1(:), tmp2(:)
      integer     :: i, j
      if (allocated(self%subsets))        deallocate(self%subsets)
      allocate(self%subsets(self%n_trees, maxval(self%cls_pops)))
      do i = 1, self%n_trees
         allocate(tmp1(self%cls_pops(i)))
         allocate(tmp2( maxval(self%cls_pops) - size(tmp1)), source = -1)
         tmp1 = pack([(j, j=1,size(labels))], labels == i)
         self%subsets(i,:) = [tmp1, tmp2]
         deallocate(tmp1, tmp2)  
      end do 
      if(allocated(self%root_array)) deallocate(self%root_array)
      allocate(self%root_array(self%n_trees))
   end subroutine set_subsets

   pure subroutine set_distmat(self, dist_mat)
      class(multi_dendro), intent(inout) :: self 
      real, intent(in)          :: dist_mat(:,:) 
      self%dist_mat = dist_mat
   end subroutine set_distmat

   subroutine build_multi_dendro(self, linkage)
      class(multi_dendro), intent(inout)   :: self
      integer, intent(in)  :: linkage
      real, allocatable    :: sub_distmat(:,:), height(:)
      integer, allocatable :: refs(:), merge_mat(:,:)
      integer  :: icls, i, j, ncls_ap, nref_sub
      type(hclust) :: hc
      ncls_ap = self%n_trees
      allocate(self%node_store(ncls_ap))
      allocate(self%heights(ncls_ap))
      do icls = 1, ncls_ap
         nref_sub = self%cls_pops(icls)
         allocate(refs(nref_sub))
         refs = self%subsets(icls, 1:nref_sub)
         allocate(sub_distmat(nref_sub, nref_sub))
         sub_distmat = 1.0
         do i=1,nref_sub
            do j=i+1,nref_sub
               sub_distmat(i,j) = self%dist_mat(refs(i), refs(j))
               sub_distmat(j,i) = sub_distmat(i,j)
            end do
         end do
         ! if cluster size is 1, done 
         if(nref_sub == 1) then 
            allocate(self%node_store(icls)%nodes(1))
            self%node_store(icls)%nodes(1)%ref_idx = refs(1)
            allocate(self%node_store(icls)%nodes(1)%subset(1))
            self%node_store(icls)%nodes(1)%subset = [refs(1)]
            self%node_store(icls)%root_idx = 1
            self%heights(icls) = 0
            deallocate(refs, sub_distmat)
            cycle 
         end if 
         ! Aggl. Clustering
         allocate(merge_mat(2, nref_sub-1))
         allocate(height(nref_sub - 1))
         call hc%new(nref_sub, sub_distmat, linkage)
         call hc%cluster(merge_mat, height)
         call hc%kill()
         ! Tracking Merges
         call merge2node(merge_mat, height, refs, self%node_store(icls)%nodes, self%node_store(icls)%root_idx)
         ! point root to rest of the tree 
         associate(r => self%node_store(icls)%nodes(self%node_store(icls)%root_idx))
            self%root_array(icls)%left  => r%left
            self%root_array(icls)%right => r%right
            self%root_array(icls)%ref_idx = r%ref_idx
            self%root_array(icls)%level   = r%level
            if (allocated(self%root_array(icls)%subset)) deallocate(self%root_array(icls)%subset)
            allocate(self%root_array(icls)%subset(size(r%subset)))
            self%root_array(icls)%subset = r%subset
         end associate
         self%heights(icls) = nref_sub - 1
         deallocate(sub_distmat, refs, merge_mat, height)
      end do 
      contains 
         subroutine merge2node(merge_mat, height, refs, nodes, root_idx)
            integer, intent(in) :: merge_mat(:,:)
            real,    intent(in) :: height(:)
            integer, intent(in) :: refs(:)
            type(s2_node), pointer, intent(inout) :: nodes(:)
            integer, intent(out) :: root_idx
            integer :: k, s, l, r, m, best, n_nodes
            real    :: best_sum, sum 
            integer, allocatable :: tmp(:)
            n_nodes = 2*nref_sub - 1
            allocate(nodes(n_nodes))
            ! allocating nodes 
            do k = 1, n_nodes
               nullify(nodes(k)%left, nodes(k)%right, nodes(k)%parent)
               nodes(k)%level      = 0
               nodes(k)%visit      = .false.
               nodes(k)%is_pop     = .false.
               nodes(k)%ref_idx    = 0
               ! nodes(k)%F_ptcl2ref = 0.0
               if (allocated(nodes(k)%subset)) deallocate(nodes(k)%subset)
               ! can set all leaves
               if (k <= nref_sub) then
                  nodes(k)%ref_idx = refs(k)
                  nodes(k)%is_pop  = .true.
                  nodes(k)%level   = 0 
                  allocate(nodes(k)%subset(1))
                  nodes(k)%subset = [refs(k)]
               end if
            end do
            ! set internal nodes / root. 
            m = size(refs)
            do s = 1, nref_sub - 1
               l = merge_mat(1, s)
               r = merge_mat(2, s)
               k = m + s
               ! assigning pointers 
               nodes(k)%left  => nodes(l)
               nodes(k)%right => nodes(r)
               nodes(l)%parent => nodes(k)
               nodes(r)%parent => nodes(k)
               nodes(k)%level      = s
               nodes(k)%is_pop     = .true.
               ! parent subset is union of children 
               allocate(tmp(size(nodes(l)%subset) + size(nodes(r)%subset)))
               tmp = [nodes(l)%subset, nodes(r)%subset]
               allocate(nodes(k)%subset(size(tmp)))
               nodes(k)%subset = tmp
               deallocate(tmp)
               ! Calculate a new medoid in merged set 
               ! probably large opportunity to memoize (can store partial sums )
               best = -1
               best_sum = huge(1.0)
               do i = 1, size(nodes(k)%subset)
                  sum = 0.0
                  do j = 1, size(nodes(k)%subset)
                     sum = sum + self%dist_mat(nodes(k)%subset(i), nodes(k)%subset(j))
                  end do
                  if (sum < best_sum) then
                     best_sum = sum 
                     best = nodes(k)%subset(i)
                  end if
               end do
               nodes(k)%ref_idx = best
            end do
            root_idx = n_nodes
         end subroutine merge2node
   end subroutine build_multi_dendro   

   recursive subroutine print_tree(root)
      type(s2_node), pointer  :: root
      ! point to root 
      if(associated(root)) then 
         print *, root%subset, '//', root%ref_idx
         call print_tree(root%right)
         call print_tree(root%left)
      end if 
   end subroutine print_tree
   ! find cluster to which an aribitrary reference belongs to 
   subroutine make_subtree_ll(self, ref_idx, tree_idx)
      type(multi_dendro), intent(inout)   :: self 
      integer, intent(in)        :: ref_idx
      integer, intent(out)       :: tree_idx
      type(s2_node), allocatable :: roots(:) 
      integer, allocatable :: refs(:), tmp(:)
      type(linked_list) :: list
      integer           :: i, nrefs, ncls, icls, iref
      allocate(roots(self%n_trees))
      nrefs = size(self%dist_mat)
      allocate(refs(nrefs))
      do icls = 1, ncls
         allocate(tmp(self%cls_pops(icls)))
         tmp = self%node_store(icls)%nodes(self%node_store(icls)%root_idx)%subset
         if(icls == 1) then
            refs(1:size(tmp)) = tmp 
         else 
            refs(self%cls_pops(icls - 1):self%cls_pops(icls - 1) + self%cls_pops(icls)) = tmp
         end if 
         deallocate(tmp)
      end do 
      do iref = 1, nrefs
         call list%push_back(iref)
      end do 
      do iref = 1, nrefs
         call list%at_int(iref, refs(iref))
      end do 
   end subroutine make_subtree_ll

   ! subroutine search_tree4ref(root, ref_indx)
   !    type(s2_node), pointer, intent(inout) :: root 
   !    integer, intent(out)  :: indxs(2)
   !    integer, intent(in)   :: ref_idx   
   !    ! found which subtree to which reference belongs to 
   !    ! dfs from bottom, maybe just follow to root node... 
   ! end subroutine search_tree4ref


   ! recursive subroutine walk_anywhere() 


   ! end subroutine walk_anywhere()

   recursive subroutine print_s2_tree(root, unit, indent, show_subset, max_subset)
      type(s2_node), intent(in) :: root
      integer,        intent(in), optional :: unit
      integer,        intent(in), optional :: indent
      logical,        intent(in), optional :: show_subset
      integer,        intent(in), optional :: max_subset

      integer :: u, ind, ms, nshow
      logical :: do_subset
      character(len=:), allocatable :: pad

      if (present(unit)) then
         u = unit
      else
         u = output_unit
      end if
   
      if (present(indent)) then
         ind = indent
      else
         ind = 0
      end if
   
      if (present(show_subset)) then
         do_subset = show_subset
      else
         do_subset = .false.
      end if
   
      if (present(max_subset)) then
         ms = max_subset
      else
         ms = 10
      end if

      pad = repeat(" ", ind)

      write(u,'(a,"node: lvl=",i0," ref=",i0," pop=",l1," visit=",l1)') &
         pad, root%level, root%ref_idx, root%is_pop, root%visit

      if (do_subset) then
         if (allocated(root%subset)) then
            if (size(root%subset) == 0) then
               write(u,'(a,"  subset(n=0): <empty>")') pad
            else
               nshow = min(size(root%subset), ms)
               write(u,'(a,"  subset(n=",i0,"): ",*(i0,1x))') pad, size(root%subset), root%subset(1:nshow)
               if (size(root%subset) > nshow) write(u,'(a,"  ... (truncated)")') pad
            end if
         else
            write(u,'(a,"  subset: <unallocated>")') pad
         end if
      end if

      if (associated(root%left)) then
         write(u,'(a,"  L:")') pad
         call print_s2_tree(root%left, unit=u, indent=ind+4, show_subset=do_subset, max_subset=ms)
      end if

      if (associated(root%right)) then
         write(u,'(a,"  R:")') pad
         call print_s2_tree(root%right, unit=u, indent=ind+4, show_subset=do_subset, max_subset=ms)
      end if
   end subroutine print_s2_tree

   subroutine print_multi_dendro(self, unit, show_subset, max_subset)
      class(multi_dendro), intent(in) :: self
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_subset
      integer, intent(in), optional :: max_subset

      integer :: u, i
      logical :: do_subset
      integer :: ms

      if (present(unit)) then
         u = unit
      else
         u = output_unit
      end if
   
      if (present(show_subset)) then
         do_subset = show_subset
      else
         do_subset = .false.
      end if
   
      if (present(max_subset)) then
         ms = max_subset
      else
         ms = 10
      end if

      write(u,'("multi_dendro: n_trees=",i0)') self%n_trees
      if (allocated(self%heights)) write(u,'("heights: ",*(i0,1x))') self%heights

      if (.not. allocated(self%root_array)) then
         write(u,'("root_array: <unallocated>")')
         return
      end if

      do i = 1, size(self%root_array)
         write(u,'(/,"=== Tree ",i0," (root ref=",i0,") ===")') i, self%root_array(i)%ref_idx
         call print_s2_tree(self%root_array(i), unit=u, indent=0, show_subset=do_subset, max_subset=ms)
      end do
   end subroutine print_multi_dendro      
end module simple_tree
