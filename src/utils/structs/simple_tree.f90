module simple_tree  
use simple_srch_sort_loc
use simple_hclust
use simple_string
use simple_hash
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
               call hpsort(tmp)
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
   subroutine get_cls_indx(self, ref_idx, tree_idx)
      type(multi_dendro), intent(inout)   :: self 
      integer, intent(in)        :: ref_idx
      integer, intent(out)       :: tree_idx 
      integer, allocatable :: refs(:), tmp(:)
      type(hash), allocatable     :: hsh(:)
      integer           :: i, nrefs, ncls, icls, iref
      character(len=KEYLEN)      :: key
      ncls  = self%n_trees
      if(ncls == 1) then 
         tree_idx = 1
         return 
      end if 
      allocate(hsh(ncls))
      do icls = 1, ncls
         allocate(tmp(self%cls_pops(icls)))
         tmp = self%node_store(icls)%nodes(self%node_store(icls)%root_idx)%subset
         nrefs = size(tmp)
         hsh(icls) = hash()
         ! set keys for each hash table
         do iref = 1, nrefs 
            write(key, '(I0)') tmp(iref)
            key = trim(key)  
            call hsh(icls)%set(key, iref)
         end do 
         deallocate(tmp)
         write(key, '(I0)') ref_idx
         key = trim(key) 
         if(hsh(icls)%isthere(trim(key))) then 
            tree_idx = icls
            return 
         end if
      end do 
   end subroutine get_cls_indx

   ! Find first instance of ref_idx in subtree
   recursive function search_tree4ref(root, ref_idx) result(final)
      type(s2_node), pointer, intent(in) :: root
      integer, intent(in)                :: ref_idx
      type(s2_node), pointer             :: final
      type(s2_node), pointer             :: cur
      integer, allocatable               :: tmp(:)
      integer                            :: j, n
      logical                            :: in_left
      ! root should be immutable 
      final => null()
      cur   => root
      do while (associated(cur))
         ! return if ref_idx found   
         if (cur%ref_idx == ref_idx) then
            final => cur
            return
         end if
         in_left = .false.
         ! heap sort + binary search to check if ref_idx is in the subset
         if (associated(cur%left)) then
            n = size(cur%left%subset)
            if (n > 0) then
            allocate(tmp(n))
            tmp = cur%left%subset
            j = locate(tmp, n, ref_idx)
            if (j >= 1 .and. j < n) then
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
   end function search_tree4ref

   ! Walk from Node
   recursive subroutine walk_from_node(cur, indxs, objs, done) 
      type(s2_node), pointer, intent(inout)  :: cur
      integer, intent(out)  :: indxs(2)
      real, intent(in)    :: objs(2)
      logical     :: done 
      cur%left%F_ptcl2ref  = objs(1)
      cur%right%F_ptcl2ref = objs(2)
      ! need to cache
      if(cur%F_ptcl2ref > objs(1) .and. cur%F_ptcl2ref > objs(2)) then 
         done = .true. 
         return 
      end if 
      if(objs(1) > objs(2)) then 
         cur => cur%left 
      else 
         cur => cur%right    
      end if
      if (.not.(associated(cur%left) .and. associated(cur%right))) then
         done =.true. 
         return
      end if
      indxs(1) = cur%left%ref_idx
      indxs(2) = cur%right%ref_idx 
   end subroutine walk_from_node

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
