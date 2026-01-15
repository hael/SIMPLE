module simple_tree  
use simple_srch_sort_loc
use iso_c_binding
use iso_fortran_env
implicit none 
#include "simple_local_flags.inc"

type :: s1_node
   type(s1_node), pointer  :: left => null(), right => null(), parent => null()
   real     :: inpl_ang
   logical  :: visited = .true. 
   integer, allocatable :: inpl_subset(:) ! array of inpl angles below s1 node
end type s1_node

type  :: s2_node 
   type(s2_node), pointer   :: left => null(), right => null()
   type(c_ptr)              :: parent = c_null_ptr
   real        :: F_ptcl2ref  
   integer     :: level
   logical     :: visit = .false.
   logical     :: is_pop   = .false.
   integer     :: ref_idx = 0
   integer, allocatable :: subset(:) ! array of refs remaining 
end type s2_node

type  :: multi_dendro 
   type(s2_node), allocatable :: root_array(:)
   real, allocatable :: dist_mat(:,:) ! full_distmat
   integer, allocatable :: cls_pops(:) 
   integer, allocatable :: subsets(:,:) 
   integer, allocatable :: subsets_avail(:,:)
   integer, allocatable :: medoids(:)
   integer, allocatable :: heights(:)
   integer  :: n_trees ! number of clusters
contains 
   ! Setters / Getters
   procedure   :: get_heights
   procedure   :: set_cls_pops  
   procedure   :: set_subsets
   procedure   :: set_distmat
   procedure   :: set_medoids
   ! Constructor
   procedure   :: build_multi_dendro
end type multi_dendro  

contains 

   ! getters / setters 
   ! rewrite this 
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
      if (allocated(self%subsets_avail))  deallocate(self%subsets_avail)
      allocate(self%subsets(self%n_trees, maxval(self%cls_pops)))
      do i = 1, self%n_trees
         allocate(tmp1(self%cls_pops(i)))
         allocate(tmp2( maxval(self%cls_pops) - size(tmp1)), source = -1)
         tmp1 = pack([(j, j=1,size(labels))], labels == i)
         self%subsets(i,:) = [tmp1, tmp2]
         deallocate(tmp1, tmp2)  
      end do 
      allocate(self%subsets_avail(size(self%cls_pops), maxval(self%cls_pops)))
      self%subsets_avail = self%subsets
   end subroutine set_subsets

   pure subroutine set_distmat(self, dist_mat)
      class(multi_dendro), intent(inout) :: self 
      real, intent(in)          :: dist_mat(:,:) 
      self%dist_mat = dist_mat
   end subroutine set_distmat

   pure subroutine set_medoids(self, centers)
      class(multi_dendro), intent(inout)  :: self
      integer, intent(inout)           :: centers(:)
      integer  :: i
      if(allocated(self%root_array)) deallocate(self%root_array)
      self%n_trees = size(centers)
      allocate(self%root_array(self%n_trees))
      do i = 1, size(centers)
         self%root_array(i)%ref_idx = centers(i)
      end do 
   end subroutine 

   subroutine build_multi_dendro (self)
      class(multi_dendro), intent(inout)   :: self
      real, allocatable    :: tmp(:)
      integer  :: i, j 
      ! Initialize each sub_tree 
      allocate(self%heights(self%n_trees))
      do i = 1, self%n_trees
         self%root_array(i)%level = 0 
         allocate(self%root_array(i)%subset(self%cls_pops(i)))
         self%root_array(i)%subset = self%subsets(i,1:self%cls_pops(i))
         self%root_array(i)%is_pop = .true.
         self%heights(i) = 0
         call clust_insert_s2_node(self%root_array(i), self%dist_mat, 0)
      end do 
      contains 
         recursive subroutine clust_insert_s2_node(root, dist_mat, level)
            type(s2_node), intent(inout), target:: root
            real, intent(in)             :: dist_mat(:,:)
            integer, intent(in)          :: level
            integer                    :: l, idx, new_med(2), closest, sec_closest, n, m
            real                       :: d, d1, d2 
            integer, allocatable       :: new_subset(:), tmp1(:)
            ! update subsets available at multi_dendro level so each node is getting assigned from updated pool 
            allocate(tmp1(maxval(self%cls_pops) - size(root%subset)), source = -1)
            self%subsets_avail(i,:) = [root%subset, tmp1]
            root%level = level
            if(level > self%heights(i)) self%heights(i) = level
            if(size(root%subset) < 2) return
            ! identify closest and 2nd closest to root
            closest     = -1 ; sec_closest = -1
            d1 = huge(1.0) ; d2 = huge(1.0)
            if (level == 0 .and. (size(root%subset) < 2)) return
            l = 1
            do while (l <= size(root%subset))
               idx = self%subsets_avail(i,l)
               l = l + 1
               if (idx == root%ref_idx) cycle
               d = dist_mat(idx, root%ref_idx)
               if (d < d1) then
                  d2 = d1
                  sec_closest = closest
                  d1 = d
                  closest = idx
               else if (d < d2) then
                  d2 = d
                  sec_closest = idx
               end if
            end do
            ! remove those from subset
            if(level > 0 ) then 
               allocate(new_subset(size(root%subset) - 2))
            else 
               allocate(new_subset(size(root%subset) - 3))
            end if 
            new_subset = pack(root%subset, &
            root%subset /= closest .and. root%subset /= sec_closest .and. root%subset /= root%ref_idx)

            ! push 2nd closest left 
            allocate(root%left)
            nullify(root%left%left, root%left%right)
            allocate(root%left%subset(size(new_subset)))
            root%left%subset = new_subset
            root%left%ref_idx = sec_closest
            root%left%level = level + 1 
            root%left%visit  = .false.
            root%left%is_pop = .true.
            root%left%parent = c_loc(root)

            ! push closest right 
            allocate(root%right)
            nullify(root%right%left, root%right%right)
            allocate(root%right%subset(size(new_subset)))
            root%right%subset = new_subset
            root%right%ref_idx = closest
            root%right%level  = level + 1 
            root%right%visit  = .false.
            root%right%is_pop = .true.
            root%right%parent = c_loc(root)
            ! cleanup
            if (size(root%subset) < 3) return
            deallocate(new_subset, tmp1)

            call clust_insert_s2_node(root%left,  dist_mat, level + 1)
            call clust_insert_s2_node(root%right, dist_mat, level + 1)
         end subroutine clust_insert_s2_node

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

   recursive subroutine walk_tree(root, indxs, objs, done)
      type(s2_node), pointer, intent(inout) :: root 
      integer, intent(out)  :: indxs(2)
      real, intent(in)    :: objs(2)
      logical     :: done   
      
      root%left%F_ptcl2ref  = objs(1)
      root%right%F_ptcl2ref = objs(2)

      if(objs(1) > objs(2)) then 
         root => root%left 
      else 
         root => root%right
      end if 
      root%visit = .true. 
      if (.not.(associated(root%left) .and. associated(root%right))) then
         done =.true. 
         return
      end if

      indxs(1) = root%left%ref_idx
      indxs(2) = root%right%ref_idx 
   end subroutine walk_tree

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
