module simple_tree  
use simple_srch_sort_loc
implicit none 
#include "simple_local_flags.inc"

type :: s1_node
   type(s1_node), pointer  :: left => null(), right => null(), parent => null()
   real     :: inpl_ang
   logical  :: visited = .true. 
   integer, allocatable :: inpl_subset(:) ! array of inpl angles below s1 node
end type s1_node

type  :: s2_node 
   type(s2_node), pointer   :: left => null(), right => null(), parent => null()
   real        :: F_ptcl2ref  
   integer     :: level
   logical     :: visit = .false.
   logical     :: is_pop   = .false.
   integer     :: ref_idx = 0
   integer, allocatable     :: subset(:) ! array of refs remaining 
end type s2_node

type  :: multi_dendro 
   type(s2_node), allocatable :: root_array(:)
   real, allocatable :: dist_mat(:,:) ! full_distmat
   integer, allocatable :: npnts_array(:) 
   integer  :: n_trees ! number of clusters
   integer, allocatable :: subsets(:,:) 
   integer, allocatable :: subsets_avail(:,:)
   integer, allocatable :: medoids(:)
   integer, allocatable :: heights(:)
contains 
   ! Setters / Getters
   procedure   :: get_heights
   procedure   :: set_npnts_array  
   procedure   :: set_subsets
   procedure   :: set_distmat
   procedure   :: set_medoids
   ! Constructor
   procedure   :: build_multi_dendro 
end type multi_dendro  

contains 

   ! getters / setters 
   pure function get_heights(self) result(k)
      class(multi_dendro ), intent(in) :: self
      integer  :: k(self%n_trees)
      k = self%heights
   end function get_heights

   pure subroutine set_npnts_array(self, labels)
      class(multi_dendro ), intent(inout)  :: self 
      integer, intent(in)           :: labels(:)     
      integer  :: npnts_array(self%n_trees), i 
      do i = 1, maxval(labels)
         npnts_array(i) = count(labels == i)
      end do 
      allocate(self%npnts_array(self%n_trees))
      self%npnts_array = npnts_array
   end subroutine set_npnts_array

   pure subroutine set_subsets(self, labels)
      class(multi_dendro ), intent(inout)  :: self 
      integer, intent(in)           :: labels(:) 
      integer, allocatable          :: tmp1(:), tmp2(:)
      integer     :: i, j
      allocate(self%subsets(size(self%npnts_array), maxval(self%npnts_array)))
      do i = 1, size(self%npnts_array)
         allocate(tmp1(self%npnts_array(i)))
         allocate(tmp2( maxval(self%npnts_array) - size(tmp1)), source = -1)
         tmp1 = pack([(j, j=1,size(labels))], labels == i)
         self%subsets(i,:) = [tmp1, tmp2]
         deallocate(tmp1, tmp2)  
      end do 
      allocate(self%subsets_avail(size(self%npnts_array), maxval(self%npnts_array)))
      self%subsets_avail = self%subsets
   end subroutine set_subsets

   pure subroutine set_distmat(self, dist_mat)
      class(multi_dendro ), intent(inout) :: self 
      real, intent(in)          :: dist_mat(:,:) 
      self%dist_mat = dist_mat
   end subroutine set_distmat

   pure subroutine set_medoids(self, centers)
      class(multi_dendro ), intent(inout)  :: self
      integer, intent(in)           :: centers(:)
      integer  :: i
      do i = 1, size(self%npnts_array)
         self%root_array(i)%ref_idx = centers(i)
      end do 
   end subroutine 

   subroutine build_multi_dendro (self)
      class(multi_dendro ), intent(inout)   :: self
      real, allocatable    :: tmp(:)
      integer  :: i, j 
      type(s2_node), pointer       :: print 
      ! Initialize each sub_tree 
      allocate(self%root_array(size(self%npnts_array)))
      do i = 1, size(self%npnts_array)
         allocate(self%root_array(i))
         self%root_array(i)%level = 0 
         allocate(self%root_array(i)%subset(self%npnts_array(i)))
         self%root_array(i)%subset = self%subsets(i,1:self%npnts_array(i))
         call clust_insert_s2_node(self%root_array(i), self%dist_mat, 0)
      end do 
      contains 
         recursive subroutine clust_insert_s2_node(root, dist_mat, level)
            type(s2_node), intent(inout) :: root
            real, intent(in)             :: dist_mat(:,:)
            integer, intent(in)          :: level
            integer                    :: l, idx, new_med(2), closest, sec_closest, n, m
            real                       :: d, d1, d2 
            integer, allocatable       :: new_subset(:), tmp1(:)
            ! update subsets available at multi_dendro level so each node is getting assigned from updated pool 
            allocate(tmp1(maxval(self%npnts_array) - size(root%subset)), source = -1)
           
            self%subsets_avail(i,:) = [root%subset, tmp1]
            root%level = level
            if(size(root%subset) < 2) return
            ! identify closest and 2nd closest to root
            l = 1
            do while (l <= size(root%subset))
               idx = self%subsets_avail(i,l) 
               if(idx == root%ref_idx) cycle
               d = dist_mat(idx, root%ref_idx)
               if(d < d1) then 
                  d2 = d1 
                  sec_closest = closest
                  d1 = d 
                  closest = idx
               else if(d < d2) then 
                  d2 = d 
                  sec_closest = idx
               end if
                l = l + 1
            end do 
            ! remove those from subset
            allocate(new_subset(size(root%subset) - 2))
            n = 0 
            do m = 1, size(root%subset)
               if (root%subset(m) /= closest .and. root%subset(m) /= sec_closest) then
                  n = n + 1
                  new_subset(n) = root%subset(m)
               end if
            end do
            ! push 2nd closest left 
            allocate(root%left)
            nullify(root%left%left, root%left%right)
            allocate(root%left%subset(size(root%subset) - 2))
            root%left%subset = new_subset
            root%left%ref_idx = sec_closest
            root%left%level = level + 1 
            root%left%visit  = .false.
            root%left%is_pop = .true.
            ! push closest right 
            allocate(root%right)
            nullify(root%right%left, root%right%right)
            allocate(root%right%subset(size(root%subset) - 2))
            root%right%subset = new_subset
            root%right%ref_idx = closest
            root%right%level  = level + 1 
            root%right%visit  = .false.
            root%right%is_pop = .true.

            ! cleanup
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

end module simple_tree
