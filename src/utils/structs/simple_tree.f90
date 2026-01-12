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
   real        :: psi, theta ! S2 Orientation (change to oris later)
   integer     :: ref_idx = 0
   integer, allocatable     :: o_subset(:) ! array of refs remaining 
end type s2_node

type  :: dendro
   type(s2_node), pointer :: root 
   type(s2_node), allocatable :: root_array(:)
   integer  :: height !height array 
   integer  :: npnts !npnts array
   real, allocatable :: dist_mat(:,:) ! full_distmat
   integer, allocatable :: npnts_array(:) 
   integer  :: n_trees ! number of clusters
   integer, allocatable :: subsets(:,:) 
   integer, allocatable :: medoids(:)
contains 
   ! Setters / Getters
   procedure   :: get_height
   procedure   :: set_npnts 
   procedure   :: set_distmat
   ! Constructors / Methods 
   procedure   :: build_dendro
end type dendro 

contains 

   ! getters / setters 
   elemental integer function get_height(self) result(k)
      class(dendro), intent(in) :: self
      k = self%height
   end function get_height

   pure subroutine set_npnts(self, npnts)
      class(dendro), intent(inout)  :: self
      integer, intent(in)        :: npnts  
      self%npnts = npnts
   end subroutine set_npnts 

   pure subroutine set_npnts_array(self, labels)
      class(dendro), intent(inout)  :: self 
      integer, intent(in)           :: labels(:)     
      integer  :: npnts_array(self%n_trees), i 
      do i = 1, maxval(labels)
         npnts_array(i) = count(labels == i)
      end do 
      allocate(self%npnts_array(self%n_trees))
      self%npnts_array = npnts_array
   end subroutine set_npnts_array

   pure subroutine set_subsets(self, labels)
      class(dendro), intent(inout)  :: self 
      integer, intent(in)           :: labels(:) 
      integer, allocatable          :: tmp1(:), tmp2(:)
      integer     :: i
      allocate(self%subsets(size(self%npnts_array), maxval(self%npnts_array)))
      do i = 1, maxval(self%npnts_array) 
         allocate(tmp1(self%npnts_array(i)))
         allocate(tmp2( maxval(self%npnts_array) - size(tmp1)))
         tmp1 = pack(labels, labels == i)
         self%subsets(i,:) = [tmp1, tmp2]
         deallocate(tmp1, tmp2)  
      end do 
   end subroutine set_subsets

   pure subroutine set_distmat(self, dist_mat)
      class(dendro), intent(inout) :: self 
      real, intent(in)          :: dist_mat(:,:) 
      self%dist_mat = dist_mat
   end subroutine set_distmat

   pure subroutine set_medoids(self, centers)
      class(dendro), intent(inout)  :: self
      integer, intent(in)           :: centers(:)
      integer  :: i
      do i = 1, size(self%npnts_array)
         self%root_array(i)%ref_idx = centers(i)
      end do 
   end subroutine 

   subroutine build_dendro(self)
      class(dendro), intent(inout)   :: self
      real, allocatable    :: tmp(:)
      integer  :: i, j 
      ! Initialize Root of Entire Tree
      allocate(self%root)
      nullify(self%root%right, self%root%left)
      self%root%level = 0      
      allocate(self%root%o_subset(self%npnts))
      allocate(tmp(self%npnts))
      self%root%o_subset = [(i, i = 1,self%npnts )] 
      do j = 1, self%npnts
         tmp(j) = sum((self%dist_mat(j,:)))
      end do  
      self%root%ref_idx = maxloc(tmp, 1)
      deallocate(tmp)
      call clust_insert_s2_node(self%root, self%dist_mat, 0)
      ! Initialize each sub_tree 
      allocate(self%root_array(size(self%npnts_array)))
      do i = 1, size(self%npnts_array)
         allocate(self%root_array(i))
         self%root_array(i)%level = 0 
         allocate(self%root_array(i)%o_subset(self%npnts_array(i)))
         self%root_array(i)%o_subset = self%subsets(i,1:self%npnts_array(i))
         call clust_insert_s2_node(self%root_array(i), self%dist_mat, 0)
      end do 
      contains 
         recursive subroutine clust_insert_s2_node(root, dist_mat, level)
            type(s2_node), intent(inout) :: root
            real, intent(in)             :: dist_mat(:,:)
            integer, intent(in)          :: level
            integer                    :: l, idx, new_med(2), closest, sec_closest, n, m 
            real                       :: d, d1, d2 
            integer, allocatable       :: new_subset(:)
            root%level = level
            if(size(root%o_subset) < 2) return
            ! identify closest and 2nd closest to root
            do l = 1, size(root%o_subset)
               idx = root%o_subset(l)
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
            end do 

            ! remove those from subset
            allocate(new_subset(size(root%o_subset) - 2))
            n = 0 
            do m = 1, size(root%o_subset)
               if (root%o_subset(m) /= closest .and. root%o_subset(m) /= sec_closest) then
                  n = n + 1
                  new_subset(n) = root%o_subset(m)
               end if
            end do

            ! push 2nd closest left 
            allocate(root%left)
            nullify(root%left%left, root%left%right)
            allocate(root%left%o_subset(size(root%o_subset) - 2))
            root%left%o_subset = new_subset
            root%left%ref_idx = sec_closest
            root%left%level = level + 1 
            root%left%visit  = .false.
            root%left%is_pop = .true.
            
            ! push closest right 
            allocate(root%right)
            nullify(root%right%left, root%right%right)
            allocate(root%right%o_subset(size(root%o_subset) - 2))
            root%right%o_subset = new_subset
            root%right%ref_idx = closest
            root%right%level  = level + 1 
            root%right%visit  = .false.
            root%right%is_pop = .true.

            ! cleanup
            ! deallocate(sub_distmat, all_indxs, left_idx, right_idx)

            call clust_insert_s2_node(root%left,  dist_mat, level + 1)
            call clust_insert_s2_node(root%right, dist_mat, level + 1)

         end subroutine clust_insert_s2_node

   end subroutine build_dendro  

   recursive subroutine print_dendro(root)
      type(s2_node), pointer  :: root
      if(associated(root)) then 
         print *, root%o_subset, '//', root%ref_idx
         call print_dendro(root%right)
         call print_dendro(root%left)
      end if 
   end subroutine print_dendro 

   recursive subroutine walk_dendro(root, indxs, objs, done)
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
   end subroutine walk_dendro

end module simple_tree
