module simple_tree  
use simple_srch_sort_loc
implicit none 
#include "simple_local_flags.inc"

type :: s1_node
   type(s1_node), pointer  :: left => null()
   type(s1_node), pointer  :: right => null()
   type(s1_node), pointer  :: parent => null()
   real     :: inpl_ang
   logical  :: visited = .true. 
   integer, allocatable :: inpl_subset(:) ! array of inpl angles below s1 node
end type s1_node

type  :: s2_node 
   type(s2_node), pointer   :: left => null()
   type(s2_node), pointer   :: right => null()
   type(s2_node), pointer   :: parent => null()
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
   integer  :: height
   integer  :: npnts
   real, allocatable :: dist_mat(:,:) 
   ! maybe store in queue/stack ... 
contains 
   ! Setters / Getters
   procedure   :: get_height
   procedure   :: set_npnts 
   procedure   :: set_distmat
   ! Constructors / Methods 
   procedure   :: build_dendro
end type dendro 

! might need backtracking 
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

   pure subroutine set_distmat(self, dist_mat)
      class(dendro), intent(inout) :: self 
      real, intent(in)          :: dist_mat(:,:) 
      self%dist_mat = dist_mat
   end subroutine set_distmat

   subroutine build_dendro(self)
      class(dendro), intent(inout)   :: self
      real, allocatable    :: tmp(:)
      integer  :: i
      ! Initialize Root of Entire Tree
      allocate(self%root)
      nullify(self%root%right, self%root%left)
      self%root%level = 0      
      allocate(self%root%o_subset(self%npnts))
      allocate(tmp(self%npnts))
      self%root%o_subset = [(i, i = 1,self%npnts )] 
      do i = 1, self%npnts
         tmp(i) = sum((self%dist_mat(i,:)))
      end do  
      self%root%ref_idx = maxloc(tmp, 1)
      ! self%root%ref_idx is just exemplar for that cluster 
      ! find closest and 2nd closest ref in dist_mat(exemplar, :), but only consider within cluster 
      ! just track with subset. 
      ! remove left and right from the subset 
      deallocate(tmp)
      ! calculate s2_node psi, theta
      ! loop over distmat and find the o^th column which maximizes sum FM 
      ! psi, theta from o
      call clust_insert_s2_node(self%root, self%dist_mat, 0)

      contains 
         recursive subroutine clust_insert_s2_node(root, dist_mat, level)
            type(s2_node), pointer, intent(inout) :: root
            real, intent(in)                   :: dist_mat(:,:)
            integer, intent(in)                :: level
            integer                    :: sub_dim, i_dim, j_dim, i
            integer                    :: l, r
            integer                    :: new_med(2)
            real, allocatable          :: sub_distmat(:,:)
            integer, allocatable       :: all_indxs(:), left_idx(:), right_idx(:), all_indxs_c(:)
            root%level = level
            sub_dim = size(root%o_subset)
            if (sub_dim <= 1) return
            allocate(sub_distmat(sub_dim, sub_dim))
            ! Could probably just use upper part
            do i_dim = 1, sub_dim
               do j_dim = 1, sub_dim
                  sub_distmat(i_dim, j_dim) = dist_mat(root%o_subset(i_dim), root%o_subset(j_dim))
               end do
               sub_distmat(i_dim, i_dim) = 0.0
            end do
            
            new_med = maxloc(sub_distmat)
   
            l = ceiling(sub_dim / 2.)
            r = floor(sub_dim / 2.)
            allocate(all_indxs(sub_dim))
            do i = 1, sub_dim
               all_indxs(i) = root%o_subset(i)
            end do
            allocate(all_indxs_c(size(all_indxs)))
            all_indxs_c = all_indxs
            call hpsort(sub_distmat(:, new_med(1)), all_indxs)

            allocate(left_idx(l), right_idx(r))
            
            do i = 1, l
               left_idx(i) = all_indxs(i)
            end do
            do i = 1, r
               right_idx(i) = all_indxs(l + i)
            end do

            allocate(root%left)
            nullify(root%left%left, root%left%right)
            allocate(root%left%o_subset(l))
            root%left%o_subset = left_idx
            root%left%ref_idx = all_indxs_c(new_med(1))
            root%left%level = level + 1 
            root%left%visit  = .false.
            root%left%is_pop = .true.

            allocate(root%right)
            nullify(root%right%left, root%right%right)
            allocate(root%right%o_subset(r))
            root%right%o_subset = right_idx
            root%right%ref_idx = all_indxs_c(new_med(2))
            root%right%level  = level + 1 
            root%right%visit  = .false.
            root%right%is_pop = .true.

            deallocate(sub_distmat, all_indxs, left_idx, right_idx)

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
