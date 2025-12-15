module simple_hierch_tree  
use simple_srch_sort_loc
implicit none 
#include "simple_local_flags.inc"
type  :: node 
   type(node), pointer   :: left => null()
   type(node), pointer   :: right => null()
   integer     :: medoid_ref_idx 
   real        :: medoid_fm
   real        :: F_ptcl_ref2med  
   real        :: F_ptcl2ref  ! Only need in end nodes. 
   integer, allocatable     :: subset(:) ! array of refs in node 
   integer     :: level
   logical     :: visit = .false.
   logical     :: is_pop   = .false.
   integer     :: indx
end type node

type  :: dendro
   type(node), pointer :: root 
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
   procedure   :: walk_dendro ! probably need a queue for this. 
end type dendro 

! might need backtracking 
contains 

   ! getters / setters 
   pure elemental integer function get_height(self) result(k)
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
      self%root%indx = 0     
      allocate(self%root%subset(self%npnts))
      self%root%subset = [(i, i = 1,self%npnts )] 
      allocate(tmp(self%npnts))
      do i = 1, self%npnts
         tmp(i) = sum(self%dist_mat(:, i)) / self%npnts
      end do  
      self%root%medoid_fm = maxval(tmp)
      deallocate(tmp)

      call clust_insert_node(self%root, self%dist_mat, 0)
      ! call print_dendro(self%root)

      contains 
         recursive subroutine clust_insert_node(root, dist_mat, level)
            type(node), pointer, intent(inout) :: root
            real, intent(in)                   :: dist_mat(:,:)
            integer, intent(in)                :: level
            integer                    :: sub_dim, i_dim, j_dim, i
            integer                    :: n1, n2
            integer                    :: new_med(2)
            real, allocatable          :: sub_distmat(:,:)
            integer, allocatable       :: all_indxs(:), left_idx(:), right_idx(:)
            real                       :: medoid_fm1, medoid_fm2
            root%level = level
            if (.not. allocated(root%subset)) return
            sub_dim = size(root%subset)
            if (sub_dim <= 1) return 
            allocate(sub_distmat(sub_dim, sub_dim))
            do i_dim = 1, sub_dim
               do j_dim = 1, sub_dim
                  sub_distmat(i_dim, j_dim) = dist_mat(root%subset(i_dim), root%subset(j_dim))
               end do
               sub_distmat(i_dim, i_dim) = 0.0
            end do
            ! change this to upper part 
            new_med = maxloc(sub_distmat)
   
            n1 = ceiling(sub_dim / 2.)
            n2 = floor(sub_dim / 2.)
            allocate(all_indxs(sub_dim))
            do i = 1, sub_dim
               all_indxs(i) = root%subset(i)
            end do
            call hpsort(sub_distmat(:, new_med(1)), all_indxs)

            allocate(left_idx(n1), right_idx(n2))
            ! need to make sure obj function values of right > left
            do i = 1, n1
               left_idx(i) = all_indxs(i)
            end do
            do i = 1, n2
               right_idx(i) = all_indxs(n1 + i)
            end do

            medoid_fm1 = sum(sub_distmat(:, new_med(1))) / real(n1)
            medoid_fm2 = sum(sub_distmat(:, new_med(2))) / real(n2)

            allocate(root%left)
            nullify(root%left%left, root%left%right)
            root%left%medoid_fm = medoid_fm1
            allocate(root%left%subset(n1))
            root%left%subset = left_idx
            root%left%indx   = 0
            root%left%visit  = .false.
            root%left%is_pop = .true.

            allocate(root%right)
            nullify(root%right%left, root%right%right)
            root%right%medoid_fm = medoid_fm2
            allocate(root%right%subset(n2))
            root%right%subset = right_idx
            root%right%indx   = 0
            root%right%visit  = .false.
            root%right%is_pop = .true.

            deallocate(sub_distmat, all_indxs, left_idx, right_idx)

            call clust_insert_node(root%left,  dist_mat, level + 1)
            call clust_insert_node(root%right, dist_mat, level + 1)

         end subroutine clust_insert_node

   end subroutine build_dendro  

   recursive subroutine print_dendro(root)
      type(node), pointer  :: root
      if(associated(root)) then 
         print *, root%subset
         call print_dendro(root%right)
         call print_dendro(root%left)
      end if 
   end subroutine print_dendro 

   subroutine walk_dendro(self, best_node)
      class(dendro), intent(inout)  :: self
      type(node), intent(out)       :: best_node 

   end subroutine walk_dendro

end module simple_hierch_tree
