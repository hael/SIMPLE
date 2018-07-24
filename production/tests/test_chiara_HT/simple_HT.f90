module Ht_project_mod
  include 'simple_lib.f08'
  use simple_image, only : image
  use simple_ctf,   only : ctf
  use canny_no_thresh_mod
  implicit none
  type ::three_vect
  real, allocatable:: a(:)
  real, allocatable:: b(:)
  real, allocatable:: rot(:)
  end type three_vect
  contains

  subroutine print_mat(matrix)
    integer, intent(in) :: matrix(:,:)
    integer             :: j, s(2)
    s = shape(matrix)
    do j = 1, s(2)
      print *, matrix(:,j)
    enddo
  end subroutine print_mat

  subroutine build_ellipse(self,center, axis, rot)
    type(image), intent(inout) :: self
    real, intent(in)           :: center(2), axis(2), rot  !USE RADIANS
    real, allocatable          :: mat(:,:,:), teta(:)
    integer                    :: ldim(3), i, j, k
   allocate(teta(360))
   do i = 1,360
     teta(i) = deg2rad(real(i))
   end do
   ldim = self%get_ldim()
   allocate(mat(ldim(1),ldim(2),1), source = self%get_rmat())
   do k = 1,size(teta)
   do i = 1,ldim(1)
   do j = 1,ldim(2)
      if(abs(i-center(1)-axis(1)*cos(teta(k))*cos(rot)+axis(2)*sin(teta(k))*sin(rot))<1 .and. &
       & abs(j-center(1)-axis(1)*cos(teta(k))*sin(rot)-axis(2)*sin(teta(k))*cos(rot))<1) then
         mat(i,j,1) = 1.
      end if
    enddo
    enddo
    enddo
    call self%set_rmat(mat)
   end subroutine build_ellipse

  function discret_param_space(ranges, num_steps) result(vectors)
   real, dimension(2,3), intent(in)  :: ranges
   integer                           :: num_steps(3)
   real                              :: s(3)
   type(three_vect)                  :: vectors
   integer                           :: i
    allocate(vectors%a(num_steps(1)+1), vectors%b(num_steps(2)+1), vectors%rot(num_steps(3)+1))
    s = 0.
    do i =1,3
       s(i) = (ranges(2,i)-ranges(1,i)) / (real(num_steps(i)))
     enddo
    do i = 1, num_steps(1)+1
      vectors%a(i) = ranges(1,1)+(i*s(1))-s(1)
    end do
    do i = 1, num_steps(2)+1
     vectors%b(i) = ranges(1,2)+(i*s(2))-s(2)
    end do
    do i = 1, num_steps(3)+1
     vectors%rot(i) = ranges(1,3)+(i*s(3))-s(3)
    end do
    end function discret_param_space
end module Ht_project_mod

  program Ht_project
  use Ht_project_mod
  use canny_no_thresh_mod
  implicit none
  integer                 :: N, Na,Nb                                !< number of points and discretization parameters
  integer                 :: i,j,k,l                                      !< to indicize loops
  real                    :: Aa, Bb                                   !< optimal parameters
  real, allocatable       :: a(:),b(:), dist(:)                     !< discretized parameters
  type(image)             :: img_e                     !< images
  integer, allocatable    :: H(:,:)
  integer, allocatable    :: x_inds(:), y_inds(:)
  integer                 :: num_steps(3), ldim(3), max_pos(2)
  real,    allocatable    :: rmat(:,:,:), mat_HT(:,:,:)
  type(three_vect)        :: vectors
  integer                 :: box, cnt, l_sig(1)
  real                    :: ranges(2,3), center(2), b1
  !STEP 0, READING AND PREPROCESSING
  box  = 256
  ldim =       [box,box,1]
  allocate(rmat(box,box,1))
  call img_e%new(ldim,1.)
  center = [real(floor(real(box/2))), real(floor(real(box/2)))]
  call build_ellipse(img_e, center , axis = [50.,60.], rot = 0.)
  rmat = img_e%get_rmat()
  N    = count(rmat > 0.5)          !number of points of interest, whites are 1
  allocate(x_inds(N), y_inds(N), source = 0)
  cnt  = 0                          !points of interest, img --> cartesian plane
  do i = 1,box
  do j = 1,box
    if( rmat(i,j,1) > 0.5 )then
      cnt  = cnt + 1
      x_inds(cnt) = i
      y_inds(cnt) = j
    endif
  enddo
  enddo
  call img_e%vis
    !STEP 1, PARAMETER SPACE DISCETIZATION
  ranges    = reshape([40.,70.,50.,80.,-pi/2.,pi/2.],[2,3])
  num_steps = [100, 100, 30]
  vectors   = discret_param_space(ranges, num_steps)
  a         = vectors%a
  b         = vectors%b
  Na        = size(a)                          !< save them just for comfort
  Nb        = size(b)
  !STEP 2, ACCUMULATOR MATRIX AND HOUGH-LOOP
  allocate(H(Na,Nb), source = 0 )
  allocate(dist(Nb), source = 0.)
  do i=1, N
 		 do j=1, Na
 					 if( 0. <= &
 										 &  (a(j)**2*(x_inds(i)*sin(0.)+y_inds(i)*cos(0.)-center(1))**2) &
 										 & /(a(j)**2-(x_inds(i)*cos(0.)-y_inds(i)*sin(0.)-center(2))**2) &
 								& .and. (a(j)**2-(x_inds(i)*cos(0.)-y_inds(i)*sin(0.)-center(2))**2) /= 0.) then
 							 b1        = &
 												 & sqrt(abs((a(j)**2*(x_inds(i)*sin(0.)+y_inds(i)*cos(0.)-center(2))**2)/ &
 												 &          (a(j)**2-(x_inds(i)*cos(0.)-y_inds(i)*sin(0.)-center(2))**2)))
 							 if(b(1) <= b1 .and. b1 <= b(Nb)) then
 										do l= 1, Nb
 											dist(l) = abs(b1-b(l));
 										enddo
 										l_sig = minloc(dist)
 										H(j,l_sig) = H(j,l_sig)+1
 							endif
 					 endif
 		 enddo
 enddo
   max_pos = maxloc(H)
   Aa = a(max_pos(1))
   Bb = b(max_pos(2))
   print *, "max position: ", max_pos
   print *, "Optimal parameters:"
   print *, "a = ", Aa
   print *, "b= ", Bb
   print *, "Max reached ", H(max_pos(1), max_pos(2))
   print *, "Number of voting points ", N
   print *, "Number of votes ",        sum(H)
   call build_ellipse(img_e,[real(floor(real(box/2))), real(floor(real(box/2)))],[Aa,Bb],rot = 0.)
   call img_e%vis
 end program Ht_project
