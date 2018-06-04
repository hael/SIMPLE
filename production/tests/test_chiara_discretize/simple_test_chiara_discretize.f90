module discret_mod
  use simple_image, only : image !Maybe I ll need more
  use simple_ctf,   only : ctf
  implicit none
  real, parameter :: pi = 3.14159265

  type ::three_vect
  real, allocatable:: a(:)
  real, allocatable:: b(:)
  real, allocatable:: rot(:)
  end type three_vect

  type :: voting_pt
  real :: x,y   !coordinates
  integer, allocatable :: arr_voted(:,:)
  end type voting_pt

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
    real, intent(in)           :: center(2), axis(2), rot
    real, allocatable          :: mat(:,:,:), teta(:)
    integer                    :: ldim(3), i, j, k
   allocate(teta(360))
   do i = 1,360
     teta(i) = i
   end do
   ldim = self%get_ldim()
   allocate(mat(ldim(1),ldim(2),1))
   mat  = 0.
   do k = 1,size(teta)
   do i = 1,ldim(1)
   do j = 1,ldim(2)
      if(abs(i-center(1)-axis(1)*cos(teta(k))*cos(rot)+axis(2)*sin(teta(k))*sin(rot))<1 .and. &
       & abs(j-center(1)-axis(1)*cos(teta(k))*sin(rot)-axis(2)*sin(teta(k))*cos(rot))<1) then
         mat(i,j,1) = 1
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
end module discret_mod

  program discret
  use discret_mod
  implicit none
  integer                 :: N, Na,Nb                                     !< number of points and discretization parameters
  integer                 :: i,j,k,d,l, p, q                              !< to indicize loops
  integer                 :: ind(1), max_pos(2)                           !< index
  real                    :: Aa, Bb, Tt                                   !< optimal parameters
  real, allocatable       :: a(:),b(:),teta(:), t(:)                      !< discretized parameters
  type(image)             :: img_edge, img4viz, img_e                     !< images
  integer, allocatable    :: H(:,:)
  integer, allocatable    :: x_inds(:), y_inds(:)
  integer, dimension(3)   :: num_steps, ldim
  real   , allocatable    :: rmat(:,:,:), mat_HT(:,:,:)
  type(three_vect)        :: vectors
  integer                 :: box, cnt, bool, si(2), pt_vt_c(3)
  real :: center(2), axis(2), rot, a1, b1, ranges(2,3), a2, b2
  real, dimension(:), allocatable :: f1, f2, g1, g2, delta
  type(voting_pt), allocatable :: pt(:)

  !STEP 0, READING AND PREPROCESSING
  box = 256
  ldim = [box,box,1]
  allocate(rmat(box,box,1))
  center = [real(floor(real(box/2))), real(floor(real(box/2)))]
  axis = [55.,55.]
  rot  = 0.
  call img_e%new(ldim,1.)
  call img_edge%new(ldim,1.)
  call img4viz%new(ldim,1.)
  call build_ellipse(img_e, center, axis, rot)
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
  call img4viz%set_rmat(rmat)
  call img4viz%vis
    !STEP 1, PARAMETER SPACE DISCETIZATION
  ranges    = reshape([50.,65.,50.,65.,0.,60.],[2,3])
  num_steps = [15, 15, 100]
  vectors   = discret_param_space(ranges, num_steps)
  a        = vectors%a
  b        = vectors%b
  Na       = size(a)                          !< save them just for comfort
  Nb       = size(b)

  img4viz = 0.                 !To "plot" the points
  allocate(pt(N))
  do i = 1,N
     call img4viz%set([x_inds(i), y_inds(i),1], 1.)
     pt(i)%x = x_inds(i)
     pt(i)%y = y_inds(i)
     allocate(pt(i)%arr_voted(Na,Nb), source = 0)
  enddo

  !STEP 2, ACCUMULATOR MATRIX AND HOUGH-LOOP
  allocate(H(Na,Nb),      source = 0 )
  allocate(t(3600),       source = 0.)        !parameter, this declaration can be better
  do i = 1, 3600      !That s very important 4 the algorithm to work, but it influences the CPUtime a lot
    t(i) = 0+i*0.1
  enddo
  allocate(delta(size(t)), f1(size(t)), g1(size(t)), f2(size(t)), g2(size(t)), source = 0.)
  allocate(mat_HT(box,box,1), source = 0.)
  mat_HT(50,50:65,1) = 1
  mat_HT(65,50:65,1) = 1
  mat_HT(50:65,50,1) = 1
  mat_HT(50:65,65,1) = 1

  !HOUGH LOOP

do d = 1, size(t)
a1 = 55+3*cos(t(d))
b1 = 55+3*sin(t(d))
a2 = 58+3*cos(t(d))
b2 = 58+3*sin(t(d))

!I am plotting 2 HT
do p = 1, box
  do q = 1, box
   if(abs(p-a1)<1 .and. abs(q-b1)<1) then
    mat_HT(p,q,1) = 1
  end if
   if(abs(p-a2)<1 .and. abs(q-b2)<1) then
   mat_HT(p,q,1) = 1
  end if
end do
end do

do k = 1, Nb-1
do l = 1, Na-1
  if ((b(k) <= b1 .and. b1 < b(k+1) .and. a(l) <= a1 .and. a1 < a(l+1))) then
    if(pt(1)%arr_voted(l,k) == 0) then !It means this point hasn't still voted this cell
    H(l,k) = H(l,k)+1
    pt(1)%arr_voted(l,k) = 1  !It means point i has voted the cell (l,k)
  end if
end if
 if ((b(k) <= b2 .and. b2 < b(k+1) .and. a(l) <= a2 .and. a2 < a(l+1))) then
  if(pt(70)%arr_voted(l,k) == 0) then !It means this point hasn't still voted this cell
    H(l,k) = H(l,k)+1
    pt(70)%arr_voted(l,k) = 1  !It means point i has voted the cell (l,k)
  end if
  end if
enddo
enddo
enddo


call img4viz%set_rmat(mat_HT)
call img4viz%vis
max_pos = maxloc(H)
Aa = a(max_pos(1)+1)   !+1 because of indexes
Bb = b(max_pos(2)+1)
print *, "max position: ", max_pos
print *, "Optimal parameters:"
print *, "a = ", Aa
print *, "b= ", Bb
print *, "Max reached ", H(max_pos(1), max_pos(2))
print *, "Voting points", N
print *, "Number of votes ", sum(H)
!call build_ellipse(img4viz,center,[Aa,Bb],rot = 0.)
!call img4viz%vis
!print *, "Points: "
!print *, pt(1)%x, pt(1)%y
!print *, pt(70)%x, pt(70)%y
!print *, "In the right position we have: ", H(5,5)
 !H(1,5) = 78
call print_mat(H)
 H(1,5) = 78
end program discret
