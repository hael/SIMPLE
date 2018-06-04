module dont_mod
include 'simple_lib.f08'
use simple_image, only : image !Maybe I ll need more
use simple_ctf,   only : ctf
implicit none

type ::three_vect
    real, allocatable:: a(:)
    real, allocatable:: b(:)
    real, allocatable:: rot(:)
end type three_vect

type :: voting_pt
    real :: x, y                             ! coordinates
    integer, allocatable :: arr_voted(:,:)   ! to know if the point has voted a cell
end type voting_pt

contains

    subroutine print_mat(matrix)
        integer, intent(in) :: matrix(:,:)
        integer             :: j,s(2)
        s = shape(matrix)
        do j = 1, s(2)
            print *, matrix(:,j)
        enddo
    end subroutine print_mat

    subroutine build_ellipse(self,center, axis, rot)
        type(image), intent(inout) :: self
        real,        intent(in)    :: center(2), axis(2), rot      !Input in degrees
        real                       :: rot_r
        real, allocatable          :: mat(:,:,:), t(:)
        integer                    :: ldim(3), i, j, k
        ! insert controls for rot etc
        rot_r = deg2rad(rot)
        allocate(t(360))
        do i = 1,360
            t(i) = deg2rad(real(i))
        end do
        ldim = self%get_ldim()
        allocate(mat(ldim(1),ldim(2),1), source = 0.)
        do k = 1,size(t)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    if(abs(i-center(1)-axis(1)*cos(t(k))*cos(rot_r)+axis(2)*sin(t(k))*sin(rot_r))<1 .and. &
                    & abs(j-center(1)-axis(1)*cos(t(k))*sin(rot_r)-axis(2)*sin(t(k))*cos(rot_r))<1) then
                        mat(i,j,1) = 1
                    end if
                enddo
            enddo
        enddo
        call self%set_rmat(mat)
    end subroutine build_ellipse

    subroutine sobel(img_in,img_out,thresh)
        type(image), intent(inout)          :: img_in,img_out
        type(image)                         :: Gr, img_p          !gradient and padded images
        real, dimension(:,:,:), allocatable :: wc, wr             !row and column masks
        integer, parameter                  :: L = 3              !dimension of the mask
        integer, dimension(3)               :: ldim               !dimension of the image
        integer                             :: i,j,m,n            !loop indexes
        real, dimension(:,:,:), allocatable :: Dc,Dr              !Column and row derivates
        real, dimension(:,:,:), allocatable :: mat_in, mat_out    !images, just the matrix
        real, dimension(:,:,:), allocatable :: grad, mat_p        !gradient and padded image, just the matrix
        real                                :: thresh             !threshold for Sobel algorithm
        ldim = img_in%get_ldim()
        allocate(mat_in(ldim(1),ldim(2),1),grad(ldim(1),ldim(2),1),Dc(ldim(1),ldim(2),1), &
              & Dr(ldim(1),ldim(2),1),wc(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1),wr(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1), source = 0.)
        wc = (1./8.)*reshape([-1,0,1,-2,0,2,-1,0,1],[3,3,1])      !Sobel masks
        wr = (1./8.)*reshape([-1,-2,-1,0,0,0,1,2,1],[3,3,1])
        mat_in = img_in%get_rmat()
        call img_p%new([ldim(1)+L-1,ldim(2)+L-1,1],1.)
        call img_in%pad(img_p)                                    !padding
        mat_p = img_p%get_rmat()
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do m = -(L-1)/2,(L-1)/2
                    do n = -(L-1)/2,(L-1)/2
                        Dc(i,j,1) = Dc(i,j,1)+mat_p(i+m+1,j+n+1,1)*wc(m,n,1)
                        Dr(i,j,1) = Dr(i,j,1)+mat_p(i+m+1,j+n+1,1)*wr(m,n,1)
                    end do
                end do
            end do
        end do
        grad = sqrt(Dc**2 + Dr**2);   !gradient matrix
        call Gr%new(ldim,1.)          !gradient image
        call Gr%set_rmat(grad)
        allocate(mat_out(ldim(1),ldim(2),1), source=0.)
        where( grad > thresh ) mat_out = 1.
        call img_out%set_rmat(mat_out)
    end subroutine sobel

    function discret_param_space(ranges, num_steps) result(vectors)
        real, dimension(2,3), intent(in)  :: ranges
        integer                           :: num_steps(3)
        real                              :: s(3)
        type(three_vect)                  :: vectors
        integer                           :: i
        allocate(vectors%a(num_steps(1)+1), vectors%b(num_steps(2)+1), vectors%rot(num_steps(3)+1))
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

end module dont_mod

program I_dont
use dont_mod
implicit none
integer                      :: N, Na,Nb                              !< number of points and discretization parameters
integer                      :: i,j,k,d,l,p,q
integer                      :: max_pos(2)
real                         :: Aa, Bb                                !< optimal parameters
real,            allocatable :: a(:),b(:)                             !< discretized parameters
type(image)                  :: img_edge, img4viz, img_e              !< images
integer,         allocatable :: H(:,:)
integer,         allocatable :: x_inds(:), y_inds(:)
integer                      :: num_steps(3), ldim(3)
real,            allocatable :: rmat(:,:,:), mat_HT(:,:,:)
type(three_vect)             :: vectors
integer                      :: box, cnt
real                         :: center(2), axis(2), rot, a1, b1, ranges(2,3)
real,            allocatable :: f1(:), f2(:), g1(:), g2(:), delta(:), t(:)
type(voting_pt), allocatable :: pt(:)

!STEP 0, READING AND PREPROCESSING
box  = 256
ldim = [box,box,1]
allocate(rmat(box,box,1))
center = [real(floor(real(box/2))), real(floor(real(box/2)))]
axis   = [55.,55.]
rot = 0.
call img_e%new(ldim,1.)
call img_edge%new(ldim,1.)
call img4viz%new(ldim,1.)
call build_ellipse(img_e, center, axis, rot)
call sobel(img_e,img_edge,thresh = 0.06) !useless, already binary
rmat = img_edge%get_rmat()
N    = count(rmat > 0.5)                 !number of points of interest
allocate(x_inds(N), y_inds(N), source = 0)
cnt  = 0                                 !points of interest, img --> cartesian plane
do i = 1, box
    do j = 1, box
        if( rmat(i,j,1) > 0.5 )then
            cnt  = cnt + 1
            x_inds(cnt) = i
            y_inds(cnt) = j
        endif
    enddo
enddo
!call img4viz%set_rmat(rmat)

!STEP 1, PARAMETER SPACE DISCETIZATION
ranges    = reshape([50.,60.,50.,60.,0.,60.],[2,3])
num_steps = [50, 50, 100]
vectors   = discret_param_space(ranges, num_steps)
allocate(a(num_steps(1)+1), b(num_steps(2)+1), source = 0.)
a         = vectors%a
b         = vectors%b
Na        = size(a)                        !< save them just for comfort
Nb        = size(b)
img4viz   = 0.
allocate(pt(N))
do i = 1,N
    call img4viz%set([x_inds(i), y_inds(i),1], 1.)    !To "plot" the points
    pt(i)%x = x_inds(i)
    pt(i)%y = y_inds(i)
    allocate(pt(i)%arr_voted(Na,Nb), source = 0)
enddo
call img4viz%vis

!STEP 2, ACCUMULATOR MATRIX AND HOUGH-LOOP
allocate(H(Na,Nb), source = 0 )
allocate(t(3600),  source = 0.)        !parameter, this declaration can be better
do i = 1, 3600      !That s very important 4 the algorithm to work, but it influences the CPUtime a lot
    t(i) = deg2rad(0+i*0.01)
enddo
allocate(delta(size(t)), f1(size(t)), g1(size(t)), f2(size(t)), g2(size(t)),mat_HT(box,box,1), source = 0.)
mat_HT(50,50:65,1) = 1
mat_HT(65,50:65,1) = 1
mat_HT(50:65,50,1) = 1
mat_HT(50:65,65,1) = 1

!HOUGH LOOP
do i = 1, N                                !fixing one point of interest
    do d = 1, size(t)                      !fixing the value of the parameter
        f1(d)    = cos(t(d))*cos(rot)
        f2(d)    = -sin(t(d))*sin(rot)
        g1(d)    = cos(t(d))*sin(rot)
        g2(d)    = sin(t(d))*cos(rot)
        delta(d) = f1(d)*g2(d)-f2(d)*g1(d)
        if(abs(delta(d)) > 0.001) then     !Don t divide per 0
            a1 = ((pt(i)%x-center(1))*g2(d)-(pt(i)%y-center(2))*f2(d))/delta(d)
            b1 = ((pt(i)%y-center(2))*f1(d)-(pt(i)%x-center(1))*g1(d))/delta(d)
            !Plotting 2 random HT
            if(i == 1 .or. i == 70) then
                do p = 1, box
                    do q = 1, box
                        if(abs(p-a1)<1 .and. abs(q-b1)<1) mat_HT(p,q,1) = 1
                    end do
                end do
            end if
            !Voting
            do k = 1, Nb-1
                do l = 1, Na-1
                    if (b(k) <= b1 .and. b1 < b(k+1) .and. a(l) <= a1 .and. a1 < a(l+1)) then  !if (a1,b1) is in cell (l,k)
                        if(pt(i)%arr_voted(l,k) == 0) then                                     !if point i hasn't voted cell (l,k) yet
                            H(l,k) = H(l,k)+1                                                  !point i votes cell (l,k)
                            pt(i)%arr_voted(l,k) = 1                                           !remember that point i has voted the cell (l,k)
                        end if
                    end if
                enddo
            enddo
        endif
    enddo
enddo
call img4viz%set_rmat(mat_HT)
call img4viz%vis
max_pos = maxloc(H)
Aa = a(max_pos(1))
Bb = b(max_pos(2))
print *, "max position: ",   max_pos
print *, "Optimal parameters:"
print *, "a = ",             Aa
print *, "b = ",             Bb
print *, "Max reached ",     H(max_pos(1), max_pos(2))
print *, "Voting points",    N
print *, "Number of votes ", sum(H)
!call build_ellipse(img4viz,center,[Aa,Bb],rot = 0.)
!call img4viz%vis
!print *, "In the right position we have: ", H(21,21)
!call print_mat(H)
end program I_dont
