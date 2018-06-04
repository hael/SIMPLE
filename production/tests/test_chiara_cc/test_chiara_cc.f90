module cc_rot_mod
include 'simple_lib.f08'
use simple_image
use simple_ctf,  only : ctf
use simple_fftw3
implicit none
contains

  subroutine build_ellipse(self,center, axis, rot)
      type(image), intent(inout) :: self
      real,        intent(in)    :: center(2), axis(2), rot
      real, allocatable          :: mat(:,:,:), teta(:)
      integer                    :: ldim(3), i, j, k
      real :: rot_r
      ! insert controls for rot etc
      rot_r = deg2rad(rot)
      allocate(teta(360))
      do i = 1,360
          teta(i) = deg2rad(real(i))
      end do
      ldim = self%get_ldim()
      allocate(mat(ldim(1),ldim(2),1), source = self%get_rmat())
      do k = 1,size(teta)
          do i = 1,ldim(1)
              do j = 1,ldim(2)
                  if(abs(i-center(1)-axis(1)*cos(teta(k))*cos(rot_r)+axis(2)*sin(teta(k))*sin(rot_r))<1 .and. &
                  & abs(j-center(1)-axis(1)*cos(teta(k))*sin(rot_r)-axis(2)*sin(teta(k))*cos(rot_r))<1) then
                      mat(i,j,1) = 1
                  end if
              enddo
          enddo
      enddo
      call self%set_rmat(mat)
  end subroutine build_ellipse

  function build_ctf(smpd, kv, cs, fraca, dfx, dfy, angast) result(img_ctf)
      integer, parameter     :: box=256
      real                   :: smpd, kv, cs, fraca, dfx, dfy, angast
      type(image)            :: ctf_image, img4viz, img_ctf
      type(ctf)              :: tfun
      call ctf_image%new([box,box,1], smpd)
      call img4viz%new  ([box,box,1], smpd)
      tfun = ctf(smpd, kv, cs, fraca)
      call tfun%ctf2img(ctf_image, dfx, dfy, angast)
      call ctf_image%ft2img('real', img4viz)
      img_ctf = img4viz
  end function build_ctf

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


  subroutine scale_img(self, new_min, new_max)
        type(image), intent(inout) :: self
        real, intent(in) :: new_min, new_max
        real :: old_min, old_max, sc
        real, allocatable :: rmat(:,:,:)
        integer ::  ldim(3)

        ldim = self%get_ldim()
        allocate(rmat(ldim(1),ldim(2),1), source = 0.)
        rmat = self%get_rmat()
        old_min = minval(rmat)
        old_max = maxval(rmat)
        sc = (new_max - new_min)/(old_max - old_min)
        rmat = sc*rmat+new_min-sc*old_min
        call self%set_rmat(rmat)
  end subroutine scale_img

  subroutine print_mat(matrix)
      real, intent(in) :: matrix(:,:,:)
      integer          :: j,s(3)
      s = shape(matrix)
      do j = 1, s(2)
          print *, matrix(:,j,1)
      enddo
  end subroutine print_mat
end module cc_rot_mod


program cc_rot
use cc_rot_mod
use simple_image
implicit none
type(image)       :: img2, cc, img_ctf                                   !original image, rotated image, cross correlation image
integer           :: box, ldim(3), m(3), window
real              :: center(2), axes1(2), axes2(2), axes3(2), rot      !parameters to build ellipses
real, allocatable :: rmat_out(:,:,:), cc_mat(:,:,:), w_mat(:,:,:), rmat(:,:,:)

box    = 256
window = 2      !anint(box/64)
ldim   = [box,box,1]
center = [real(floor(real(box/2))), real(floor(real(box/2)))]
axes1  = [55.,67.]
axes2  = [30.,70.]
axes3  = [45.,70.]
rot    = 0.

call img2%new    ( ldim,1. )
call cc%new      ( ldim,1. )
call img_ctf%new ( ldim,1. )
allocate(rmat(ldim(1), ldim(2),1), rmat_out(ldim(1), ldim(2),1), cc_mat(ldim(1), ldim(2),1), w_mat(window+1, window+1, 1), source = 0.)

!FIRST IMAGE
img2    = 0.
cc      = 0.
img_ctf = build_ctf(smpd=1.0, kv=300., cs=2.7, fraca=0.1, dfx=0., dfy=0., angast=0.)
call img_ctf%vis
call img_ctf%rtsq_serial( 90., 0., 0., rmat_out )
call img2%set_rmat(rmat_out)
!call img2%vis
cc   = img_ctf%ccf(img2)
call cc%vis
rmat = img_ctf%get_rmat()
call scale_img(cc, 10.,20.)
cc_mat       = cc%get_rmat()
w_mat(:,:,1) = cc_mat( box/2-window/2:box/2+window/2, box/2-window/2:box/2+window/2, 1 )
print *, "mean window:", (sum(w_mat)/size(w_mat))
print *, "mean tot:", sum(cc_mat)/size(cc_mat)
print *, "Difference: ", ((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))!/sum(cc_mat)
if(((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))>1.) then
  print *, "keep"                                                 !To make it indipendent from 0.1 maybe I should normalize
else
  print *, "trash"
endif
!SECOND IMAGE
img_ctf = 0.
img2    = 0.
cc      = 0.
img_ctf = build_ctf(smpd=1.0, kv=300., cs=2.7, fraca=0.1, dfx=3., dfy=1., angast=30.)
call img_ctf%vis
call img_ctf%rtsq_serial( 90., 0., 0., rmat_out )
call img2%set_rmat(rmat_out)
!call img2%vis
cc   = img_ctf%ccf(img2)
call cc%vis
rmat = img_ctf%get_rmat()
call scale_img(cc, 10.,20.)
cc_mat       = cc%get_rmat()
w_mat(:,:,1) = cc_mat( box/2-window/2:box/2+window/2, box/2-window/2:box/2+window/2, 1 )
print *, "mean window:", (sum(w_mat)/size(w_mat))
print *, "mean tot:", sum(cc_mat)/size(cc_mat)
print *, "Difference: ", ((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))!/sum(cc_mat)
if(((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))>1.) then
  print *, "keep"                                                 !To make it indipendent from 0.1 maybe I should normalize
else
  print *, "trash"
endif
!THIRD IMAGE
img_ctf = 0.
img2    = 0.
cc      = 0.
img_ctf = build_ctf(smpd=.5, kv=270., cs=2.7, fraca=0.1, dfx=0., dfy=0., angast=1.)
call img_ctf%vis
call img_ctf%rtsq_serial( 90., 0., 0., rmat_out )
call img2%set_rmat(rmat_out)
!call img2%vis
cc   = img_ctf%ccf(img2)
call cc%vis
rmat = img_ctf%get_rmat()
call scale_img(cc, 10.,20.)
cc_mat       = cc%get_rmat()
w_mat(:,:,1) = cc_mat( box/2-window/2:box/2+window/2, box/2-window/2:box/2+window/2, 1 )
print *, "mean window:", (sum(w_mat)/size(w_mat))
print *, "mean tot:", sum(cc_mat)/size(cc_mat)
print *, "Difference: ", ((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))!/sum(cc_mat)
if(((sum(w_mat)/size(w_mat))-(sum(cc_mat)/size(cc_mat)))>1.) then
  print *, "keep"                                                 !To make it indipendent from 0.1 maybe I should normalize
else
  print *, "trash"
endif
end program cc_rot
