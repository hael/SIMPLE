module cc_rot_mod
include 'simple_lib.f08'
use simple_image
use simple_ctf,  only : ctf
use simple_fftw3
implicit none
contains

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
cc   = img_ctf%ccf(img2)
rmat = img_ctf%get_rmat()
call cc%scale_pixels(10.,20.)
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
cc   = img_ctf%ccf(img2)
call cc%vis
rmat = img_ctf%get_rmat()
call cc%scale_pixels(10.,20.)
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
cc   = img_ctf%ccf(img2)
call cc%vis
rmat = img_ctf%get_rmat()
call cc%scale_pixels(10.,20.)
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
