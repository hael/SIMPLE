program simple_test_mask
include 'simple_lib.f08'
use simple_image,    only: image
use simple_binimage, only: binimage
implicit none
real,    parameter :: SMPD=0.356, MSK=56.
integer, parameter :: EDGE=6, BOX=160
type(image)        :: spher_msk
type(binimage)     :: spher_msk_bin
integer            :: npix
! make a spherical mask
call spher_msk%disc([BOX,BOX,BOX], SMPD, MSK, npix)
! transfer to binimg instance
call spher_msk_bin%transfer2bimg(spher_msk)
! apply soft cosine edge
call spher_msk_bin%cos_edge(EDGE)
! write
call spher_msk_bin%write_bimg('spherical_mask.mrc')
end program simple_test_mask
