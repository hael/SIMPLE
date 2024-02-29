program simple_test_radial_cc
include 'simple_lib.f08'
use simple_image,        only: image
implicit none
#include "simple_local_flags.inc"
real,    parameter  :: smpd=0.358 
type(image)                 :: img1, img2
real,    allocatable        :: rad_corrs(:), rad_dists(:), filt(:)
character(len=256)          :: fn_img1, fn_img2
integer                     :: ldim_refs(3), ifoo, ldim1(3)
if( command_argument_count() /= 2 )then
    write(logfhandle,'(a)')  'Usage: simple_test_radial_cc img1.mrc img2.mrc', NEW_LINE('a')
    stop
else
    call get_command_argument(1, fn_img1)
    call get_command_argument(2, fn_img2)
endif
call find_ldim_nptcls(fn_img1, ldim1, ifoo)
ldim_refs  = [ldim1(1), ldim1(2), ldim1(3)]
call img1%new(ldim_refs, smpd)
call img2%new(ldim_refs, smpd)
call img1%read(fn_img1)
call img2%read(fn_img2)
if( img1 .eqdims. img2 ) then
    continue
else
    THROW_HARD('Nonconforming dimensions of images')
endif
if( .not. img1%same_smpd(img2) ) THROW_HARD('Nonconforming smpd of images')
call img1%radial_cc(img2,rad_corrs, rad_dists)
! create output filter for mask
!   !> \brief  converts the FSC to the optimal low-pass filter
!    function fsc2optlp( corrs ) result( filt )
filt = fsc2optlp( rad_corrs )
! create a mask for the different 3D shells
!> \brief mask  is for spherical masking
    !! \param mskrad mask radius in pixels
    !! \param which mask type
    !! \param inner include cosine edge material
    !! \param width width of inner patch
!    subroutine mask( self, mskrad, which, inner, width, backgr )

! apply mask to the 3D volume img%mask
! create a mask using SIMPLE 
call img1%kill()
call img2%kill()
end program simple_test_radial_cc


