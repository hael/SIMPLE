! Order 2D class average images by resolution
program simple_test_order2D_by_resol
include 'simple_lib.f08'
use simple_image,            only: image
use simple_srch_sort_loc,    only: hpsort
implicit none
#include "simple_local_flags.inc"
character(*),parameter           :: fn_cavgs="cavgs.mrc",fn_cavgs_even="cavgs_even.mrc", fn_cavgs_odd="cavgs_odd.mrc"
real, parameter                  :: smpd=0.358
real                             :: res_fsc05, res_fsc0143
type(image)                      :: even, odd
real,allocatable                 :: res(:), corrs(:)
integer                          :: ldim(3), nptcls, nptcls1, nptcls2, iptcl, nyq
call find_ldim_nptcls( fn_cavgs, ldim, nptcls )
!write(logfhandle,'(a,1x,i6)') "Number of 2D images in the stack", nptcls
call find_ldim_nptcls( fn_cavgs_even, ldim, nptcls1 )
!write(logfhandle,'(a,1x,3(i6))') "ldim", ldim(:)
call find_ldim_nptcls( fn_cavgs_odd, ldim, nptcls2 )
!write(logfhandle,'(a,1x,3(i6))') "ldim", ldim(:)
!write(logfhandle,'(a,1x,i6)') "Number of 2D images in the even stack", nptcls1
!write(logfhandle,'(a,1x,i6)') "Number of 2D images in the odd stack", nptcls2
if( nptcls1 /= nptcls2 ) THROW_HARD("different number of images in the Stacks") 
! Construct 2D stack images
ldim(3) = 1
!write(logfhandle,'(a,1x,3(i6))') "ldim", ldim(:)
call even%new(ldim, smpd) 
call odd%new(ldim, smpd) 
do iptcl = 1, nptcls1
      !write(logfhandle,*) 'Particle # ', iptcl
      call even%read( fn_cavgs_even, iptcl )
      call  odd%read( fn_cavgs_odd, iptcl )
      ! forward FT
      call even%fft()
      call odd%fft()
      ! calculate FSC
      res = even%get_res()
      nyq = even%get_filtsz()        
      !allocate(corrs(nyq))
      if( .not. allocated(corrs) ) allocate(corrs(nyq))
      call even%fsc(odd, corrs)
      call even%ifft
      call odd%ifft
      call get_resolution(corrs, res, res_fsc05, res_fsc0143)
      !write(*, *) 'Comparing clean particle vs clean particle...'
      !write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
      !write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143  
      write(logfhandle,*) 'Particle # ', iptcl, res_fsc05, res_fsc0143
end do
!call hpsort(res_fsc0143, order)
!call reverse(order)
!do i = 1, nptcls
!      write(logfhandle,*) 'Particle # ', order(i), res_fsc05(order(i)), res_fsc0143(order(i))
!end do
call even%kill() 
call odd%kill() 
end program simple_test_order2D_by_resol
