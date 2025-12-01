program simple_test_nano_mask
include 'simple_lib.f08'
use simple_image
use simple_image_msk, only: automask2D
implicit none
! constants
character(len=*), parameter :: STK='selected.spi'
real,             parameter :: SMPD=0.358
integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
! variables
type(image),    allocatable :: imgs(:)
real,           allocatable :: diams(:), shifts(:,:)
integer                     ::  n, i, ldim(3)
! read images
call find_ldim_nptcls(string(STK), ldim, n)
allocate(imgs(n))
do i = 1, n
    call imgs(i)%new(ldim, SMPD)
    call imgs(i)%read(string(STK), i)
end do
! mask
call automask2D(imgs, NGROW, WINSZ, EDGE, diams, shifts)
end program simple_test_nano_mask
