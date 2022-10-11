program simple_test_nano_mask
include 'simple_lib.f08'
use simple_image
use simple_masker, only: automask2D
implicit none
! constants
character(len=*), parameter :: STK='selected.spi'
real,             parameter :: SMPD=0.358
integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
! variables
type(image), allocatable :: imgs(:)
real,        allocatable :: diams(:)
integer                  ::  n, i, ldim(3)
! read images
call find_ldim_nptcls(STK, ldim, n)
allocate(imgs(n))
do i = 1, n
    call imgs(i)%new(ldim, SMPD)
    call imgs(i)%read(STK, i)
end do
! mask
call automask2D(imgs, NGROW, WINSZ, EDGE, diams)
end program simple_test_nano_mask
