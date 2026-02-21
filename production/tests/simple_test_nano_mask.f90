program simple_test_nano_mask
use simple_core_module_api
use simple_image,      only: image
use simple_image_msk,  only: automask2D
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
! constants
character(len=*), parameter :: STK='selected.spi'
real,             parameter :: SMPD=0.358
integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
! variables
type(cmdline)               :: cline
type(parameters)            :: params
type(image),    allocatable :: imgs(:)
real,           allocatable :: diams(:), shifts(:,:)
integer                     ::  n, i, ldim(3)
! setup parameters
call cline%set('smpd',    SMPD)
call cline%set('amsklp',  20.)
call cline%set('msk',     100.)
call cline%set('automsk', 'no')
call cline%set('ngrow',   real(NGROW))
call cline%set('winsz',   real(WINSZ))
call cline%set('edge',    real(EDGE))
call cline%set('part',    1.)
call params%new(cline)
! read images
call find_ldim_nptcls(string(STK), ldim, n)
allocate(imgs(n))
do i = 1, n
    call imgs(i)%new(ldim, SMPD)
    call imgs(i)%read(string(STK), i)
end do
! mask
call automask2D(params, imgs, NGROW, WINSZ, EDGE, diams, shifts)
end program simple_test_nano_mask
