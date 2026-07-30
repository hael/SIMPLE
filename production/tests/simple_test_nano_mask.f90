program simple_test_nano_mask
use simple_core_module_api
use simple_image,      only: image
use simple_image_msk,  only: automask2D
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
! constants
character(len=*), parameter :: DEFAULT_STK='selected.spi'
real,             parameter :: DEFAULT_SMPD=0.358, DEFAULT_MSKDIAM=100.
integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
! variables
type(cmdline)               :: cline
type(parameters)            :: params
type(image),    allocatable :: imgs(:)
real,           allocatable :: diams(:), shifts(:,:)
integer                     ::  n, i, ldim(3)
! setup parameters
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_nano_mask stk=<images> smpd=<pixel size> mskdiam=<diameter>'
    write(logfhandle,'(a)') 'No input keywords provided; running the default test with selected.spi'
    call cline%set('stk',      DEFAULT_STK)
    call cline%set('smpd',     DEFAULT_SMPD)
    call cline%set('mskdiam',  DEFAULT_MSKDIAM)
else
    call cline%parse_oldschool
endif
call cline%set('amsklp',  20.)
call cline%set('automsk', 'no')
call cline%set('ngrow',   real(NGROW))
call cline%set('winsz',   real(WINSZ))
call cline%set('edge',    real(EDGE))
call cline%set('part',    1.)
call cline%checkvar('stk',     1)
call cline%checkvar('smpd',    2)
call cline%checkvar('mskdiam', 3)
call cline%check
call params%new(cline)
! read images
call find_ldim_nptcls(params%stk, ldim, n)
allocate(imgs(n))
do i = 1, n
    call imgs(i)%new(ldim, params%smpd)
    call imgs(i)%read(params%stk, i)
end do
! mask
call automask2D(params, imgs, NGROW, WINSZ, EDGE, diams, shifts)
end program simple_test_nano_mask
