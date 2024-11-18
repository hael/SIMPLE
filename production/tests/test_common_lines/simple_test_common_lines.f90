program simple_test_common_lines
use simple_stack_io,   only: stack_io
use simple_image,      only: image
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
implicit none
type(stack_io)   :: stkio_r1, stkio_r2
type(image)      :: img1, img2
type(builder)    :: build
type(parameters) :: params
type(cmdline)    :: cline
character(len=*), parameter :: stkname1 = 'reprojs_fine.mrcs'
character(len=*), parameter :: stkname2 = 'reprojs.mrcs'
real,             parameter :: smpd     = 1.72
integer :: nptcls1, nptcls2, iptcl, ldim(3)
call stkio_r1%open(stkname1, smpd, 'read', bufsz=100)
call stkio_r2%open(stkname2, smpd, 'read', bufsz=100)
nptcls1 = stkio_r1%get_nptcls()
nptcls2 = stkio_r2%get_nptcls()
ldim    = stkio_r1%get_ldim()
call img1%new(ldim, smpd)
call img2%new(ldim, smpd)
! read
call stkio_r1%read(iptcl, img1)
do iptcl = 1, nptcls2
    call stkio_r2%read(iptcl, img2)
    ! investigating img1 vs all img2 in stkio_r2
end do
call stkio_r1%close
call stkio_r2%close
end program simple_test_common_lines