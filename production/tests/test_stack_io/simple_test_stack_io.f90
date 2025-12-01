program simple_test_stack_io
include 'simple_lib.f08'
use simple_stack_io, only: stack_io
use simple_image,    only: image
implicit none
type(stack_io)              :: stkio_r, stkio_w
type(image)                 :: img
character(len=*), parameter :: stkname = 'cavgs_iter030_ranked.mrc'
real,             parameter :: smpd    = 1.3
integer                     :: nptcls, iptcl, ldim(3)
call stkio_r%open(string(stkname), smpd, 'read', bufsz=100)
call stkio_w%open(string('outstk_written.mrc'), smpd, 'write', box=256, is_ft=.false.)
nptcls = stkio_r%get_nptcls()
ldim   = stkio_r%get_ldim()
call img%new(ldim, smpd)
! read
do iptcl = 1, nptcls
    call stkio_r%read(iptcl, img)
    call img%write(string('outstk_read.mrc'), iptcl)
end do
call stkio_r%close
! write
do iptcl = 1, nptcls
    call img%read(string(stkname), iptcl)
    call stkio_w%write(iptcl, img)
end do
call stkio_w%close
! readwrite
call stkio_r%open(string(stkname), smpd, 'read', bufsz=100)
call stkio_w%open(string('outstk_read_written.mrc'), smpd, 'write', box=256, is_ft=.false.)
nptcls = stkio_r%get_nptcls()
ldim   = stkio_r%get_ldim()
call img%new(ldim, smpd)
do iptcl = 1, nptcls
    call stkio_r%read(iptcl, img)
    call stkio_w%write(iptcl, img)
end do
call stkio_r%close
call stkio_w%close
end program simple_test_stack_io
