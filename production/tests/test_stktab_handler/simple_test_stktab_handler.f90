program simple_test_stktab_handler
use simple_stktab_handler
use simple_image,          only: image
use simple_defs
implicit none

character(len=STDLEN)         :: ftab='filetable.txt'
type(image)                   :: img
type(stktab_handler)          :: stkhandle
character(len=:), allocatable :: stkname
integer                       :: nmics, nptcls, ldim(3), iptcl, ind

! call test_stktab_handler
! stop

call stkhandle%new(ftab)
nmics  = stkhandle%get_nmics()
nptcls = stkhandle%get_nptcls()
ldim   = stkhandle%get_ldim()
print *, 'nmics:  ', nmics
print *, 'nptcls: ', nptcls
print *, 'ldim:   ', ldim
call img%new(ldim, 1.0)
do iptcl=1,nptcls
	call stkhandle%get_stkname_and_ind(iptcl, stkname, ind)
	call img%read(stkname, ind)
	call img%write('stk_from_test_pimg_reader.mrc', iptcl)
	deallocate(stkname)
end do

end program simple_test_stktab_handler
