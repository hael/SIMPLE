program simple_test_mrc_validation
include 'simple_lib.f08' 
use simple_atoms, only: atoms
use simple_image
implicit none
character(len=STDLEN)         :: vol_file
character(len=:), allocatable :: smpd_char
type(image)                   :: vol 
real                          :: smpd
integer                       :: ldim(3), ifoo, slen
if( command_argument_count() /= 2 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mrc_validation vol.mrc smpd'
    write(logfhandle,'(a)') 'vol.mrc : volume' 
    write(logfhandle,'(a)') 'smpd    : SMPD value in Angstrom per voxel ' 
else
    call get_command_argument(1, vol_file)
    call get_command_argument(2, length=slen)
    allocate(character(slen) :: smpd_char)
    call get_command_argument(2, smpd_char)
    read(smpd_char, *) smpd
endif
call find_ldim_nptcls(string(trim(vol_file)),  ldim, ifoo)
print *, trim(vol_file), ldim
call vol%new(ldim, smpd)
call vol%read(string(trim(vol_file)))
call vol%write(string('vol_simple.mrc'))
call vol%kill
end program simple_test_mrc_validation
