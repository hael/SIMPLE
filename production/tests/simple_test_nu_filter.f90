program simple_test_nu_filter
use simple_core_module_api
use simple_image,     only: image
use simple_nu_filter, only: setup_dmats
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN)         :: vol_file
character(len=:), allocatable :: smpd_char
type(image)                   :: vol
integer                       :: ldim(3), ifoo, slen
real                          :: smpd, t_start, t_end
if( command_argument_count() /= 2 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_nu_filter vol.mrc smpd'
    write(logfhandle,'(a)') 'vol.mrc : input volume (.mrc/.spi)'
    write(logfhandle,'(a)') 'smpd    : sampling distance in Angstrom per voxel'
    stop
endif
call get_command_argument(1, vol_file)
call get_command_argument(2, length=slen)
allocate(character(slen) :: smpd_char)
call get_command_argument(2, smpd_char)
read(smpd_char, *) smpd
call find_ldim_nptcls(string(trim(vol_file)), ldim, ifoo)
if( ldim(3) <= 1 )then
    THROW_HARD('simple_test_nu_filter expects a 3D input volume')
endif
call vol%new(ldim, smpd)
call vol%read(string(trim(vol_file)))
call vol%set_smpd(smpd)
call cpu_time(t_start)
call setup_dmats(vol, vol)
call cpu_time(t_end)
write(logfhandle,'(a,f10.3,a)') 'setup_dmats elapsed CPU time: ', t_end - t_start, ' s'
call vol%kill
end program simple_test_nu_filter