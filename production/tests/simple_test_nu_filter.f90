program simple_test_nu_filter
use simple_core_module_api
use simple_image,     only: image
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN)         :: even_file, odd_file, out_even_file, out_odd_file
character(len=:), allocatable :: smpd_char
type(image)                   :: vol_even, vol_odd, vol_even_nu, vol_odd_nu
integer                       :: ldim_even(3), ldim_odd(3), ifoo, slen
real                          :: smpd, t_start, t_end
if( command_argument_count() /= 3 .and. command_argument_count() /= 5 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_nu_filter even.mrc odd.mrc smpd [out_even.mrc out_odd.mrc]'
    write(logfhandle,'(a)') 'even.mrc : input even volume (.mrc/.spi)'
    write(logfhandle,'(a)') 'odd.mrc  : input odd volume (.mrc/.spi)'
    write(logfhandle,'(a)') 'smpd    : sampling distance in Angstrom per voxel'
    write(logfhandle,'(a)') '[out_even.mrc out_odd.mrc] : optional output volume paths'
    stop
endif
call get_command_argument(1, even_file)
call get_command_argument(2, odd_file)
call get_command_argument(3, length=slen)
allocate(character(slen) :: smpd_char)
call get_command_argument(3, smpd_char)
read(smpd_char, *) smpd
if( command_argument_count() == 5 )then
    call get_command_argument(4, out_even_file)
    call get_command_argument(5, out_odd_file)
else
    out_even_file = 'nu_filtered_even.mrc'
    out_odd_file  = 'nu_filtered_odd.mrc'
endif
call find_ldim_nptcls(string(trim(even_file)), ldim_even, ifoo)
call find_ldim_nptcls(string(trim(odd_file )), ldim_odd,  ifoo)
if( ldim_even(3) <= 1 .or. ldim_odd(3) <= 1 )then
    THROW_HARD('simple_test_nu_filter expects 3D input volumes')
endif
if( any(ldim_even /= ldim_odd) )then
    THROW_HARD('simple_test_nu_filter expects even/odd volumes with identical dimensions')
endif
!$ nthr_glob = omp_get_max_threads()
!$ call omp_set_num_threads(nthr_glob)
call vol_even%new(ldim_even, smpd)
call vol_odd%new(ldim_odd,  smpd)
call vol_even%read(string(trim(even_file)))
call vol_odd%read(string(trim(odd_file)))
call vol_even%set_smpd(smpd)
call vol_odd%set_smpd(smpd)
call cpu_time(t_start)
call setup_nu_dmats(vol_even, vol_odd)
call optimize_nu_cutoff_finds()
call nu_filter_vols(vol_even_nu, vol_odd_nu)
call cpu_time(t_end)
call vol_even_nu%write(string(trim(out_even_file)))
call vol_odd_nu%write(string(trim(out_odd_file)))
write(logfhandle,'(a,f10.3,a)') 'nonuniform filtering elapsed CPU time: ', t_end - t_start, ' s'
write(logfhandle,'(a,a)') 'wrote even output volume: ', trim(out_even_file)
write(logfhandle,'(a,a)') 'wrote odd output volume: ', trim(out_odd_file)
call vol_even_nu%kill
call vol_odd_nu%kill
call vol_even%kill
call vol_odd%kill
call cleanup_nu_filter()
end program simple_test_nu_filter