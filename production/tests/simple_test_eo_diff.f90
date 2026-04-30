!@descr: tests randomization of phases below noise power
program simple_test_eo_diff
use simple_core_module_api
use simple_image, only: image
use simple_refine3D_fnames, only: refine3D_state_halfvol_fname, refine3D_state_vol_fname
implicit none
#include "simple_local_flags.inc"
integer,          parameter :: box           = 300
integer,          parameter :: ldim(3)       = [box,box,box]
real,             parameter :: smpd          = 1.2156
type(string)       :: name_vol, name_vol_even, name_vol_odd
type(image)       :: vol, vol_even, vol_odd, vol_noise
integer           :: filtsz
real, allocatable :: res(:), corrs(:)
name_vol      = refine3D_state_vol_fname(1)
name_vol_even = refine3D_state_halfvol_fname(1, 'even')
name_vol_odd  = refine3D_state_halfvol_fname(1, 'odd')
call vol%new(ldim, smpd)
call vol_even%new(ldim, smpd)
call vol_odd%new(ldim, smpd)
call vol%read(name_vol)
call vol_even%read(name_vol_even)
call vol_odd%read(name_vol_odd)
call vol_noise%copy(vol_even)
call vol_noise%subtr(vol_odd)
call vol_noise%write(string('noisevol_state01.mrc'))
res = vol_even%get_res()
call vol_even%fft
call vol_odd%fft
call vol_noise%fft
filtsz = vol_even%get_filtsz()
allocate(corrs(filtsz))
call vol_even%ran_phases_below_noise_power(vol_odd)
call vol_even%ifft
call vol_even%write(string('new_impl.mrc'))
end program simple_test_eo_diff
