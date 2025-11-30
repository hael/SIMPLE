program simple_test_eo_diff
include 'simple_lib.f08'
use simple_image, only: image
implicit none
#include "simple_local_flags.inc"
character(len=*), parameter :: name_vol      = 'recvol_state01.mrc'
character(len=*), parameter :: name_vol_even = 'recvol_state01_even.mrc'
character(len=*), parameter :: name_vol_odd  = 'recvol_state01_odd.mrc'
integer,          parameter :: box           = 300
integer,          parameter :: ldim(3)       = [box,box,box]
real,             parameter :: smpd          = 1.2156
type(image)       :: vol, vol_even, vol_odd, vol_noise
integer           :: filtsz
real, allocatable :: res(:), corrs(:)
call vol%new(ldim, smpd)
call vol_even%new(ldim, smpd)
call vol_odd%new(ldim, smpd)
call vol%read(string(name_vol))
call vol_even%read(string(name_vol_even))
call vol_odd%read(string(name_vol_odd))
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
