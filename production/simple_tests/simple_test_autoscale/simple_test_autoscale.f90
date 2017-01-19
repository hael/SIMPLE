program simple_test_autoscale
use simple_math, only: autoscale
implicit none
real,    parameter :: smpd_in = 2.43, msk_in = 36.0
integer, parameter :: box_in = 128
integer            :: box_new
real               :: smpd_new, scale, msk_new
call autoscale( box_in, msk_in, smpd_in, box_new, msk_new, smpd_new, scale )
print *, 'box_new: ',  box_new
print *, 'msk_new: ',  msk_new
print *, 'smpd_new: ', smpd_new
print *, 'scale: ', scale
end program simple_test_autoscale
