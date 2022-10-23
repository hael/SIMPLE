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
integer           :: k, filtsz
real              :: res_fsc05, res_fsc0143
real, allocatable :: noise_power(:), even_power(:), power(:), res(:), corrs(:)
call vol%new(ldim, smpd)
call vol_even%new(ldim, smpd)
call vol_odd%new(ldim, smpd)
call vol%read(name_vol)
call vol_even%read(name_vol_even)
call vol_odd%read(name_vol_odd)
call vol_noise%copy(vol_even)
call vol_noise%subtr(vol_odd)
call vol_noise%write('noisevol_state01.mrc')
! call vol_noise%spectrum('power', noise_power)
! call vol_even%spectrum('power',  even_power)
! call vol%spectrum('power',       power)
res = vol_even%get_res()
! do k=1,size(res)
!     print *, res(k), power(k), noise_power(k)
! end do
call vol_even%fft
call vol_odd%fft
call vol_noise%fft
filtsz = vol_even%get_filtsz()
allocate(corrs(filtsz))
! call vol_even%fsc(vol_odd, corrs)
! write(logfhandle,'(a)') 'ORIGINAL RESOLUTION ESTIMATION'
! do k=1,size(res)
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(k), '>>> CORRELATION:', corrs(k)
! end do
! call get_resolution(corrs, res, res_fsc05, res_fsc0143)
! write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
! write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
call vol_even%ran_phases_below_noise_power(vol_odd)
! call vol_odd%fcomps_below_noise_power_stats(vol_noise)
! call vol_even%fsc(vol_odd, corrs)
! write(logfhandle,'(a)') 'RESOLUTION ESTIMATION AFTER NOISE FILTERING'
! do k=1,size(res)
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(k), '>>> CORRELATION:', corrs(k)
! end do
! call get_resolution(corrs, res, res_fsc05, res_fsc0143)
! write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
! write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
call vol_even%ifft
call vol_even%write('new_impl.mrc')
end program simple_test_eo_diff
