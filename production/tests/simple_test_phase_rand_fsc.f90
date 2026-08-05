program simple_test_phase_rand_fsc
use simple_core_module_api
use simple_image, only: image
use simple_image_bin, only: image_bin
use simple_fsc,   only: phase_rand_fsc
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"
integer, parameter :: BOX = 32
real,    parameter :: SMPD = 1.0
type(image) :: even, odd, envmask
type(image_bin) :: binary_mask
real :: voxels(BOX,BOX,BOX), expected, denom
real, allocatable :: fsc(:), fsc_masked(:), fsc_random_masked(:), fsc_unmasked(:)
integer :: i, j, k, n, k_rand, npix

do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            voxels(i,j,k) = sin(0.37 * real(i) + 0.19 * real(j) + 0.11 * real(k)) + &
                &0.25 * cos(1.31 * real(i) - 0.73 * real(j) + 0.47 * real(k))
        enddo
    enddo
enddo
call even%new([BOX,BOX,BOX], SMPD)
call odd%new([BOX,BOX,BOX], SMPD)
call even%set_rmat(voxels, .false.)
call odd%set_rmat(-voxels, .false.)
call envmask%disc([BOX,BOX,BOX], SMPD, 10., npix)
call binary_mask%transfer2bimg(envmask)
call binary_mask%cos_edge(3, envmask)
n = even%get_filtsz()
call phase_rand_fsc(even, odd, envmask, 1, n, fsc, fsc_masked, fsc_random_masked, fsc_unmasked)

if( any(.not.ieee_is_finite(fsc)) .or. any(.not.ieee_is_finite(fsc_masked)) .or. &
    &any(.not.ieee_is_finite(fsc_random_masked)) .or. any(.not.ieee_is_finite(fsc_unmasked)) )then
    THROW_HARD('phase-randomized FSC produced a non-finite value')
endif
k_rand = 0
do k = 2,n
    if( fsc_unmasked(k) < 0.8 )then
        k_rand = k
        exit
    endif
enddo
if( k_rand == 0 .or. k_rand > n-2 ) THROW_HARD('phase-randomized FSC test did not reach correction branch')
do k = 1,n
    if( k < k_rand + 2 )then
        expected = fsc_masked(k)
    else
        denom = 1. - fsc_random_masked(k)
        if( fsc_random_masked(k) > fsc_masked(k) .or. denom <= TINY )then
            expected = 0.
        else
            expected = max(0., min(1., (fsc_masked(k) - fsc_random_masked(k)) / denom))
        endif
    endif
    if( abs(fsc(k) - expected) > 1.e-5 ) THROW_HARD('phase-randomized FSC correction mismatch')
enddo

call even%kill
call odd%kill
call envmask%kill
call binary_mask%kill_bimg
deallocate(fsc, fsc_masked, fsc_random_masked, fsc_unmasked)
write(logfhandle,'(A)') 'Phase-randomized FSC solvent-correction test passed'

end program simple_test_phase_rand_fsc
