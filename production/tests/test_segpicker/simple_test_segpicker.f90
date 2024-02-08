program simple_test_segpicker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,    only: image
use simple_tvfilter, only : tvfilter
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: mic       = '/home/elmlundho/cache/NanoX/NP_Minyoung_dry/intgs/FoilHole_6388490_Data_6383078_31_20240207_130359_EER_intg.mrc'
character(len=*), parameter :: pcontrast = 'black'
real,             parameter :: SHRINK    = 4.
real,             parameter :: lp        = 10.
real,             parameter :: lambda    = 3.
integer :: ldim(3), nframes, ldim_shrink(3)
real    :: smpd, smpd_shrink
type(image)    :: img_mic, img_mic_shrink
type(tvfilter) :: tvf

call find_ldim_nptcls(mic, ldim, nframes, smpd)
print *, 'ldim(1)    ', ldim(1)
print *, 'ldim(2)    ', ldim(2)
print *, 'ldim(3)    ', ldim(3)
print *, 'nframes    ', nframes
print *, 'smpd_found ', smpd

! shrink micrograph
ldim_shrink(1) = round2even(real(ldim(1))/SHRINK)
ldim_shrink(2) = round2even(real(ldim(2))/SHRINK)
ldim_shrink(3) = 1
smpd_shrink    = smpd * SHRINK
call img_mic%new(ldim, smpd)
call img_mic%read(mic)
call img_mic%mul(real(product(ldim))) ! to prevent numerical underflow when performing FFT
call img_mic%fft
call img_mic_shrink%new(ldim_shrink, smpd_shrink)
call img_mic_shrink%set_ft(.true.)
call img_mic%clip(img_mic_shrink)
select case(trim(pcontrast))
    case('black')
        ! flip contrast (assuming black particle contrast on input)
        call img_mic_shrink%mul(-1.)
    case('white')
        ! nothing to do
    case DEFAULT
        THROW_HARD('uknown pcontrast parameter, use (black|white)')
end select
call img_mic_shrink%bp(0., lp)
call img_mic_shrink%ifft
call img_mic_shrink%write('mic_shrink_lp.mrc')
! TV denoising
call tvf%new()
call tvf%apply_filter(img_mic_shrink, lambda)
call tvf%kill
call img_mic_shrink%write('mic_shrink_lp_tv.mrc')


end program simple_test_segpicker