program simple_test_ctf
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
use simple_ctf,   only: ctf
implicit none
#include "simple_local_flags.inc"

integer, parameter :: LDIM(3) = [256,256,1]
real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.0, ANGAST = 0., KV = 300., CS = 2.0, AC = 0.1, PHSH = 0.
type(image)        :: img, img_spec
type(ctf)          :: tfun
real, allocatable  :: tvals(:,:)
real    :: maxSpaFreqSq
integer :: logi_lims(3,2), fdims(3)

call img%new(LDIM, SMPD)
call img_spec%new(LDIM, SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call tfun%ctf2img(img, DFX, DFY, ANGAST, PHSH)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg.mrc')

call tfun%ctf_1stzero2img(img, DFX, DFY, ANGAST, PHSH)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stzero.mrc')

call tfun%ctf_1stpeak2img(img, DFX, DFY, ANGAST, PHSH)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stpeak.mrc')


logi_lims = img%loop_lims(2)
fdims     = img%get_array_shape()
allocate(tvals(fdims(1),fdims(2)))
call img%zero_and_flag_ft

img = cmplx(1.0,0.0)
call tfun%eval_and_apply(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, CTFLIMFLAG_PI, maxSpaFreqSq)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stzero_pi.mrc')
call tfun%eval_and_apply_before1stpeak(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, maxSpaFreqsq, CTFLIMFLAG_PI)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stzero_pi_recovered.mrc')

img = cmplx(1.0,0.0)
call tfun%eval_and_apply(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, CTFLIMFLAG_EL, maxSpaFreqSq)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stzero_el.mrc')
call tfun%eval_and_apply_before1stpeak(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, maxSpaFreqsq, CTFLIMFLAG_EL)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stzero_el_recovered.mrc')

img = cmplx(1.0,0.0)
call tfun%eval_and_apply(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, CTFLIMFLAG_PIO2, maxSpaFreqSq)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stpeak_pio2.mrc')
call tfun%eval_and_apply_before1stpeak(img, CTFFLAG_YES, logi_lims, fdims(1:2), tvals, DFX, DFY, ANGAST, PHSH, maxSpaFreqsq, CTFLIMFLAG_PIO2)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stpeak_pio2_recovered.mrc')

end program simple_test_ctf
