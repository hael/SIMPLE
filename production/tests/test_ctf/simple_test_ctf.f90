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

call img%new(LDIM, SMPD)
call img_spec%new(LDIM, SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call tfun%ctf2img(img, DFX, DFY, ANGAST, PHSH)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg.mrc')

call tfun%ctf_1stpeak2img(img, DFX, DFY, ANGAST, PHSH)
call img%ft2img('real', img_spec)
call img_spec%write('ctfimg_1stpeak.mrc')

end program simple_test_ctf
