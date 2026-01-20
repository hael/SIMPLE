program simple_test_ctf
!$ use omp_lib
!$ use omp_lib_kinds
use simple_core_module_api
use simple_image, only: image
use simple_ctf,   only: ctf, memoize4ctf_apply, unmemoize4ctf_apply
implicit none
#include "simple_local_flags.inc"
integer, parameter :: LDIM(3) = [256,256,1]
real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.0, ANGAST = 0., KV = 300., CS = 2.0, AC = 0.1
type(image)        :: img, img_spec
type(ctf)          :: tfun
call img%new(LDIM, SMPD)
call img_spec%new(LDIM, SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call memoize4ctf_apply(img)
call tfun%ctf2img(img, DFX, DFY, ANGAST)
call img%ft2img('real', img_spec)
call img_spec%write(string('ctfimg.mrc'))
call unmemoize4ctf_apply
end program simple_test_ctf
