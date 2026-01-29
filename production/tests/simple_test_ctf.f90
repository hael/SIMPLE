program simple_test_ctf
use simple_core_module_api
use simple_image,           only: image
use simple_ctf,             only: ctf
use simple_memoize_ft_maps, only: memoize_ft_maps
implicit none
#include "simple_local_flags.inc"
integer, parameter :: LDIM(3) = [256,256,1]
real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.0, ANGAST = 0., KV = 300., CS = 2.0, AC = 0.1
type(image)        :: img, img_spec
type(ctf)          :: tfun
call img%new(LDIM, SMPD)
call img_spec%new(LDIM, SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call memoize_ft_maps(LDIM, SMPD)
call img%ctf2img(tfun, DFX, DFY, ANGAST)
call img%ft2img('real', img_spec)
call img_spec%write(string('ctfimg.mrc'))
end program simple_test_ctf
