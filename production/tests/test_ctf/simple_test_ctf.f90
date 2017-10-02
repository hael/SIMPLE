program simple_test_ctf
use simple_ctf,   only: ctf
use simple_image, only: image
implicit none
type(ctf)   :: tfun
type(image) :: img, img_msk
integer, parameter :: BOX=256
real,    parameter :: SMPD=1.26, KV=300., CS=2.0, AC=0.1, DFX=2.5, DFY=2.1, ANGAST=30., HPLIM=20.0, LPLIM=5.0
call img%new([BOX,BOX,1], SMPD)
call img_msk%new([BOX,BOX,1], SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call tfun%ctf2pspecimg(img, DFX, DFY, ANGAST)
call img%vis
call img_msk%resmsk(HPLIM, LPLIM)
call img%mul(img_msk)
call img%vis
end program simple_test_ctf
