program simple_test_ctf
#include "simple_lib.f08"
use simple_ctf,   only: ctf
use simple_image, only: image
use simple_rnd,   only: ran3
use simple_ctf_estimate
use simple_timer
implicit none
! erravg(microns):    3.46422195E-04
! errmax(microns):    1.91450119E-04
! time(s)        :    22.322969166000000
type(ctf)               :: tfun
type(image)             :: img, img_msk
real                    :: dfx_found, dfy_found, angast_found, cc, err
real                    :: errmax, erravg, dfx_ran, phshift_found, ctfres
integer, parameter      :: BOX=512, NTST=5
real,    parameter      :: SMPD=1.26, KV=300., CS=2.0, AC=0.1, DFX=2.23
real,    parameter      :: DFY=2.21, ANGAST=30., HPLIM=20.0, LPLIM=5.0
integer                 :: itst
integer(timer_int_kind) :: tctf
call img%new([BOX,BOX,1], SMPD)
call img_msk%new([BOX,BOX,1], SMPD)
tfun   = ctf(SMPD, KV, CS, AC)
tctf   = tic()
errmax = 0.
erravg = 0.
! do itst=1,NTST
!     print *, 'test: ', itst
!     dfx_ran = 0.5 + ran3() * 4.5
!     call tfun%ctf2pspecimg(img, dfx_ran, dfx_ran, 0.)
!     call ctf_estimate_init(img, SMPD, KV, CS, AC, [0.5,5.0], [HPLIM,LPLIM], 'no' )
!     call ctf_estimate_srch( dfx_found, dfy_found, angast_found, phshift_found, cc, 'diag'//int2str(itst)//'.mrc' )
!     print *, 'ctfres: ', ctfres
!     err = abs(dfx_found - dfx_ran)
!     if( err > errmax ) errmax = err
!     erravg = erravg + err
! end do
print *, 'erravg(microns): ', erravg
print *, 'errmax(microns): ', errmax
print *, 'time(s)        : ', toc(tctf)
end program simple_test_ctf
