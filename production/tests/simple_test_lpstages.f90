program simple_test_lpstages
use simple_core_module_api
implicit none
#include "simple_local_flags.inc"
integer, parameter :: LDIM(3)=[256,256,1], BOX=LDIM(1), FILTSZ=BOX/2, NSTAGES=10
real,    parameter :: SMPD=1.3, LPSTART_DEFAULT=20., LPSTART_LB=10., LPFINAL=6., TOL=1.e-5
real               :: frc(FILTSZ) = 1.
type(lp_crop_inf)  :: lpinfo(NSTAGES)
type(lp_crop_inf)  :: lpinfo_single(1)
call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, lpinfo, l_cavgs=.false.)
if( abs(lpinfo(NSTAGES)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages final stage is not lpfinal')
call lpstages_fast(BOX, 1, SMPD, LPSTART_DEFAULT, LPFINAL, lpinfo_single)
if( abs(lpinfo_single(1)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages_fast single stage is not lpstop')
call lpstages_setlims(BOX, 1, SMPD, LPSTART_DEFAULT, LPFINAL, lpinfo_single)
if( abs(lpinfo_single(1)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages_setlims single stage is not lpstop')
end program simple_test_lpstages
