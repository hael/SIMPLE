program simple_test_lpstages
use simple_core_module_api
implicit none
#include "simple_local_flags.inc"
integer, parameter :: LDIM(3)=[256,256,1], BOX=LDIM(1), FILTSZ=BOX/2, NSTAGES=10
real,    parameter :: SMPD=1.3, LPSTART_DEFAULT=20., LPSTART_LB=10., LPFINAL=6., TOL=1.e-5
real               :: frc(FILTSZ) = 1.
type(lp_crop_inf)  :: lpinfo(NSTAGES)
type(lp_crop_inf)  :: lpinfo_single(1)
integer            :: i
call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, lpinfo, l_cavgs=.false.)
if( abs(lpinfo(NSTAGES)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages final stage is not lpfinal')
call lpstages_fast(BOX, 1, SMPD, LPSTART_LB, LPFINAL, lpinfo_single)
if( abs(lpinfo_single(1)%lp - LPSTART_DEFAULT) > TOL ) THROW_HARD('lpstages_fast single stage ignores start floor')
call lpstages_setlims(BOX, 1, SMPD, LPSTART_DEFAULT, LPFINAL, lpinfo_single)
if( abs(lpinfo_single(1)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages_setlims single stage is not lpstop')
call lpstages_fast(BOX, NSTAGES, SMPD, LPSTART_LB, LPFINAL, lpinfo)
if( abs(lpinfo(1)%lp - LPSTART_DEFAULT) > TOL ) THROW_HARD('lpstages_fast does not enforce start floor')
if( abs(lpinfo(NSTAGES)%lp - LPFINAL) > TOL ) THROW_HARD('lpstages_fast does not preserve lpstop')
if( any(lpinfo(2:NSTAGES)%lp >= lpinfo(1:NSTAGES-1)%lp) ) THROW_HARD('lpstages_fast does not march toward lpstop')
frc = [(1.0 - 0.02 * real(i - 1), i=1,FILTSZ)]
call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_LB, LPFINAL, lpinfo, &
    &l_cavgs=.false.)
if( abs(lpinfo(1)%lp - LPSTART_DEFAULT) > TOL ) THROW_HARD('lpstages does not enforce lpstart floor')
if( any(lpinfo(2:NSTAGES)%lp > lpinfo(1:NSTAGES-1)%lp) ) THROW_HARD('lpstages does not march toward lpfinal')
end program simple_test_lpstages
