program simple_test_lpstages
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"
integer, parameter :: LDIM(3)=[256,256,1], BOX=LDIM(1), FILTSZ=BOX/2, NSTAGES=10
real,    parameter :: SMPD=1.3, LPSTART_DEFAULT=20., LPSTART_LB=10., LPFINAL=6.
real               :: frc(FILTSZ) = 1.
type(lp_crop_inf)  :: lpinfo(NSTAGES)
call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, lpinfo, l_cavgs=.false.)
end program simple_test_lpstages
