program simple_test_lpstages
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"

integer, parameter :: LDIM(3)=[256,256,1], BOX=LDIM(1), FILTSZ=BOX/2, NSTAGES=10
real,    parameter :: SMPD=1.3, LPLIMS(2)=[10.,6.], LPSTART_DEFAULT=20.
real               :: frc(FILTSZ) = 0.
type(lp_crop_inf)  :: lpinfo(NSTAGES)

call lpstages(BOX, NSTAGES, frc, SMPD, LPLIMS, LPSTART_DEFAULT, lpinfo, verbose=.true.)

end program simple_test_lpstages