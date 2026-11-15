!@descr: shared low-level definitions for the FLEX heterogeneity workflow
module simple_defs_flex
use simple_defs, only: dp
implicit none
private

public :: FLEX_MAX_BW_GROW, FLEX_ACCUM_BYTE_BUDGET, FLEX_FSC_SIGNAL_THRESHOLD
public :: COV_LP_OVER_NYQUIST, COV_LP_DEFAULT, COV_MINBOX, FLEX_PCG_MAXITS_DEFAULT, FLEX_PCG_RTOL_DEFAULT

integer,  parameter :: FLEX_MAX_BW_GROW = 4
real(dp), parameter :: FLEX_ACCUM_BYTE_BUDGET = 8.0d9
real(dp), parameter :: FLEX_FSC_SIGNAL_THRESHOLD = 0.143_dp

! flex_pca working sampling. 2*smpd_crop is the working Nyquist; the 1.25 safety factor keeps the
! covariance band clear of it
real,     parameter :: COV_LP_OVER_NYQUIST = 2.5
real,     parameter :: COV_LP_DEFAULT      = 16.0 !< default lp (A); smpd_target = lp / COV_LP_OVER_NYQUIST; also the flex_pca UI default
integer,  parameter :: COV_MINBOX          = 64   !< smallest working box

!> flex_pca rec_backend=pcg defaults (also the flex_pca UI defaults): fixed 4-iteration budget of the
!! warm-started basis M-step (rtol=0 disarms the rtol and FLEX_PCG_XTOL stops). State maps are
!! reconstructed by the reconstruction service (reconstruct3D) with the run's maxits_pcg.
integer,  parameter :: FLEX_PCG_MAXITS_DEFAULT = 4
real,     parameter :: FLEX_PCG_RTOL_DEFAULT   = 0.0

end module simple_defs_flex
