!@descr: shared low-level definitions for the FLEX heterogeneity workflow
module simple_defs_flex
use simple_defs, only: dp
implicit none
private

public :: FLEX_MAX_BW_GROW, FLEX_ACCUM_BYTE_BUDGET, FLEX_FSC_SIGNAL_THRESHOLD

integer,  parameter :: FLEX_MAX_BW_GROW = 4
real(dp), parameter :: FLEX_ACCUM_BYTE_BUDGET = 8.0d9
real(dp), parameter :: FLEX_FSC_SIGNAL_THRESHOLD = 0.143_dp

end module simple_defs_flex
