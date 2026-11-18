!@descr: settings of the dose-series workflow (dose_filter3D)
module simple_defs_dose
use simple_defs, only: dp
implicit none

real(dp), parameter :: DOSE_PLOT_RES(5) = [20._dp, 10._dp, 8._dp, 6._dp, 4._dp] !< resolutions (A) of the weights and step tables and figures

end module simple_defs_dose
