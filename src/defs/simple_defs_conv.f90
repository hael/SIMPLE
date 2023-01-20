module simple_defs_conv
! overlap limits
real, parameter :: OVERLAP_2D          = 0.80   ! based on Susan's observations on the small stuff, early stopping reduces overfitting
real, parameter :: OVERLAP_2D_FRAC     = 0.99   ! this limit is for when fractional update is on
real, parameter :: OVERLAP_2D_NANO     = 0.998
real, parameter :: OVERLAP_3D          = 0.99
real, parameter :: OVERLAP_STATE       = 0.98
real, parameter :: OVERLAP_STATE_JOINT = 0.96
real, parameter :: OVERLAP_STATE_HET   = 0.98
! fraction of search space scanned limits
real, parameter :: FRACSRCHSPACE_2D    = 90.0  ! based on Susan's observations on the small stuff, early stopping reduces overfitting
real, parameter :: FRACSRCHSPACE_3D    = 99.0
real, parameter :: FRACSRCHSPACE_FRAC  = 99.0  ! this limit is for when fractional update is on
real, parameter :: FRACSRCHSPACE_HET   = 99.0
! other limits
real, parameter :: MSK_FRAC            = 0.07
real, parameter :: MINSHIFT            = 10.0
real, parameter :: MAXSHIFT            = 12.0
end module simple_defs_conv
