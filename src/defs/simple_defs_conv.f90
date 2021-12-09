module simple_defs_conv
real, parameter :: MI_CLASS_LIM_3D      = 0.996
real, parameter :: MI_CLASS_LIM_2D      = 0.80  ! based on Susan's observations on the small stuff, early stopping reduces overfitting
real, parameter :: MI_CLASS_LIM_2D_FRAC = 0.99  ! this limit is for when fractional update is on (does it happen in cluster2D???)
real, parameter :: MI_STATE_LIM         = 0.98
real, parameter :: MI_STATE_JOINT_LIM   = 0.96
real, parameter :: FRAC_LIM             = 90.0  ! based on Susan's observations on the small stuff, early stopping reduces overfitting
real, parameter :: FRAC_LIM_FRAC        = 99.0
real, parameter :: MSK_FRAC             = 0.06
real, parameter :: MINSHIFT             = 5.0
real, parameter :: MAXSHIFT             = 6.0
real, parameter :: HET_MI_STATE_LIM     = 0.98
real, parameter :: HET_FRAC_LIM         = 99.0
real, parameter :: MI_CLASS_LIM_2D_NANO = 0.998
end module simple_defs_conv
