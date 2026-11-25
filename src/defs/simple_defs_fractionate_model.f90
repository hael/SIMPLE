!@descr: fractionate_model total-volume output and density-envelope settings
module simple_defs_fractionate_model
implicit none
private
public :: FRACTIONATE_DISTANCE_FBODY
public :: FRACTIONATE_TOTAL_DIR, FRACTIONATE_MASK_LP, FRACTIONATE_MASK_BINWIDTH, FRACTIONATE_MASK_EDGE

character(len=*), parameter :: FRACTIONATE_DISTANCE_FBODY = 'distance_to_surface'
character(len=*), parameter :: FRACTIONATE_TOTAL_DIR      = 'rec_frames_total'
real,             parameter :: FRACTIONATE_MASK_LP        = 15.
integer,          parameter :: FRACTIONATE_MASK_BINWIDTH  = 1, FRACTIONATE_MASK_EDGE = 1
end module simple_defs_fractionate_model
