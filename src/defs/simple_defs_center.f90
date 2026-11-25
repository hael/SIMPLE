!@descr: settings and file names of the particle-centering workflow (center_ptcls)
module simple_defs_center
implicit none

integer,          parameter :: CENTER_NSAMPLE     = 10000         !< the references' sample: enough particles for ten low-resolution classes
integer,          parameter :: CENTER_NCLS        = 10            !< ten references centred twice as many trpv1 picks within 3 A as the rotational average (Phase 1)
integer,          parameter :: CENTER_REF_MAXITS  = 15            !< the references' iterations: refine2D stops earlier on the sample (7 to 13), more changed nothing
real,             parameter :: CENTER_LP          = 20.           !< the references' band (A)
real,             parameter :: CENTER_LPSTOP      = 12.           !< the pass's band (A): finer than the references', which centred more particles
real,             parameter :: CENTER_SHIFT_SIGMA = 10.           !< the shift prior (A): a 3-sigma disk of 30 A covers the picking error's 95th percentile
character(len=*), parameter :: CENTER_SAMPLE_PROJ = 'center_sample.simple'
character(len=*), parameter :: CENTER_PASS_PROJ   = 'center_pass.simple'
character(len=*), parameter :: CENTER_PASS_DIR    = 'center_pass' !< the pass runs as a first iteration (every particle not yet centred) in its own folder

end module simple_defs_center
