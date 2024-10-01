program test_dock_vols
use simple_dock_vols,  only: dock_vols
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
implicit none

character(len=*), parameter :: VREF = 'vol1.mrc', VTARG = 'vol2.mrc', VTARG_DOCKED = 'vol2_docked.mrc'
real,             parameter :: SMPD = 1.2156, HP = 100., LP = 8.0, MSKDIAM = 175.
integer,          parameter :: NTHR = 0
type(dock_vols)  :: dvols
type(parameters) :: params
type(cmdline)    :: cline

call cline%set('nthr', NTHR)
call params%new(cline)
call dvols%new(VREF, VTARG, SMPD, HP, LP, MSKDIAM)
call dvols%srch_rots
call dvols%rotate_target(VTARG, VTARG_DOCKED)



end program test_dock_vols
