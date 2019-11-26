! executes the parallel (or distributed workflows) of SIMPLE
program single_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface, list_single_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_commander_base, only: execute_commander
use simple_commander_sim,  only: simulate_atoms_commander
use simple_commander_cluster2D
use simple_commander_tseries
use simple_spproj_hlev
implicit none
#include "simple_local_flags.inc"
type(tseries_import_commander)               :: xtseries_import
type(tseries_import_particles_commander)     :: xtseries_import_particles
type(tseries_ctf_estimate_commander)         :: xtseries_ctf_estimate
type(tseries_make_pickavg_commander)         :: xtseries_make_pickavg
type(tseries_motion_correct_commander_distr) :: xmcorr_distr
type(tseries_track_commander_distr)          :: xtrack_distr
type(center2D_nano_commander_distr)          :: xcenter2D_distr
type(cluster2D_nano_commander_hlev)          :: xcluster2D_distr
type(estimate_diam_commander)                :: xestimate_diam
type(simulate_atoms_commander)               :: xsimulate_atoms
type(refine3D_nano_commander_distr)          :: xrefine3D_distr
! OTHER DECLARATIONS
character(len=STDLEN) :: args, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
! parse command line
call get_command_argument(1, args, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(args, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, args, pos )
prg = args(pos+1:) ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_single_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! set global defaults
if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
select case(prg)
    ! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
    case( 'tseries_import' )
       call xtseries_import%execute(cline)
    case( 'tseries_import_particles' )
       call xtseries_import_particles%execute(cline)
   case( 'tseries_make_pickavg')
       call xtseries_make_pickavg%execute(cline)
   case( 'tseries_ctf_estimate' )
       call xtseries_ctf_estimate%execute(cline)
    case( 'motion_correct, tseries_motion_correct' )
        call xmcorr_distr%execute( cline )
    case( 'track', 'tseries_track' )
        call xtrack_distr%execute( cline )
    case( 'center2D', 'center2D_nano' )
        call xcenter2D_distr%execute(cline)
    case( 'cluster2D', 'cluster2D_nano' )
        call xcluster2D_distr%execute(cline)
    case( 'estimate_diam')
        call xestimate_diam%execute(cline)
    case( 'simulate_atoms' )
        call xsimulate_atoms%execute(cline)
    case( 'refine3D', 'refine3D_nano')
        call execute_commander(xrefine3D_distr, cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
end program single_exec
