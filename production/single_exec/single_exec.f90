! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
program single_exec
include 'simple_lib.f08'
use simple_user_interface,       only: make_user_interface, list_single_prgs_in_ui
use simple_cmdline,              only: cmdline, cmdline_err
use simple_commander_sim,        only: simulate_atoms_commander
use simple_commander_preprocess, only: map_cavgs_selection_commander
use simple_commander_imgproc,    only: estimate_diam_commander, pspec_int_rank_commander
use simple_commander_rec,        only: random_rec_commander_distr
use simple_exec_helpers
use simple_commander_project
use simple_commander_cluster2D
use simple_commander_tseries
use simple_commander_oris
use simple_commander_atoms
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(new_project_commander)                   :: xnew_project
type(update_project_commander)                :: xupdate_project
type(print_project_info_commander)            :: xprint_project_info
type(print_project_field_commander)           :: xprint_project_field
type(tseries_import_commander)                :: xtseries_import
type(import_particles_commander)              :: ximport_particles
type(tseries_import_particles_commander)      :: xtseries_import_particles
type(prune_project_commander_distr)           :: xprune_project

! TIME-SERIES PRE-PROCESSING PROGRAMS
type(tseries_make_pickavg_commander)          :: xtseries_make_pickavg
type(tseries_motion_correct_commander_distr)  :: xmcorr
type(tseries_track_particles_commander_distr) :: xtrack
type(graphene_subtr_commander)                :: xgraphene_subtr

! PARTICLE 3D RECONSTRUCTION PROGRAMS
type(analysis2D_nano_commander)               :: xanalysis2D_nano
type(center2D_nano_commander)                 :: xcenter2D
type(cluster2D_nano_commander)                :: xcluster2D
type(map_cavgs_selection_commander)           :: xmap_cavgs_selection
type(estimate_diam_commander)                 :: xestimate_diam
type(simulate_atoms_commander)                :: xsimulate_atoms
type(refine3D_nano_commander)                 :: xrefine3D_nano
type(autorefine3D_nano_commander)             :: xautorefine3D_nano
type(tseries_reconstruct3D_distr)             :: xtseries_reconstruct3D
type(tseries_swap_stack_commander)            :: xtseries_swap_stack

! VALIDATION PROGRAMS
type(vizoris_commander)                       :: xvizoris

! MODEL BUILDING/ANALYSIS PROGRAMS
type(detect_atoms_commander)                  :: xdetect_atoms
type(atoms_stats_commander)                   :: xatoms_stats
type(tseries_atoms_analysis_commander)        :: xtseries_atoms_analysis

! OTHER DECLARATIONS
character(len=STDLEN) :: args, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
integer(timer_int_kind)                     :: t0
real(timer_int_kind)                        :: rt_exec

! start timer
t0 = tic()
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
! generate script for queue submission?
call script_exec(cline, trim(prg), 'single_exec')
! set global defaults
if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')

select case(prg)

    ! PROJECT MANAGEMENT PROGRAMS
    case( 'new_project' )
        call xnew_project%execute(cline)
    case( 'update_project' )
        call xupdate_project%execute(cline)
    case( 'print_project_info' )
        call xprint_project_info%execute(cline)
    case( 'print_project_field' )
        call xprint_project_field%execute(cline)
    case( 'tseries_import' )
        call xtseries_import%execute(cline)
    case( 'import_particles')
        call ximport_particles%execute(cline)
    case( 'tseries_import_particles' )
        call xtseries_import_particles%execute(cline)
    case( 'prune_project' )
        call xprune_project%execute( cline )

    ! TIME-SERIES PRE-PROCESSING PROGRAMS
    case( 'tseries_make_pickavg')
        call xtseries_make_pickavg%execute(cline)
    case( 'tseries_motion_correct' )
        call xmcorr%execute( cline )
    case( 'tseries_track_particles' )
        call xtrack%execute( cline )
    case( 'graphene_subtr' )
        call cline%set('mkdir', 'no')
        call xgraphene_subtr%execute( cline )

    ! PARTICLE 3D RECONSTRUCTION PROGRAMS
    case( 'analysis2D_nano' )
        call xanalysis2D_nano%execute(cline)
    case( 'center2D_nano' )
        call xcenter2D%execute(cline)
    case( 'cluster2D_nano' )
        call xcluster2D%execute(cline)
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case( 'estimate_diam')
        call cline%set('mkdir', 'no')
        call xestimate_diam%execute(cline)
    case( 'simulate_atoms' )
        call cline%set('mkdir', 'no')
        call xsimulate_atoms%execute(cline)
    case( 'refine3D_nano')
        call xrefine3D_nano%execute(cline)
    case( 'autorefine3D_nano')
        call xautorefine3D_nano%execute(cline)
    case( 'tseries_reconstruct3D')
        call xtseries_reconstruct3D%execute(cline)
    case( 'tseries_swap_stack')
        call xtseries_swap_stack%execute(cline)

    ! VALIDATION PROGRAMS
    case( 'vizoris' )
        call xvizoris%execute(cline)

    ! MODEL BUILDING/ANALYSIS PROGRAMS
    case( 'detect_atoms' )
        call cline%set('mkdir', 'no')
        call xdetect_atoms%execute(cline)
    case( 'atoms_stats' )
        call cline%set('mkdir', 'yes')
        call xatoms_stats%execute(cline)
    case( 'tseries_atoms_analysis' )
        call xtseries_atoms_analysis%execute(cline)

    ! UNSUPPORTED
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
call simple_print_git_version('4d0e7a6')
! end timer and print
rt_exec = toc(t0)
call simple_print_timer(rt_exec)
end program single_exec
