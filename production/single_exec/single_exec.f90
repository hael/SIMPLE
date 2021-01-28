! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
program single_exec
include 'simple_lib.f08'
use simple_user_interface,       only: make_user_interface, list_single_prgs_in_ui
use simple_cmdline,              only: cmdline, cmdline_err
use simple_commander_base,       only: execute_commander
use simple_commander_sim,        only: simulate_atoms_commander
use simple_commander_preprocess, only: map_cavgs_selection_commander
use simple_commander_imgproc,    only: estimate_diam_commander, pspec_int_rank_commander
use simple_commander_rec,        only: random_rec_commander_distr
use simple_commander_project
use simple_commander_cluster2D
use simple_commander_tseries
use simple_commander_oris
use simple_spproj_hlev
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

! RECONSTRUCTION PROGRAMS
type(tseries_ctf_estimate_commander)          :: xtseries_ctf_estimate
type(tseries_make_pickavg_commander)          :: xtseries_make_pickavg
type(tseries_motion_correct_commander_distr)  :: xmcorr_distr
type(tseries_track_particles_commander_distr) :: xtrack_distr
type(motion_refine_nano_commander)            :: xmotion_refine
type(graphene_subtr_commander)                :: xgraphene_subtr
type(center2D_nano_commander_distr)           :: xcenter2D_distr
type(cluster2D_nano_commander_hlev)           :: xcluster2D_distr
type(map_cavgs_selection_commander)           :: xmap_cavgs_selection
type(estimate_diam_commander)                 :: xestimate_diam
type(simulate_atoms_commander)                :: xsimulate_atoms
type(random_rec_commander_distr)              :: xrndrec
type(refine3D_nano_commander_distr)           :: xrefine3D_distr
type(initial_3Dmodel_nano_commander_distr)    :: xinitial_3Dmodel_nano_distr

! VALIDATION PROGRAMS
type(vizoris_commander)                       :: xvizoris
type(validate_nano_commander)                 :: xvalidate_nano

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
    case( 'tseries_import_particles' )
        call xtseries_import_particles%execute(cline)
    case( 'import_particles')
        call ximport_particles%execute(cline)
    case( 'prune_project' )
        call xprune_project%execute( cline )

    ! RECONSTRUCTION PROGRAMS
    case( 'tseries_make_pickavg')
        call xtseries_make_pickavg%execute(cline)
    case( 'tseries_ctf_estimate' )
        call xtseries_ctf_estimate%execute(cline)
    case( 'tseries_motion_correct' )
        call xmcorr_distr%execute( cline )
    case( 'tseries_track_particles' )
        call xtrack_distr%execute( cline )
    case( 'motion_refine_nano' )
        call xmotion_refine%execute( cline )
    case( 'graphene_subtr' )
        call xgraphene_subtr%execute( cline )
    case( 'center2D_nano' )
        call xcenter2D_distr%execute(cline)
    case( 'cluster2D_nano' )
        call xcluster2D_distr%execute(cline)
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case( 'estimate_diam')
        call cline%set('mkdir', 'no')
        call xestimate_diam%execute(cline)
    case( 'simulate_atoms' )
         call cline%set('mkdir', 'no')
        call xsimulate_atoms%execute(cline)
    case( 'random_rec')
        call xrndrec%execute(cline)
    case( 'refine3D_nano')
        call execute_commander(xrefine3D_distr, cline)
    case( 'initial_3Dmodel_nano')
        call xinitial_3Dmodel_nano_distr%execute(cline)

    ! VALIDATION PROGRAMS
    case( 'vizoris' )
        call xvizoris%execute(cline)
    case( 'validate_nano')
        call xvalidate_nano%execute(cline)

    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
end program single_exec
