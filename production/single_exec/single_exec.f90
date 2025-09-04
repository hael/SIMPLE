! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
program single_exec
include 'simple_lib.f08'
use simple_user_interface,       only: make_user_interface, list_single_prgs_in_ui
use simple_cmdline,              only: cmdline, cmdline_err
use simple_commander_sim,        only: simulate_atoms_commander
use simple_commander_preprocess, only: map_cavgs_selection_commander
use simple_commander_imgproc,    only: estimate_diam_commander
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
type(denoise_trajectory_commander)            :: xden_traj

! PARTICLE 3D RECONSTRUCTION PROGRAMS
type(analysis2D_nano_commander)               :: xanalysis2D_nano
type(center2D_nano_commander)                 :: xcenter2D
type(cluster2D_nano_commander)                :: xcluster2D
type(map_cavgs_selection_commander)           :: xmap_cavgs_selection
type(ppca_denoise_classes_commander)          :: xppca_denoise_classes
type(estimate_diam_commander)                 :: xestimate_diam
type(simulate_atoms_commander)                :: xsimulate_atoms
type(refine3D_nano_commander)                 :: xrefine3D_nano
type(extract_substk_commander)                :: xextract_substk
type(extract_subproj_commander)               :: xextract_subproj
type(autorefine3D_nano_commander)             :: xautorefine3D_nano
type(tseries_reconstruct3D_distr)             :: xtseries_reconstruct3D
type(tseries_core_finder_commander)           :: xtseries_core_finder
type(tseries_swap_stack_commander)            :: xtseries_swap_stack

! VALIDATION PROGRAMS
type(vizoris_commander)                       :: xvizoris
type(cavgsproc_nano_commander)                :: xcavgsproc
type(cavgseoproc_nano_commander)              :: xcavgseoproc
type(map2model_fsc_commander)                 :: xmap2model_fsc
type(map_validation_commander)                :: xmap_validation
type(model_validation_commander)              :: xmodel_validation
type(model_validation_eo_commander)           :: xmodel_validation_eo
type(ptclsproc_nano_commander)                :: xptclsproc

! MODEL BUILDING/ANALYSIS PROGRAMS
type(pdb2mrc_commander)                       :: xpdb2mrc
type(detect_atoms_commander)                  :: xdetect_atoms
type(conv_atom_denoise_commander)             :: xconv_atom_denoise
type(atoms_stats_commander)                   :: xatoms_stats
type(tseries_atoms_rmsd_commander)            :: xtseries_atoms_rmsd
type(tseries_core_atoms_analysis_commander)   :: xtseries_core_atoms_analysis
type(tseries_make_projavgs_commander)         :: xtseries_make_projavgs

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
    case( 'tseries_make_pickavg' )
        call xtseries_make_pickavg%execute(cline)
    case( 'tseries_motion_correct' )
        call xmcorr%execute( cline )
    case( 'tseries_track_particles' )
        call xtrack%execute( cline )
    case( 'graphene_subtr' )
        call cline%set('mkdir', 'no')
        call xgraphene_subtr%execute( cline )
    case( 'denoise_trajectory' )
        call cline%set('mkdir', 'no')
        call xden_traj%execute( cline )

    ! PARTICLE 3D RECONSTRUCTION PROGRAMS
    case( 'analysis2D_nano' )
        call xanalysis2D_nano%execute(cline)
    case( 'center2D_nano' )
        call xcenter2D%execute(cline)
    case( 'cluster2D_nano' )
        call xcluster2D%execute(cline)
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case( 'ppca_denoise_classes' )
        call xppca_denoise_classes%execute(cline)
    case( 'estimate_diam' )
        call cline%set('mkdir', 'no')
        call xestimate_diam%execute(cline)
    case( 'simulate_atoms' )
        call cline%set('mkdir', 'no')
        call xsimulate_atoms%execute(cline)
    case( 'refine3D_nano' )
        call xrefine3D_nano%execute(cline)
    case( 'extract_substk' )
        call xextract_substk%execute(cline)
    case( 'extract_subproj' )
        call xextract_subproj%execute(cline)
    case( 'autorefine3D_nano' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'autorefine3D_nano', 'single_exec')
        else
            call xautorefine3D_nano%execute(cline)
        endif
    case( 'tseries_reconstruct3D' )
        call xtseries_reconstruct3D%execute(cline)
    case( 'tseries_core_finder' )
        call xtseries_core_finder%execute(cline)
    case( 'tseries_swap_stack' ) 
        call xtseries_swap_stack%execute(cline)

    ! VALIDATION PROGRAMS
    case( 'vizoris' )
        call xvizoris%execute(cline)
    case( 'cavgsproc_nano' )
        call xcavgsproc%execute(cline)
    case( 'cavgseoproc_nano' )
        call xcavgseoproc%execute(cline)
    case( 'map2model_fsc' )
        call xmap2model_fsc%execute(cline)
    case( 'map_validation' )
        call xmap_validation%execute(cline)
    case( 'model_validation' )
        call xmodel_validation%execute(cline)
    case( 'model_validation_eo' )
        call xmodel_validation_eo%execute(cline)
    case( 'ptclsproc_nano' )
        call xptclsproc%execute(cline)

    ! MODEL BUILDING/ANALYSIS PROGRAMS
    case( 'pdb2mrc' )
        call cline%set('mkdir', 'no')
        call xpdb2mrc%execute(cline)
    case( 'conv_atom_denoise')
        call xconv_atom_denoise%execute(cline)
    case( 'detect_atoms' )
        call cline%set('mkdir', 'no')
        call xdetect_atoms%execute(cline)
    case( 'atoms_stats' )
        call cline%set('mkdir', 'yes')
        call xatoms_stats%execute(cline)
    case( 'tseries_atoms_rmsd' )
        call xtseries_atoms_rmsd%execute(cline)
    case( 'tseries_core_atoms_analysis' )
        call xtseries_core_atoms_analysis%execute(cline)
    case( 'tseries_make_projavgs' )
        call xtseries_make_projavgs%execute(cline)

    ! UNSUPPORTED
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
call simple_print_git_version('dc6dd592')
! end timer and print
rt_exec = toc(t0)
call simple_print_timer(rt_exec)
end program single_exec

