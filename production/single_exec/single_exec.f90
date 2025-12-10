! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
program single_exec
include 'simple_lib.f08'
use simple_user_interface,        only: make_user_interface, list_single_prgs_in_ui
use simple_cmdline,               only: cmdline, cmdline_err
use simple_commanders_sim,        only: commander_simulate_atoms
use simple_commanders_preprocess, only: commander_map_cavgs_selection
use simple_commanders_imgproc,    only: commander_estimate_diam
use simple_exec_helpers
use simple_commanders_project
use simple_commanders_cluster2D
use simple_commanders_tseries
use simple_commanders_oris
use simple_commanders_atoms
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(commander_new_project)                   :: xnew_project
type(commander_update_project)                :: xupdate_project
type(commander_print_project_info)            :: xprint_project_info
type(commander_print_project_field)           :: xprint_project_field
type(commander_tseries_import)                :: xtseries_import
type(commander_import_particles)              :: ximport_particles
type(commander_tseries_import_particles)      :: xtseries_import_particles
type(commander_prune_project_distr)           :: xprune_project

! TIME-SERIES PRE-PROCESSING PROGRAMS
type(commander_tseries_make_pickavg)          :: xtseries_make_pickavg
type(commander_tseries_motion_correct_distr)  :: xmcorr
type(commander_tseries_track_particles_distr) :: xtrack
type(commander_graphene_subtr)                :: xgraphene_subtr
type(commander_denoise_trajectory)            :: xden_traj

! PARTICLE 3D RECONSTRUCTION PROGRAMS
type(commander_analysis2D_nano)               :: xanalysis2D_nano
type(commander_center2D_nano)                 :: xcenter2D
type(commander_cluster2D_nano)                :: xcluster2D
type(commander_map_cavgs_selection)           :: xmap_cavgs_selection
type(commander_ppca_denoise_classes)          :: xppca_denoise_classes
type(commander_estimate_diam)                 :: xestimate_diam
type(commander_simulate_atoms)                :: xsimulate_atoms
type(commander_refine3D_nano)                 :: xrefine3D_nano
type(commander_extract_substk)                :: xextract_substk
type(commander_extract_subproj)               :: xextract_subproj
type(commander_autorefine3D_nano)             :: xautorefine3D_nano
type(commander_tseries_reconstruct3D_distr)   :: xtseries_reconstruct3D
type(commander_tseries_core_finder)           :: xtseries_core_finder
type(commander_tseries_swap_stack)            :: xtseries_swap_stack

! VALIDATION PROGRAMS
type(commander_vizoris)                       :: xvizoris
type(commander_cavgsproc_nano)                :: xcavgsproc
type(commander_cavgseoproc_nano)              :: xcavgseoproc
type(commander_map2model_fsc)                 :: xmap2model_fsc
type(commander_map_validation)                :: xmap_validation
type(commander_model_validation)              :: xmodel_validation
type(commander_model_validation_eo)           :: xmodel_validation_eo
type(commander_ptclsproc_nano)                :: xptclsproc

! MODEL BUILDING/ANALYSIS PROGRAMS
type(commander_pdb2mrc)                       :: xpdb2mrc
type(commander_detect_atoms)                  :: xdetect_atoms
type(commander_conv_atom_denoise)             :: xconv_atom_denoise
type(commander_atoms_stats)                   :: xatoms_stats
type(commander_atoms_register)                :: xatoms_register
type(commander_crys_score)                    :: xcrys_score
type(commander_tseries_atoms_rmsd)            :: xtseries_atoms_rmsd
type(commander_tseries_core_atoms_analysis)   :: xtseries_core_atoms_analysis
type(commander_tseries_make_projavgs)         :: xtseries_make_projavgs

! OTHER DECLARATIONS
character(len=STDLEN)      :: args, prg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent

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
call script_exec(cline, string(trim(prg)), string('single_exec'))
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
        l_silent = .true.
    case( 'print_project_field' )
        call xprint_project_field%execute(cline)
        l_silent = .true.
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
            call restarted_exec(cline, string('autorefine3D_nano'), string('single_exec'))
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
    case( 'atoms_register' )
        call cline%set('mkdir', 'no')
        call xatoms_register%execute(cline)
    case( 'crys_score' )
        call cline%set('mkdir', 'no')
        call xcrys_score%execute(cline)
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
if( .not. l_silent )then
    call simple_print_git_version('94ac3e79')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program single_exec

