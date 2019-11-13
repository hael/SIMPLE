! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface,list_shmem_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_spproj_hlev
use simple_commander_project
use simple_commander_checks
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_refine3D
use simple_commander_rec
use simple_commander_relion
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
use simple_commander_resolest
use simple_projection_frcs
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(new_project_commander)           :: xnew_project
type(update_project_commander)        :: xupdate_project
type(print_project_info_commander)    :: xprint_project_info
type(print_project_field_commander)   :: xprint_project_field
type(import_movies_commander)         :: ximport_movies
type(import_boxes_commander)          :: ximport_boxes
type(import_particles_commander)      :: ximport_particles
type(import_cavgs_commander)          :: ximport_cavgs
type(merge_stream_projects_commander) :: xmerge_stream_projects
type(replace_project_field_commander) :: xreplace_project_field
type(selection_commander)             :: xselection
type(export_relion_commander)         :: xexport_relion

! SINGLE-PARTICLE WORKFLOW PROGRAMS
type(map_cavgs_selection_commander) :: xmap_cavgs_selection
type(cluster_cavgs_commander)       :: xcluster_cavgs
type(write_classes_commander)       :: xwrite_classes
type(symaxis_search_commander)      :: xsymsrch
type(symmetry_test_commander)       :: xsymtst
type(radial_sym_test_commander)     :: xradsymtst
type(symmetrize_map_commander)      :: xsymmetrize_map
type(dock_volpair_commander)        :: xdock_volpair
type(postprocess_commander)         :: xpostprocess

! IMAGE PROCESSING PROGRAMS
type(pspec_stats_commander)   :: xpspecstats
type(mask_commander)          :: xmask
type(fsc_commander)           :: xfsc
type(local_res_commander)     :: xlocal_res
type(centervol_commander)     :: xcenter
type(reproject_commander)     :: xreproject
type(volops_commander)        :: xvolops
type(convert_commander)       :: xconvert
type(ctfops_commander)        :: xctfops
type(filter_commander)        :: xfilter
type(normalize_commander)     :: xnormalize
type(scale_commander)         :: xscale
type(stack_commander)         :: xstack
type(stackops_commander)      :: xstackops
type(shift_commander)         :: xshift
type(estimate_diam_commander) :: xestimate_diam

! ORIENTATION PROCESSING PROGRAMS
type(make_oris_commander) :: xmake_oris
type(orisops_commander)   :: xorisops
type(oristats_commander)  :: xoristats
type(vizoris_commander)   :: xvizoris

! PRINT INFO PROGRAMS
type(info_image_commander)        :: xinfo_image
type(info_stktab_commander)       :: xinfo_stktab
type(print_fsc_commander)         :: xprint_fsc
type(print_magic_boxes_commander) :: xprint_magic_boxes

! SIMULATOR PROGRAMS
type(simulate_noise_commander)       :: xsimulate_noise
type(simulate_particles_commander)   :: xsimulate_particles
type(simulate_movie_commander)       :: xsimulate_movie
type(simulate_subtomogram_commander) :: xsimulate_subtomogram
type(simulate_atoms_commander)       :: xsimulate_atoms

! TIME-SERIES (NANO-PARTICLE) PROGRAMS
type(tseries_import_commander)           :: xtseries_import
type(tseries_import_particles_commander) :: xtseries_import_particles
type(tseries_ctf_estimate_commander)     :: xtseries_ctf_estimate
type(detect_atoms_commander)             :: xdetect_atoms
type(atoms_rmsd_commander)               :: xatoms_rmsd
type(radial_dependent_stats_commander)   :: xradial_dependent_stats
type(atom_cluster_analysis_commander)    :: xatom_cluster_analysis
type(nano_softmask_commander)            :: xnano_softmask

! SYSTEM INTERACTION PROGRAMS
type(mkdir_commander) :: xmkdir

! OTHER DECLARATIONS
character(len=STDLEN) :: xarg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_shmem_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse

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
    case( 'import_movies' )
        call ximport_movies%execute(cline)
    case( 'import_boxes' )
        call ximport_boxes%execute(cline)
    case( 'import_particles' )
        call ximport_particles%execute(cline)
    case( 'import_cavgs' )
        call ximport_cavgs%execute(cline)
    case( 'merge_stream_projects' )
        call xmerge_stream_projects%execute(cline)
    case( 'replace_project_field' )
        call xreplace_project_field%execute(cline)
    case( 'selection', 'report_selection' )
        call xselection%execute(cline)
    case( 'export_relion' )
        call xexport_relion%execute(cline)

    ! SINGLE-PARTICLE WORKFLOW PROGRAMS
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case('cluster_cavgs')
        call xcluster_cavgs%execute(cline)
    case('write_classes')
        call xwrite_classes%execute(cline)
    case( 'symaxis_search' )
        call xsymsrch%execute( cline )
    case( 'symmetry_test' )
        call xsymtst%execute( cline )
      case( 'radial_sym_test' )
          call xradsymtst%execute( cline )
    case( 'symmetrize_map' )
        call xsymmetrize_map%execute(cline)
    case( 'dock_volpair' )
        call xdock_volpair%execute(cline)
    case( 'postprocess' )
        call xpostprocess%execute(cline)

    ! IMAGE PROCESSING PROGRAMS
    case('pspec_stats')
        call xpspecstats%execute(cline)
    case( 'mask' )
        call xmask%execute(cline)
    case( 'fsc' )
        call xfsc%execute(cline)
    case( 'local_resolution' )
        call xlocal_res%execute(cline)
    case( 'center' )
        call xcenter%execute(cline)
    case( 'reproject' )
        call xreproject%execute(cline)
    case( 'volops' )
        call xvolops%execute(cline)
    case( 'convert' )
        call xconvert%execute(cline)
    case( 'ctfops' )
        call xctfops%execute(cline)
    case( 'filter' )
        call xfilter%execute(cline)
    case( 'normalize' )
        call xnormalize%execute(cline)
    case( 'scale' )
        call xscale%execute(cline)
    case( 'stack' )
        call xstack%execute(cline)
    case( 'stackops' )
        call xstackops%execute(cline)
    case( 'shift' )
        call xshift%execute(cline)
    case( 'estimate_diam')
        call xestimate_diam%execute(cline)

    ! ORIENTATION PROCESSING PROGRAMS
    case( 'make_oris' )
        call xmake_oris%execute(cline)
    case( 'orisops' )
        call xorisops%execute(cline)
    case( 'oristats' )
        call xoristats%execute(cline)
    case( 'vizoris' )
        call xvizoris%execute(cline)

    ! PRINT INFO PROGRAMS
    case( 'info_image' )
        call xinfo_image%execute(cline)
    case( 'info_stktab' )
        call xinfo_stktab%execute(cline)
    case( 'print_fsc' )
        call xprint_fsc%execute(cline)
    case( 'print_magic_boxes' )
        call xprint_magic_boxes%execute(cline)

    ! SIMULATOR PROGRAMS
    case( 'simulate_noise' )
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        call xsimulate_particles%execute(cline)
    case( 'simulate_movie' )
        call xsimulate_movie%execute(cline)
    case( 'simulate_subtomogram' )
        call xsimulate_subtomogram%execute(cline)
    case( 'simulate_atoms' )
        call xsimulate_atoms%execute(cline)

    ! TIME-SERIES (NANO-PARTICLE) PROGRAMS
    case( 'tseries_import' )
        call xtseries_import%execute(cline)
    case( 'tseries_import_particles' )
        call xtseries_import_particles%execute(cline)
    case( 'tseries_ctf_estimate' )
        call xtseries_ctf_estimate%execute(cline)
    case('detect_atoms')
        call xdetect_atoms%execute(cline)
    case('radial_dependent_stats')
        call xradial_dependent_stats%execute(cline)
    case('atom_cluster_analysis')
        call xatom_cluster_analysis%execute(cline)
    case('nano_softmask')
        call xnano_softmask%execute(cline)
    case('atoms_rmsd')
        call xatoms_rmsd%execute(cline)


    ! SYSTEM INTERACTION PROGRAMS
    case( 'mkdir' )
        call xmkdir%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
end program simple_exec
