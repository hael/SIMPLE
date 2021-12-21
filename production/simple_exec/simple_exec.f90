! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface,list_simple_prgs_in_ui
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
use simple_commander_cluster2D_stream
use simple_commander_abinitio
use simple_commander_refine3D
use simple_commander_cluster3D
use simple_commander_rec
use simple_commander_relion
use simple_commander_sim
use simple_commander_volops
use simple_commander_resolest
use simple_spproj_hlev
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(new_project_commander)                 :: xnew_project
type(update_project_commander)              :: xupdate_project
type(print_project_info_commander)          :: xprint_project_info
type(print_project_field_commander)         :: xprint_project_field
type(import_movies_commander)               :: ximport_movies
type(import_boxes_commander)                :: ximport_boxes
type(import_particles_commander)            :: ximport_particles
type(import_cavgs_commander)                :: ximport_cavgs
type(merge_stream_projects_commander)       :: xmerge_stream_projects
type(replace_project_field_commander)       :: xreplace_project_field
type(selection_commander)                   :: xselection
type(export_relion_commander)               :: xexport_relion

! PRE-PROCESSING WORKFLOWS
type(preprocess_commander_distr)            :: xpreprocess
type(preprocess_commander_stream)           :: xpreprocess_stream
type(extract_commander_distr)               :: xextract_distr
type(reextract_commander_distr)             :: xreextract_distr
type(motion_correct_commander_distr)        :: xmotion_correct_distr
type(gen_pspecs_and_thumbs_commander_distr) :: xgen_pspecs_and_thumbs
type(motion_correct_tomo_commander_distr)   :: xmotion_correct_tomo_distr
type(ctf_estimate_commander_distr)          :: xctf_estimate_distr
type(pick_commander_distr)                  :: xpick_distr

! CLUSTER2D WORKFLOWS
type(make_cavgs_commander_distr)            :: xmake_cavgs_distr
type(cluster2D_autoscale_commander_hlev)    :: xcluster2D_hlev
type(cluster2D_commander_stream)            :: xcluster2D_stream
type(cleanup2D_commander_hlev)              :: xcleanup2D_distr

! AB INITIO 3D RECONSTRUCTION WORKFLOW
type(initial_3Dmodel_commander_hlev)        :: xinitial_3Dmodel

! REFINE3D WORKFLOWS
type(calc_pspec_commander_distr)            :: xcalc_pspec_distr
type(refine3D_commander_distr)              :: xrefine3D_distr
type(reconstruct3D_commander_distr)         :: xreconstruct3D_distr

! CLUSTER3D WORKFLOWS
type(cluster3D_commander)                   :: xcluster3D
type(cluster3D_refine_commander)            :: xcluster3D_refine

! OTHER SINGLE-PARTICLE WORKFLOW PROGRAMS
type(map_cavgs_selection_commander)         :: xmap_cavgs_selection
type(cluster_cavgs_commander)               :: xcluster_cavgs
type(write_classes_commander)               :: xwrite_classes
type(symaxis_search_commander)              :: xsymsrch
type(symmetry_test_commander)               :: xsymtst
type(symmetrize_map_commander)              :: xsymmetrize_map
type(dock_volpair_commander)                :: xdock_volpair
type(postprocess_commander)                 :: xpostprocess
type(remoc_commander)                       :: xremoc

! IMAGE PROCESSING PROGRAMS
type(mask_commander)                        :: xmask
type(fsc_commander)                         :: xfsc
type(nonuniform_filter_commander)           :: xnonuniform_filter
type(centervol_commander)                   :: xcenter
type(reproject_commander)                   :: xreproject
type(volops_commander)                      :: xvolops
type(convert_commander)                     :: xconvert
type(ctfops_commander)                      :: xctfops
type(filter_commander)                      :: xfilter
type(normalize_commander)                   :: xnormalize
type(scale_commander)                       :: xscale
type(stack_commander)                       :: xstack
type(stackops_commander)                    :: xstackops

! ORIENTATION PROCESSING PROGRAMS
type(make_oris_commander)                   :: xmake_oris
type(orisops_commander)                     :: xorisops
type(oristats_commander)                    :: xoristats
type(vizoris_commander)                     :: xvizoris

! PRINT INFO PROGRAMS
type(info_image_commander)                  :: xinfo_image
type(info_stktab_commander)                 :: xinfo_stktab
type(print_fsc_commander)                   :: xprint_fsc
type(print_magic_boxes_commander)           :: xprint_magic_boxes

! SIMULATOR PROGRAMS
type(simulate_noise_commander)              :: xsimulate_noise
type(simulate_particles_commander)          :: xsimulate_particles
type(simulate_movie_commander)              :: xsimulate_movie
type(simulate_subtomogram_commander)        :: xsimulate_subtomogram

! MISCELLANEOUS WORKFLOWS
type(scale_project_commander_distr)         :: xscale_project
type(prune_project_commander_distr)         :: xprune_project

! SYSTEM INTERACTION PROGRAMS
type(mkdir_commander)                       :: xmkdir

! PARALLEL PROCESSING PROGRAMS
type(split_commander)                       :: xsplit

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
    call list_simple_prgs_in_ui
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

    ! PRE-PROCESSING WORKFLOWS
    case( 'preprocess' )
        call xpreprocess%execute(cline)
    case( 'preprocess_stream' )
        call xpreprocess_stream%execute(cline)
    case( 'extract' )
        call xextract_distr%execute(cline)
    case( 'reextract' )
        call xreextract_distr%execute(cline)
    case( 'motion_correct' )
        call xmotion_correct_distr%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'motion_correct_tomo' )
        call xmotion_correct_tomo_distr%execute(cline)
    case( 'ctf_estimate' )
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        call xpick_distr%execute(cline)

    ! CLUSTER2D WORKFLOWS
    case( 'make_cavgs' )
        call xmake_cavgs_distr%execute(cline)
    case( 'cleanup2D' )
        call xcleanup2D_distr%execute(cline)
    case( 'cluster2D' )
        call xcluster2D_hlev%execute(cline)
    case( 'cluster2D_stream' )
        call xcluster2D_stream%execute(cline)

    ! AB INITIO 3D RECONSTRUCTION WORKFLOW
    case( 'initial_3Dmodel' )
        call xinitial_3Dmodel%execute(cline)

    ! REFINE3D WORKFLOWS
    case( 'calc_pspec' )
        call xcalc_pspec_distr%execute(cline)
    case( 'refine3D' )
        call xrefine3D_distr%execute(cline)
    case( 'reconstruct3D' )
        call xreconstruct3D_distr%execute( cline )

    ! CLUSTER3D WORKFLOWS
    case( 'cluster3D' )
        call xcluster3D%execute( cline )
    case( 'cluster3D_refine' )
        call xcluster3D_refine%execute( cline )

    ! OTHER SINGLE-PARTICLE WORKFLOW PROGRAMS
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
    case( 'symmetrize_map' )
        call xsymmetrize_map%execute(cline)
    case( 'dock_volpair' )
        call xdock_volpair%execute(cline)
    case( 'postprocess' )
        call xpostprocess%execute(cline)
    case( 'remoc' )
        call xremoc%execute(cline)

    ! IMAGE PROCESSING PROGRAMS
    case( 'mask' )
        call xmask%execute(cline)
    case( 'fsc' )
        call xfsc%execute(cline)
    case( 'nonuniform_filter' )
        call xnonuniform_filter%execute(cline)
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

    ! MISCELLANEOUS WORKFLOWS
    case( 'scale_project' )
        call xscale_project%execute( cline )
    case( 'prune_project' )
        call xprune_project%execute( cline )

    ! SYSTEM INTERACTION PROGRAMS
    case( 'mkdir' )
        call xmkdir%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS
    case( 'split' )
        call xsplit%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
end program simple_exec
