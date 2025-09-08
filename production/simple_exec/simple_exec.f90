! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface,list_simple_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_exec_helpers
use simple_commander_project
use simple_commander_starproject
use simple_commander_checks
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_pspec
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_cluster2D_stream
use simple_commander_cavgs
use simple_commander_abinitio
use simple_commander_abinitio2D
use simple_commander_refine3D
use simple_commander_rec
use simple_commander_relion
use simple_commander_sim
use simple_commander_volops
use simple_commander_resolest
use simple_commander_euclid
use simple_commander_atoms
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(new_project_commander)                 :: xnew_project
type(update_project_commander)              :: xupdate_project
type(print_project_info_commander)          :: xprint_project_info
type(print_project_field_commander)         :: xprint_project_field
type(zero_project_shifts_commander)         :: xzero_project_shifts
type(import_movies_commander)               :: ximport_movies
type(import_boxes_commander)                :: ximport_boxes
type(import_particles_commander)            :: ximport_particles
type(import_cavgs_commander)                :: ximport_cavgs
type(replace_project_field_commander)       :: xreplace_project_field
type(selection_commander)                   :: xselection
type(export_relion_commander)               :: xexport_relion
type(import_starproject_commander)          :: ximport_starproject
type(export_starproject_commander)          :: xexport_starproject
type(assign_optics_groups_commander)        :: xassign_optics_groups
type(merge_projects_commander)              :: xmerge_projects

! PRE-PROCESSING WORKFLOWS
type(preprocess_commander_distr)            :: xpreprocess
type(extract_commander_distr)               :: xextract_distr
type(reextract_commander_distr)             :: xreextract_distr
type(motion_correct_commander_distr)        :: xmotion_correct_distr
type(gen_pspecs_and_thumbs_commander_distr) :: xgen_pspecs_and_thumbs
type(analyze_pspecs_commander)              :: xpow_anal
type(ctf_estimate_commander_distr)          :: xctf_estimate_distr
type(pick_commander_distr)                  :: xpick_distr

! CLUSTER2D WORKFLOWS
type(make_cavgs_commander_distr)            :: xmake_cavgs_distr
type(abinitio2D_commander)                  :: xabinitio2D
type(cluster2D_autoscale_commander)         :: xcluster2D_hlev
type(cluster2D_commander_subsets)           :: xcluster2D_subsets
type(cleanup2D_commander_hlev)              :: xcleanup2D_distr
type(map_cavgs_selection_commander)         :: xmap_cavgs_selection
type(map_cavgs_states_commander)            :: xmap_cavgs_states
type(sample_classes_commander)              :: xsample_classes
type(cluster_cavgs_commander)               :: xcluster_cavgs
type(select_clusters_commander)             :: xsel_clusts
type(match_cavgs_commander)                 :: xmatch_cavgs
type(match_cavgs2afm_commander)            :: xmatch_cavgs2afm
type(score_ptcls_commander)                 :: xscore_ptcls
type(write_classes_commander)               :: xwrite_classes
type(consolidate_chunks_commander)          :: xconsolidate_chunks

! AB INITIO 3D RECONSTRUCTION WORKFLOW
type(estimate_lpstages_commander)           :: xestimate_lpstages
type(noisevol_commander)                    :: xnoisevol
type(abinitio3D_cavgs_commander)            :: xabinitio3D_cavgs
type(abinitio3D_cavgs_fast_commander)       :: xabinitio3D_cavgs_fast
type(abinitio3D_commander)                  :: xabinitio3D
type(multivol_assign_commander)             :: xmultivol_assign
type(abinitio3D_parts_commander)            :: xabinitio3D_parts

! REFINE3D WORKFLOWS
type(calc_pspec_commander_distr)            :: xcalc_pspec_distr
type(refine3D_distr_commander)              :: xrefine3D_distr
type(refine3D_auto_commander)               :: xrefine3D_auto
type(reconstruct3D_commander_distr)         :: xreconstruct3D

! OTHER SINGLE-PARTICLE WORKFLOW PROGRAMS
type(symaxis_search_commander)              :: xsymsrch
type(symmetry_test_commander)               :: xsymtst
type(symmetrize_map_commander)              :: xsymmetrize_map
type(dock_volpair_commander)                :: xdock_volpair
type(postprocess_commander)                 :: xpostprocess
type(automask_commander)                    :: xautomask
type(auto_spher_mask_commander)             :: xauto_spher_mask
type(fractionate_movies_commander_distr)    :: xfractionate_movies
type(comparemc_commander)                   :: xcomparemc

! IMAGE PROCESSING PROGRAMS
type(binarize_commander)                    :: xbinarize
type(mask_commander)                        :: xmask
type(automask2D_commander)                  :: xautomask2D
type(fsc_commander)                         :: xfsc
type(clin_fsc_commander)                    :: xclin_fsc
type(nununiform_filter3D_commander)         :: xnununiform_filter3D
type(centervol_commander)                   :: xcenter
type(reproject_commander)                   :: xreproject
type(volanalyze_commander)                  :: xvolanalyze
type(volops_commander)                      :: xvolops
type(convert_commander)                     :: xconvert
type(ctfops_commander)                      :: xctfops
type(ctf_phaseflip_commander)               :: xctf_phaseflip
type(filter_commander)                      :: xfilter
type(normalize_commander)                   :: xnormalize
type(ppca_denoise_commander)                :: xppca_denoise
type(ppca_denoise_classes_commander)        :: xppca_denoise_classes
type(ppca_volvar_commander)                 :: xppca_volvar
type(scale_commander)                       :: xscale
type(stack_commander)                       :: xstack
type(stackops_commander)                    :: xstackops
type(uniform_filter2D_commander)            :: xuniform_filter2D
type(uniform_filter3D_commander)            :: xuniform_filter3D
type(icm2D_commander)                       :: xicm2D
type(icm3D_commander)                       :: xicm3D
type(make_pickrefs_commander)               :: xmake_pickrefs

! ORIENTATION PROCESSING PROGRAMS
type(check_states_commander)                :: xcheck_states
type(make_oris_commander)                   :: xmake_oris
type(orisops_commander)                     :: xorisops
type(oristats_commander)                    :: xoristats
type(oriconsensus_commander)                :: xoriconsensus
type(vizoris_commander)                     :: xvizoris

! PRINT INFO PROGRAMS
type(info_image_commander)                  :: xinfo_image
type(info_stktab_commander)                 :: xinfo_stktab
type(print_fsc_commander)                   :: xprint_fsc
type(print_magic_boxes_commander)           :: xprint_magic_boxes
type(print_dose_weights_commander)          :: xprint_dose_weights

! SIMULATOR PROGRAMS
type(simulate_noise_commander)              :: xsimulate_noise
type(simulate_particles_commander)          :: xsimulate_particles
type(simulate_movie_commander)              :: xsimulate_movie
type(simulate_subtomogram_commander)        :: xsimulate_subtomogram

! MISCELLANEOUS WORKFLOWS
type(afm_commander)                         :: xafm
type(scale_project_commander_distr)         :: xscale_project
type(map2model_fsc_commander)               :: xmap2model_fsc
type(map_validation_commander)              :: xmap_validation
type(model_validation_commander)            :: xmodel_validation
type(model_validation_eo_commander)         :: xmodel_validation_eo
type(projops_commander)                     :: xprojops
type(prune_project_commander_distr)         :: xprune_project
type(pdb2mrc_commander)                     :: xpdb2mrc   
type(sharpvol_commander)                    :: xsharpvol  

! SYSTEM INTERACTION PROGRAMS
type(mkdir_commander)                       :: xmkdir

! PARALLEL PROCESSING PROGRAMS
type(split_commander)                       :: xsplit

! OTHER DECLARATIONS
character(len=STDLEN)                       :: xarg, prg, entire_line
type(cmdline)                               :: cline
integer                                     :: cmdstat, cmdlen, pos
integer(timer_int_kind)                     :: t0
real(timer_int_kind)                        :: rt_exec

! start timer
t0 = tic()
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
! generate script for queue submission?
call script_exec(cline, trim(prg), 'simple_exec')

select case(trim(prg))

    ! PROJECT MANAGEMENT PROGRAMS
    case( 'new_project' )
        call xnew_project%execute(cline)
    case( 'update_project' )
        call xupdate_project%execute(cline)
    case( 'print_project_info' )
        call xprint_project_info%execute(cline)
    case( 'print_project_field' )
        call xprint_project_field%execute(cline)
    case( 'zero_project_shifts' )
        call xzero_project_shifts%execute(cline)
    case( 'import_movies' )
        call ximport_movies%execute(cline)
    case( 'import_boxes' )
        call ximport_boxes%execute(cline)
    case( 'import_particles' )
        call ximport_particles%execute(cline)
    case( 'import_cavgs' )
        call ximport_cavgs%execute(cline)
    case( 'replace_project_field' )
        call xreplace_project_field%execute(cline)
    case( 'selection', 'report_selection' )
        call xselection%execute(cline)
    case( 'export_relion' )
        call xexport_relion%execute(cline)
    case( 'import_starproject' )
        call ximport_starproject%execute(cline)
    case( 'export_starproject' )
        call xexport_starproject%execute(cline)
    case( 'assign_optics_groups' )
        call xassign_optics_groups%execute(cline)
    case( 'merge_projects' )
        call xmerge_projects%execute(cline)

    ! PRE-PROCESSING WORKFLOWS
    case( 'preprocess' )
        call xpreprocess%execute(cline)
    case( 'extract' )
        call xextract_distr%execute(cline)
    case( 'reextract' )
        call xreextract_distr%execute(cline)
    case( 'motion_correct' )
        call xmotion_correct_distr%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'analyze_pspecs' )
        call xpow_anal%execute(cline)
    case( 'ctf_estimate' )
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        call xpick_distr%execute(cline)

    ! CLUSTER2D WORKFLOWS
    case( 'make_cavgs' )
        call xmake_cavgs_distr%execute(cline)
    case( 'abinitio2D' )
        call xabinitio2D%execute(cline)
    case( 'cleanup2D' )
        call xcleanup2D_distr%execute(cline)
    case( 'cluster2D' )
        call xcluster2D_hlev%execute(cline)
    case( 'cluster2D_subsets' )
        call xcluster2D_subsets%execute(cline)
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case( 'map_cavgs_states' )
        call xmap_cavgs_states%execute(cline)
    case('sample_classes')
        call xsample_classes%execute(cline)
    case( 'cluster_cavgs' )
        call xcluster_cavgs%execute(cline)
    case('select_clusters')
        call xsel_clusts%execute(cline)
    case( 'match_cavgs' )
        call xmatch_cavgs%execute(cline)
    case ( 'match_cavgs2afm' )
        call xmatch_cavgs2afm%execute(cline)
    case( 'score_ptcls' )
        call xscore_ptcls%execute(cline)
    case( 'write_classes' )
        call xwrite_classes%execute(cline)
    case('consolidate_chunks')
        call xconsolidate_chunks%execute(cline)

    ! AB INITIO 3D RECONSTRUCTION WORKFLOW
    case('estimate_lpstages')
        call xestimate_lpstages%execute(cline)
    case( 'noisevol' )
        call xnoisevol%execute(cline)
    case( 'abinitio3D_cavgs' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'abinitio3D_cavgs', 'simple_exec')
        else
            call xabinitio3D_cavgs%execute(cline)
        endif
    case( 'abinitio3D_cavgs_fast' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'abinitio3D_cavgs_fast', 'simple_exec')
        else
            call xabinitio3D_cavgs_fast%execute(cline)
        endif
    case( 'abinitio3D' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'abinitio3D', 'simple_exec')
        else
            call xabinitio3D%execute(cline)
        endif
    case('multivol_assign')
        call xmultivol_assign%execute(cline)
    case( 'abinitio3D_parts' )
        call xabinitio3D_parts%execute(cline)

    ! REFINE3D WORKFLOWS
    case( 'calc_pspec' )
        call xcalc_pspec_distr%execute(cline)
    case( 'refine3D' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'refine3D', 'simple_exec')
        else
            call xrefine3D_distr%execute(cline)
        endif
    case( 'refine3D_auto' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, 'refine3D_auto', 'simple_exec')
        else
            call xrefine3D_auto%execute(cline)
        endif
    case( 'reconstruct3D' )
        call xreconstruct3D%execute( cline )

    ! OTHER SINGLE-PARTICLE WORKFLOW PROGRAMS
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
    case( 'automask' )
        call xautomask%execute(cline)
    case( 'auto_spher_mask' )
        call xauto_spher_mask%execute(cline)
    case( 'fractionate_movies' )
        call xfractionate_movies%execute(cline)
    case( 'comparemc' )
        call xcomparemc%execute(cline)

    ! IMAGE PROCESSING PROGRAMS
    case( 'binarize' )
        call xbinarize%execute(cline)
    case( 'mask' )
        call xmask%execute(cline)
    case( 'automask2D' )
        call xautomask2D%execute(cline)
    case( 'fsc' )
        call xfsc%execute(cline)
    case( 'clin_fsc' )
        call xclin_fsc%execute(cline)
    case( 'nununiform_filter3D' )
        call xnununiform_filter3D%execute(cline)
    case( 'center' )
        call xcenter%execute(cline)
    case( 'reproject' )
        call xreproject%execute(cline)
    case( 'volops' )
        call xvolops%execute(cline)
    case( 'volanalyze' )
        call xvolanalyze%execute(cline)
    case( 'convert' )
        call xconvert%execute(cline)
    case( 'ctfops' )
        call xctfops%execute(cline)
    case( 'ctf_phaseflip' )
        call xctf_phaseflip%execute(cline)
    case( 'filter' )
        call xfilter%execute(cline)
    case( 'normalize' )
        call xnormalize%execute(cline)
    case( 'ppca_denoise' )
        call xppca_denoise%execute(cline)
    case( 'ppca_denoise_classes' )
        call xppca_denoise_classes%execute(cline)
    case( 'ppca_volvar' )
        call xppca_volvar%execute(cline)
    case( 'scale' )
        call xscale%execute(cline)
    case( 'stack' )
        call xstack%execute(cline)
    case( 'stackops' )
        call xstackops%execute(cline)
    case( 'uniform_filter2D' )
        call xuniform_filter2D%execute(cline)
    case( 'uniform_filter3D' )
        call xuniform_filter3D%execute(cline)
    case( 'icm2D' )
        call xicm2D%execute(cline)
    case( 'icm3D' )
        call xicm3D%execute(cline)
    case( 'make_pickrefs' )
        call xmake_pickrefs%execute(cline)

    ! ORIENTATION PROCESSING PROGRAMS
    case( 'check_states' )
        call xcheck_states%execute(cline)
    case( 'make_oris' )
        call xmake_oris%execute(cline)
    case( 'orisops' )
        call xorisops%execute(cline)
    case( 'oristats' )
        call xoristats%execute(cline)
    case( 'oriconsensus' )
        call xoriconsensus%execute(cline)
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
    case( 'print_dose_weights' )
        call xprint_dose_weights%execute(cline)

    ! SIMULATOR PROGRAMS
    case( 'simulate_noise' )
        call xsimulate_noise%execute(cline)
    case( 'simulate_particles' )
        call xsimulate_particles%execute(cline)
    case( 'simulate_movie' )
        call xsimulate_movie%execute(cline)
    case( 'simulate_subtomogram' )
        call xsimulate_subtomogram%execute(cline)

    ! VALIDATION PROGRAMS
    case( 'map2model_fsc' )
        call xmap2model_fsc%execute(cline)
    case( 'map_validation' )
        call xmap_validation%execute(cline)
    case( 'model_validation' )
        call xmodel_validation%execute(cline)
    case( 'model_validation_eo' )
        call xmodel_validation_eo%execute(cline)

    ! MISCELLANEOUS WORKFLOWS
    case( 'afm' )
        call xafm%execute( cline )
    case( 'pdb2mrc' )
        call xpdb2mrc%execute( cline )
    case( 'projops' )
        call xprojops%execute( cline )   
    case( 'prune_project' )
        call xprune_project%execute( cline )
    case( 'scale_project' )
        call xscale_project%execute( cline )
    case( 'sharpvol' )
        call xsharpvol%execute( cline )


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
call simple_print_git_version('dc6dd592')
! end timer and print
rt_exec = toc(t0)
call simple_print_timer(rt_exec)
end program simple_exec
