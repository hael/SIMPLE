! executes the shared-memory parallelised programs in SIMPLE
program simple_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_user_interface, only: make_user_interface,list_simple_prgs_in_ui
use simple_commanders_abinitio
use simple_commanders_abinitio2D
use simple_commanders_atoms
use simple_commanders_cavgs
use simple_commanders_checks
use simple_commanders_cluster2D
use simple_commanders_distr
use simple_commanders_euclid
use simple_commanders_imgproc
use simple_commanders_mask
use simple_commanders_misc
use simple_commanders_oris
use simple_commanders_preprocess
use simple_commanders_project
use simple_commanders_rec
use simple_commanders_refine3D
use simple_commanders_relion
use simple_commanders_resolest
use simple_commanders_sim
use simple_commanders_starproject
use simple_commanders_validate
use simple_commanders_volops
use simple_stream_cluster2D_subsets
use simple_exec_helpers
implicit none
#include "simple_local_flags.inc"

! PROJECT MANAGEMENT PROGRAMS
type(commander_new_project)                 :: xnew_project
type(commander_update_project)              :: xupdate_project
type(commander_print_project_info)          :: xprint_project_info
type(commander_print_project_field)         :: xprint_project_field
type(commander_zero_project_shifts)         :: xzero_project_shifts
type(commander_import_movies)               :: ximport_movies
type(commander_import_boxes)                :: ximport_boxes
type(commander_import_particles)            :: ximport_particles
type(commander_import_cavgs)                :: ximport_cavgs
type(commander_replace_project_field)       :: xreplace_project_field
type(commander_selection)                   :: xselection
type(commander_export_relion)               :: xexport_relion
type(commander_import_starproject)          :: ximport_starproject
type(commander_export_starproject)          :: xexport_starproject
type(commander_assign_optics_groups)        :: xassign_optics_groups
type(commander_merge_projects)              :: xmerge_projects
type(commander_extract_subproj)             :: xextract_subproj

! PRE-PROCESSING WORKFLOWS
type(commander_preprocess_distr)            :: xpreprocess
type(commander_extract_distr)               :: xextract_distr
type(commander_extract_distr)               :: xreextract_distr
type(commander_motion_correct_distr)        :: xmotion_correct_distr
type(commander_gen_pspecs_and_thumbs_distr) :: xgen_pspecs_and_thumbs
type(commander_ctf_estimate_distr)          :: xctf_estimate_distr
type(commander_pick_distr)                  :: xpick_distr

! CLUSTER2D WORKFLOWS
type(commander_make_cavgs_distr)            :: xmake_cavgs_distr
type(commander_abinitio2D)                  :: xabinitio2D
type(commander_cluster2D_autoscale)         :: xcluster2D_hlev
type(stream_cluster2D_subsets)              :: xcluster2D_subsets
type(commander_cleanup2D_hlev)              :: xcleanup2D_distr
type(commander_map_cavgs_selection)         :: xmap_cavgs_selection
type(commander_map_cavgs_states)            :: xmap_cavgs_states
type(commander_sample_classes)              :: xsample_classes
type(commander_cluster_cavgs)               :: xcluster_cavgs
type(commander_cluster_stack)               :: xcluster_stack
type(commander_select_clusters)             :: xsel_clusts
type(commander_match_cavgs)                 :: xmatch_cavgs
type(commander_match_stacks)                :: xmatch_stacks
type(commander_score_ptcls)                 :: xscore_ptcls
type(commander_write_classes)               :: xwrite_classes
type(commander_write_mic_filetab)           :: xwrite_mic_filetab

! AB INITIO 3D RECONSTRUCTION WORKFLOW
type(commander_estimate_lpstages)           :: xestimate_lpstages
type(commander_noisevol)                    :: xnoisevol
type(commander_abinitio3D_cavgs)            :: xabinitio3D_cavgs
type(commander_abinitio3D_cavgs_fast)       :: xabinitio3D_cavgs_fast
type(commander_abinitio3D)                  :: xabinitio3D
type(commander_multivol_assign)             :: xmultivol_assign

! REFINE3D WORKFLOWS
type(commander_calc_pspec_distr)            :: xcalc_pspec_distr
type(commander_refine3D_distr)              :: xrefine3D_distr
type(commander_refine3D_auto)               :: xrefine3D_auto
type(commander_reconstruct3D_distr)         :: xreconstruct3D

! OTHER SINGLE-PARTICLE WORKFLOW PROGRAMS
type(commander_symaxis_search)              :: xsymsrch
type(commander_symmetry_test)               :: xsymtst
type(commander_symmetrize_map)              :: xsymmetrize_map
type(commander_dock_volpair)                :: xdock_volpair
type(commander_postprocess)                 :: xpostprocess
type(commander_automask)                    :: xautomask
type(commander_auto_spher_mask)             :: xauto_spher_mask
type(commander_fractionate_movies_distr)    :: xfractionate_movies
type(commander_comparemc)                   :: xcomparemc

! VALIDATION WORKFLOWS
type(commander_mini_stream)                 :: xmini_stream
type(commander_check_refpick)               :: xcheck_refpick

! IMAGE PROCESSING PROGRAMS
type(commander_binarize)                    :: xbinarize
type(commander_mask)                        :: xmask
type(commander_automask2D)                  :: xautomask2D
type(commander_fsc)                         :: xfsc
type(commander_clin_fsc)                    :: xclin_fsc
type(commander_nununiform_filter3D)         :: xnununiform_filter3D
type(commander_centervol)                   :: xcenter
type(commander_reproject)                   :: xreproject
type(commander_volanalyze)                  :: xvolanalyze
type(commander_volops)                      :: xvolops
type(commander_convert)                     :: xconvert
type(commander_ctfops)                      :: xctfops
type(commander_ctf_phaseflip)               :: xctf_phaseflip
type(commander_filter)                      :: xfilter
type(commander_normalize)                   :: xnormalize
type(commander_ppca_denoise)                :: xppca_denoise
type(commander_ppca_denoise_classes)        :: xppca_denoise_classes
type(commander_ppca_volvar)                 :: xppca_volvar
type(commander_scale)                       :: xscale
type(commander_stack)                       :: xstack
type(commander_stackops)                    :: xstackops
type(commander_uniform_filter2D)            :: xuniform_filter2D
type(commander_uniform_filter3D)            :: xuniform_filter3D
type(commander_icm2D)                       :: xicm2D
type(commander_icm3D)                       :: xicm3D

! ORIENTATION PROCESSING PROGRAMS
type(commander_check_states)                :: xcheck_states
type(commander_make_oris)                   :: xmake_oris
type(commander_orisops)                     :: xorisops
type(commander_oristats)                    :: xoristats
type(commander_oriconsensus)                :: xoriconsensus
type(commander_vizoris)                     :: xvizoris

! PRINT INFO PROGRAMS
type(commander_info_image)                  :: xinfo_image
type(commander_info_stktab)                 :: xinfo_stktab
type(commander_print_fsc)                   :: xprint_fsc
type(commander_print_magic_boxes)           :: xprint_magic_boxes
type(commander_print_dose_weights)          :: xprint_dose_weights

! SIMULATOR PROGRAMS
type(commander_simulate_noise)              :: xsimulate_noise
type(commander_simulate_particles)          :: xsimulate_particles
type(commander_simulate_movie)              :: xsimulate_movie
type(commander_simulate_subtomogram)        :: xsimulate_subtomogram

! MISCELLANEOUS WORKFLOWS
type(commander_afm)                         :: xafm
type(commander_scale_project_distr)         :: xscale_project
type(commander_map2model_fsc)               :: xmap2model_fsc
type(commander_map_validation)              :: xmap_validation
type(commander_model_validation)            :: xmodel_validation
type(commander_model_validation_eo)         :: xmodel_validation_eo
type(commander_projops)                     :: xprojops
type(commander_prune_project_distr)         :: xprune_project
type(commander_pdb2mrc)                     :: xpdb2mrc   
type(commander_sharpvol)                    :: xsharpvol  

! SYSTEM INTERACTION PROGRAMS
type(commander_mkdir)                       :: xmkdir

! PARALLEL PROCESSING PROGRAMS
type(commander_split)                       :: xsplit

! OTHER DECLARATIONS
character(len=STDLEN)                       :: xarg, prg
character(len=XLONGSTRLEN)                  :: entire_line
type(cmdline)                               :: cline
integer                                     :: cmdstat, cmdlen, pos
integer(timer_int_kind)                     :: t0
real(timer_int_kind)                        :: rt_exec
logical                                     :: l_silent

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
call script_exec(cline, string(trim(prg)), string('simple_exec'))
l_silent = .false.
select case(trim(prg))

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
    case( 'extract_subproj' )
        call xextract_subproj%execute(cline)

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
    case( 'ctf_estimate' )
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        call xpick_distr%execute(cline)

    ! CLUSTER2D WORKFLOWS
    case( 'make_cavgs' )
        call xmake_cavgs_distr%execute(cline)
    case( 'abinitio2D' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('abinitio2D'), string('simple_exec'))
        else
            call xabinitio2D%execute(cline)
        endif
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
    case( 'cluster_stack' )
        call xcluster_stack%execute(cline)
    case('select_clusters')
        call xsel_clusts%execute(cline)
    case( 'match_cavgs' )
        call xmatch_cavgs%execute(cline)
    case( 'match_stacks' )
        call xmatch_stacks%execute(cline)
    case( 'score_ptcls' )
        call xscore_ptcls%execute(cline)
    case( 'write_classes' )
        call xwrite_classes%execute(cline)
    case( 'write_mic_filetab' )
        call xwrite_mic_filetab%execute(cline)

    ! AB INITIO 3D RECONSTRUCTION WORKFLOW
    case('estimate_lpstages')
        call xestimate_lpstages%execute(cline)
    case( 'noisevol' )
        call xnoisevol%execute(cline)
    case( 'abinitio3D_cavgs' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('abinitio3D_cavgs'), string('simple_exec'))
        else
            call xabinitio3D_cavgs%execute(cline)
        endif
    case( 'abinitio3D_cavgs_fast' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('abinitio3D_cavgs_fast'), string('simple_exec'))
        else
            call xabinitio3D_cavgs_fast%execute(cline)
        endif
    case( 'abinitio3D' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('abinitio3D'), string('simple_exec'))
        else
            call xabinitio3D%execute(cline)
        endif
    case('multivol_assign')
        call xmultivol_assign%execute(cline)

    ! REFINE3D WORKFLOWS
    case( 'calc_pspec' )
        call xcalc_pspec_distr%execute(cline)
    case( 'refine3D' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('refine3D'), string('simple_exec'))
        else
            call xrefine3D_distr%execute(cline)
        endif
    case( 'refine3D_auto' )
        if( cline%defined('nrestarts') )then
            call restarted_exec(cline, string('refine3D_auto'), string('simple_exec'))
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

    ! VALIDATION WORKFLOWS
    case( 'mini_stream' )
        call xmini_stream%execute(cline)
    case( 'check_refpick' )
        call xcheck_refpick%execute(cline)

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
        l_silent = .true.
    case( 'print_magic_boxes' )
        call xprint_magic_boxes%execute(cline)
        l_silent = .true.
    case( 'print_dose_weights' )
        call xprint_dose_weights%execute(cline)
        l_silent = .true.

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
if( .not. l_silent )then
    call simple_print_git_version('94ac3e79')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_exec
