! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_user_interface, only: make_user_interface, write_ui_json, print_ui_latex
use simple_symanalyzer,    only: print_subgroups
use simple_commander_project
use simple_commander_checks
use simple_commander_distr
use simple_commander_misc
use simple_commander_imgproc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_refine3D
use simple_commander_rec
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
implicit none
#include "simple_local_flags.inc"

! PRE-PROCESSING PROGRAMS
type(preprocess_commander)            :: xpreprocess
type(extract_commander)               :: xextract
type(reextract_commander)             :: xreextract
type(motion_correct_commander)        :: xmotion_correct
type(gen_pspecs_and_thumbs_commander) :: xgen_pspecs_and_thumbs
type(ctf_estimate_commander)          :: xctf_estimate
type(map_cavgs_selection_commander)   :: xmap_cavgs_selection
type(pick_extract_commander)          :: xpick_extract
type(pick_commander)                  :: xpick
type(make_pickrefs_commander)         :: xmake_pickrefs
type(pick_commander_chiara)           :: xpickchiara

! CLUSTER2D PROGRAMS
type(make_cavgs_commander)            :: xmake_cavgs
type(cluster2D_commander)             :: xcluster2D
type(cavgassemble_commander)          :: xcavgassemble
type(check_2Dconv_commander)          :: xcheck_2Dconv
type(rank_cavgs_commander)            :: xrank_cavgs
type(export_cavgs_commander)          :: xexport_cavgs

! REFINE3D PROGRAMS
type(nspace_commander)                :: xnspace
type(refine3D_commander)              :: xprime3D
type(check_3Dconv_commander)          :: xcheck_3Dconv

! RECONSTRUCTION PROGRAMS
type(volassemble_eo_commander)        :: xvolassemble_eo
type(reconstruct3D_commander)         :: xreconstruct3D

! CHECKER PROGRAMS
type(check_box_commander)             :: xcheck_box
type(check_nptcls_commander)          :: xcheck_nptcls

! VOLOPS PROGRAMS
type(postprocess_commander)           :: xpostprocess
type(volume_smat_commander)           :: xvolume_smat
type(automask_commander)              :: xautomask

! GENERAL IMAGE PROCESSING PROGRAMS
type(scale_commander)                 :: xscale
type(binarise_commander)              :: xbinarise
type(edge_detector_commander)         :: xdetector

! MISCELLANOUS PROGRAMS
type(masscen_commander)               :: xmasscen
type(cluster_smat_commander)          :: xcluster_smat
type(print_dose_weights_commander)    :: xprint_dose_weights
type(res_commander)                   :: xres
type(stk_corr_commander)              :: xstk_corr
type(kstest_commander)                :: xkstst

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(rotmats2oris_commander)         :: xrotmats2oris
type(txt2project_commander)          :: xtxt2project
type(project2txt_commander)          :: xproject2txt
type(print_project_vals_commander)   :: xprint_project_vals
type(multivariate_zscore_commander)  :: xmultizscore
type(o_peaksstats_commander)         :: xo_peaksstats

! TIME-SERIES ANALYSIS PROGRAMS
type(tseries_track_commander)        :: xtseries_track
type(tseries_split_commander)        :: xtseries_split

! PARALLEL PROCESSING PROGRAMS
type(merge_nnmat_commander)          :: xmerge_nnmat
type(merge_similarities_commander)   :: xmerge_similarities
type(split_pairs_commander)          :: xsplit_pairs
type(split_commander)                :: xsplit

! OTHER DECLARATIONS
character(len=STDLEN) :: xarg, prg
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
! this parses all key=value pairs on the command line
call cline%parse_private

select case(prg)

    ! PRIVATE UTILITY PROGRAMS
    case( 'write_ui_json' )
        call write_ui_json
    case( 'print_ui_latex' )
        call print_ui_latex
    case( 'print_sym_subgroups' )
        call print_subgroups

    ! PRE-PROCESSING PROGRAMS
    case( 'preprocess' )
        call xpreprocess%execute(cline)
    case( 'extract' )
        call xextract%execute(cline)
    case( 'reextract' )
        call xreextract%execute(cline)
    case( 'motion_correct' )
        call xmotion_correct%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'ctf_estimate' )
        call xctf_estimate%execute(cline)
    case( 'map_cavgs_selection' ) ! LACKS DESCRIPTION
        call xmap_cavgs_selection%execute(cline)
    case( 'pick_extract' )        ! LACKS DESCRIPTION
        call xpick_extract%execute(cline)
    case( 'pick' )
        call xpick%execute(cline)
    case( 'make_pickrefs' )       ! LACKS DESCRIPTION
        call xmake_pickrefs%execute(cline)
    case ('pick_chiara')          ! LACKS DESCRIPTION
        call xpickchiara%execute(cline)

    ! CLUSTER2D PROGRAMS
    case( 'make_cavgs' )
        call xmake_cavgs%execute(cline)
    case( 'cluster2D' )
        call xcluster2D%execute(cline)
    case( 'cavgassemble' )       ! LACKS DESCRIPTION
        call xcavgassemble%execute(cline)
    case( 'check_2Dconv' )       ! LACKS DESCRIPTION
        call xcheck_2Dconv%execute(cline)
    case( 'rank_cavgs' )         ! LACKS DESCRIPTION
        call xrank_cavgs%execute(cline)
    case( 'export_cavgs' )       ! LACKS DESCRIPTION
        call xexport_cavgs%execute(cline)

    ! REFINE3D PROGRAMS
    case( 'nspace' )
        call xnspace%execute(cline)
    case( 'refine3D' )
        call xprime3D%execute(cline)
    case( 'check_3Dconv' )      ! LACKS DESCRIPTION
        call xcheck_3Dconv%execute(cline)

    ! RECONSTRUCTION PROGRAMS
    case( 'reconstruct3D' )
        call xreconstruct3D%execute(cline)
    case( 'volassemble_eo' )    ! LACKS DESCRIPTION
        call xvolassemble_eo%execute(cline)

    ! CHECKER PROGRAMS
    case( 'check_box' )         ! LACKS DESCRIPTION
        call xcheck_box%execute(cline)
    case( 'check_nptcls' )      ! LACKS DESCRIPTION
        call xcheck_nptcls%execute(cline)

    ! VOLOPS PROGRAMS
    case( 'postprocess' )
        call xpostprocess%execute(cline)
    case( 'volume_smat' )       ! LACKS DESCRIPTION
        call xvolume_smat%execute(cline)
    case( 'automask' )          ! LACKS DESCRIPTION
        call xautomask%execute(cline)

    ! GENERAL IMAGE PROCESSING PROGRAMS
    case( 'scale' )
        call xscale%execute(cline)
    case( 'binarise' )         ! LACKS DESCRIPTION
        call xbinarise%execute(cline)
    case('edge_detect')        ! LACKS DESCRIPTION
        call xdetector%execute(cline)

    ! MISCELLANOUS PROGRAMS
    case( 'masscen' )          ! LACKS DESCRIPTION
        call xmasscen%execute(cline)
    case( 'cluster_smat' )     ! LACKS DESCRIPTION
        call xcluster_smat%execute(cline)
    case( 'print_dose_weights' ) ! DOESN'T DO ANYTHING
        call xprint_dose_weights%execute(cline)
    case( 'res' )              ! LACKS DESCRIPTION
        call xres%execute(cline)
    case( 'stk_corr' )         ! LACKS DESCRIPTION
        call xstk_corr%execute(cline)
    case( 'kstest' )           ! LACKS DESCRIPTION
        call xkstst%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS
    case( 'rotmats2oris' )        ! LACKS DESCRIPTION
        call xrotmats2oris%execute(cline)
    case( 'txt2project' )         ! LACKS DESCRIPTION
        call xtxt2project%execute(cline)
    case( 'project2txt' )         ! LACKS DESCRIPTION
        call xproject2txt%execute(cline)
    case( 'print_project_vals' )  ! LACKS DESCRIPTION
        call xprint_project_vals%execute(cline)
    case( 'multivariate_zscore' ) ! LACKS DESCRIPTION
        call xmultizscore%execute(cline)
    case( 'o_peaksstats')         ! LACKS DESCRIPTION
        call xo_peaksstats%execute(cline)

    ! TIME-SERIES ANALYSIS PROGRAMS
    case( 'tseries_track' )
        call xtseries_track%execute(cline)
    case( 'tseries_split' )      ! LACKS DESCRIPTION
        call xtseries_split%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS
    case( 'merge_nnmat' )        ! LACKS DESCRIPTION
        call xmerge_nnmat%execute(cline)
    case( 'merge_similarities' ) ! LACKS DESCRIPTION
        call xmerge_similarities%execute(cline)
    case( 'split_pairs' )        ! LACKS DESCRIPTION
        call xsplit_pairs%execute(cline)
    case( 'split' )              ! LACKS DESCRIPTION
        call xsplit%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
    end select

end program simple_private_exec
