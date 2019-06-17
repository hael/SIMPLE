! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_user_interface, only: write_ui_json, print_ui_latex
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
type(postprocess_commander)           :: xpostprocess   ! DUPLICATED
type(reproject_commander)             :: xreproject     ! DUPLICATED
type(volume_smat_commander)           :: xvolume_smat
type(automask_commander)              :: xautomask

! GENERAL IMAGE PROCESSING PROGRAMS
type(scale_commander)                 :: xscale         ! DUPLICATED
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
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: xarg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos

! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
call cline%parse_oldschool

select case(prg)

    ! not bona fide simple programs
    case( 'write_ui_json')
        call write_ui_json
    case( 'print_ui_latex')
        call print_ui_latex
    case( 'print_sym_subgroups')
        call print_subgroups

    ! PRE-PROCESSING PROGRAMS

    case( 'preprocess' )
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',              5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',           8.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate', 30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',  5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',         20.)
        if( .not. cline%defined('outfile')         ) call cline%set('outfile', 'simple_unidoc'//METADATA_EXT)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',   'black')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        call xpreprocess%execute(cline)
    case( 'extract' )
        ! set defaults
        call cline%set('nthr',1.)
        if( .not. cline%defined('outside')         ) call cline%set('outside',   'no')
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream', 'no')
        if( cline%defined('stream') )then
            if( cline%get_carg('stream').eq.'no' .and. .not.cline%defined('outfile') )then
                THROW_HARD('OUTFILE must be defined with STREAM=NO')
            endif
        endif
        call xextract%execute(cline)
    case( 'reextract' )
        ! set defaults
        call cline%set('nthr',1.)
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        call xreextract%execute(cline)
    case( 'motion_correct' )
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'simple_unidoc'//METADATA_EXT)
        call xmotion_correct%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'ctf_estimate' )
        ! set defaults
        if( .not. cline%defined('hp') ) call cline%set('hp', 30.)
        if( .not. cline%defined('lp') ) call cline%set('lp',  5.)
        call xctf_estimate%execute(cline)
    case( 'map_cavgs_selection' )
        call xmap_cavgs_selection%execute(cline)
    case( 'pick_extract' )
        call xpick_extract%execute(cline)
    case( 'pick' )
        call xpick%execute(cline)
    case( 'make_pickrefs' )
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        if( .not.cline%defined('refs') .and. .not.cline%defined('vol1') )then
            THROW_HARD('One of REFS, VOL1 & PROJFILE must be informed!')
        endif
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast','black')
        call xmake_pickrefs%execute(cline)
    case ('pick_chiara')
        call xpickchiara%execute(cline)

    ! CLUSTER2D PROGRAMS

    case( 'make_cavgs' )
        call xmake_cavgs%execute(cline)
    case( 'cluster2D' )
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',     30.)
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',    30.)
        call xcluster2D%execute(cline)
    case( 'cavgassemble' )
        call xcavgassemble%execute(cline)
    case( 'check_2Dconv' )
        ! set defaults
        call cline%set('oritype', 'ptcl2D')
        call xcheck_2Dconv%execute(cline)
    case( 'rank_cavgs' )
        ! set defaults
        call cline%set('oritype', 'cls2D')
        call xrank_cavgs%execute(cline)
    case( 'export_cavgs' )
        ! set defaults
        call cline%set('oritype', 'cls2D')
        call xexport_cavgs%execute(cline)

    ! REFINE3D PROGRAMS

    case( 'nspace' )
        call xnspace%execute(cline)
    case( 'refine3D' )
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        if( .not. cline%defined('refine') )then
            call cline%set('refine',  'single')
        else
            if( cline%get_carg('refine').eq.'multi' .and. .not. cline%defined('nstates') )then
                THROW_HARD('refine=MULTI requires specification of NSTATES')
            endif
        endif
        call xprime3D%execute(cline)
    case( 'check_3Dconv' )
        ! set defaults
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! execute
        call xcheck_3Dconv%execute(cline)

    ! RECONSTRUCTION PROGRAMS

    case( 'reconstruct3D' )
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 5.) ! to assure that shifts are being used
        call xreconstruct3D%execute(cline)
    case( 'volassemble_eo' )
        call xvolassemble_eo%execute(cline)

    ! CHECKER PROGRAMS

    case( 'check_box' )
        call xcheck_box%execute(cline)
    case( 'check_nptcls' )
        call xcheck_nptcls%execute(cline)

    ! VOLOPS PROGRAMS

    case( 'postprocess' )
        call xpostprocess%execute(cline)
    case( 'reproject' )
        ! set defaults
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        call xreproject%execute(cline)
    case( 'volume_smat' )
        call xvolume_smat%execute(cline)
    case( 'automask' )
        call xautomask%execute(cline)

    ! GENERAL IMAGE PROCESSING PROGRAMS

    case( 'scale' )
        call xscale%execute(cline)
    case( 'binarise' )
        call xbinarise%execute(cline)
    case('edge_detect')
        call xdetector%execute(cline)

    ! MISCELLANOUS PROGRAMS

    case( 'masscen' )
        call xmasscen%execute(cline)
    case( 'cluster_smat' )
        call xcluster_smat%execute(cline)
    case( 'print_dose_weights' )
        call xprint_dose_weights%execute(cline)
    case( 'res' )
        call xres%execute(cline)
    case( 'stk_corr' )
        call xstk_corr%execute(cline)
    case( 'kstest' )
        call xkstst%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS

    case( 'rotmats2oris' )
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call xrotmats2oris%execute(cline)
    case( 'txt2project' )
        call xtxt2project%execute(cline)
    case( 'project2txt' )
        call xproject2txt%execute(cline)
    case( 'print_project_vals' )
        call xprint_project_vals%execute(cline)
    case( 'multivariate_zscore' )
        call xmultizscore%execute(cline)
    case( 'o_peaksstats')
        call xo_peaksstats%execute(cline)

    ! TIME-SERIES ANALYSIS PROGRAMS

    case( 'tseries_track' )
        ! set defaults
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        call xtseries_track%execute(cline)
    case( 'tseries_split' )
        call xtseries_split%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS

    case( 'merge_nnmat' )
        call xmerge_nnmat%execute(cline)
    case( 'merge_similarities' )
        call xmerge_similarities%execute(cline)
    case( 'split_pairs' )
        call xsplit_pairs%execute(cline)
    case( 'split' )
        call xsplit%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
    end select

end program simple_private_exec
