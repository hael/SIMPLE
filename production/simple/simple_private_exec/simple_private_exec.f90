! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
include 'simple_lib.f08'
use simple_cmdline, only: cmdline, cmdline_err
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
type(preprocess_commander)           :: xpreprocess
type(powerspecs_commander)           :: xpowerspecs
type(motion_correct_commander)       :: xmotion_correct
type(ctf_estimate_commander)         :: xctf_estimate
type(map_cavgs_selection_commander)  :: xmap_cavgs_selection
type(pick_extract_commander)         :: xpick_extract
type(pick_commander)                 :: xpick

! CLUSTER2D PROGRAMS
type(make_cavgs_commander)           :: xmake_cavgs
type(cluster2D_commander)            :: xcluster2D
type(cavgassemble_commander)         :: xcavgassemble
type(check_2Dconv_commander)         :: xcheck_2Dconv
type(rank_cavgs_commander)           :: xrank_cavgs

! REFINE3D PROGRAMS
type(nspace_commander)               :: xnspace
type(refine3D_init_commander)        :: xrefine3D_init
type(refine3D_commander)             :: xprime3D
type(check_3Dconv_commander)         :: xcheck_3Dconv

! RECONSTRUCTION PROGRAMS
type(volassemble_eo_commander)       :: xvolassemble_eo
type(reconstruct3D_commander)        :: xreconstruct3D
type(volassemble_commander)          :: xvolassemble

! CHECKER PROGRAMS
type(check_box_commander)            :: xcheck_box
type(check_nptcls_commander)         :: xcheck_nptcls

! VOLOPS PROGRAMS
type(postprocess_commander)          :: xpostprocess   ! DUPLICATED
type(reproject_commander)            :: xreproject     ! DUPLICATED
type(volume_smat_commander)          :: xvolume_smat
type(dock_volpair_commander)         :: xdock_volpair
type(symmetrize_map_commander)       :: xsymmetrize_map
type(automask_commander)             :: xautomask

! GENERAL IMAGE PROCESSING PROGRAMS
type(scale_commander)                :: xscale         ! DUPLICATED
type(binarise_commander)             :: xbinarise
type(edge_detector_commander)        :: xdetector

! MISCELLANOUS PROGRAMS
type(masscen_commander)              :: xmasscen
type(cluster_smat_commander)         :: xcluster_smat
type(intgpeaks_commander)            :: xintgpeaks
type(print_dose_weights_commander)   :: xprint_dose_weights
type(res_commander)                  :: xres
type(stk_corr_commander)             :: xstk_corr

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(rotmats2oris_commander)              :: xrotmats2oris
type(txt2project_commander)               :: xtxt2project
type(project2txt_commander)               :: xproject2txt
type(print_project_header_commander)      :: xprint_project_header
type(print_project_vals_commander)        :: xprint_project_vals
type(update_project_stateflags_commander) :: xupdate_project_stateflags
type(multivariate_zscore_commander)       :: xmultizscore

! TIME-SERIES ANALYSIS PROGRAMS
type(tseries_extract_commander)      :: xtseries_extract
type(tseries_track_commander)        :: xtseries_track
type(tseries_backgr_subtr_commander) :: xtseries_backgr_subtr
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
select case(prg)

    ! PRE-PROCESSING PROGRAMS

    case( 'preprocess' )
        ! executes motion_correct, ctf_estimate and pick in sequence
        keys_required(1)   = 'projfile'
        ! set optional keys
        keys_optional(1)   = 'nthr'
        keys_optional(2)   = 'refs'
        keys_optional(3)   = 'fbody'
        keys_optional(4)   = 'dose_rate'
        keys_optional(5)   = 'exp_time'
        keys_optional(6)   = 'lpstart'
        keys_optional(7)   = 'lpstop'
        keys_optional(8)   = 'trs'
        keys_optional(9)   = 'pspecsz'
        keys_optional(10)  = 'numlen'
        keys_optional(11)  = 'startit'
        keys_optional(12)  = 'scale'
        keys_optional(13)  = 'nframesgrp'
        keys_optional(14)  = 'fromf'
        keys_optional(15)  = 'tof'
        keys_optional(16)  = 'hp_ctf_estimate'
        keys_optional(17)  = 'lp_ctf_estimate'
        keys_optional(18)  = 'lp_pick'
        keys_optional(19)  = 'dfmin'
        keys_optional(20)  = 'dfmax'
        keys_optional(21)  = 'dfstep'
        keys_optional(22)  = 'astigtol'
        keys_optional(23)  = 'phaseplate'
        keys_optional(24)  = 'thres'
        keys_optional(25)  = 'rm_outliers'
        keys_optional(26)  = 'nsig'
        keys_optional(27)  = 'dopick'
        keys_optional(28)  = 'ndev'
        keys_optional(29)  = 'pcontrast'
        keys_optional(30)  = 'ctfreslim'
        call cline%parse_oldschool(keys_required(:1), keys_optional(:30))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',              5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',           8.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',        512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate', 30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',  5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',         20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        call xpreprocess%execute(cline)
    case( 'powerspecs' )
        ! for generating powerspectra from a stack or filetable
        keys_required(1) = 'smpd'
        keys_required(2) = 'fbody'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'stk'
        keys_optional(3) = 'filetab'
        keys_optional(4) = 'pspecsz'
        keys_optional(5) = 'speckind'
        keys_optional(6) = 'startit'
        keys_optional(7) = 'lp'
        keys_optional(8) = 'clip'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        call xpowerspecs%execute(cline)
    case( 'motion_correct' )
        ! for movie alignment
        keys_required(1)  = 'projfile'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'fbody'
        keys_optional(3)  = 'dose_rate'
        keys_optional(4)  = 'exp_time'
        keys_optional(5)  = 'lpstart'
        keys_optional(6)  = 'lpstop'
        keys_optional(7)  = 'trs'
        keys_optional(8)  = 'pspecsz'
        keys_optional(9)  = 'dir'
        keys_optional(10) = 'startit'
        keys_optional(11) = 'scale'
        keys_optional(12) = 'nframesgrp'
        keys_optional(13) = 'tomo'
        keys_optional(14) = 'fromf'
        keys_optional(15) = 'tof'
        keys_optional(16) = 'nsig'
        call cline%parse_oldschool(keys_required(:1), keys_optional(:16))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        call xmotion_correct%execute(cline)
    case( 'ctf_estimate' )
        ! for fitting the CTF
        keys_required(1)  = 'projfile'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'pspecsz'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'lp'
        keys_optional(5) = 'dfmin'
        keys_optional(6) = 'dfmax'
        keys_optional(7) = 'astigtol'
        keys_optional(8) = 'phaseplate'
        keys_optional(9) = 'dir'
        call cline%parse_oldschool(keys_required(:1), keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        call xctf_estimate%execute(cline)
    case( 'map_cavgs_selection' )
        ! for mapping class average selection to project
        keys_required(1)  = 'stk'
        keys_required(2)  = 'stk2'
        keys_required(3)  = 'projfile'
        call cline%parse_oldschool(keys_required(:3))
        call xmap_cavgs_selection%execute(cline)
    case( 'motion_correct_ctf_estimate' )
        ! pipelined motion_correct + ctf_estimate
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'kv'
        keys_required(4)  = 'cs'
        keys_required(5)  = 'fraca'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'fbody'
        keys_optional(3)  = 'dose_rate'
        keys_optional(4)  = 'exp_time'
        keys_optional(5)  = 'lpstart'
        keys_optional(6)  = 'lpstop'
        keys_optional(7)  = 'trs'
        keys_optional(8)  = 'pspecsz'
        keys_optional(9)  = 'numlen'
        keys_optional(10) = 'startit'
        keys_optional(11) = 'scale'
        keys_optional(12) = 'nframesgrp'
        keys_optional(13) = 'fromf'
        keys_optional(14) = 'tof'
        keys_optional(15) = 'nsig'
        keys_optional(16) = 'outfile'
        keys_optional(17) = 'hp'
        keys_optional(18) = 'lp'
        keys_optional(19) = 'dfmin'
        keys_optional(20) = 'dfmax'
        keys_optional(21) = 'dfstep'
        keys_optional(22) = 'astigtol'
        keys_optional(23) = 'phaseplate'
        call cline%parse_oldschool(keys_required(:5), keys_optional(:23))
        ! set defaults
        call cline%set('dopick', 'no')
        call cline%set('prg', 'preprocess')
        if( .not. cline%defined('trs')             ) call cline%set('trs',              5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',           8.)
        if( .not. cline%defined('pscpecsz')        ) call cline%set('pscpecsz',       512.)
        if( .not. cline%defined('hp_ctfestimate')  ) call cline%set('hp_ctfestimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',  5.)
        call xpreprocess%execute(cline)
    case( 'pick_extract' )
        ! for template-based particle picking
        keys_required(1) = 'projfile'
        keys_required(2) = 'refs'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'ndev'
        keys_optional(5) = 'box_extract'
        keys_optional(6) = 'pcontrast'
        keys_optional(7) = 'outside'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:7))
        if( .not. cline%defined('pcontrast') )call cline%set('pcontrast', 'black')
        call xpick_extract%execute(cline)
    case( 'pick' )
        ! for template-based particle picking
        keys_required(1) = 'projfile'
        keys_required(2) = 'refs'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'ndev'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:4))
        call xpick%execute(cline)

    ! CLUSTER2D PROGRAMS

    case( 'make_cavgs' )
        ! for producing class averages or initial random references for cluster2D execution
        keys_required(1)  = 'projfile'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncls'
        keys_optional(3)  = 'mul'
        keys_optional(4)  = 'tseries'
        keys_optional(5)  = 'outfile'
        keys_optional(6)  = 'refs'
        keys_optional(7)  = 'remap_cls'
        call cline%parse_oldschool(keys_required(:1), keys_optional(:7))
        call xmake_cavgs%execute(cline)
    case( 'cluster2D' )
        ! is a reference-free 2D alignment/clustering algorithm adopted from the PRIME
        ! probabilistic ab initio 3D reconstruction algorithm
        keys_required(1)  = 'projfile'
        keys_required(2)  = 'msk'
        keys_required(3)  = 'ncls'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'objfun'
        keys_optional(3)  = 'phaseplate'
        keys_optional(4)  = 'refs'
        keys_optional(5)  = 'hp'
        keys_optional(6)  = 'lp'
        keys_optional(7)  = 'lpstart'
        keys_optional(8)  = 'lpstop'
        keys_optional(9)  = 'cenlp'
        keys_optional(10) = 'trs'
        keys_optional(11) = 'inner'
        keys_optional(12) = 'width'
        keys_optional(13) = 'startit'
        keys_optional(14) = 'maxits'
        keys_optional(15) = 'center'
        keys_optional(16) = 'refine'
        keys_optional(17) = 'match_filt'
        keys_optional(18) = 'shellw'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:18))
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',     30.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',       'no')
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',    30.)
        if( .not. cline%defined('refine')    ) call cline%set('refine',   'snhc')
        call xcluster2D%execute(cline)
    case( 'cavgassemble' )
        ! for assembling class averages when the clustering
        ! program (cluster2D) has been executed in distributed mode
        keys_required(1) = 'projfile'
        keys_required(2) = 'nparts'
        keys_required(3) = 'ncls'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'refs'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:2))
        call xcavgassemble%execute(cline)
    case( 'check_2Dconv' )
        ! for convergence checking and run-time stats printing (3D)
        keys_required(1) = 'projfile'
        call cline%parse_oldschool(keys_required(:1))
        ! set defaults
        call cline%set('oritype', 'ptcl2D')
        call xcheck_2Dconv%execute(cline)
    case( 'rank_cavgs' )
        ! for ranking class averages
        keys_required(1) = 'projfile'
        keys_required(2) = 'stk'
        ! set optional keys
        keys_optional(1) = 'outstk'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        ! set defaults
        call cline%set('oritype', 'cls2D')
        call xrank_cavgs%execute(cline)

    ! REFINE3D PROGRAMS

    case( 'nspace' )
        ! for calculating the expected resolution obtainable with different values of nspace
        ! (number of discrete projection directions used for discrete search)
        keys_required(1)  = 'moldiam'
        call cline%parse_oldschool(keys_required=keys_required(:1))
        call xnspace%execute(cline)
    case( 'refine3D_init' )
        ! initialization of 3D refinement
        keys_required(1) = 'msk'
        keys_required(2) = 'pgrp'
        keys_required(3) = 'projfile'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'inner'
        keys_optional(3) = 'nspace'
        keys_optional(4) = 'nran'
        keys_optional(5) = 'shellw'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:5))
        ! set defaults
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        call xrefine3D_init%execute(cline)
    case( 'refine3D' )
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'projfile'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'vol2'
        keys_optional(3)  = 'trs'
        keys_optional(4)  = 'hp'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'cenlp'
        keys_optional(7)  = 'focusmsk'
        keys_optional(8) = 'objfun'
        keys_optional(9) = 'lpstop'
        keys_optional(10) = 'lplim_crit'
        keys_optional(11) = 'eo'
        keys_optional(12) = 'refine'
        keys_optional(13) = 'frac'
        keys_optional(14) = 'mskfile'
        keys_optional(15) = 'inner'
        keys_optional(16) = 'width'
        keys_optional(17) = 'nspace'
        keys_optional(18) = 'nstates'
        keys_optional(19) = 'startit'
        keys_optional(20) = 'maxits'
        keys_optional(21) = 'shbarrier'
        keys_optional(22) = 'noise'
        keys_optional(23) = 'nnn'
        keys_optional(24) = 'rrate'
        keys_optional(25) = 'update_frac'
        keys_optional(26) = 'shellw'
        keys_optional(27) = 'clsfrcs'
        call cline%parse_oldschool(keys_required(:4), keys_optional(:27))
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        if( .not. cline%defined('refine') )then
            call cline%set('refine',  'single')
        else
            if( cline%get_carg('refine').eq.'multi' .and. .not. cline%defined('nstates') )then
                THROW_HARD('refine=MULTI requires specification of NSTATES')
            endif
        endif
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        call xprime3D%execute(cline)
    case( 'check_3Dconv' )
        ! for convergence checking and run-time stats printing (3D)
        keys_required(1) = 'projfile'
        keys_required(2) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'nstates'
        keys_optional(3) = 'eo'
        keys_optional(4) = 'nspace'
        keys_optional(5) = 'refine'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:5))
        ! set defaults
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! execute
        call xcheck_3Dconv%execute(cline)

    ! RECONSTRUCTION PROGRAMS

    case( 'reconstruct3D' )
        keys_required(1)  = 'projfile'
        keys_required(2)  = 'pgrp'
        keys_required(3)  = 'msk'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'eo'
        keys_optional(3)  = 'frac'
        keys_optional(4)  = 'mskfile'
        keys_optional(5)  = 'shellw'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:5))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        call xreconstruct3D%execute(cline)
    case( 'volassemble_eo' )
        keys_required(1) = 'nparts'
        keys_required(2) = 'projfile'
        keys_required(3) = 'msk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'state'
        keys_optional(3) = 'nstates'
        keys_optional(4) = 'mskfile'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:4))
        call xvolassemble_eo%execute(cline)
    case( 'volassemble' )
        keys_required(1) = 'nparts'
        keys_required(2) = 'projfile'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'state'
        keys_optional(3) = 'nstates'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:3))
        call xvolassemble%execute(cline)

    ! CHECKER PROGRAMS

    case( 'check_box' )
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        call cline%parse_oldschool( keys_optional=keys_optional(:2))
        call xcheck_box%execute(cline)
    case( 'check_nptcls' )
        keys_required(1) = 'stk'
        call cline%parse_oldschool(keys_required(:1))
        call xcheck_nptcls%execute(cline)

    ! VOLOPS PROGRAMS

    case( 'postprocess' )
        ! for post-processing of volumes
        keys_required(1)  = 'msk'
        keys_required(2)  = 'projfile'
        ! set optional keys
        keys_optional(1)  = 'fsc'
        keys_optional(2)  = 'lp'
        keys_optional(3)  = 'mw'
        keys_optional(4)  = 'bfac'
        keys_optional(5)  = 'automsk'
        keys_optional(6)  = 'amsklp'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'binwidth'
        keys_optional(9)  = 'thres'
        keys_optional(10) = 'mskfile'
        keys_optional(11) = 'vol_filt'
        keys_optional(12) = 'inner'
        keys_optional(13) = 'mirr'
        ! set defaults
        call cline%parse_oldschool(keys_required(:2), keys_optional(:13))
        call xpostprocess%execute(cline)
    case( 'reproject' )
        ! for re-projecting a volume using interpolation in Fourier space
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nspace'
        keys_optional(2)  = 'outstk'
        keys_optional(3)  = 'oritab'
        keys_optional(4)  = 'nthr'
        keys_optional(5)  = 'rnd'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'pgrp'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'top'
        keys_optional(10) = 'msk'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:10))
        ! set defaults
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        call xreproject%execute(cline)
    case( 'volume_smat' )
        ! for creating a similarity matrix based on volume2volume correlation
        keys_required(1) = 'vollist'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        keys_required(4) = 'msk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'hp'
        call cline%parse_oldschool(keys_required(:4), keys_optional(:2))
        call xvolume_smat%execute(cline)
    case( 'dock_volpair' )
        ! for docking a pair of volumes. vol1 is reference and vol2 target
        keys_required(1) = 'vol1'
        keys_required(2) = 'vol2'
        keys_required(3) = 'smpd'
        keys_required(4) = 'lpstart'
        keys_required(5) = 'lpstop'
        keys_required(6) = 'msk'
        keys_required(7) = 'trs'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'outvol'
        keys_optional(3) = 'dockmode'
        call cline%parse_oldschool(keys_required(:7), keys_optional(:3))
        ! set defaults
        call xdock_volpair%execute(cline)
    case( 'symmetrize_map' )
        ! for finding the symmetry axis and average over the symmetry-related rotations
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        keys_required(4) = 'msk'
        keys_required(5) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'outvol'
        call cline%parse_oldschool(keys_required(:5), keys_optional(:2))
        call xsymmetrize_map%execute(cline)
    case( 'automask' )
        ! for volumetric envelope masking
        keys_required(1) = 'msk'
        keys_required(2) = 'amsklp'
        keys_required(3) = 'mw'
        keys_required(4) = 'thres'
        keys_required(5) = 'vol1'
        keys_required(6) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'edge'
        keys_optional(2) = 'binwidth'
        keys_optional(3) = 'nthr'
        call cline%parse_oldschool(keys_required(:6), keys_optional(:3))
        call xautomask%execute(cline)

    ! GENERAL IMAGE PROCESSING PROGRAMS

    case( 'scale' )
        !  provides re-scaling and clipping routines for MRC or SPIDER stacks and volumes
        keys_required(1)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'stk'
        keys_optional(2)  = 'vol1'
        keys_optional(3)  = 'filetab'
        keys_optional(4)  = 'msk'
        keys_optional(5)  = 'scale'
        keys_optional(6)  = 'scale2'
        keys_optional(7)  = 'newbox'
        keys_optional(8)  = 'clip'
        keys_optional(9)  = 'outvol'
        keys_optional(10) = 'outstk'
        keys_optional(11) = 'outstk2'
        call cline%parse_oldschool(keys_required(:1),keys_optional(:11))
        call xscale%execute(cline)
    case( 'binarise' )
        ! for binarisation of stacks and volumes
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'stk'
        keys_optional(3)  = 'vol1'
        keys_optional(4)  = 'thres'
        keys_optional(5)  = 'npix'
        keys_optional(6)  = 'grow'
        keys_optional(7)  = 'edge'
        keys_optional(8)  = 'neg'
        keys_optional(9)  = 'outvol'
        keys_optional(10) = 'outstk'
        call cline%parse_oldschool(keys_optional=keys_optional(:10))
        call xbinarise%execute(cline)
    case('edge_detect')
        keys_required(1) = 'detector'
        keys_required(2) = 'stk'
        keys_required(3) = 'automatic'
        keys_optional(1) = 'outstk'
        keys_optional(2) = 'thres'
        keys_optional(3) = 'npix'
        keys_optional(4) = 'thres_low'
        keys_optional(5) = 'thres_up'        
        call cline%parse_oldschool(keys_required(:3),keys_optional(:5))
        call xdetector%execute(cline)

    ! MISCELLANOUS PROGRAMS

    case( 'masscen' )
        ! for centering images acccording to their centre of mass
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'neg'
        keys_optional(3) = 'outstk'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        call xmasscen%execute(cline)
    case( 'cluster_smat' )
        ! for clustering a similarity matrix and use an combined cluster validation
        ! index to assess the quality of the clustering
        keys_required(1) = 'nptcls'
        keys_required(2) = 'fname'
        keys_required(3) = 'ncls'
        keys_required(4) = 'label'
        ! set optional keys
        keys_optional(1) = 'nthr'
        call cline%parse_oldschool(keys_required(:4), keys_optional(:1))
        call xcluster_smat%execute(cline)
    case( 'intgpeaks' )
        keys_required(1) = 'vol1'
        keys_required(2) = 'pdbfile'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'outfile'
        keys_optional(2) = 'msk'
        keys_optional(3) = 'inner'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        call xintgpeaks%execute(cline)
    case( 'print_dose_weights' )
        ! for printing the dose weights applied to individual frames
        keys_required(1) = 'nframes'
        keys_required(2) = 'exp_time'
        keys_required(3) = 'dose_rate'
        keys_required(4) = 'box'
        keys_required(5) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'kv'
        call cline%parse_oldschool(keys_required(:5),keys_optional(:1))
        call xprint_dose_weights%execute(cline)
    case( 'res' )
        ! for checking the low-pass resolution limit for a given Fourier index
        keys_required(1) = 'smpd'
        keys_required(2) = 'find'
        keys_required(3) = 'box'
        call cline%parse_oldschool(keys_required(:3))
        call xres%execute(cline)
    case( 'stk_corr' )
        ! for checking the low-pass resolution limit for a given Fourier index
        keys_required(1) = 'stk'
        keys_required(2) = 'stk2'
        keys_required(3) = 'smpd'
        keys_required(4) = 'msk'
        keys_optional(1) = 'lp'
        call cline%parse_oldschool(keys_required(:4), keys_optional(:1))
        call xstk_corr%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS

    case( 'rotmats2oris' )
        ! converts a text file (9 records per line) describing rotation matrices into a SIMPLE oritab
        keys_required(1)  = 'infile'
        ! set optional keys
        keys_optional(1)  = 'outfile'
        keys_optional(2)  = 'oritype'
        call cline%parse_oldschool( keys_required(:1), keys_optional(:2) )
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call xrotmats2oris%execute(cline)
    case( 'txt2project' )
        ! adds or replaces a text oritab in a binary *.simple project file
        keys_required(1) = 'oritab'
        keys_required(2) = 'projfile'
        keys_required(3) = 'oritype'
        call cline%parse_oldschool(keys_required(:3))
        call xtxt2project%execute(cline)
    case( 'project2txt' )
        ! converts a binary *.simple project file to a text oritab
        keys_required(1) = 'projfile'
        keys_required(2) = 'oritype'
        ! set optional keys
        keys_optional(1)  = 'outfile'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        call xproject2txt%execute(cline)
    case( 'print_project_header' )
        ! converts a binary *.simple project file to a text oritab<project2txt/end>
        keys_required(1) = 'projfile'
        call cline%parse_oldschool(keys_required(:1))
        call xprint_project_header%execute(cline)
    case( 'print_project_vals' )
        keys_required(1) = 'projfile'
        keys_required(2) = 'keys'
        keys_required(3) = 'oritype'
        call cline%parse_oldschool(keys_required(:3))
        call xprint_project_vals%execute(cline)
    case( 'update_project_stateflags' )
        keys_required(1) = 'projfile'
        keys_required(2) = 'infile'
        keys_required(3) = 'oritype'
        call cline%parse_oldschool(keys_required(:3))
        call xupdate_project_stateflags%execute(cline)
    case( 'multivariate_zscore' )
        keys_required(1) = 'keys'
        keys_required(2) = 'projfile'
        ! set optional keys
        keys_optional(2) = 'oritype'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        call xmultizscore%execute(cline)

    ! TIME-SERIES ANALYSIS PROGRAMS

     case( 'tseries_extract' )
        ! for creating overlapping chunks of nframesgrp frames from time-series data
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'nframesgrp'
        call cline%parse_oldschool(keys_required(:3))
        call xtseries_extract%execute(cline)
    case( 'tseries_track' )
        ! for particle tracking in time-series data
        keys_required(1) = 'filetab'
        keys_required(2) = 'fbody'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'boxfile'
        keys_optional(3) = 'xcoord'
        keys_optional(4) = 'ycoord'
        keys_optional(5) = 'offset'
        keys_optional(6) = 'box'
        keys_optional(7) = 'neg'
        keys_optional(8) = 'cenlp'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        call xtseries_track%execute(cline)
    case('tseries_backgr_subtr')
        ! for background subtraction in time-series data. The goal is to subtract the two graphene
        ! peaks @ 2.14 A and @ 1.23 A. This is done by band-pass filtering the background image,
        ! recommended (and default settings) are hp=5.0 lp=1.1 and width=5.0.
        keys_required(1) = 'stk'
        keys_required(2) = 'stk_backgr'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'hp'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'width'
        keys_optional(5) = 'deftab'
        keys_optional(6) = 'outstk'
        call cline%parse_oldschool(keys_required(:3), keys_optional(:6))
        ! set defaults
        if( .not. cline%defined('hp')    ) call cline%set('hp',    5.0)
        if( .not. cline%defined('lp')    ) call cline%set('lp',    1.1)
        if( .not. cline%defined('width') ) call cline%set('width', 5.0)
        call xtseries_backgr_subtr%execute(cline)
    case( 'tseries_split' )
        ! for splitting a time-series stack and its associated orientations
        keys_required(1) = 'stk'
        keys_required(2) = 'oritab'
        keys_required(3) = 'smpd'
        keys_required(4) = 'chunksz'
        keys_required(5) = 'stepsz'
        call cline%parse_oldschool(keys_required(:5))
        call xtseries_split%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS

    case( 'merge_nnmat' )
        ! for merging partial nearest neighbour matrices calculated in distributed mode
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        keys_required(3) = 'nnn'
        call cline%parse_oldschool( keys_required(:3) )
        call xmerge_nnmat%execute(cline)
    case( 'merge_similarities' )
        ! for merging similarities calculated between pairs of objects into a
        ! similarity matrix that can be inputted to cluster_smat
        keys_required(1) = 'nptcls'
        ! set optional keys
        keys_optional(1) = 'nparts'
        call cline%parse_oldschool(keys_required(:1), keys_optional(:1))
        call xmerge_similarities%execute(cline)
    case( 'split_pairs' )
        ! for splitting calculations between pairs of objects into balanced partitions
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        call cline%parse_oldschool(keys_required(:2))
        call xsplit_pairs%execute(cline)
    case( 'split' )
        ! for splitting of image stacks into partitions for parallel execution
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'nparts'
        call cline%parse_oldschool(keys_required=keys_required(:3))
        call xsplit%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
    end select

end program simple_private_exec
