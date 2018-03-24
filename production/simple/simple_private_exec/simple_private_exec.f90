! shared-memory parallelised programs executed by distributed commanders
program simple_private_exec
#include "simple_lib.f08"
use simple_cmdline, only: cmdline, cmdline_err
use simple_gen_doc
use simple_commander_checks
use simple_commander_comlin
use simple_commander_distr
use simple_commander_imgproc
use simple_commander_mask
use simple_commander_misc
use simple_commander_oris
use simple_commander_preprocess
use simple_commander_cluster2D
use simple_commander_refine3D
use simple_commander_project
use simple_commander_rec
use simple_commander_sim
use simple_commander_volops
use simple_commander_tseries
implicit none

! PRE-PROCESSING PROGRAMS
type(preprocess_commander)           :: xpreprocess
type(select_frames_commander)        :: xselect_frames
type(boxconvs_commander)             :: xboxconvs
type(powerspecs_commander)           :: xpowerspecs
type(motion_correct_commander)       :: xmotion_correct
type(ctf_estimate_commander)         :: xctf_estimate
type(pick_commander)                 :: xpick

! CLUSTER2D PROGRAMS
type(make_cavgs_commander)           :: xmake_cavgs
type(cluster2D_commander)            :: xcluster2D
type(cavgassemble_commander)         :: xcavgassemble
type(check2D_conv_commander)         :: xcheck2D_conv
type(rank_cavgs_commander)           :: xrank_cavgs

! PRIME3D PROGRAMS
type(npeaks_commander)               :: xnpeaks
type(nspace_commander)               :: xnspace
type(refine3D_init_commander)         :: xrefine3D_init
type(rec_test_commander)             :: xrec_test
type(multiptcl_init_commander)       :: xmultiptcl_init
type(prime3D_commander)              :: xprime3D
type(check3D_conv_commander)         :: xcheck3D_conv

! COMMON-LINES PROGRAMS
type(symsrch_commander)              :: xsymsrch

! SYMMETRY PROGRAMS
type(sym_aggregate_commander)        :: xsym_aggregate
type(dsymsrch_commander)             :: xdsymsrch

! MASK PROGRAMS
type(resmask_commander)              :: xresmask

! RECONSTRUCTION PROGRAMS
type(volassemble_eo_commander)       :: xvolassemble_eo
type(reconstruct3D_commander)        :: xreconstruct3D
type(volassemble_commander)          :: xvolassemble

! CHECKER PROGRAMS
type(check_box_commander)            :: xcheck_box
type(check_nptcls_commander)         :: xcheck_nptcls

! VOLOPS PROGRAMS
type(postprocess_commander)          :: xpostprocess ! DUPLICATED
type(project_commander)              :: xproject     ! DUPLICATED
type(volaverager_commander)          :: xvolaverager
type(volume_smat_commander)          :: xvolume_smat
type(dock_volpair_commander)         :: xdock_volpair

! GENERAL IMAGE PROCESSING PROGRAMS
type(scale_commander)                :: xscale       ! DUPLICATED
type(binarise_commander)             :: xbinarise
type(corrcompare_commander)          :: xcorrcompare
type(image_diff_commander)           :: ximage_diff
type(image_smat_commander)           :: ximage_smat

! MISCELLANOUS PROGRAMS
type(masscen_commander)              :: xmasscen
type(cluster_smat_commander)         :: xcluster_smat
type(intgpeaks_commander)            :: xintgpeaks
type(print_dose_weights_commander)   :: xprint_dose_weights
type(res_commander)                  :: xres

! ORIENTATION DATA MANAGEMENT PROGRAMS
type(map2ptcls_commander)            :: xmap2ptcls
type(cluster_oris_commander)         :: xcluster_oris
type(rotmats2oris_commander)         :: xrotmats2oris
type(txt2project_commander)          :: xtxt2project
type(project2txt_commander)          :: xproject2txt
type(manage_project_commander)       :: xmanage_project
type(print_project_info_commander)   :: xprint_project_info

! TIME-SERIES ANALYSIS PROGRAMS
type(tseries_extract_commander)      :: xtseries_extract
type(tseries_track_commander)        :: xtseries_track
type(tseries_backgr_subtr_commander) :: xtseries_backgr_subtr
type(tseries_split_commander)        :: xtseries_split

! PARALLEL PROCESSING PROGRAMS
type(merge_algndocs_commander)       :: xmerge_algndocs
type(merge_nnmat_commander)          :: xmerge_nnmat
type(merge_similarities_commander)   :: xmerge_similarities
type(split_pairs_commander)          :: xsplit_pairs
! type(split_commander)                :: xsplit

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
if( str_has_substr(prg, 'simple') ) stop 'giving program names with simple* prefix is depreciated'
select case(prg)

    ! PRE-PROCESSING PROGRAMS

    case( 'preprocess' )
        !==Program preprocess
        !
        ! <preprocess/begin>is a program that executes motion_correct, ctf_estimate and pick in sequence
        ! <preprocess/end>
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
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:30))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',              5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',           8.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',        512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate', 30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',  5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',         20.)
        if( .not. cline%defined('outfile')         ) call cline%set('outfile', 'simple_unidoc'//METADATA_EXT)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('opt')             ) call cline%set('opt',        'simplex')
        call xpreprocess%execute(cline)
    case( 'select_frames' )
        !==Program select_frames
        !
        ! <select_frames/begin>is a program for selecting contiguous segments of frames from DDD movies
        ! <select_frames/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'fbody'
        keys_required(3) = 'fromf'
        keys_required(4) = 'tof'
        keys_required(5) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'startit'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5), keys_optional(:1))
        ! execute
        call xselect_frames%execute(cline)
    case( 'boxconvs' )
        !==Program boxconvs
        !
        ! <boxconvs/begin>is a program for averaging overlapping boxes across a micrograph
        ! in order to check if gain correction was appropriately done<boxconvs/end>
        !
        ! set required keys
        keys_required(1) = 'fbody'
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'filetab'
        keys_optional(3) = 'boxconvsz'
        keys_optional(4) = 'startit'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('boxconvsz') ) call cline%set('boxconvsz', 512.)
        ! execute
        call xboxconvs%execute(cline)
    case( 'powerspecs' )
        !==Program powerspecs
        !
        ! <powerspecs/begin>is a program for generating powerspectra from a stack or filetable<powerspecs/end>
        !
        ! set required keys
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        ! execute
        call xpowerspecs%execute(cline)
    case( 'motion_correct' )
        !==Program motion_correct
        !
        ! <motion_correct/begin>is a program for movie alignment or motion_correctring based the same principal strategy as
        ! Grigorieffs program (hence the name). There are two important differences: automatic weighting of
        ! the frames using a correlation-based M-estimator and continuous optimisation of the shift parameters.
        ! Input is a textfile with absolute paths to movie files in addition to a few input parameters, some
        ! of which deserve a comment. If dose_rate and exp_time are given the individual frames will be
        ! low-pass filtered accordingly (dose-weighting strategy). If scale is given, the movie will be Fourier
        ! cropped according to the down-scaling factor (for super-resolution movies). If nframesgrp is given
        ! the frames will be pre-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given,
        ! a contiguous subset of frames will be averaged without any dose-weighting applied.
        ! <motion_correct/end>
        !

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
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:16))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'simple_unidoc'//METADATA_EXT)
        if( .not. cline%defined('opt')     ) call cline%set('opt', 'simplex')
        ! execute
        call xmotion_correct%execute(cline)
    case( 'ctf_estimate' )
        !==Program ctf_estimate
        !
        ! <ctf_estimate/begin>is a SIMPLE program for fitting the CTF<ctf_estimate/end>
        !
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
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        ! execute
        call xctf_estimate%execute(cline)
    case( 'motion_correct_ctf_estimate' )
        !==Program motion_correct_ctf_estimate
        !
        ! <motion_correct_ctf_estimate/begin>is a pipelined motion_correct + ctf_estimate program<motion_correct_ctf_estimate/end>
        !
        ! set required keys
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
        ! parse command line
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
        if( .not. cline%defined('outfile')         ) call cline%set('outfile', 'simple_unidoc'//METADATA_EXT)
        if( .not. cline%defined('opt')             ) call cline%set('opt', 'simplex')
        ! execute
        call xpreprocess%execute(cline)
    case( 'pick' )
        !==Program pick
        !
        ! <pick/begin>is a template-based picker program<pick/end>
        !
        ! set required keys
        keys_required(1) = 'refs'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'ndev'
        keys_optional(5) = 'dir'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1), keys_optional(:5))
        ! execute
        call xpick%execute(cline)

    ! CLUSTER2D PROGRAMS

    case( 'make_cavgs' )
        !==Program make_cavgs
        !
        ! <make_cavgs/begin>is used  to produce class averages or initial random references
        ! for cluster2D execution. <make_cavgs/end>
        !
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncls'
        keys_optional(3)  = 'filwidth'
        keys_optional(4)  = 'mul'
        keys_optional(5)  = 'tseries'
        keys_optional(6)  = 'outfile'
        keys_optional(7)  = 'refs'
        keys_optional(8)  = 'remap_cls'
        keys_optional(9)  = 'weights2D'
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        ! execute
        call xmake_cavgs%execute(cline)
    case( 'cluster2D' )
        !==Program cluster2D
        !
        ! <cluster2D/begin>is a reference-free 2D alignment/clustering algorithm adopted from the prime3D
        ! probabilistic ab initio 3D reconstruction algorithm<cluster2D/end>
        !
        ! set required keys
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
        keys_optional(16) = 'weights2D'
        keys_optional(17) = 'refine'
        keys_optional(18) = 'match_filt'
        keys_optional(19) = 'dyncls'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:19))
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',     30.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',       'no')
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',    30.)
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D','no')
        ! execute
        call xcluster2D%execute(cline)
    case( 'cavgassemble' )
        !==Program cavgassemble
        !
        ! <cavgassemble/begin>is a program that assembles class averages when the clustering
        ! program (cluster2D) has been executed in distributed mode<cavgassemble/end>
        !
        ! set required keys
        keys_required(1) = 'projfile'
        keys_required(2) = 'nparts'
        keys_required(3) = 'ncls'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'refs'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:2))
        ! execute
        call xcavgassemble%execute(cline)
    case( 'check2D_conv' )
        !==Program check2D_conv
        !
        ! <check2D_conv/begin>is a program for checking if a cluster2D run has converged.
        ! The statistics outputted include (1) the overlap between the distribution of parameters
        ! for succesive runs. (2) The percentage of search space scanned, i.e. how many reference
        ! images are evaluated on average. (3) The average correlation between the images and
        ! their corresponding best matching reference section. If convergence to a local optimum
        ! is achieved, the fraction increases. Convergence is achieved if the parameter distribution
        ! overlap is larger than 0.95 and more than 99% of the reference sections need to be
        ! searched to find an improving solution<check2D_conv/end>
        !
        ! set required keys
        keys_required(1) = 'projfile'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1))
        ! execute
        call xcheck2D_conv%execute(cline)
    case( 'rank_cavgs' )
        !==Program rank_cavgs
        !
        ! <rank_cavgs/begin>is a program for ranking class averages by decreasing population, given the
        ! stack of class averages (stk argument) and the 2D orientations document (oritab)
        ! generated by cluster2D<rank_cavgs/end>
        !
        ! set required keys
        keys_required(1) = 'projfile'
        keys_required(2) = 'stk'
        ! set optional keys
        keys_optional(1) = 'outstk'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        ! set defaults
        call cline%set('oritype', 'cls2D')
        ! execute
        call xrank_cavgs%execute(cline)

    ! PRIME3D PROGRAMS

    case( 'npeaks' )
        !==Program npeaks
        !
        ! <npeaks/begin>is a program for checking the number of nonzero orientation weights (number of correlation peaks
        ! included in the weighted reconstruction)<npeaks/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'box'
        keys_required(3) = 'lp'
        ! set optional keys
        keys_optional(1) = 'nspace'
        keys_optional(2) = 'moldiam'
        keys_optional(3) = 'pgrp'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        ! set defaults
        if( .not. cline%defined('lp') )     call cline%set('lp',       20.)
        ! execute
        call xnpeaks%execute(cline)
    case( 'nspace' )
        !==Program nspace
        !
        ! <nspace/begin>is a program for calculating the expected resolution obtainable with different values of nspace
        ! (number of discrete projection directions used for discrete search)<nspace/end>
        !
        ! set required keys
        keys_required(1)  = 'moldiam'
        ! parse command line
        call cline%parse_oldschool(keys_required=keys_required(:1))
        ! execute
        call xnspace%execute(cline)
    case( 'refine3D_init' )
        !==Program refine3D_init
        !
        ! <refine3D_init/begin>is a program for generating a random initial model for initialisation of PRIME3D.
        ! If the data set is large (>5000 images), generating a random model can be slow. To speedup, set
        ! nran to some smaller number, resulting in nran images selected randomly for
        ! reconstruction<refine3D_init/end>
        !
        ! set required keys

        keys_required(1) = 'msk'
        keys_required(2) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'inner'
        keys_optional(3) = 'nspace'
        keys_optional(4) = 'nran'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        ! execute
        call xrefine3D_init%execute(cline)
    case( 'multiptcl_init' )
        !==Program multiptcl_init
        !
        ! <multiptcl_init/begin>is a program for generating random initial models for initialisation of PRIME3D
        ! when run in multiparticle mode<multiptcl_init/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'ctf'
        keys_required(3)  = 'pgrp'
        keys_required(4)  = 'nstates'
        keys_required(5)  = 'msk'
        keys_required(6)  = 'oritab'
        ! set optionnal keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'inner'
        keys_optional(4)  = 'width'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'eo'
        keys_optional(7)  = 'frac'
        keys_optional(8)  = 'state2split'
        keys_optional(9)  = 'norec'
        keys_optional(10) = 'mul'
        keys_optional(11) = 'zero'
        keys_optional(12) = 'tseries'
        keys_optional(13) = 'center'
        keys_optional(14) = 'stk'
        keys_optional(15) = 'stktab'
        keys_optional(16) = 'phaseplate'
        ! parse command line
        call cline%parse_oldschool(keys_required(:6), keys_optional(:16))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 3.) ! to assure that shifts are being used
        !execute
        call xmultiptcl_init%execute(cline)
    case( 'refine3D' )
        !==Program refine3D
        !
        ! <refine3D/begin>is an ab inito reconstruction/refinement program based on probabilistic
        ! projection matching. There are a daunting number of options in refine3D. If you
        ! are processing class averages we recommend that you instead use the simple_distr_exec prg=
        ! initial_3Dmodel route.<refine3D/end>
        !
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'vol2'
        keys_optional(3)  = 'oritab'
        keys_optional(4)  = 'deftab'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'hp'
        keys_optional(7)  = 'lp'
        keys_optional(8)  = 'cenlp'
        keys_optional(9) = 'focusmsk'
        keys_optional(10) = 'objfun'
        keys_optional(11) = 'lpstop'
        keys_optional(12) = 'lplim_crit'
        keys_optional(13) = 'eo'
        keys_optional(14) = 'refine'
        keys_optional(15) = 'frac'
        keys_optional(16) = 'mskfile'
        keys_optional(17) = 'inner'
        keys_optional(18) = 'width'
        keys_optional(19) = 'nspace'
        keys_optional(20) = 'nstates'
        keys_optional(21) = 'startit'
        keys_optional(22) = 'maxits'
        keys_optional(23) = 'shbarrier'
        keys_optional(24) = 'noise'
        keys_optional(25) = 'nnn'
        keys_optional(26) = 'rrate'
        keys_optional(27) = 'update_frac'
        keys_optional(28) = 'stk'
        keys_optional(29) = 'stktab'
        keys_optional(30) = 'phaseplate'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5), keys_optional(:30))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine') )then
            call cline%set('refine',  'single')
        else
            if( cline%get_carg('refine').eq.'multi' )then
                if( .not. cline%defined('nstates') ) stop 'refine=MULTI requires specification of NSTATES'
                if( .not. cline%defined('oritab')  ) stop 'refine=MULTI requires ORITAB input'
            endif
        endif
        ! execute
        call xprime3D%execute(cline)
    case( 'rec_test' )
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'msk'
        keys_required(3) = 'ctf'
        keys_required(4) = 'pgrp'
        keys_required(5) = 'oritab'
        keys_required(6) = 'stk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        call cline%parse_oldschool(keys_required(:6), keys_optional(:1))
        ! execute
        call xrec_test%execute(cline)
    case( 'check3D_conv' )
        !==Program check3D_conv
        !
        ! <check3D_conv/begin>is a program for checking if a PRIME3D run has converged. The statistics
        ! outputted include (1) angle of feasible region, which is proportional to the angular
        ! resolution of the set of discrete projection directions being searched. (2) The average angular
        ! distance between orientations in the present and previous iteration. In the early iterations,
        ! the distance is large because a diverse set of orientations is explored. If convergence to a
        ! local optimum is achieved, the distance decreases. (3) The percentage of search space scanned,
        ! i.e. how many reference images are evaluated on average. (4) The average correlation between
        ! the images and their corresponding best matching reference sections. (5) The average standard
        ! deviation of the Euler angles. Convergence is achieved if the angular distance between the
        ! orientations in successive iterations falls significantly below the angular resolution of the
        ! search space and more than 99% of the reference sections need to be matched on average
        ! <check3D_conv/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'box'
        keys_required(3) = 'oritab'
        keys_required(4) = 'nptcls'
        keys_required(5) = 'pgrp'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'nstates'
        keys_optional(3) = 'eo'
        keys_optional(4) = 'nspace'
        keys_optional(5) = 'find'
        keys_optional(6) = 'refine'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5), keys_optional(:6))
        ! set defaults
        if( .not. cline%defined('lp') .and. .not.cline%defined('lp') )call cline%set('lp', 20.)
        ! execute
        call xcheck3D_conv%execute(cline)

    ! COMMON-LINES PROGRAMS

    case( 'symsrch' )
        !==Program symsrch
        !
        ! <symsrch/begin>is a program for searching for the principal symmetry axis of a volume
        ! reconstructed without assuming any point-group symmetry. The program takes as input an
        ! asymmetrical 3D reconstruction. The alignment document for all the particle images
        ! that have gone into the 3D reconstruction and the desired point-group symmetry needs to
        ! be inputted. The 3D reconstruction is then projected in 50 (default option) even directions,
        ! common lines-based optimisation is used to identify the principal symmetry axis, the rotational
        ! transformation is applied to the inputted orientations, and a new alignment document is produced.
        ! Input this document to reconstruct3D together with the images and the point-group symmetry to generate a
        ! symmetrised map.<symsrch/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'pgrp'
        keys_required(5) = 'oritab'
        keys_required(6) = 'outfile'
        keys_required(7) = 'lp'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'cenlp'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'nspace'
        keys_optional(5) = 'center'
        ! parse command line
        call cline%parse_oldschool(keys_required(:7), keys_optional(:5))
        ! set defaults
        if( .not. cline%defined('nspace') )then
            call cline%set('nptcls', 50.) ! 50 projections 4 symsrch
            call cline%set('nspace', 50.) ! 50 projections 4 symsrch
        else
            call cline%set('nptcls', cline%get_rarg('nspace'))
        endif
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp', 30.)
        call cline%set('compare', 'no')
        ! execute
        call xsymsrch%execute(cline)

    ! SYMMETRY PROGRAMs

    case( 'sym_aggregate' )
        !==Program symsrch
        !
        ! <sym_aggregate/begin>is a program for robust identifiaction of the symmetry axis
        ! of a map using image-to-volume simiarity validation of the axis<sym_aggregate/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'pgrp'
        keys_required(5) = 'oritab'
        keys_required(6) = 'oritab2'
        keys_required(7) = 'outfile'
        keys_required(8) = 'lp'
        keys_required(9) = 'stk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'cenlp'
        keys_optional(3) = 'hp'
        ! parse command line
        call cline%parse_oldschool(keys_required(:9), keys_optional(:3))
        ! set defaults
        call cline%set('eo','no')
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        ! execute
        call xsym_aggregate%execute(cline)
    case( 'dsymsrch' )
        !==Program symsrch
        !
        ! <dsymsrch/begin>is a program for identifying rotational symmetries in class averages of
        ! D-symmetric molecules and generating a cylinder that matches the shape.<dsymsrch/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'msk'
        keys_required(3) = 'pgrp'
        keys_required(4) = 'stk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'cenlp'
        keys_optional(3) = 'outfile'
        keys_optional(4) = 'outvol'
        ! parse command line
        call cline%parse_oldschool(keys_required(:4), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 30.)
        ! execute
        call xdsymsrch%execute(cline)

    ! MASK PROGRAMS

    case( 'resmask' )
        !==Program resmask
        !
        ! <resmask/begin>is a program for 3D envelope masking for resolution estimation<resmask/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'msk'
        keys_required(3) = 'mskfile'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3))
        ! execute
        call xresmask%execute(cline)

    ! RECONSTRUCTION PROGRAMS

    case( 'reconstruct3D' )
        ! set required keys
        keys_required(1)  = 'projfile'
        keys_required(2)  = 'pgrp'
        keys_required(3)  = 'msk'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'eo'
        keys_optional(3)  = 'frac'
        keys_optional(4)  = 'mskfile'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:4))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        ! execute
        call xreconstruct3D%execute(cline)
    case( 'volassemble_eo' )
        !==Program volassemble_eo
        !
        ! <volassemble_eo/begin>is a program that assembles volume(s) when the reconstruction
        ! program (reconstruct3D with eo=yes) has been executed in distributed mode. inner applies a soft-edged
        ! inner mask. An inner mask is used for icosahedral virus reconstruction, because the
        ! DNA or RNA core is often unordered and  if not removed it may negatively impact the
        ! alignment. The width parameter controls the fall-off of the edge of the
        ! inner mask<volassemble_eo/end>
        !
        ! set required keys
        keys_required(1) = 'nparts'
        keys_required(2) = 'projfile'
        keys_required(3) = 'msk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'state'
        keys_optional(3) = 'nstates'
        keys_optional(4) = 'mskfile'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:4))
        ! execute
        call xvolassemble_eo%execute(cline)
    case( 'volassemble' )
        !==Program volassemble
        !
        ! <volassemble/begin>is a program that assembles volume(s) when the reconstruction program
        ! (reconstruct3D) has been executed in distributed mode. odd is used to assemble the odd reconstruction,
        ! even is used to assemble the even reconstruction, eo is used to assemble both the even and the
        ! odd reconstruction and state is used to assemble the inputted state. Normally, you do not fiddle with
        ! these parameters. They are used internally<volassemble/end>
        !
        ! set required keys
        keys_required(1) = 'nparts'
        keys_required(2) = 'projfile'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'state'
        keys_optional(3) = 'nstates'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:3))
        ! execute
        call xvolassemble%execute(cline)

    ! CHECKER PROGRAMS

    case( 'check_box' )
        !==Program check_box
        !
        ! <check_box/begin>is a program for checking the image dimensions of MRC and SPIDER
        !  stacks and volumes<check_box/end>
        !
        ! set optional keys
        keys_optional(1) = 'stk'
        keys_optional(2) = 'vol1'
        ! parse command line
        call cline%parse_oldschool( keys_optional=keys_optional(:2))
        ! execute
        call xcheck_box%execute(cline)
    case( 'check_nptcls' )
        !==Program check_nptcls
        !
        ! <check_nptcls/begin>is a program for checking the number of images in MRC and SPIDER
        ! stacks<check_nptcls/end>
        !
        ! set optional keys
        keys_required(1) = 'stk'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1))
        ! execute
        call xcheck_nptcls%execute(cline)

    ! VOLOPS PROGRAMS

    case( 'postprocess' )
        !==Program postprocess
        !
        ! <postprocess/begin>is a program for post-processing of volumes<postprocess/end>
        !
        ! set required keys
        keys_required(1)  = 'vol1'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:13))
        ! execute
        call xpostprocess%execute(cline)
    case( 'project' )
        !==Program project
        !
        ! <project/begin>is a program for projecting a volume using interpolation in Fourier space. Input is a SPIDER or
        ! MRC volume. Output is a stack of projection images of the same format as the inputted volume. Projections
        ! are generated by extraction of central sections from the Fourier volume and back transformation of the 2D FTs.
        ! nspace controls the number of projection images generated with quasi-even projection directions. The
        ! oritab parameter allows you to input the orientations that you wish to have your volume projected in. If
        ! rnd=yes, random rather than quasi-even projections are generated, trs then controls the halfwidth of
        ! the random origin shift. Less commonly used parameters are pgrp, which controls the point-group symmetry
        ! c (rotational), d (dihedral), t (tetrahedral), o (octahedral) or i (icosahedral). The point-group symmetry is
        ! used to restrict the set of projections to within the asymmetric unit.
        ! neg inverts the contrast of the projections. <project/end>
        !
        ! set required keys
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:10))
        ! set defaults
        if( .not. cline%defined('wfun')  ) call cline%set('wfun', 'kb')
        if( .not. cline%defined('winsz') ) call cline%set('winsz', 1.5)
        if( .not. cline%defined('alpha') ) call cline%set('alpha', 2.)
        ! execute
        call xproject%execute(cline)
    case( 'volaverager' )
        !==Program volaverager
        !
        ! <volaverager/begin>is a program for averaging volumes according to state label in oritab
        ! <volaverager/end>
        !
        ! set required keys
        keys_required(1) = 'vollist'
        keys_required(2) = 'oritab'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        ! execute
        call xvolaverager%execute(cline)
    case( 'volume_smat' )
        !==Program volume_smat
        !
        ! <volume_smat/begin>is a program for creating a similarity matrix based on volume2volume
        ! correlation<volume_smat/end>
        !
        ! set required keys
        keys_required(1) = 'vollist'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        keys_required(4) = 'msk'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'hp'
        ! parse command line
        call cline%parse_oldschool(keys_required(:4), keys_optional(:2))
        ! execute
        call xvolume_smat%execute(cline)
    case( 'dock_volpair' )
        !==Program dock_volpair
        !
        ! <dock_volpair/begin>is a program for docking a pair of volumes. vol1 is reference and vol2 target.
        ! <dock_volpair/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'vol2'
        keys_required(3) = 'smpd'
        keys_required(4) = 'lp'
        keys_required(5) = 'msk'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'dockmode'
        keys_optional(3) = 'outvol'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5), keys_optional(:3))
        ! execute
        call xdock_volpair%execute(cline)

    ! GENERAL IMAGE PROCESSING PROGRAMS

    case( 'scale' )
        !==Program scale
        !
        ! <scale/begin>is a program that provides re-scaling and clipping routines for MRC or SPIDER stacks
        ! and volumes<scale/end>
        !
        ! set required keys
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:1),keys_optional(:11))
        ! execute
        call xscale%execute(cline)
    case( 'binarise' )
        !==Program binarise
        !
        ! <binarise/begin>is a program for binarisation of stacks and volumes<binarise/end>
        !
        ! set optional keys
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
        ! parse command line
        call cline%parse_oldschool(keys_optional=keys_optional(:10))
        ! execute
        call xbinarise%execute(cline)
    case( 'corrcompare' )
        !==Program corrcompare
        !
        ! <corrcompare/begin>is a program for comparing stacked images using real-space and Fourier-based approaches
        ! <corrcompare/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'stk2'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'stats'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'smpd'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:4))
        ! execute
        call xcorrcompare%execute(cline)
     case( 'image_diff' )
        !==Program corrcompare
        !
        ! <image_diff/begin>is a program for comparing stacked images using differences
        ! <image_diff/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'stk2'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'stats'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'smpd'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:4))
        ! execute
        call ximage_diff%execute(cline)
    case( 'image_smat' )
        !==Program image_smat
        !
        ! <image_smat/begin>is a program for creating a similarity matrix based on common line correlation. The idea
        ! being that it should be possible to cluster images based on their 3D similarity witout having a 3D model
        ! by only operating on class averages and find averages that fit well together in 3D<image_smat/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'msk'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'nthr'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:4))
        ! execute
        call ximage_smat%execute(cline)

    ! MISCELLANOUS PROGRAMS

    case( 'masscen' )
        !==Program masscen
        !
        ! <masscen/begin>is a program for centering images acccording to their
        ! centre of mass<masscen/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        ! set optional keys
        keys_optional(1) = 'msk'
        keys_optional(2) = 'neg'
        keys_optional(3) = 'outstk'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        ! execute
        call xmasscen%execute(cline)
    case( 'cluster_smat' )
        !==Program cluster_smat
        !
        ! <cluster_smat/begin>is a program for clustering a similarity matrix and use
        ! an combined cluster validation index to assess the quality of the clustering
        ! based on the number of clusters<cluster_smat/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'fname'
        keys_required(3) = 'ncls'
        keys_required(4) = 'label'
        ! set optional keys
        keys_optional(1) = 'nthr'
        ! parse command line
        call cline%parse_oldschool(keys_required(:4), keys_optional(:1))
        ! execute
        call xcluster_smat%execute(cline)
    case( 'intgpeaks' )
        !==Program intgpeaks
        !
        ! <intgpeaks/begin>is a program for <intgpeaks/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'pdbfile'
        keys_required(3) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'outfile'
        keys_optional(2) = 'msk'
        keys_optional(3) = 'inner'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        ! execute
        call xintgpeaks%execute(cline)
    case( 'print_dose_weights' )
        !==Program print_dose_weights
        !
        ! <print_dose_weights/begin>is a program for printing the dose weights applied to individual frames<print_dose_weights/end>
        !
        ! set required keys
        keys_required(1) = 'nframes'
        keys_required(2) = 'exp_time'
        keys_required(3) = 'dose_rate'
        keys_required(4) = 'box'
        keys_required(5) = 'smpd'
        ! set optional keys
        keys_optional(1) = 'kv'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5),keys_optional(:1))
        ! execute
        call xprint_dose_weights%execute(cline)
    case( 'res' )
        !==Program res
        !
        ! <res/begin>is a program for checking the low-pass resolution limit for a given Fourier index<res/end>
        !
        !set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'find'
        keys_required(3) = 'box'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3))
        !execute
        call xres%execute(cline)

    ! ORIENTATION DATA MANAGEMENT PROGRAMS

    case( 'cluster_oris' )
        !==Program cluster_oris
        !
        ! <cluster_oris/begin>is a program for clustering orientations based on geodesic distance<cluster_oris/end>
        !
        ! Required keys
        keys_required(1) = 'oritab'
        keys_required(2) = 'ncls'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'oritype'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2), keys_optional(:2))
        ! execute
        call xcluster_oris%execute(cline)
    case( 'map2ptcls_doc' )
        !==Program map2ptcls
       !
       ! <map2ptcls/begin>is a program for mapping parameters that have been obtained using class averages to
       ! the individual particle images<map2ptcls/end>
       !
       ! set required keys
       keys_required(1) = 'oritab'
       ! set optional keys
       keys_optional(1) = 'nthr'
       keys_optional(2) = 'oritab3D'
       keys_optional(3) = 'deftab'
       keys_optional(4) = 'outfile'
       keys_optional(5) = 'mul'
       ! parse command line
       call cline%parse_oldschool(keys_required(:1), keys_optional(:5))
       ! set defaults
       if( .not. cline%defined('outfile') ) call cline%set('outfile', 'mapped_ptcls_params.txt')
       ! execute
       call xmap2ptcls%execute(cline)
    case( 'rotmats2oris' )
        !==Program rotmats2oris
        !
        ! <rotmats2oris/begin>converts a text file (9 records per line) describing
        ! rotation matrices into a SIMPLE oritab<rotmats2oris/end>
        !
        ! Required keys
        keys_required(1)  = 'infile'
        ! set optional keys
        keys_optional(1)  = 'outfile'
        keys_optional(2)  = 'oritype'
        ! parse command line
        call cline%parse_oldschool( keys_required(:1), keys_optional(:2) )
        ! set defaults
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        ! execute
        call xrotmats2oris%execute(cline)
    case( 'txt2project' )
        !==Program txt2project
        !
        ! <txt2project/begin>adds or replaces a text oritab in a binary *.simple project file<txt2project/end>
        !
        ! Required keys
        keys_required(1) = 'oritab'
        keys_required(2) = 'projfile'
        keys_required(3) = 'oritype'
        call cline%parse_oldschool(keys_required(:3))
        ! execute
        call xtxt2project%execute(cline)
    case( 'project2txt' )
        !==Program project2txt
        !
        ! <project2txt/begin>converts a binary *.simple project file to a text oritab<project2txt/end>
        !
        ! Required keys
        keys_required(1) = 'projfile'
        keys_required(2) = 'oritype'
        ! set optional keys
        keys_optional(1)  = 'outfile'
        call cline%parse_oldschool(keys_required(:2), keys_optional(:1))
        ! execute
        call xproject2txt%execute(cline)
    case( 'manage_project' )
        !==Program manageproject
        !
        ! </begin></end>
        !
        ! set required keys
        keys_required(1) = 'ctf'
        ! set optional keys
        keys_optional(1)  = 'smpd'
        keys_optional(2)  = 'cs'
        keys_optional(3)  = 'kv'
        keys_optional(4)  = 'fraca'
        keys_optional(5)  = 'projfile'
        keys_optional(6)  = 'projname'
        keys_optional(7)  = 'phaseplate'
        keys_optional(8)  = 'user_email'
        keys_optional(9)  = 'time_per_image'
        keys_optional(10) = 'user_account'
        keys_optional(11) = 'user_project'
        keys_optional(12) = 'qsys_partition'
        keys_optional(13) = 'qsys_qos'
        keys_optional(14) = 'qsys_reservation'
        keys_optional(15) = 'job_memory_per_task'
        keys_optional(16) = 'stk'
        keys_optional(17) = 'stktab'
        keys_optional(18) = 'plaintexttab'
        keys_optional(19) = 'oritab'
        keys_optional(20) = 'deftab'
        keys_optional(21) = 'dfunit'
        keys_optional(22) = 'angastunit'
        keys_optional(23) = 'phshiftunit'
        keys_optional(24) = 'oritype'
        keys_optional(25) = 'filetab'
        ! parse command line
        call cline%parse_oldschool(keys_required, keys_optional(:25))
        ! execute
        call xmanage_project%execute(cline)
    case( 'print_project_info' )
        !==Program print_project_info
        !
        ! <print_project_info/begin>prints information abourt a *.simple project file<print_project_info/end>
        !
        ! Required keys
        keys_required(1) = 'projfile'
        call cline%parse_oldschool(keys_required(:1))
        ! execute
        call xprint_project_info%execute(cline)

    ! TIME-SERIES ANALYSIS PROGRAMS

     case( 'tseries_extract' )
        !==Program tseries_extract
        !
        ! <tseries_extract/begin>is a program for creating overlapping chunks of nframesgrp frames from time-series data
        ! <tseries_extract/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'nframesgrp'
        ! parse command line
        call cline%parse_oldschool(keys_required(:3))
        ! execute
        call xtseries_extract%execute(cline)
    case( 'tseries_track' )
        !==Program tseries_track
        !
        ! <tseries_track/begin>is a program for particle tracking in time-series data
        ! <tseries_track/end>
        !
        ! set required keys
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        ! execute
        call xtseries_track%execute(cline)
    case('tseries_backgr_subtr')
        !==Program tseries_backgr_subtr
        !
        ! <tseries_backgr_subtr/begin>is a program for background subtraction in time-series data.
        ! The goal is to subtract the two graphene peaks @ 2.14 A and @ 1.23 A. This is done by
        ! band-pass filtering the background image, recommended (and default settings) are hp=5.0
        ! lp=1.1 and width=5.0. <tseries_backgr_subtr/end>
        !
        ! set required keys
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
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:6))
        ! set defaults
        if( .not. cline%defined('hp')    ) call cline%set('hp',    5.0)
        if( .not. cline%defined('lp')    ) call cline%set('lp',    1.1)
        if( .not. cline%defined('width') ) call cline%set('width', 5.0)
        ! execute
        call xtseries_backgr_subtr%execute(cline)
    case( 'tseries_split' )
        !==Program tseries_split
        !
        ! <tseries_split/begin>is a program for splitting a time-series stack and its associated orientations
        ! <tseries_split/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'oritab'
        keys_required(3) = 'smpd'
        keys_required(4) = 'chunksz'
        keys_required(5) = 'stepsz'
        ! parse command line
        call cline%parse_oldschool(keys_required(:5))
        ! execute
        call xtseries_split%execute(cline)

    ! PARALLEL PROCESSING PROGRAMS

    case( 'merge_algndocs' )
        !==Program merge_algndocs
        !
        ! <merge_algndocs/begin>is a program for merging alignment documents from SIMPLE
        ! runs in distributed mode<merge_algndocs/end>
        !
        ! set required keys
        keys_required(1) = 'fbody'
        keys_required(2) = 'nptcls'
        keys_required(3) = 'ndocs'
        ! set optional keys
        keys_optional(1) = 'projfile'
        keys_optional(2) = 'numlen'
        keys_optional(3) = 'oritype' ! needs to be required when we move to *.simple format
        ! parse command line
        call cline%parse_oldschool(keys_required(:3), keys_optional(:3))
        ! execute
        call xmerge_algndocs%execute(cline)
    case( 'merge_nnmat' )
        !==Program merge_nnmat
        !
        ! <merge_nnmat/begin>is a program for merging partial nearest neighbour matrices calculated
        ! in distributed mode<merge_nnmat/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        keys_required(3) = 'nnn'
        ! parse command line
        call cline%parse_oldschool( keys_required(:3) )
        ! execute
        call xmerge_nnmat%execute(cline)
    case( 'merge_similarities' )
        !==Program merge_similarities
        !
        ! <merge_similarities/begin>is a program for merging similarities calculated between pairs of objects
        ! into a similarity matrix that can be inputted to cluster_smat<merge_similarities/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        ! set optional keys
        keys_optional(1) = 'nparts'
        ! parse command line
        call cline%parse_oldschool(keys_required(:1), keys_optional(:1))
        ! execute
        call xmerge_similarities%execute(cline)
    case( 'split_pairs' )
        !==Program split_pairs
        !
        ! <split_pairs/begin>is a program for splitting calculations between pairs of objects
        ! into balanced partitions<split_pairs/end>
        !
        ! set required keys
        keys_required(1) = 'nptcls'
        keys_required(2) = 'nparts'
        ! parse command line
        call cline%parse_oldschool(keys_required(:2))
        ! execute
        call xsplit_pairs%execute(cline)
    ! case( 'split' )
    !     !==Program split
    !     !
    !     ! <split/begin>is a program for splitting of image stacks into partitions for parallel execution.
    !     ! This is done to reduce I/O latency<split/end>
    !     !
    !     ! set required keys
    !     keys_required(1) = 'smpd'
    !     keys_required(2) = 'stk'
    !     ! set optional keys
    !     keys_optional(1) = 'nparts'
    !     keys_optional(2) = 'neg'
    !     ! parse command line
    !     call cline%parse_oldschool(keys_required(:2), keys_optional(:2))
    !
    !     ! execute
    !     call xsplit%execute(cline)
    ! case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_private_exec
