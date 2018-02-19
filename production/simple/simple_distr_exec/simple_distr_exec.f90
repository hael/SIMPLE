! executes the parallel (or distributed workflows) of SIMPLE
program simple_distr_exec
use simple_defs
use simple_gen_doc
use simple_cmdline, only: cmdline,cmdline_err
use simple_strings, only: str_has_substr
use simple_fileio,  only: extract_abspath
use simple_commander_stream_wflows
use simple_commander_distr_wflows
use simple_commander_hlev_wflows
implicit none

! PRE-PROCESSING
type(preprocess_stream_commander)            :: xpreprocess_stream
type(motion_correct_ctffind_distr_commander) :: xmotion_correct_ctffind_distr
type(motion_correct_distr_commander)         :: xmotion_correct_distr
type(motion_correct_tomo_distr_commander)    :: xmotion_correct_tomo_distr
type(powerspecs_distr_commander)             :: xpowerspecs_distr
type(ctffind_distr_commander)                :: xctffind_distr
type(ctf_estimate_distr_commander)           :: xctf_estimate_distr
type(pick_distr_commander)                   :: xpick_distr

! PRIME2D
type(make_cavgs_distr_commander)             :: xmake_cavgs_distr
type(prime2D_autoscale_commander)            :: xprime2D_distr
type(prime2D_stream_distr_commander)         :: xprime2D_stream_distr

! 3D SIMILARITY MATRIX GENERATION WITH COMMON LINES
type(comlin_smat_distr_commander)            :: xcomlin_smat_distr

! PRIME3D
type(prime3D_init_distr_commander)           :: xprime3D_init_distr
type(prime3D_distr_commander)                :: xprime3D_distr
type(reconstruct3D_distr_commander)          :: xreconstruct3D_distr
type(symsrch_distr_commander)                :: xsymsrch_distr

! TIME-SERIES WORKFLOWS
type(tseries_track_distr_commander)          :: xtseries_track_distr

! HIGH-LEVEL WORKFLOWS
type(initial_3Dmodel_commander)              :: xinitial_3Dmodel
type(cluster3D_commander)                    :: xcluster3D
type(cluster3D_refine_commander)             :: xcluster3D_refine

! SUPORTING DISTRIBUTED WORKFLOWS
type(scale_stk_parts_commander)              :: xscale_stk_parts

! OTHER DECLARATIONS
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: arg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
logical               :: describe

! parse command line
call get_command_argument(1, arg, cmdlen, cmdstat)
call get_command(entire_line)
if( str_has_substr(entire_line, 'prg=list') ) call list_all_simple_distr_programs
describe = str_has_substr(entire_line, 'describe=yes')
pos = index(arg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, arg, pos )
prg = arg(pos+1:) ! this is the program name
if( str_has_substr(prg, 'simple_') ) stop 'giving program names with simple_* prefix is depreciated'

select case(prg)

    ! PRE-PROCESSING STREAM, LINKING MOTION_CORRECT + CTFFIND + PICK

    case( 'preprocess' )
        !==Program preprocess
        !
        ! <preprocess/begin>is a distributed workflow that executes motion_correct, ctffind and pick in sequence
        ! and in streaming mode as the microscope collects the data <preprocess/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'kv'
        keys_required(3)  = 'cs'
        keys_required(4)  = 'fraca'
        keys_required(5)  = 'dir_movies'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'refs'
        keys_optional(3)  = 'fbody'
        keys_optional(4)  = 'dose_rate'
        keys_optional(5)  = 'exp_time'
        keys_optional(6)  = 'lpstart'
        keys_optional(7)  = 'lpstop'
        keys_optional(8)  = 'trs'
        keys_optional(9)  = 'pspecsz_motion_correct'
        keys_optional(10) = 'pspecsz_ctffind'
        keys_optional(11) = 'numlen'
        keys_optional(12) = 'startit'
        keys_optional(13) = 'scale'
        keys_optional(14) = 'nframesgrp'
        keys_optional(15) = 'fromf'
        keys_optional(16) = 'tof'
        keys_optional(17) = 'hp_ctffind'
        keys_optional(18) = 'lp_ctffind'
        keys_optional(19) = 'lp_pick'
        keys_optional(20) = 'dfmin'
        keys_optional(21) = 'dfmax'
        keys_optional(22) = 'dfstep'
        keys_optional(23) = 'astigtol'
        keys_optional(24) = 'phaseplate'
        keys_optional(25) = 'thres'
        keys_optional(26) = 'rm_outliers'
        keys_optional(27) = 'nsig'
        keys_optional(28) = 'dopick'
        keys_optional(29) = 'ndev'
        keys_optional(30) = 'box_extract'
        keys_optional(31) = 'pcontrast'
        ! parse command line
        if( describe ) call print_doc_preprocess
        call cline%parse(keys_required(:6), keys_optional(:31))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            8.)
        if( .not. cline%defined('pspecsz_motion_correct')  ) call cline%set('pspecsz_motion_correct',  512.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz_ctffind', 512.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp_ctffind',       30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp_ctffind',        5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('stream')          ) call cline%set('stream',  'yes')
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('opt')             ) call cline%set('opt', 'simplex')
        call xpreprocess_stream%execute(cline)

    ! PIPELINED MOTION_CORRECT + CTFFIND

    case( 'motion_correct_ctffind' )
        !==Program motion_correct_ctffind
        !
        ! <motion_correct_ctffind/begin>is a pipelined distributed workflow: motion_correct + ctffind program<motion_correct_ctffind/end>
        !
        ! set required keys
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'kv'
        keys_required(4)  = 'cs'
        keys_required(5)  = 'fraca'
        keys_required(6)  = 'nparts'
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
        if( describe ) call print_doc_motion_correct_ctffind
        call cline%parse(keys_required(:6), keys_optional(:23))
        ! set defaults
        call cline%set('dopick', 'no'     )
        call cline%set('prg',    'preprocess')
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            8.)
        if( .not. cline%defined('pspecsz_motion_correct')  ) call cline%set('pspecsz_motion_correct',  512.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz_ctffind', 512.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp_ctffind',       30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp_ctffind',        5.)
        if( .not. cline%defined('opt')             ) call cline%set('opt', 'simplex')
        ! execute
        call xmotion_correct_ctffind_distr%execute(cline)

    ! MOTION_CORRECT_MOVIES

    case( 'motion_correct' )
        !==Program motion_correct
        !
        ! <motion_correct/begin>is a distributed workflow for movie alignment or motion_correctring based the same
        ! principal strategy as Grigorieffs program (hence the name). There are two important
        ! differences: automatic weighting of the frames using a correlation-based M-estimator and
        ! continuous optimisation of the shift parameters. Input is a textfile with absolute paths
        ! to movie files in addition to a few input parameters, some of which deserve a comment. If
        ! dose_rate and exp_time are given the individual frames will be low-pass filtered accordingly
        ! (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to
        ! the down-scaling factor (for super-resolution movies). If nframesgrp is given the frames will
        ! be pre-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given, a
        ! contiguous subset of frames will be averaged without any dose-weighting applied.
        ! <motion_correct/end>
        !
        ! set required keys
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'fbody'
        keys_optional(3)  = 'dose_rate'
        keys_optional(4)  = 'exp_time'
        keys_optional(5)  = 'lpstart'
        keys_optional(6)  = 'lpstop'
        keys_optional(7)  = 'trs'
        keys_optional(8)  = 'kv'
        keys_optional(9)  = 'pspecsz'
        keys_optional(10) = 'numlen'
        keys_optional(11) = 'startit'
        keys_optional(12) = 'scale'
        keys_optional(13) = 'nframesgrp'
        keys_optional(14) = 'fromf'
        keys_optional(15) = 'tof'
        keys_optional(16) = 'nsig'
        ! parse command line
        if( describe ) call print_doc_motion_correct
        call cline%parse(keys_required(:3), keys_optional(:16))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('opt')     ) call cline%set('opt', 'simplex')
        ! execute
        call xmotion_correct_distr%execute(cline)
    case( 'motion_correct_tomo' )
        !==Program motion_correct_tomo
        !
        ! <motion_correct_tomo/begin>is a distributed workflow for movie alignment or motion_correctring of tomographic movies.
        ! Input is a textfile with absolute paths to movie files in addition to a few input parameters, some
        ! of which deserve a comment. The exp_doc document should contain per line exp_time=X and dose_rate=Y.
        ! It is asssumed that the input list of movies (one per tilt) are ordered temporally. This is necessary
        ! for correct dose-weighting of tomographic tilt series. If scale is given, the movie will be Fourier
        ! cropped according to the down-scaling factor (for super-resolution movies). If nframesgrp is given
        ! the frames will be pre-averaged in the given chunk size (Falcon 3 movies). <motion_correct_tomo/end>
        !
        ! set required keys
        keys_required(1)  = 'tomoseries'
        keys_required(2)  = 'exp_doc'
        keys_required(3)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'lpstart'
        keys_optional(3)  = 'lpstop'
        keys_optional(4)  = 'trs'
        keys_optional(5)  = 'kv'
        keys_optional(6)  = 'pspecsz'
        keys_optional(7)  = 'numlen'
        keys_optional(8)  = 'startit'
        keys_optional(9)  = 'scale'
        keys_optional(10) = 'nframesgrp'
        keys_optional(11) = 'nsig'
        ! parse command line
        if( describe ) call print_doc_motion_correct_tomo
        call cline%parse(keys_required(:3), keys_optional(:11))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',       5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',  15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',    8.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',   'yes')
        if( .not. cline%defined('opt')     ) call cline%set('opt', 'simplex')
        ! execute
        call xmotion_correct_tomo_distr%execute(cline)

    ! GENERATE POWERSPECTRA

    case( 'powerspecs' )
        !==Program powerspecs
        !
        ! <powerspecs/begin>is a program for generating powerspectra from a stack or filetable<powerspecs/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'fbody'
        keys_required(3) = 'filetab'
        keys_required(4) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'pspecsz'
        keys_optional(2) = 'speckind'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'clip'
        ! parse command line
        if( describe ) call print_doc_powerspecs
        call cline%parse(keys_required(:4), keys_optional(:4))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        ! execute
        call xpowerspecs_distr%execute(cline)

    ! CTFFIND

    case( 'ctffind' )
        !==Program ctffind
        !
        ! <ctffind/begin>is a distributed workflow that wraps CTFFIND4 (Grigorieff lab)<ctffind/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'kv'
        keys_required(4) = 'cs'
        keys_required(5) = 'fraca'
        keys_required(6) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'pspecsz'
        keys_optional(2) = 'hp'
        keys_optional(3) = 'lp'
        keys_optional(4) = 'dfmin'
        keys_optional(5) = 'dfmax'
        keys_optional(6) = 'dfstep'
        keys_optional(7) = 'astigtol'
        keys_optional(8) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_ctffind
        call cline%parse(keys_required(:6), keys_optional(:8))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        ! execute
        call xctffind_distr%execute(cline)
    case( 'ctf_estimate' )
        !==Program ctf_estimate
        !
        ! <ctf_estimate/begin>is a distributed SIMPLE workflow for fitting the CTF<ctf_estimate/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'smpd'
        keys_required(3) = 'kv'
        keys_required(4) = 'cs'
        keys_required(5) = 'fraca'
        keys_required(6) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'pspecsz'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'lp'
        keys_optional(5) = 'dfmin'
        keys_optional(6) = 'dfmax'
        keys_optional(7) = 'astigtol'
        keys_optional(8) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_ctf_estimate
        call cline%parse(keys_required(:6), keys_optional(:8))
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz',   512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',         30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',          5.)
        ! execute
        call xctf_estimate_distr%execute(cline)

    ! PARTICLE PICKER

    case( 'pick' )
        !==Program pick
        !
        ! <pick/begin>is a distributed workflow for template-based particle picking<pick/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'refs'
        keys_required(3) = 'smpd'
        keys_required(4) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'lp'
        keys_optional(3) = 'thres'
        keys_optional(4) = 'ndev'
        ! parse command line
        if( describe ) call print_doc_pick
        call cline%parse(keys_required(:4), keys_optional(:4))
        ! execute
        call xpick_distr%execute(cline)

    ! PRIME2D

    case( 'make_cavgs' )
        !==Program make_cavgs
        !
        ! <make_cavgs/begin>is a distributed workflow used for producing class averages or
        ! initial random references for cluster2D execution. <make_cavgs/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'ctf'
        keys_required(3)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncls'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'oritab'
        keys_optional(5)  = 'filwidth'
        keys_optional(6)  = 'mul'
        keys_optional(7)  = 'outfile'
        keys_optional(8)  = 'refs'
        keys_optional(9)  = 'remap_classes'
        keys_optional(10) = 'weights2D'
        keys_optional(11) = 'balance'
        keys_optional(12) = 'stk'
        keys_optional(13) = 'stktab'
        keys_optional(14) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_make_cavgs
        call cline%parse(keys_required(:3), keys_optional(:14))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        ! execute
        call xmake_cavgs_distr%execute(cline)
    case( 'cluster2D' )
        !==Program cluster2D
        !
        ! <cluster2D/begin>is a distributed workflow implementing a reference-free 2D alignment/clustering
        ! algorithm adopted from the prime3D probabilistic ab initio 3D reconstruction algorithm<cluster2D/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'msk'
        keys_required(3)  = 'ncls'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'refs'
        keys_optional(4)  = 'oritab'
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
        keys_optional(16) = 'autoscale'
        keys_optional(17) = 'weights2D'
        keys_optional(18) = 'balance'
        keys_optional(19) = 'match_filt'
        keys_optional(20) = 'stk'
        keys_optional(21) = 'stktab'
        keys_optional(22) = 'dyncls'
        keys_optional(23) = 'phaseplate'
        keys_optional(24) = 'opt'
        keys_optional(25) = 'objfun'
        ! documentation
        ! if( describe ) call print_doc_cluster2D
        call cline%parse( keys_required(:5), keys_optional(:25) )
        ! sanity checks
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'no' )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     50. )
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no' )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        ! execute
        call xprime2D_distr%execute(cline)
    case( 'cluster2D_stream' )
        !==Program cluster2D_stream
        !
        ! <cluster2D_stream/begin>is a distributed workflow implementing reference-free 2D alignment/clustering
        ! algorithm adopted from the prime3D probabilistic ab initio 3D reconstruction algorithm<cluster2D_stream/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'msk'
        keys_required(3)  = 'ctf'
        keys_required(4)  = 'nparts'
        keys_required(5)  = 'dir_ptcls'
        keys_required(6)  = 'ncls_start'
        keys_required(7)  = 'nptcls_per_cls'
        ! set optional keys
        keys_optional(1)  = 'match_filt'
        keys_optional(2)  = 'nthr'
        keys_optional(3)  = 'hp'
        keys_optional(4)  = 'lp'
        keys_optional(5)  = 'cenlp'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'inner'
        keys_optional(8)  = 'width'
        keys_optional(9)  = 'filwidth'
        keys_optional(10) = 'center'
        keys_optional(11) = 'phaseplate'
        keys_optional(12) = 'autoscale'
        keys_optional(13) = 'weights2D'
        keys_optional(14) = 'opt'
        keys_optional(15) = 'objfun'
        ! documentation
        ! if( describe ) call print_doc_cluster2D_stream
        call cline%parse( keys_required(:7), keys_optional(:15) )
        ! set defaults
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',       8.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',         'no')
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D',  'no')
        call xprime2D_stream_distr%execute(cline)

    ! 3D SIMILARITY MATRIX GENERATION WITH COMMON LINES

    case( 'comlin_smat' )
        !==Program comlin_smat
        !
        ! <comlin_smat/begin>is a distributed workflow for creating a similarity matrix based on
        ! common line correlation. The idea being that it should be possible to cluster images based
        ! on their 3D similarity witout having a 3D model by only operating on class averages
        ! and find averages that fit well together in 3D<comlin_smat/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'lp'
        keys_required(4) = 'msk'
        keys_required(5) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'hp'
        keys_optional(2) = 'trs'
        ! parse command line
        if( describe ) call print_doc_comlin_smat
        call cline%parse(keys_required(:5), keys_optional(:2))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('trs') ) call cline%set('trs', 3.0)
        ! execute
        call xcomlin_smat_distr%execute(cline)

    ! PRIME3D

    case('prime3D_init')
        !==Program prime3D_init
        !
        ! <prime3D_init/begin>is a distributed workflow for generating a random initial model for
        ! initialisation of PRIME3D. If the data set is large (>5000 images), generating a random
        ! model can be slow. To speedup, set nran to some smaller number, resulting in nran images
        ! selected randomly for reconstruction<prime3D_init/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'msk'
        keys_required(3)  = 'ctf'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'lp'
        keys_optional(4)  = 'inner'
        keys_optional(5)  = 'width'
        keys_optional(6)  = 'nspace'
        keys_optional(7)  = 'nran'
        keys_optional(8)  = 'npeaks'
        keys_optional(9)  = 'stk'
        keys_optional(10) = 'stktab'
        keys_optional(11) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_prime3D_init
        call cline%parse(keys_required(:5), keys_optional(:11))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xprime3D_init_distr%execute( cline )
    case( 'refine3D' )
        !==Program refine3D
        !
        ! <refine3D/begin>is a distributed workflow for ab inito reconstruction/refinement based on
        ! probabilistic projection matching. PRIME is short for PRobabilistic Initial 3D Model
        ! generation for Single-particle cryo-Electron microscopy. There are a daunting number of
        ! options in refine3D. If you are processing class averages we recommend that you instead
        ! use the simple_distr_exec prg=initial_3Dmodel route for executing refine3D. Automated
        ! workflows for single- and multi-particle refinement using refine3D are planned for the
        ! next release (3.0)<refine3D/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'msk'
        keys_required(3)  = 'ctf'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'vol1'
        keys_optional(4)  = 'oritab'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'hp'
        keys_optional(7)  = 'lp'
        keys_optional(8)  = 'cenlp'
        keys_optional(9) = 'dynlp'
        keys_optional(10) = 'lpstart'
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
        keys_optional(21) = 'npeaks'
        keys_optional(22) = 'startit'
        keys_optional(23) = 'maxits'
        keys_optional(24) = 'shbarrier'
        keys_optional(25) = 'noise'
        keys_optional(26) = 'nnn'
        keys_optional(27) = 'norec'
        keys_optional(28) = 'balance'
        keys_optional(29) = 'center'
        keys_optional(30) = 'pproc'
        keys_optional(31) = 'stk'
        keys_optional(32) = 'stktab'
        keys_optional(33) = 'weights3D'
        keys_optional(34) = 'phaseplate'
        keys_optional(35) = 'opt'
        keys_optional(36) = 'update_frac'
        keys_optional(37) = 'focusmsk'
        keys_optional(38) = 'objfun'
        ! documentation
        ! if( describe ) call print_doc_refine3D
        ! parse command line
        call cline%parse( keys_required(:5), keys_optional(:38) )
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( cline%defined('lp') .or. cline%defined('find') ) call cline%set('dynlp',   'no')
        if( .not. cline%defined('cenlp')                   ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine')                  ) call cline%set('refine',  'no')
        if( .not. cline%defined('pproc')                   ) call cline%set('pproc',  'yes')
        if( .not. cline%defined('eo') )then
            call cline%set('eo', 'no')
        else
            if( cline%get_carg('eo').ne.'no' )call cline%set('dynlp','no')
        endif
        ! execute
        call xprime3D_distr%execute(cline)
    case( 'reconstruct3D' )
        !==Program reconstruct3D
        !
        ! <reconstruct3D/begin>is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,
        ! given input orientations and state assignments. The algorithm is based on direct Fourier inversion
        ! with a Kaiser-Bessel (KB) interpolation kernel. This window function reduces the real-space ripple
        ! artifacts associated with direct moving windowed-sinc interpolation. The feature sought when
        ! implementing this algorithm was to enable quick, reliable reconstruction from aligned individual
        ! particle images. mul is used to scale the origin shifts if down-sampled
        ! were used for alignment and the original images are used for reconstruction. ctf=yes or ctf=flip
        ! turns on the Wiener restoration. If the images were phase-flipped set ctf=flip. amsklp, mw, and edge
        ! control the solvent mask: the low-pass limit used to generate the envelope; the molecular weight of
        ! the molecule (protein assumed but it works reasonably well also for RNA; slight modification of mw
        ! might be needed). The inner parameter controls the radius of the soft-edged mask used to remove
        ! the unordered DNA/RNA core of spherical icosahedral viruses<reconstruct3D/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'oritab'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'eo'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'frac'
        keys_optional(5)  = 'mskfile'
        keys_optional(6)  = 'mul'
        keys_optional(7)  = 'balance'
        keys_optional(8)  = 'stk'
        keys_optional(9)  = 'stktab'
        keys_optional(10) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_reconstruct3D
        call cline%parse(keys_required(:6), keys_optional(:10))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs',  5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        ! execute
        call xreconstruct3D_distr%execute( cline )
    case( 'symsrch' )
        !==Program symsrch
        !
        ! <symsrch/begin>is a distributed workflow for searching for the principal symmetry axis of a volume
        ! reconstructed without assuming any point-group symmetry. The program takes as input an
        ! asymmetrical 3D reconstruction. The alignment document for all the particle images
        ! that have gone into the 3D reconstruction and the desired point-group symmetry needs to
        ! be inputted. The 3D reconstruction is then projected in 50 (default option) even directions,
        ! common lines-based optimisation is used to identify the principal symmetry axis, the rotational
        ! transformation is applied to the inputted orientations, and a new alignment document is produced.
        ! Input this document to reconstruct3D together with the images and the point-group symmetry to generate a
        ! symmetrised map. If you are unsure about the point-group, you should use the compare=yes mode and
        ! input the highest conceviable point-group. The program then calculates probabilities for all lower
        ! groups inclusive.<symsrch/end>
        !
        ! set required keys
        keys_required(1) = 'vol1'
        keys_required(2) = 'smpd'
        keys_required(3) = 'msk'
        keys_required(4) = 'pgrp'
        keys_required(5) = 'oritab'
        keys_required(6) = 'outfile'
        keys_required(7) = 'lp'
        keys_required(8) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'cenlp'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'nspace'
        keys_optional(5) = 'center'
        ! parse command line
        if( describe ) call print_doc_symsrch
        call cline%parse(keys_required(:8), keys_optional(:5))
        ! sanity check
        if(cline%defined('compare'))stop 'Distributed execution of SYMSRCH does not support the COMPARE argument'
        ! set defaults
        if( .not. cline%defined('nspace') )then
            call cline%set('nptcls', 50.) ! 50 projections 4 symsrch
            call cline%set('nspace', 50.) ! 50 projections 4 symsrch
        else
            call cline%set('nptcls', cline%get_rarg('nspace'))
        endif
        if( .not.cline%defined('cenlp')   ) call cline%set('cenlp', 30.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        ! execute
        call xsymsrch_distr%execute( cline )

    ! TIME-SERIES DISTRIBUTED WORKFLOWS

    case( 'tseries_track' )
        !==Program tseries_track
        !
        ! <tseries_track/begin>is a distributed workflow for particle tracking
        ! in time-series data <tseries_track/end>
        !
        ! set required keys
        keys_required(1) = 'filetab'
        keys_required(2) = 'fbody'
        keys_required(3) = 'smpd'
        keys_required(4) = 'boxfile'
        keys_required(5) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'offset'
        keys_optional(3) = 'cenlp'
        ! parse command line
        if( describe ) call print_doc_tseries_track
        call cline%parse(keys_required(:5), keys_optional(:3))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        ! execute
        call xtseries_track_distr%execute( cline )

    ! HIGH-LEVEL DISTRIBUTED WORKFLOWS

    case( 'initial_3Dmodel' )
        !==Program initial_3Dmodel
        !
        ! <initial_3Dmodel/begin>is a distributed workflow for generating an initial
        ! 3D model from class averages obtained with cluster2D<initial_3Dmodel/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'hp'
        keys_optional(3)  = 'lpstart'
        keys_optional(4)  = 'lpstop'
        keys_optional(5)  = 'frac'
        keys_optional(6)  = 'inner'
        keys_optional(7)  = 'width'
        keys_optional(8)  = 'nspace'
        keys_optional(9)  = 'autoscale'
        keys_optional(10) = 'pgrp_known'
        keys_optional(11) = 'center'
        keys_optional(12) = 'opt'
        keys_optional(13) = 'update_frac'
        keys_optional(14) = 'objfun'
        ! parse command line
        if( describe ) call print_doc_initial_3Dmodel
        call cline%parse(keys_required(:5), keys_optional(:14))
        ! set defaults
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        ! execute
        call xinitial_3Dmodel%execute( cline )
    case( 'cluster3D' )
        !==Program cluster3D
        !
        ! <cluster3D/begin>is a distributed workflow for heterogeneity analysis by 3D clustering
        ! <cluster3D/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'oritab'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'nstates'
        keys_required(7)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'lp'
        keys_optional(3)  = 'lpstop'
        keys_optional(4)  = 'eo'
        keys_optional(5)  = 'frac'
        keys_optional(6)  = 'inner'
        keys_optional(7)  = 'width'
        keys_optional(8)  = 'nspace'
        keys_optional(9)  = 'balance'
        keys_optional(10) = 'stk'
        keys_optional(11) = 'stktab'
        keys_optional(12) = 'phaseplate'
        keys_optional(13) = 'opt'
        keys_optional(14) = 'oritab2'
        keys_optional(15) = 'mskfile'
        keys_optional(16) = 'startit'
        keys_optional(17) = 'objfun'
        ! parse command line
        ! if( describe ) call print_doc_cluster3D
        call cline%parse(keys_required(:7), keys_optional(:17))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('eo') .and. .not. cline%defined('lp') ) call cline%set('eo','yes')
        if( cline%defined('lp') ) call cline%set('eo','no')
        ! execute
        call xcluster3D%execute( cline )
    case( 'cluster3D_refine' )
        !==Program cluster3D_refine
        !
        ! <cluster3D_refine/begin>is a distributed workflow for refinement of heterogeneity analysis
        ! by cluster3D<cluster3D_refine/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'oritab'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'ctf'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'lp'
        keys_optional(3)  = 'lpstop'
        keys_optional(4)  = 'eo'
        keys_optional(5)  = 'frac'
        keys_optional(6)  = 'inner'
        keys_optional(7)  = 'width'
        keys_optional(8)  = 'nspace'
        keys_optional(9)  = 'balance'
        keys_optional(10) = 'stk'
        keys_optional(11) = 'stktab'
        keys_optional(12) = 'phaseplate'
        keys_optional(13) = 'opt'
        keys_optional(14) = 'msklist'
        keys_optional(15) = 'vollist'
        keys_optional(16) = 'state'
        keys_optional(17) = 'trs'
        keys_optional(18) = 'objfun'
        keys_optional(19) = 'update_frac'
        ! parse command line
        ! if( describe ) call print_doc_cluster3D_refine
        call cline%parse(keys_required(:6), keys_optional(:19))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('eo')       ) call cline%set('eo',       'no')
        ! execute
        call xcluster3D_refine%execute( cline )

    ! SUPPORTING DISTRIBUTED WORKFLOWS

    case( 'scale_stk_parts' )
        !==Program scale_stk_parts
        !
        ! <scale_stk_parts/begin>is a distributed workflow for scaling
        ! partial stacks<scale_stk_parts/end>
        !
        ! set required keys
        keys_required(1) = 'smpd'
        keys_required(2) = 'newbox'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'stktab'
        keys_optional(3) = 'nparts'
        ! parse command line
        if( describe ) call print_doc_scale_stk_parts
        call cline%parse(keys_required(:2), keys_optional(:3))
        ! sanity check
        if( cline%defined('nparts') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'nparts or stktab need to be part of command line!'
        endif
        ! execute
        call xscale_stk_parts%execute( cline )
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_distr_exec
