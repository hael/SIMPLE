! executes the parallel (or distributed workflows) of SIMPLE
program simple_distr_exec
use simple_defs
use simple_cmdline,      only: cmdline,cmdline_err
use simple_strings,      only: str_has_substr
use simple_filehandling, only: extract_abspath
use simple_gen_doc
use simple_commander_stream_wflows
use simple_commander_distr_wflows
use simple_commander_hlev_wflows
implicit none

! PRE-PROCESSING
type(preproc_stream_commander)           :: xpreproc_stream
type(unblur_ctffind_distr_commander)     :: xunblur_ctffind_distr
type(unblur_distr_commander)             :: xunblur_distr
type(unblur_tomo_movies_distr_commander) :: xunblur_tomo_distr
type(ctffind_distr_commander)            :: xctffind_distr
type(pick_distr_commander)               :: xpick_distr
! PRIME2D
type(makecavgs_distr_commander)          :: xmakecavgs_distr
type(prime2D_autoscale_commander)        :: xprime2D_distr
! 3D SIMILARITY MATRIX GENERATION WITH COMMON LINES
type(comlin_smat_distr_commander)        :: xcomlin_smat_distr
! PRIME3D
type(prime3D_init_distr_commander)       :: xprime3D_init_distr
type(prime3D_distr_commander)            :: xprime3D_distr
type(cont3D_distr_commander)             :: xcont3D_distr
type(recvol_distr_commander)             :: xrecvol_distr
type(symsrch_distr_commander)            :: xsymsrch_distr
! TIME-SERIES WORKFLOWS
type(tseries_track_distr_commander)      :: xtseries_track_distr
! HIGH-LEVEL WORKFLOWS
type(ini3D_from_cavgs_commander)         :: xini3D_from_cavgs
type(het_ensemble_commander)             :: xhet_ensemble

! OTHER DECLARATIONS
integer, parameter    :: MAXNKEYS=100, KEYLEN=32
character(len=KEYLEN) :: keys_required(MAXNKEYS)='', keys_optional(MAXNKEYS)=''
character(len=STDLEN) :: arg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
logical               :: describe
call get_command_argument(1, arg, cmdlen, cmdstat)
call get_command(entire_line)
if( str_has_substr(entire_line, 'prg=list') ) call list_all_simple_distr_programs
describe = str_has_substr(entire_line, 'describe=yes')
pos = index(arg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, arg, pos )
prg = arg(pos+1:) ! this is the program name
if( str_has_substr(prg, 'simple_') ) stop 'giving program names with simple_* prefix is depreciated'

select case(prg)

    ! PRE-PROCESSING STREAM, LINKING UNBLUR + CTFFIND + PICK

    case( 'preproc' )
        !==Program preproc
        !
        ! <preproc/begin>is a distributed workflow that executes unblur, ctffind and pick in sequence
        ! and in streaming mode as the microscope collects the data <preproc/end>
        !
        ! set required keys
        keys_required(1)  = 'smpd'
        keys_required(2)  = 'kv'
        keys_required(3)  = 'cs'
        keys_required(4)  = 'fraca'
        keys_required(5)  = 'dir_movies'
        keys_required(6)  = 'dir_target'
        keys_required(7)  = 'ncunits'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'refs'
        keys_optional(3)  = 'fbody'
        keys_optional(4)  = 'dose_rate'
        keys_optional(5)  = 'exp_time'
        keys_optional(6)  = 'lpstart'
        keys_optional(7)  = 'lpstop'
        keys_optional(8)  = 'trs'
        keys_optional(9)  = 'pspecsz_unblur'
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
        keys_optional(29) = 'fromm'
        ! parse command line
        if( describe ) call print_doc_preproc
        call cline%parse(keys_required(:7), keys_optional(:29))
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',                5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',             8.)
        if( .not. cline%defined('pspecsz_unblur')  ) call cline%set('pspecsz_unblur',   256.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz_ctffind', 1024.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp_ctffind',        30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp_ctffind',         5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',           20.)
        if( .not. cline%defined('stream')          ) call cline%set('stream',  'yes')
        call xpreproc_stream%execute(cline)

    ! PIPELINED UNBLUR + CTFFIND

    case( 'unblur_ctffind' )
        !==Program unblur_ctffind
        !
        ! <unblur_ctffind/begin>is a pipelined distributed workflow: unblur + ctffind program<unblur_ctffind/end> 
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
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'fbody'
        keys_optional(4)  = 'dose_rate'
        keys_optional(5)  = 'exp_time'
        keys_optional(6)  = 'lpstart'
        keys_optional(7)  = 'lpstop'
        keys_optional(8)  = 'trs'
        keys_optional(9)  = 'pspecsz'
        keys_optional(10) = 'numlen'
        keys_optional(11) = 'startit'
        keys_optional(12) = 'scale'
        keys_optional(13) = 'nframesgrp'
        keys_optional(14) = 'fromf'
        keys_optional(15) = 'tof'
        keys_optional(16) = 'nsig'
        keys_optional(17) = 'outfile'
        keys_optional(18) = 'hp'
        keys_optional(19) = 'lp'
        keys_optional(20) = 'dfmin'
        keys_optional(21) = 'dfmax'
        keys_optional(22) = 'dfstep'
        keys_optional(23) = 'astigtol'
        keys_optional(24) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_unblur_ctffind
        call cline%parse(keys_required(:6), keys_optional(:24))
        ! set defaults
        call cline%set('dopick', 'no'     )
        call cline%set('prg',    'preproc')
        if( .not. cline%defined('trs')             ) call cline%set('trs',                5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',             8.)
        if( .not. cline%defined('pspecsz_unblur')  ) call cline%set('pspecsz_unblur',   512.)
        if( .not. cline%defined('pspecsz_ctffind') ) call cline%set('pspecsz_ctffind', 1024.)
        if( .not. cline%defined('hp_ctffind')      ) call cline%set('hp_ctffind',        30.)
        if( .not. cline%defined('lp_ctffind')      ) call cline%set('lp_ctffind',         5.)
        ! execute
        call xunblur_ctffind_distr%execute(cline)

    ! UNBLUR_MOVIES

    case( 'unblur' )
        !==Program unblur
        !
        ! <unblur/begin>is a distributed workflow for movie alignment or unblurring based the same 
        ! principal strategy as Grigorieffs program (hence the name). There are two important 
        ! differences: automatic weighting of the frames using a correlation-based M-estimator and 
        ! continuous optimisation of the shift parameters. Input is a textfile with absolute paths 
        ! to movie files in addition to a few input parameters, some of which deserve a comment. If 
        ! dose_rate and exp_time are given the individual frames will be low-pass filtered accordingly
        ! (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to 
        ! the down-scaling factor (for super-resolution movies). If nframesgrp is given the frames will 
        ! be pre-averaged in the given chunk size (Falcon 3 movies). If fromf/tof are given, a 
        ! contiguous subset of frames will be averaged without any dose-weighting applied. 
        ! <unblur/end>
        !
        ! set required keys
        keys_required(1)  = 'filetab'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'fbody'
        keys_optional(4)  = 'dose_rate'
        keys_optional(5)  = 'exp_time'
        keys_optional(6)  = 'lpstart'
        keys_optional(7)  = 'lpstop'
        keys_optional(8)  = 'trs'
        keys_optional(9)  = 'kv'
        keys_optional(10) = 'pspecsz'
        keys_optional(11) = 'numlen'
        keys_optional(12) = 'startit'
        keys_optional(13) = 'scale'
        keys_optional(14) = 'nframesgrp'
        keys_optional(15) = 'fromf'
        keys_optional(16) = 'tof'
        keys_optional(17) = 'nsig'
        ! parse command line
        if( describe ) call print_doc_unblur
        call cline%parse(keys_required(:3), keys_optional(:17))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',      5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart', 15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',   8.)
        ! execute
        call xunblur_distr%execute(cline)
    case( 'unblur_tomo' )
        !==Program unblur_tomo
        !
        ! <unblur_tomo/begin>is a distributed workflow for movie alignment or unblurring of tomographic movies.
        ! Input is a textfile with absolute paths to movie files in addition to a few input parameters, some
        ! of which deserve a comment. The exp_doc document should contain per line exp_time=X and dose_rate=Y.
        ! It is asssumed that the input list of movies (one per tilt) are ordered temporally. This is necessary
        ! for correct dose-weighting of tomographic tilt series. If scale is given, the movie will be Fourier 
        ! cropped according to the down-scaling factor (for super-resolution movies). If nframesgrp is given 
        ! the frames will be pre-averaged in the given chunk size (Falcon 3 movies). <unblur_tomo/end>
        !
        ! set required keys
        keys_required(1)  = 'tomoseries'
        keys_required(2)  = 'exp_doc'
        keys_required(3)  = 'smpd'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'lpstart'
        keys_optional(4)  = 'lpstop'
        keys_optional(5)  = 'trs'
        keys_optional(6)  = 'kv'
        keys_optional(7)  = 'pspecsz'
        keys_optional(8)  = 'numlen'
        keys_optional(9)  = 'startit'
        keys_optional(10) = 'scale'
        keys_optional(11) = 'nframesgrp'
        keys_optional(12) = 'nsig'
        ! parse command line
        if( describe ) call print_doc_unblur_tomo
        call cline%parse(keys_required(:3), keys_optional(:12))
        ! set defaults
        if( .not. cline%defined('trs')     ) call cline%set('trs',       5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',  15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',    8.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',   'yes')
        ! execute
        call xunblur_tomo_distr%execute(cline)

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
        keys_optional(1) = 'ncunits'
        keys_optional(2) = 'pspecsz'
        keys_optional(3) = 'hp'
        keys_optional(4) = 'lp'
        keys_optional(5) = 'dfmin'
        keys_optional(6) = 'dfmax'
        keys_optional(7) = 'dfstep'
        keys_optional(8) = 'astigtol'
        keys_optional(9) = 'phaseplate'
        ! parse command line
        if( describe ) call print_doc_ctffind
        call cline%parse(keys_required(:6), keys_optional(:9))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 1024.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         5.)
        ! execute
        call xctffind_distr%execute(cline)

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
        keys_optional(4) = 'rm_outliers'
        ! parse command line
        if( describe ) call print_doc_pick
        call cline%parse(keys_required(:4), keys_optional(:4))
        ! execute
        call xpick_distr%execute(cline)

    ! PRIME2D

    case( 'makecavgs' )
        !==Program makecavgs
        !
        ! <makecavgs/begin>is a distributed workflowused for producing class averages or 
        ! initial random references for prime2D execution. <makecavgs/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'ctf'
        keys_required(4)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'ncls'
        keys_optional(4)  = 'deftab'
        keys_optional(5)  = 'oritab'
        keys_optional(6)  = 'filwidth'
        keys_optional(7)  = 'mul'
        keys_optional(8)  = 'outfile'
        keys_optional(9)  = 'refs'
        keys_optional(10) = 'remap_classes'
        keys_optional(11) = 'weights2D'
        keys_optional(12) = 'balance'
        ! parse command line
        if( describe ) call print_doc_makecavgs
        call cline%parse(keys_required(:4), keys_optional(:12))
        ! set defaults
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        ! execute
        call xmakecavgs_distr%execute(cline)
    case( 'prime2D' )
        !==Program prime2D
        !
        ! <prime2D/begin>is a distributed workflow implementing reference-free 2D alignment/clustering 
        ! algorithm adopted from the prime3D probabilistic ab initio 3D reconstruction algorithm<prime2D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ncls'
        keys_required(5)  = 'ctf'
        ! set optional keys
        keys_optional(1)  = 'nparts'
        keys_optional(2)  = 'chunksz'
        keys_optional(3)  = 'nthr'
        keys_optional(4)  = 'ncunits'
        keys_optional(5)  = 'deftab'
        keys_optional(6)  = 'refs'
        keys_optional(7)  = 'oritab'
        keys_optional(8)  = 'hp'
        keys_optional(9)  = 'lp'
        keys_optional(10) = 'lpstart'
        keys_optional(11) = 'lpstop'
        keys_optional(12) = 'cenlp'
        keys_optional(13) = 'trs'
        keys_optional(14) = 'automsk'
        keys_optional(15) = 'amsklp'
        keys_optional(16) = 'edge'
        keys_optional(17) = 'inner'
        keys_optional(18) = 'width'
        keys_optional(19) = 'startit'
        keys_optional(20) = 'maxits'
        keys_optional(21) = 'filwidth'
        keys_optional(22) = 'center'
        keys_optional(23) = 'autoscale'
        keys_optional(24) = 'oritab3D'
        keys_optional(25) = 'weights2D'
        keys_optional(26) = 'refine'
        keys_optional(27) = 'balance'
        ! documentation
        if( describe ) call print_doc_prime2D
        call cline%parse( keys_required(:5), keys_optional(:27) )
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15.)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',         'no')
        if( .not. cline%defined('amsklp')    ) call cline%set('amsklp',     20.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30.)
        if( .not. cline%defined('edge')      ) call cline%set('edge',       10.)
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30.)
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        if( cline%defined('nparts') .and. cline%defined('chunksz') )then
            stop 'nparts and chunksz cannot simultaneously be part of command line'
        else if(cline%defined('nparts') )then
            ! ok
        else if( cline%defined('chunksz') )then
            ! ok
        else
            stop 'eiter nparts or chunksz need to be part of command line'
        endif
        call xprime2D_distr%execute(cline)

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
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp' 
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'lp'
        keys_optional(5)  = 'inner'
        keys_optional(6)  = 'width'
        keys_optional(7)  = 'nspace'
        keys_optional(8)  = 'nran'
        keys_optional(9)  = 'npeaks'      
        ! parse command line
        if( describe ) call print_doc_prime3D_init
        call cline%parse(keys_required(:6), keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        ! execute
        call xprime3D_init_distr%execute( cline )
    case( 'prime3D' )
        !==Program prime3D
        !
        ! <prime3D/begin>is a distributed workflow for ab inito reconstruction/refinement based on
        ! probabilistic projection matching. PRIME is short for PRobabilistic Initial 3D Model
        ! generation for Single-particle cryo-Electron microscopy. There are a daunting number of
        ! options in PRIME3D. If you are processing class averages we recommend that you instead
        ! use the simple_distr_exec prg=ini3D_from_cavgs route for executing PRIME3D. Automated
        ! workflows for single- and multi-particle refinement using prime3D are planned for the
        ! next release (3.0)<prime3D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'deftab'
        keys_optional(4)  = 'vol1'
        keys_optional(5)  = 'oritab'
        keys_optional(6)  = 'trs'
        keys_optional(7)  = 'hp'
        keys_optional(8)  = 'lp'
        keys_optional(9)  = 'cenlp'
        keys_optional(10) = 'dynlp'
        keys_optional(11) = 'lpstart'
        keys_optional(12) = 'lpstop'
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
        keys_optional(27) = 'rrate'
        keys_optional(28) = 'norec'
        keys_optional(29) = 'nsub'
        keys_optional(30) = 'lp_grid'
        keys_optional(31) = 'balance'
        keys_optional(32) = 'center'
        ! documentation
        if( describe ) call print_doc_prime3D
        ! parse command line
        call cline%parse( keys_required(:6), keys_optional(:32) )
        ! set defaults
        if( .not. cline%defined('nspace')                  ) call cline%set('nspace', 1000.)
        if( cline%defined('lp') .or. cline%defined('find') ) call cline%set('dynlp',   'no')
        if( .not. cline%defined('cenlp')                   ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine')                  ) call cline%set('refine',  'no')
        if( .not. cline%defined('eo') )then
            call cline%set('eo', 'no')
        else
            if( cline%get_carg('eo').eq.'yes' )call cline%set('dynlp','no')
        endif
        ! execute
        call xprime3D_distr%execute(cline)
    case( 'cont3D' )
        !==Program cont3D
        !
        ! <cont3D/begin><cont3D/end>
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'ctf'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'oritab'
        keys_required(7)  = 'trs'
        keys_required(8)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'deftab'
        keys_optional(3)  = 'vol1'
        keys_optional(4)  = 'hp'
        keys_optional(5)  = 'lp'
        keys_optional(6)  = 'lpstop'
        keys_optional(7)  = 'frac'
        keys_optional(8)  = 'mskfile'
        keys_optional(9)  = 'inner'
        keys_optional(10) = 'width'
        keys_optional(11) = 'startit'
        keys_optional(12) = 'maxits'
        keys_optional(13) = 'refine'
        keys_optional(14) = 'eo'
        keys_optional(15) = 'athres'
        ! documentation
        if( describe ) call print_doc_cont3D
        ! parse command line
        call cline%parse( keys_required(:8), keys_optional(:15) )
        ! set defaults
        if( cline%defined('eo') )then
            if( cline%get_carg('eo').eq.'yes')then
                ! alles klar
            else
                if( .not.cline%defined('lp'))stop 'Low-pass must be defined with EO=NO'
            endif
        else
            call cline%set('eo','no')
        endif
        call cline%set('dynlp', 'no')
        if(.not.cline%defined('nspace'))call cline%set('nspace', 1000.)
        if(.not.cline%defined('refine'))call cline%set('refine', 'yes')
        ! execute
        call xcont3D_distr%execute(cline)
    case( 'recvol' )
        !==Program recvol
        !
        ! <recvol/begin>is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,
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
        ! the unordered DNA/RNA core of spherical icosahedral viruses<recvol/end>
        !
        ! set required keys
        keys_required(1) = 'stk'
        keys_required(2) = 'smpd'
        keys_required(3) = 'oritab'
        keys_required(4) = 'msk'
        keys_required(5) = 'ctf'
        keys_required(6) = 'pgrp'
        keys_required(7) = 'nparts'
        ! set optional keys
        keys_optional(1) = 'nthr'
        keys_optional(2) = 'ncunits'
        keys_optional(3) = 'eo'
        keys_optional(4) = 'deftab'
        keys_optional(5) = 'frac'
        keys_optional(6) = 'mskfile'
        keys_optional(7) = 'mul'
        keys_optional(8) = 'state'
        keys_optional(9) = 'balance'
        ! parse command line
        if( describe ) call print_doc_recvol
        call cline%parse(keys_required(:7), keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('trs') ) call cline%set('trs',  5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')  ) call cline%set('eo', 'no')
        ! execute
        call xrecvol_distr%execute( cline )
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
        ! Input this document to recvol together with the images and the point-group symmetry to generate a
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
        ! parse command line
        if( describe ) call print_doc_symsrch
        call cline%parse(keys_required(:8), keys_optional(:4))
        ! set defaults
        if(cline%defined('compare'))stop 'Distributed execution of SYMSRCH does not support the COMPARE argument'
        ! set defaults
        if( .not. cline%defined('nspace') )then
            call cline%set('nptcls', 50.) ! 50 projections 4 symsrch
            call cline%set('nspace', 50.) ! 50 projections 4 symsrch
        else
            call cline%set('nptcls', cline%get_rarg('nspace'))
        endif
        if(.not.cline%defined('cenlp')) call cline%set('cenlp', 30.)
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
        keys_required(5) = 'ncunits'
        ! set optional keys
        keys_optional(1) = 'lp'
        keys_optional(2) = 'offset'
        keys_optional(3) = 'cenlp'
        ! parse command line
        if( describe ) call print_doc_tseries_track
        call cline%parse(keys_required(:5), keys_optional(:2))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        ! execute
        call xtseries_track_distr%execute( cline )

    ! HIGH-LEVEL DISTRIBUTED WORKFLOWS

    case( 'ini3D_from_cavgs' )
        !==Program ini3D_from_cavgs
        !
        ! <ini3D_from_cavgs/begin>is a distributed workflow for generating an initial 
        ! 3D model from class averages obtained with prime2D<ini3D_from_cavgs/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'msk'
        keys_required(4)  = 'pgrp'
        keys_required(5)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'ncunits'
        keys_optional(3)  = 'hp'
        keys_optional(4)  = 'lp'
        keys_optional(5)  = 'lpstop' 
        keys_optional(6)  = 'frac'
        keys_optional(7)  = 'automsk'
        keys_optional(8)  = 'amsklp'
        keys_optional(9)  = 'edge'
        keys_optional(10) = 'binwidth'
        keys_optional(11) = 'inner'
        keys_optional(12) = 'width'
        keys_optional(13) = 'nspace'
        keys_optional(14) = 'shbarrier'
        keys_optional(15) = 'autoscale'
        keys_optional(16) = 'pgrp_known'
        keys_optional(17) = 'center'
        ! parse command line
        if( describe ) call print_doc_ini3D_from_cavgs
        call cline%parse(keys_required(:5), keys_optional(:17))
        ! set defaults
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 20.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',   10.)
        ! execute
        call xini3D_from_cavgs%execute( cline )
    case( 'het_ensemble' )
        !==Program het_ensemble
        !
        ! <het_ensemble/begin>is a distributed workflow for heterogeneity analysis 
        ! based on ensemble learning<het_ensemble/end> 
        !
        ! set required keys
        keys_required(1)  = 'stk'
        keys_required(2)  = 'smpd'
        keys_required(3)  = 'oritab'
        keys_required(4)  = 'msk'
        keys_required(5)  = 'pgrp'
        keys_required(6)  = 'ctf'
        keys_required(7)  = 'nstates'
        keys_required(8)  = 'lp'
        keys_required(9)  = 'nparts'
        ! set optional keys
        keys_optional(1)  = 'nthr'
        keys_optional(2)  = 'frac'
        keys_optional(3)  = 'automsk'
        keys_optional(4)  = 'amsklp'
        keys_optional(5)  = 'edge'
        keys_optional(6)  = 'binwidth'
        keys_optional(7)  = 'inner'
        keys_optional(8)  = 'width'
        keys_optional(9)  = 'nspace'
        keys_optional(10) = 'balance'
        ! parse command line
        if( describe ) call print_doc_het_ensemble
        call cline%parse(keys_required(:9), keys_optional(:9))
        ! set defaults
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 20.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',   10.)
        ! execute
        call xhet_ensemble%execute( cline )
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select
end program simple_distr_exec
