! executes the parallel (or distributed workflows) of SIMPLE
program simple_distr_exec
use simple_defs
use simple_gen_doc
use simple_user_interface, only: make_user_interface
use simple_cmdline,        only: cmdline, cmdline_err
use simple_strings,        only: str_has_substr
use simple_fileio,         only: extract_abspath
use simple_commander_stream_wflows
use simple_commander_distr_wflows
use simple_commander_hlev_wflows
implicit none

! PRE-PROCESSING
type(preprocess_distr_commander)                  :: xpreprocess
type(preprocess_stream_commander)                 :: xpreprocess_stream
type(motion_correct_distr_commander)              :: xmotion_correct_distr
type(motion_correct_tomo_distr_commander)         :: xmotion_correct_tomo_distr
type(powerspecs_distr_commander)                  :: xpowerspecs_distr
type(ctf_estimate_distr_commander)                :: xctf_estimate_distr
type(pick_distr_commander)                        :: xpick_distr

! CLUSTER2D
type(make_cavgs_distr_commander)                  :: xmake_cavgs_distr
type(cluster2D_autoscale_commander)               :: xcluster2D_distr
type(cluster2D_stream_distr_commander)            :: xcluster2D_stream_distr

! PRIME3D
type(refine3D_init_distr_commander)               :: xrefine3D_init_distr
type(prime3D_distr_commander)                     :: xprime3D_distr
type(reconstruct3D_distr_commander)               :: xreconstruct3D_distr
type(symsrch_distr_commander)                     :: xsymsrch_distr

! TIME-SERIES WORKFLOWS
type(tseries_track_distr_commander)               :: xtseries_track_distr

! HIGH-LEVEL WORKFLOWS
type(initial_3Dmodel_commander)                   :: xinitial_3Dmodel
type(cluster3D_commander)                         :: xcluster3D
type(cluster3D_refine_commander)                  :: xcluster3D_refine

! SUPORTING DISTRIBUTED WORKFLOWS
type(scale_project_distr_commander)               :: xscale_project

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

call make_user_interface

select case(prg)

    ! PRE-PROCESSING STREAM, LINKING MOTION_CORRECT + CTF_ESTIMATE + PICK

    case( 'preprocess' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            8.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( .not. cline%defined('opt')             ) call cline%set('opt',        'simplex')
        call xpreprocess%execute(cline)
    case( 'preprocess_stream' )
        call cline%parse()
        ! set defaults
        call cline%set('stream','yes')
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            8.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',         'yes')
        if( .not. cline%defined('opt')             ) call cline%set('opt',        'simplex')
        call xpreprocess_stream%execute(cline)
    ! MOTION_CORRECT_MOVIES

    case( 'motion_correct' )
        call cline%parse()
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
        keys_optional(12) = 'dir'
        ! parse command line
        if( describe ) call print_doc_motion_correct_tomo
        call cline%parse_oldschool(keys_required(:3), keys_optional(:12))
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
        call cline%parse_oldschool(keys_required(:4), keys_optional(:4))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        ! execute
        call xpowerspecs_distr%execute(cline)

    ! CTF_ESTIMATE

    case( 'ctf_estimate' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz',   512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',         30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',          5.)
        ! execute
        call xctf_estimate_distr%execute(cline)

    ! PARTICLE PICKER

    case( 'pick' )
        call cline%parse()
        ! execute
        call xpick_distr%execute(cline)

    ! CLUSTER2D

    case( 'make_cavgs' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        ! execute
        call xmake_cavgs_distr%execute(cline)
    case( 'cluster2D' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'no' )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     50. )
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no' )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        ! execute
        call xcluster2D_distr%execute(cline)
    case( 'cluster2D_stream' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',         'no')
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D',  'no')
        call xcluster2D_stream_distr%execute(cline)

    ! PRIME3D

    case( 'refine3D_init' )
        call cline%parse()
        call xrefine3D_init_distr%execute( cline )
    case( 'refine3D' )
        call cline%parse()
        ! set defaults
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('refine') ) call cline%set('refine',  'single')
        if( .not. cline%defined('eo')     ) call cline%set('eo', 'no')
        ! execute
        call xprime3D_distr%execute(cline)
    case( 'reconstruct3D' )
        call cline%parse()
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
        call cline%parse_oldschool(keys_required(:8), keys_optional(:5))
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
        call cline%parse_oldschool(keys_required(:5), keys_optional(:3))
        ! set defaults
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')   ) call cline%set('neg', 'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',    2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp', 5.0)
        ! execute
        call xtseries_track_distr%execute( cline )

    ! HIGH-LEVEL DISTRIBUTED WORKFLOWS

    case( 'initial_3Dmodel' )
        call cline%parse()
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
        keys_optional(9)  = 'stk'
        keys_optional(10) = 'stktab'
        keys_optional(11) = 'phaseplate'
        keys_optional(12) = 'oritab2'
        keys_optional(13) = 'mskfile'
        keys_optional(14) = 'startit'
        keys_optional(15) = 'objfun'
        keys_optional(16) = 'refine'
        ! parse command line
        ! if( describe ) call print_doc_cluster3D
        call cline%parse_oldschool(keys_required(:7), keys_optional(:16))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('refine') ) call cline%set('refine','cluster')
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
        keys_optional(9)  = 'stk'
        keys_optional(10) = 'stktab'
        keys_optional(11) = 'phaseplate'
        keys_optional(12) = 'msklist'
        keys_optional(13) = 'vollist'
        keys_optional(14) = 'state'
        keys_optional(15) = 'trs'
        keys_optional(16) = 'objfun'
        keys_optional(17) = 'update_frac'
        ! parse command line
        ! if( describe ) call print_doc_cluster3D_refine
        call cline%parse_oldschool(keys_required(:6), keys_optional(:17))
        ! sanity check
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            ! all ok
        else
            stop 'stk or stktab need to be part of command line!'
        endif
        ! set defaults
        if( .not. cline%defined('eo') )call cline%set('eo', 'no')
        ! execute
        call xcluster3D_refine%execute( cline )

    ! SUPPORTING DISTRIBUTED WORKFLOWS

    case( 'scale_project' )
        call cline%parse()
        call xscale_project%execute(cline )
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_distr_exec
