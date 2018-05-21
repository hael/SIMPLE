! executes the parallel (or distributed workflows) of SIMPLE
program simple_distr_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface, list_distr_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_commander_base, only: execute_commander
use simple_commander_distr_wflows
use simple_commander_stream_wflows
use simple_commander_hlev_wflows
implicit none

! PRE-PROCESSING
type(preprocess_distr_commander)          :: xpreprocess
type(preprocess_stream_commander)         :: xpreprocess_stream
type(motion_correct_distr_commander)      :: xmotion_correct_distr
type(motion_correct_tomo_distr_commander) :: xmotion_correct_tomo_distr
type(powerspecs_distr_commander)          :: xpowerspecs_distr
type(ctf_estimate_distr_commander)        :: xctf_estimate_distr
type(pick_distr_commander)                :: xpick_distr

! CLUSTER2D
type(make_cavgs_distr_commander)          :: xmake_cavgs_distr
type(cluster2D_autoscale_commander)       :: xcluster2D_distr
type(cluster2D_stream_distr_commander)    :: xcluster2D_stream_distr

! REFINE3D
type(refine3D_init_distr_commander)       :: xrefine3D_init_distr
type(refine3D_distr_commander)             :: xprime3D_distr
type(reconstruct3D_distr_commander)       :: xreconstruct3D_distr
type(symsrch_distr_commander)             :: xsymsrch_distr

! HIGH-LEVEL WORKFLOWS
type(initial_3Dmodel_commander)           :: xinitial_3Dmodel
type(cluster3D_commander)                 :: xcluster3D
type(cluster3D_refine_commander)          :: xcluster3D_refine

! TIME-SERIES WORKFLOWS
type(tseries_track_distr_commander)       :: xtseries_track_distr

! SUPORTING DISTRIBUTED WORKFLOWS
type(scale_project_distr_commander)       :: xscale_project

! OTHER DECLARATIONS
character(len=STDLEN) :: args, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos

! parse command line
call get_command_argument(1, args, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(args, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, args, pos )
prg = args(pos+1:) ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') ) call list_distr_prgs_in_ui

select case(prg)

    ! PRE-PROCESSING

    case( 'preprocess' )
        call cline%parse()
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          15.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            8.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        call xpreprocess%execute(cline)
    case( 'preprocess_stream' )
        call cline%parse()
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
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        call xpreprocess_stream%execute(cline)
    case( 'motion_correct' )
        call cline%parse()
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        call xmotion_correct_distr%execute(cline)
    case( 'motion_correct_tomo' )
        call cline%parse()
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   15.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     8.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',    'yes')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        call xmotion_correct_tomo_distr%execute(cline)
    case( 'powerspecs' )
        call cline%parse()
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('clip')    ) call cline%set('clip',    256.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        6.)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        call xpowerspecs_distr%execute(cline)
    case( 'ctf_estimate' )
        call cline%parse()
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz',   512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',         30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',          5.)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',    'yes')
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        call cline%parse()
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',    'yes')
        call xpick_distr%execute(cline)

    ! CLUSTER2D

    case( 'make_cavgs' )
        call cline%parse()
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',    'yes')
        call xmake_cavgs_distr%execute(cline)
    case( 'cluster2D' )
        call cline%parse()
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'no' )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     50. )
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D', 'no' )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        call xcluster2D_distr%execute(cline)
    case( 'cluster2D_stream' )
        call cline%parse()
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',         'no')
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('weights2D') ) call cline%set('weights2D',  'no')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        call xcluster2D_stream_distr%execute(cline)

    ! REFINE3D

    case( 'refine3D_init' )
        call cline%parse()
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xrefine3D_init_distr%execute( cline )
    case( 'refine3D' )
        call cline%parse()
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('refine') ) call cline%set('refine', 'single')
        if( .not. cline%defined('eo')     ) call cline%set('eo',         'no')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',     'yes')
        call xprime3D_distr%execute(cline)
    case( 'reconstruct3D' )
        call cline%parse()
        if( .not. cline%defined('trs')  ) call cline%set('trs',      5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')   ) call cline%set('eo',     'no')
        if( .not. cline%defined('mkdir')) call cline%set('mkdir', 'yes')
        call xreconstruct3D_distr%execute( cline )
    case( 'symsrch' )
        call cline%parse()
        if( .not. cline%defined('nspace') )then
            call cline%set('nspace', 50.) ! 50 projections 4 symsrch
        endif
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',    30.)
        if( .not. cline%defined('center') ) call cline%set('center', 'yes')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'yes')
        call xsymsrch_distr%execute( cline )

    ! HIGH-LEVEL DISTRIBUTED WORKFLOWS

    case( 'initial_3Dmodel' )
        call cline%parse()
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        call execute_commander(xinitial_3Dmodel, cline)
    case( 'cluster3D' )
        call cline%parse()
        if( .not. cline%defined('refine') )  call cline%set('refine', 'cluster')
        if( .not. cline%defined('eo') .and. .not. cline%defined('lp') ) call cline%set('eo', 'yes')
        if( cline%defined('lp') )            call cline%set('eo','no')
        if( .not. cline%defined('mkdir')  )  call cline%set('mkdir', 'yes')
        call xcluster3D%execute( cline )
    case( 'cluster3D_refine' )
        call cline%parse()
        if( .not. cline%defined('eo')    ) call cline%set('eo',     'no')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xcluster3D_refine%execute( cline )

    ! TIME-SERIES DISTRIBUTED WORKFLOWS

    case( 'tseries_track' )
        call cline%parse()
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')   ) call cline%set('neg',   'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',      2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',   5.0)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call xtseries_track_distr%execute( cline )

    ! SUPPORTING DISTRIBUTED WORKFLOWS

    case( 'scale_project' )
        call cline%parse()
        call xscale_project%execute(cline )
    case DEFAULT
        write(*,'(a,a)') 'program key (prg) is: ', trim(prg)
        stop 'unsupported program'
    end select

end program simple_distr_exec
