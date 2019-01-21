! executes the parallel (or distributed workflows) of SIMPLE
program simple_distr_exec
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface, list_distr_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_commander_base, only: execute_commander
use simple_spproj_hlev
use simple_commander_distr_wflows
use simple_commander_stream_wflows
use simple_commander_hlev_wflows
implicit none
#include "simple_local_flags.inc"

! PRE-PROCESSING WORKFLOWS
type(preprocess_distr_commander)            :: xpreprocess
type(preprocess_stream_commander)           :: xpreprocess_stream
type(motion_correct_distr_commander)        :: xmotion_correct_distr
type(gen_pspecs_and_thumbs_distr_commander) :: xgen_pspecs_and_thumbs
type(motion_correct_tomo_distr_commander)   :: xmotion_correct_tomo_distr
type(ctf_estimate_distr_commander)          :: xctf_estimate_distr
type(pick_distr_commander)                  :: xpick_distr
type(pick_extract_stream_distr_commander)   :: xpick_extract_stream_distr

! CLUSTER2D WORKFLOWS
type(make_cavgs_distr_commander)            :: xmake_cavgs_distr
type(cluster2D_autoscale_commander)         :: xcluster2D_distr
type(cluster2D_stream_distr_commander)      :: xcluster2D_stream_distr
type(cleanup2D_commander)                   :: xcleanup2D_distr

! AB INITIO 3D RECONSTRUCTION WORKFLOW
type(initial_3Dmodel_commander)             :: xinitial_3Dmodel

! REFINE3D WORKFLOWS
type(refine3D_init_distr_commander)         :: xrefine3D_init_distr
type(refine3D_distr_commander)              :: xprime3D_distr
type(reconstruct3D_distr_commander)         :: xreconstruct3D_distr

! CLUSTER3D WORKFLOWS
type(cluster3D_commander)                   :: xcluster3D
type(cluster3D_refine_commander)            :: xcluster3D_refine

! TIME-SERIES WORKFLOWS
type(tseries_track_distr_commander)         :: xtseries_track_distr

! MISCELLANEOUS WORKFLOWS
type(scale_project_distr_commander)         :: xscale_project

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
if( str_has_substr(entire_line, 'prg=list') )then
    call list_distr_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! set global defaults
if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')

select case(prg)

    ! PRE-PROCESSING WORKFLOWS

    case( 'preprocess' )
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          20.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            6.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        call xpreprocess%execute(cline)
    case( 'preprocess_stream' )
        call cline%set('stream','yes')
        if( .not. cline%defined('trs')             ) call cline%set('trs',               5.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',          20.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            6.)
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')          ) call cline%set('stream',         'yes')
        call xpreprocess_stream%execute(cline)
    case( 'motion_correct' )
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   20.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     6.)
        call xmotion_correct_distr%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'motion_correct_tomo' )
        if( .not. cline%defined('trs')     ) call cline%set('trs',        5.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',   20.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',     6.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',    'yes')
        call xmotion_correct_tomo_distr%execute(cline)
    case( 'ctf_estimate' )
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        call xpick_distr%execute(cline)
    case( 'pick_extract_stream' )
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        call xpick_extract_stream_distr%execute(cline)

    ! CLUSTER2D WORKFLOWS

    case( 'make_cavgs' )
        call xmake_cavgs_distr%execute(cline)
    case( 'cleanup2D' )
        if( .not. cline%defined('lp')        ) call cline%set('lp',         15. )
        if( .not. cline%defined('ncls')      ) call cline%set('ncls',      200. )
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'no' )
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     30. )
        if( .not. cline%defined('center')    ) call cline%set('center',    'no' )
        call xcleanup2D_distr%execute(cline)
    case( 'cluster2D' )
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',    15. )
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',      8. )
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'yes')
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      30. )
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',     50. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        call execute_commander(xcluster2D_distr, cline)
    case( 'cluster2D_stream' )
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'yes')
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('lpthresh')  ) call cline%set('lpthresh',    30.)
        if( .not. cline%defined('ndev')      ) call cline%set('ndev',        1.5)
        call xcluster2D_stream_distr%execute(cline)

    ! AB INITIO 3D RECONSTRUCTION WORKFLOW

    case( 'initial_3Dmodel' )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('eo')        ) call cline%set('eo',        'yes')
        call execute_commander(xinitial_3Dmodel, cline)

    ! REFINE3D WORKFLOWS

    case( 'refine3D_init' )
        call xrefine3D_init_distr%execute( cline )
    case( 'refine3D' )
        if( .not. cline%defined('cenlp')  ) call cline%set('cenlp',       30.)
        if( .not. cline%defined('refine') ) call cline%set('refine', 'single')
        if( .not. cline%defined('eo')     ) call cline%set('eo',         'no')
        if( cline%get_carg('eo').eq.'no' .and. .not.cline%defined('lp') )then
            THROW_HARD('The resolution limit should be set with LP=XX or EO=YES')
        endif
        call execute_commander(xprime3D_distr, cline)
    case( 'reconstruct3D' )
        if( .not. cline%defined('trs')  ) call cline%set('trs',      5.) ! to assure that shifts are being used
        if( .not. cline%defined('eo')   ) call cline%set('eo',     'no')
        call xreconstruct3D_distr%execute( cline )

    ! CLUSTER3D WORKFLOWS

    case( 'cluster3D' )
        if( .not. cline%defined('refine') )  call cline%set('refine', 'cluster')
        if( .not. cline%defined('eo') .and. .not. cline%defined('lp') ) call cline%set('eo', 'yes')
        if( cline%defined('lp') )            call cline%set('eo','no')
        call xcluster3D%execute( cline )
    case( 'cluster3D_refine' )
        if( .not. cline%defined('eo') ) call cline%set('eo', 'no')
        call xcluster3D_refine%execute( cline )

    ! TIME-SERIES WORKFLOWS

    case( 'tseries_track' )
        call cline%set('nthr', 1.0)
        if( .not. cline%defined('neg')   ) call cline%set('neg',   'yes')
        if( .not. cline%defined('lp')    ) call cline%set('lp',      2.0)
        if( .not. cline%defined('cenlp') ) call cline%set('cenlp',   5.0)
        call xtseries_track_distr%execute( cline )

    ! SUPPORTING WORKFLOWS

    case( 'scale_project' )
        call xscale_project%execute(cline )
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
end program simple_distr_exec
