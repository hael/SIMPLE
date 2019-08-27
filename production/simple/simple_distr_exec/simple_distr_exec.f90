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
type(extract_distr_commander)               :: xextract_distr
type(reextract_distr_commander)             :: xreextract_distr
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
type(refine3D_distr_commander)              :: xrefine3D_distr
type(reconstruct3D_distr_commander)         :: xreconstruct3D_distr

! CLUSTER3D WORKFLOWS
type(cluster3D_commander)                   :: xcluster3D
type(cluster3D_refine_commander)            :: xcluster3D_refine

! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
type(tseries_track_distr_commander)         :: xtseries_track_distr
type(cleanup2D_nano_commander)              :: xcleanup2D_nano
type(refine3D_nano_distr_commander)         :: xrefine3D_nano_distr

! MISCELLANEOUS WORKFLOWS
type(scale_project_distr_commander)         :: xscale_project
type(prune_project_distr_commander)         :: xprune_project

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
        call xpreprocess%execute(cline)
    case( 'preprocess_stream' )
        call xpreprocess_stream%execute(cline)
    case( 'extract' )
        call xextract_distr%execute(cline)
    case( 'reextract' )
        call xreextract_distr%execute(cline)
    case( 'motion_correct' )
        call xmotion_correct_distr%execute(cline)
    case( 'gen_pspecs_and_thumbs' )
        call xgen_pspecs_and_thumbs%execute(cline)
    case( 'motion_correct_tomo' )
        call xmotion_correct_tomo_distr%execute(cline)
    case( 'ctf_estimate' )
        call xctf_estimate_distr%execute(cline)
    case( 'pick' )
        call xpick_distr%execute(cline)
    case( 'pick_extract_stream' )
        call xpick_extract_stream_distr%execute(cline)

    ! CLUSTER2D WORKFLOWS
    case( 'make_cavgs' )
        call xmake_cavgs_distr%execute(cline)
    case( 'cleanup2D' )
        call xcleanup2D_distr%execute(cline)
    case( 'cluster2D' )
        call execute_commander(xcluster2D_distr, cline)
    case( 'cluster2D_stream' )
        call xcluster2D_stream_distr%execute(cline)

    ! AB INITIO 3D RECONSTRUCTION WORKFLOW
    case( 'initial_3Dmodel' )
        call execute_commander(xinitial_3Dmodel, cline)

    ! REFINE3D WORKFLOWS
    case( 'refine3D' )
        call execute_commander(xrefine3D_distr, cline)
    case( 'reconstruct3D' )
        call xreconstruct3D_distr%execute( cline )

    ! CLUSTER3D WORKFLOWS
    case( 'cluster3D' )
        call xcluster3D%execute( cline )
    case( 'cluster3D_refine' )
        call xcluster3D_refine%execute( cline )

    ! TIME-SERIES (NANO-PARTICLE) WORKFLOWS
    case( 'tseries_track' )
        call xtseries_track_distr%execute( cline )
    case( 'cleanup2D_nano' )
        call xcleanup2D_nano%execute(cline)
    case( 'refine3D_nano')
        call execute_commander(xrefine3D_nano_distr, cline)

    ! SUPPORTING WORKFLOWS
    case( 'scale_project' )
        call xscale_project%execute( cline )
    case( 'prune_project' )
        call xprune_project%execute( cline )
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
end program simple_distr_exec
