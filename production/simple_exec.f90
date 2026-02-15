!@descr: executes SIMPLE workflows
program simple_exec
use simple_exec_api
use simple_exec_project,    only: exec_project_commander
use simple_exec_preproc,    only: exec_preproc_commander
use simple_exec_cluster2D,  only: exec_cluster2D_commander
use simple_exec_cavgproc,   only: exec_cavgproc_commander
use simple_exec_abinitio3D, only: exec_abinitio3D_commander
use simple_exec_refine3D,   only: exec_refine3D_commander
use simple_exec_denoise,    only: exec_denoise_commander
use simple_exec_filter,     only: exec_filter_commander
use simple_exec_image,      only: exec_image_commander
use simple_exec_mask,       only: exec_mask_commander
use simple_exec_ori,        only: exec_ori_commander
use simple_exec_print,      only: exec_print_commander
use simple_exec_sim,        only: exec_sim_commander
use simple_exec_volume,     only: exec_volume_commander
use simple_exec_res,        only: exec_res_commander
use simple_exec_dock,       only: exec_dock_commander
implicit none
#include "simple_local_flags.inc"

! PARALLEL UTILITIES
type(commander_split)                       :: xsplit

! STREAM VALIDATION
type(commander_check_refpick)               :: xcheck_refpick
type(commander_mini_stream)                 :: xmini_stream

! MODEL ANALYSIS
type(commander_model_validation)            :: xmodel_validation

! OTHER DECLARATIONS
character(len=STDLEN)      :: xarg, prg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent, l_did_execute

! start timer
t0 = tic()
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_simple_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('simple_exec'))
l_silent      = .false.
l_did_execute = .false.

call exec_project_commander(trim(prg),    cline, l_silent, l_did_execute)
call exec_preproc_commander(trim(prg),    cline, l_silent, l_did_execute)
call exec_cluster2D_commander(trim(prg),  cline, l_silent, l_did_execute)
call exec_cavgproc_commander(trim(prg),   cline, l_silent, l_did_execute)
call exec_abinitio3D_commander(trim(prg), cline, l_silent, l_did_execute)
call exec_refine3D_commander(trim(prg),   cline, l_silent, l_did_execute)
call exec_denoise_commander(trim(prg),    cline, l_silent, l_did_execute)
call exec_filter_commander(trim(prg),     cline, l_silent, l_did_execute)
call exec_image_commander(trim(prg),      cline, l_silent, l_did_execute)
call exec_mask_commander(trim(prg),       cline, l_silent, l_did_execute)
call exec_ori_commander(trim(prg),        cline, l_silent, l_did_execute)
call exec_print_commander(trim(prg),      cline, l_silent, l_did_execute)
call exec_sim_commander(trim(prg),        cline, l_silent, l_did_execute)
call exec_volume_commander(trim(prg),     cline, l_silent, l_did_execute)
call exec_res_commander(trim(prg),        cline, l_silent, l_did_execute)
call exec_dock_commander(trim(prg),       cline, l_silent, l_did_execute)

select case(trim(prg))

    !====================================================================
    ! PARALLEL UTILITIES
    !====================================================================
    case( 'split' )
        call xsplit%execute(cline)

    

    !====================================================================
    ! STREAM VALIDATION
    !====================================================================
    case( 'check_refpick' )
        call xcheck_refpick%execute(cline)
    case( 'mini_stream' )
        call xmini_stream%execute(cline)

    !====================================================================
    ! MODEL ANALYSIS
    !====================================================================
    case( 'model_validation' )
        call xmodel_validation%execute(cline)

    case default
        THROW_HARD('prg='//trim(prg)//' is unsupported')

end select

call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('0b2bcd6b')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_exec
