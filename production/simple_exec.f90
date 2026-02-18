!@descr: executes SIMPLE workflows
program simple_exec
use simple_exec_api
implicit none
#include "simple_local_flags.inc"
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
call make_ui
if( str_has_substr(entire_line, 'prg=list') )then
    call list_simple_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('simple_exec'))
l_silent      = .false.
l_did_execute = .false. ! will be set to true if one program was executed
call exec_project_commander(   trim(prg), cline, l_silent, l_did_execute)
call exec_preproc_commander(   trim(prg), cline, l_silent, l_did_execute)
call exec_cluster2D_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_cavgproc_commander(  trim(prg), cline, l_silent, l_did_execute)
call exec_abinitio3D_commander(trim(prg), cline, l_silent, l_did_execute)
call exec_refine3D_commander(  trim(prg), cline, l_silent, l_did_execute)
call exec_denoise_commander(   trim(prg), cline, l_silent, l_did_execute)
call exec_filter_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_image_commander(     trim(prg), cline, l_silent, l_did_execute)
call exec_mask_commander(      trim(prg), cline, l_silent, l_did_execute)
call exec_ori_commander(       trim(prg), cline, l_silent, l_did_execute)
call exec_print_commander(     trim(prg), cline, l_silent, l_did_execute)
call exec_sim_commander(       trim(prg), cline, l_silent, l_did_execute)
call exec_volume_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_res_commander(       trim(prg), cline, l_silent, l_did_execute)
call exec_dock_commander(      trim(prg), cline, l_silent, l_did_execute)
call exec_validate_commander(  trim(prg), cline, l_silent, l_did_execute)
call exec_other_commander(     trim(prg), cline, l_silent, l_did_execute)
if( .not. l_did_execute )then
    THROW_HARD('Program "'//trim(prg)//'" not recognized. Use prg=list to see available programs.')
endif
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('68efd3e9')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_exec
