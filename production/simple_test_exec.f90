!@descr: executes SIMPLE tests workflows
program simple_test_exec
use simple_test_exec_api
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN)             :: xarg, prg
character(len=XLONGSTRLEN)        :: entire_line
type(cmdline)                     :: cline
integer                           :: cmdstat, cmdlen, pos
integer(timer_int_kind)           :: t0
real(timer_int_kind)              :: rt_exec
logical                           :: l_silent, l_did_execute
! start timer
t0 = tic()
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err(cmdstat, cmdlen, xarg, pos)
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_test_ui
if( str_has_substr(entire_line, 'test=list') )then
    call list_simple_test_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('simple_exec'))
l_silent      = .false.
l_did_execute = .false. ! will be set to true if one program was executed
call exec_test_class_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_test_fft_commander(      trim(prg), cline, l_silent, l_did_execute)
call exec_test_geometry_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_test_highlevel_commander(trim(prg), cline, l_silent, l_did_execute)
call exec_test_io_commander(       trim(prg), cline, l_silent, l_did_execute)
call exec_test_masks_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_test_network_commander(  trim(prg), cline, l_silent, l_did_execute)
call exec_test_numerics_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_test_optimize_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_test_parallel_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_test_stats_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_test_utils_commander(    trim(prg), cline, l_silent, l_did_execute)
if( .not. l_did_execute )then
    THROW_HARD('Program test "'//trim(prg)//'" not recognized. Use prg=list to see available programs.')
endif
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('935ca398')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_test_exec
