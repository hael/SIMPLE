!@descr: executes the SINGLE (Structure Identification of Nanoparticles with Liquid-cell Em) workflows
program single_exec
use single_exec_api
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN)      :: args, prg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent, l_did_execute
! start timer
t0 = tic()
! parse command line
call get_command_argument(1, args, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(args, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, args, pos )
prg = args(pos+1:) ! this is the program name
! make UI
call make_ui
if( str_has_substr(entire_line, 'prg=list') )then
    call list_single_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('single_exec'))
! set global defaults
if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
l_silent      = .false.
l_did_execute = .false. ! will be set to true if one program was executed
call exec_tseries_commander(    trim(prg), cline, l_silent, l_did_execute)
call exec_trajectory_commander( trim(prg), cline, l_silent, l_did_execute)
call exec_nano2D_commander(     trim(prg), cline, l_silent, l_did_execute)
call exec_nano3D_commander(     trim(prg), cline, l_silent, l_did_execute)
call exec_map_commander(        trim(prg), cline, l_silent, l_did_execute)
call exec_atom_commander(       trim(prg), cline, l_silent, l_did_execute)
call exec_validate_commander(   trim(prg), cline, l_silent, l_did_execute)
if( .not. l_did_execute )then
    THROW_HARD('Program "'//trim(prg)//'" not recognized. Use prg=list to see available programs.')
endif
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('d3fa2ff0')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program single_exec

