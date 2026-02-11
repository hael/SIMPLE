!@descr: executes SIMPLE tests workflows
program simple_test_exec
use simple_exec_api
implicit none
#include "simple_local_flags.inc"

type(commander_test_sim_workflow) :: xtest_sim_workflow 

! OTHER DECLARATIONS
character(len=STDLEN)      :: xarg, prg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent

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
if( str_has_substr(entire_line, 'test=list') )then
    call list_simple_test_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse_test
! generate script for queue submission?
call script_exec(cline, string(trim(prg)), string('simple_test_exec'))
l_silent = .false.
select case(trim(prg))

    case( 'sim_workflow' )
        call xtest_sim_workflow%execute(cline)

    case default
        THROW_HARD('test='//trim(prg)//' is unsupported')

end select

call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version('7b8f21ce')
    ! end timer and print
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_test_exec
