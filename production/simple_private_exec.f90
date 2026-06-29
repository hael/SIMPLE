!@descr: executes shared-memory parallelized programs executed by distributed commanders
program simple_private_exec
use simple_cmdline,             only: cmdline_err
use simple_defs,                only: STDLEN
use simple_private_exec_driver, only: run_private_exec_from_command_line, run_coarray_direct
implicit none
#include "simple_local_flags.inc"

character(len=STDLEN) :: xarg
integer               :: cmdstat, cmdlen

call get_command_argument(1, xarg, cmdlen, cmdstat)
if( cmdstat /= 0 ) call cmdline_err(cmdstat, cmdlen, xarg, 0)
if( trim(xarg) == '--coarray' )then
    call run_coarray_direct
else
    call run_private_exec_from_command_line
endif

end program simple_private_exec
