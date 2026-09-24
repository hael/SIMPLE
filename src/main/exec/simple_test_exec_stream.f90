!@descr: execution of the stream workflow test program (preproc)
module simple_test_exec_stream
use simple_cmdline,                only: cmdline
use simple_commanders_test_stream, only: commander_test_preproc
implicit none

public :: exec_test_stream_commander
private

type(commander_test_preproc) :: xpreproc

contains

    subroutine exec_test_stream_commander( which, cline, l_silent, l_did_execute )
        character(len=*), intent(in)    :: which
        class(cmdline),   intent(inout) :: cline
        logical,          intent(out)   :: l_silent
        logical,          intent(inout) :: l_did_execute

        if( l_did_execute ) return
        l_silent      = .false.
        l_did_execute = .true.
        select case( trim(which) )
            case( 'preproc' )
                call xpreproc%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_stream_commander

end module simple_test_exec_stream
