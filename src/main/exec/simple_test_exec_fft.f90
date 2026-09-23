!@descr: execution of test fft processing commanders
module simple_test_exec_fft
use simple_cmdline,             only: cmdline
use simple_commanders_test_fft, only: commander_test_gencorrs_fft
implicit none

public :: exec_test_fft_commander
private

type(commander_test_gencorrs_fft) :: xgencorrs_fft

contains

    subroutine exec_test_fft_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'gencorrs_fft' )
                call xgencorrs_fft%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_fft_commander

end module simple_test_exec_fft
