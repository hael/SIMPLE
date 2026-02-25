!@descr: aggregates SIMPLE stream ui program constructors
module simple_ui_stream_group
use simple_ui_hash,   only: ui_hash
use simple_ui_stream, only: construct_stream_programs, print_stream_programs
implicit none

public :: add_stream_programs, print_stream_programs_group
private

contains

    subroutine add_stream_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call construct_stream_programs(prgtab)
    end subroutine add_stream_programs

    subroutine print_stream_programs_group( logfhandle )
        integer, intent(in) :: logfhandle
        call print_stream_programs(logfhandle)
    end subroutine print_stream_programs_group
    
end module simple_ui_stream_group
