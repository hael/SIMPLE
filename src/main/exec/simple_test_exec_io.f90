!@descr: execution of test input/output processing commanders
module simple_test_exec_io
use simple_cmdline,            only: cmdline
use simple_commanders_test_io, only: commander_test_imgfile, commander_test_inside_write, &
                                     commander_test_io, commander_test_io_parallel, &
                                     commander_test_mrc2jpeg, commander_test_mrc_validation, &
                                     commander_test_stack_io, commander_test_star_export, &
                                     commander_test_starfile_test
implicit none

public :: exec_test_io_commander
private

type(commander_test_imgfile)        :: ximgfile
type(commander_test_inside_write)   :: xinside_write
type(commander_test_io)             :: xio
type(commander_test_io_parallel)    :: xio_parallel
type(commander_test_mrc2jpeg)       :: xmrc2jpeg
type(commander_test_mrc_validation) :: xmrc_validation
type(commander_test_stack_io)       :: xstack_io
type(commander_test_star_export)    :: xstar_export
type(commander_test_starfile_test)  :: xstarfile_test

contains

    subroutine exec_test_io_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'imgfile' )
                call ximgfile%execute(cline)
            case( 'inside_write' )
                call xinside_write%execute(cline)
            case( 'io' )
                call xio%execute(cline)
            case( 'io_parallel' )
                call xio_parallel%execute(cline)
            case( 'mrc2jpeg' )
                call xmrc2jpeg%execute(cline)
            case( 'mrc_validation' )
                call xmrc_validation%execute(cline)
            case( 'stack_io' )
                call xstack_io%execute(cline)
            case( 'star_export' )
                call xstar_export%execute(cline)
            case( 'starfile_test' )
                call xstarfile_test%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_io_commander

end module simple_test_exec_io
