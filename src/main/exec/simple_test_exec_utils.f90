!@descr: execution of test utils processing commanders
module simple_test_exec_utils
use simple_cmdline,               only: cmdline
use simple_commanders_test_utils, only: commander_test_ansi_colors, commander_test_binoris_test, &
                                        commander_test_binoris_io_test, commander_test_cmdline, &
                                        commander_test_install, commander_test_nice, &
                                        commander_test_serialize, commander_test_stringmatch

implicit none

public :: exec_test_utils_commander
private

type(commander_test_ansi_colors)     :: xansi_colors
type(commander_test_binoris_test)    :: xbinoris_test
type(commander_test_binoris_io_test) :: xbinoris_io_test
type(commander_test_cmdline)         :: xcmdline
type(commander_test_install)         :: xinstall
type(commander_test_nice)            :: xnice
type(commander_test_serialize)       :: xserialize
type(commander_test_stringmatch)     :: xstringmatch

contains

    subroutine exec_test_utils_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'ansi_colors' )
                call xansi_colors%execute(cline)
            case( 'binoris_test' )
                call xbinoris_test%execute(cline)
            case( 'binoris_io_test' )
                call xbinoris_io_test%execute(cline)
            case( 'cmdline' )
                call xcmdline%execute(cline)
            case( 'install' )
                call xinstall%execute(cline)
            case( 'nice' )
                call xnice%execute(cline)
            case( 'serialize' )
                call xserialize%execute(cline)
            case( 'stringmatch' )
                call xstringmatch%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_utils_commander

end module simple_test_exec_utils
