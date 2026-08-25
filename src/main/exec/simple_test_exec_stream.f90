!@descr: execution of stream test programs
module simple_test_exec_stream
use simple_cmdline,                only: cmdline
use simple_commanders_test_stream, only: commander_test_abinitio2D_stream, commander_test_assign_optics, &
    &commander_test_gen_pickrefs, commander_test_master, commander_test_pick_extract, commander_test_preproc, &
    &commander_test_sieve_cavgs
implicit none

public :: exec_test_stream_commander
private

type(commander_test_abinitio2D_stream) :: xabinitio2D_stream
type(commander_test_assign_optics)     :: xassign_optics
type(commander_test_gen_pickrefs)      :: xgen_pickrefs
type(commander_test_master)            :: xmaster
type(commander_test_pick_extract)      :: xpick_extract
type(commander_test_preproc)           :: xpreproc
type(commander_test_sieve_cavgs)       :: xsieve_cavgs

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
            case( 'master' )
                call xmaster%execute(cline)
            case( 'preproc' )
                call xpreproc%execute(cline)
            case( 'assign_optics' )
                call xassign_optics%execute(cline)
            case( 'gen_pickrefs' )
                call xgen_pickrefs%execute(cline)
            case( 'pick_extract' )
                call xpick_extract%execute(cline)
            case( 'sieve_cavgs' )
                call xsieve_cavgs%execute(cline)
            case( 'abinitio2D_stream' )
                call xabinitio2D_stream%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_stream_commander

end module simple_test_exec_stream
