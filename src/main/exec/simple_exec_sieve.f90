!@descr: execution of project management commanders
module simple_exec_sieve
use simple_cmdline,                 only: cmdline
use simple_commanders_sieve,        only: commander_sieve_ptcls
implicit none

public :: exec_sieve_commander
private

type(commander_sieve_ptcls) :: xsieve_ptcls

contains

    subroutine exec_sieve_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(adjustl(which)))
            case( 'particle_sieving' )
                call xsieve_ptcls%execute(cline)
        end select
    end subroutine exec_sieve_commander

end module simple_exec_sieve
