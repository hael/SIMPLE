!@descr: execution of simulation commanders
module simple_exec_sim
use simple_cmdline,          only: cmdline
use simple_commanders_sim,   only: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie
use simple_commanders_atoms, only: commander_pdb2mrc
implicit none

public :: exec_sim_commander
private

type(commander_pdb2mrc)            :: xpdb2mrc
type(commander_simulate_movie)     :: xsimulate_movie
type(commander_simulate_noise)     :: xsimulate_noise
type(commander_simulate_particles) :: xsimulate_particles

contains

    subroutine exec_sim_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'pdb2mrc' )
                call xpdb2mrc%execute(cline)
            case( 'simulate_movie' )
                call xsimulate_movie%execute(cline)
            case( 'simulate_noise' )
                call xsimulate_noise%execute(cline)
            case( 'simulate_particles' )
                call xsimulate_particles%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_sim_commander


end module simple_exec_sim
