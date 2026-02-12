!@descr: "sim" UI api (concrete implementation)
module simple_ui_api_sim
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: pdb2mrc
type(ui_program), target :: simulate_movie
type(ui_program), target :: simulate_noise
type(ui_program), target :: simulate_particles

contains

    subroutine register_ui_sim(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('pdb2mrc',            pdb2mrc,            prgtab)
        call add_ui_program('simulate_movie',     simulate_movie,     prgtab)
        call add_ui_program('simulate_noise',     simulate_noise,     prgtab)
        call add_ui_program('simulate_particles', simulate_particles, prgtab)
    end subroutine register_ui_sim

end module simple_ui_api_sim
