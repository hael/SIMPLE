!@descr: "single_ui_api_trajectory" UI api (concrete implementation)
module single_ui_api_trajectory
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: extract_substk
type(ui_program), target :: graphene_subtr
type(ui_program), target :: trajectory_denoise
type(ui_program), target :: trajectory_import_particles
type(ui_program), target :: trajectory_make_projavgs
type(ui_program), target :: trajectory_reconstruct3D
type(ui_program), target :: trajectory_swap_stack

contains

    subroutine register_single_ui_trajectory(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('extract_substk', extract_substk, prgtab)
        call add_ui_program('graphene_subtr', graphene_subtr, prgtab)
        call add_ui_program('trajectory_denoise', trajectory_denoise, prgtab)
        call add_ui_program('trajectory_import_particles', trajectory_import_particles, prgtab)
        call add_ui_program('trajectory_make_projavgs', trajectory_make_projavgs, prgtab)
        call add_ui_program('trajectory_reconstruct3D', trajectory_reconstruct3D, prgtab)
        call add_ui_program('trajectory_swap_stack', trajectory_swap_stack, prgtab)
    end subroutine register_single_ui_trajectory

end module single_ui_api_trajectory
