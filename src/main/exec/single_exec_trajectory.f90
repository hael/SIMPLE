module single_exec_trajectory
use simple_cmdline,                 only: cmdline
use single_commanders_trajectory,   only: commander_import_trajectory, commander_graphene_subtr,&
&commander_trajectory_denoise, commander_extract_substk, commander_trajectory_swap_stack
use single_commanders_nano3D,       only: commander_trajectory_reconstruct3D_distr
use single_commanders_experimental, only: commander_trajectory_make_projavgs
implicit none

public :: exec_trajectory_commander
private

type(commander_extract_substk)                 :: xextract_substk
type(commander_graphene_subtr)                 :: xgraphene_subtr
type(commander_import_trajectory)              :: ximport_trajectory
type(commander_trajectory_denoise)             :: xden_traj
type(commander_trajectory_make_projavgs)       :: xtrajectory_make_projavgs
type(commander_trajectory_reconstruct3D_distr) :: xtrajectory_reconstruct3D
type(commander_trajectory_swap_stack)          :: xtrajectory_swap_stack

contains

    subroutine exec_trajectory_commander(which, cline, l_silent, l_did_execute)
        character(len=*), intent(in)    :: which
        class(cmdline),   intent(inout) :: cline
        logical,          intent(inout) :: l_did_execute
        logical,          intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'extract_substk' )
                call xextract_substk%execute(cline)
            case( 'graphene_subtr' )
                call cline%set('mkdir', 'no')
                call xgraphene_subtr%execute( cline )
            case( 'import_trajectory' )
                call ximport_trajectory%execute(cline)
            case( 'trajectory_denoise' )
                call cline%set('mkdir', 'no')
                call xden_traj%execute( cline )
            case( 'trajectory_make_projavgs' )
                call xtrajectory_make_projavgs%execute(cline)
            case( 'trajectory_reconstruct3D' )
                call xtrajectory_reconstruct3D%execute(cline)
            case( 'trajectory_swap_stack' )
                call xtrajectory_swap_stack%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_trajectory_commander

end module single_exec_trajectory