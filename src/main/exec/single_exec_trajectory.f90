!@descr: execution of SINGLE trajectory commanders
module single_exec_trajectory
use simple_cmdline,                 only: cmdline
use single_commanders_trajectory,   only: commander_trajectory_import, commander_trajectory_graphene_subtr,&
&commander_trajectory_denoise, commander_trajectory_extract_substk, commander_trajectory_swap_stack
use single_commanders_nano3D,       only: commander_trajectory_reconstruct3D_distr
use single_commanders_trajectory_solve3D, only: commander_trajectory_solve3D
use single_commanders_experimental, only: commander_trajectory_make_projavgs
implicit none

public :: exec_trajectory_commander
private

type(commander_trajectory_extract_substk)                 :: xtrajectory_extract_substk
type(commander_trajectory_graphene_subtr)                 :: xtrajectory_graphene_subtr
type(commander_trajectory_import)              :: xtrajectory_import
type(commander_trajectory_denoise)             :: xden_traj
type(commander_trajectory_make_projavgs)       :: xtrajectory_make_projavgs
type(commander_trajectory_reconstruct3D_distr) :: xtrajectory_reconstruct3D
type(commander_trajectory_swap_stack)          :: xtrajectory_swap_stack
type(commander_trajectory_solve3D)             :: xtrajectory_solve3D

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
            case( 'trajectory_extract_substk' )
                call xtrajectory_extract_substk%execute(cline)
            case( 'trajectory_graphene_subtr' )
                call cline%set('mkdir', 'no')
                call xtrajectory_graphene_subtr%execute( cline )
            case( 'trajectory_import' )
                call xtrajectory_import%execute(cline)
            case( 'trajectory_denoise' )
                call cline%set('mkdir', 'no')
                call xden_traj%execute( cline )
            case( 'trajectory_make_projavgs' )
                call xtrajectory_make_projavgs%execute(cline)
            case( 'trajectory_reconstruct3D' )
                call xtrajectory_reconstruct3D%execute(cline)
            case( 'trajectory_swap_stack' )
                call xtrajectory_swap_stack%execute(cline)
            case( 'trajectory_solve3D' )
                call xtrajectory_solve3D%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_trajectory_commander

end module single_exec_trajectory