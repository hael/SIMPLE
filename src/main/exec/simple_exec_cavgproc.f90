!@descr: execution of cavgs processing commanders
module simple_exec_cavgproc
use simple_cmdline,           only: cmdline
use simple_commanders_cavgs,  only: commander_cluster_cavgs, commander_cluster_cavgs_selection,&
commander_select_clusters, commander_match_cavgs
use simple_commanders_stkops, only: commander_cluster_stack, commander_match_stacks
implicit none

public :: exec_cavgproc_commander
private

type(commander_cluster_cavgs)           :: xcluster_cavgs
type(commander_cluster_cavgs_selection) :: xcluster_cavgs_selection
type(commander_cluster_stack)           :: xcluster_stack
type(commander_match_cavgs)             :: xmatch_cavgs
type(commander_match_stacks)            :: xmatch_stacks
type(commander_select_clusters)         :: xsel_clusts

contains

    subroutine exec_cavgproc_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'cluster_cavgs' )
                call xcluster_cavgs%execute(cline)
            case( 'cluster_cavgs_selection' )
                call xcluster_cavgs_selection%execute(cline)
            case( 'cluster_stack' )
                call xcluster_stack%execute(cline)
            case( 'match_cavgs' )
                call xmatch_cavgs%execute(cline)
            case( 'match_stacks' )
                call xmatch_stacks%execute(cline)
            case( 'select_clusters' )
                 call xsel_clusts%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_cavgproc_commander

end module simple_exec_cavgproc
