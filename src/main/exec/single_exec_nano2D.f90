module single_exec_nano2D
use simple_cmdline,            only: cmdline
use single_commanders_nano2D,  only: commander_analysis2D_nano, commander_center2D_nano, commander_cluster2D_nano
use simple_commanders_imgproc, only: commander_estimate_diam
implicit none

public :: exec_nano2D_commander
private

type(commander_analysis2D_nano) :: xanalysis2D_nano
type(commander_center2D_nano)   :: xcenter2D
type(commander_cluster2D_nano)  :: xcluster2D
type(commander_estimate_diam)   :: xestimate_diam

contains

    subroutine exec_nano2D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'analysis2D_nano' )
                call xanalysis2D_nano%execute(cline)
            case( 'center2D_nano' )
                call xcenter2D%execute(cline)
            case( 'cluster2D_nano' )
                call xcluster2D%execute(cline)
            case( 'estimate_diam' )
                call cline%set('mkdir', 'no')
                call xestimate_diam%execute(cline)
            case default    
                l_did_execute = .false.
        end select
    end subroutine exec_nano2D_commander

end module single_exec_nano2D