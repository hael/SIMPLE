module simple_exec_abinitio3D
use simple_cmdline,             only: cmdline
use simple_string,              only: string
use simple_exec_helpers,        only: restarted_exec
use simple_commanders_abinitio, only: commander_abinitio3D_cavgs, commander_abinitio3D, commander_multivol_assign
use simple_commanders_volops,   only: commander_noisevol
use simple_commanders_resolest, only: commander_estimate_lpstages
implicit none

public :: exec_abinitio3D_commander
private

type(commander_abinitio3D)        :: xabinitio3D
type(commander_abinitio3D_cavgs)  :: xabinitio3D_cavgs
type(commander_estimate_lpstages) :: xestimate_lpstages
type(commander_multivol_assign)   :: xmultivol_assign
type(commander_noisevol)          :: xnoisevol

contains

    subroutine exec_abinitio3D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'abinitio3D' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('abinitio3D'), string('simple_exec'))
                else
                    call xabinitio3D%execute(cline)
                endif
            case( 'abinitio3D_cavgs' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('abinitio3D_cavgs'), string('simple_exec'))
                else
                    call xabinitio3D_cavgs%execute(cline)
                endif
            case( 'estimate_lpstages' )
                call xestimate_lpstages%execute(cline)
            case( 'multivol_assign' )
                call xmultivol_assign%execute(cline)
            case( 'noisevol' )
                call xnoisevol%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_abinitio3D_commander

end module simple_exec_abinitio3D
