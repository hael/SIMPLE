!@descr: execution of solve3D commanders
module simple_exec_solve3D
use simple_cmdline,             only: cmdline
use simple_string,              only: string
use simple_exec_helpers,        only: restarted_exec, exec_screen
use simple_commanders_solve3D, only: commander_solve3D_cavgs, commander_solve3D,&
                                    & commander_solve3D_cavgs_conditional_restarts, commander_solve3D_addon
use simple_commanders_volops,   only: commander_noisevol
use simple_commanders_resolest, only: commander_estimate_lpstages
implicit none

public :: exec_solve3D_commander
private

type(commander_solve3D)                               :: xsolve3D
type(commander_solve3D_cavgs)                         :: xsolve3D_cavgs
type(commander_solve3D_addon)                         :: xsolve3D_addon
type(commander_solve3D_cavgs_conditional_restarts) :: xsolve3D_cavgs_conditional_restarts
type(commander_estimate_lpstages)                     :: xestimate_lpstages
type(commander_noisevol)                              :: xnoisevol

contains

    subroutine exec_solve3D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'solve3D' )
                if( cline%defined('screen') )then
                    call exec_screen(cline, string('solve3D'), string('simple_exec'))
                else if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('solve3D'), string('simple_exec'))
                else
                    call xsolve3D%execute(cline)
                endif
            case( 'solve3D_cavgs' )
                if( cline%defined('nrestarts_collapse') )then
                    call xsolve3D_cavgs_conditional_restarts%execute(cline)
                else if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('solve3D_cavgs'), string('simple_exec'))
                else
                    call xsolve3D_cavgs%execute(cline)
                endif
            case( 'solve3D_addon' )
                call xsolve3D_addon%execute(cline)
            case( 'estimate_lpstages' )
                call xestimate_lpstages%execute(cline)
            case( 'noisevol' )
                call xnoisevol%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_solve3D_commander

end module simple_exec_solve3D
