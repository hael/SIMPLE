module simple_exec_refine3D
use simple_cmdline,             only: cmdline
use simple_string,              only: string
use simple_exec_helpers,        only: restarted_exec
use simple_commanders_mask,     only: commander_automask
use simple_commanders_volops,   only: commander_postprocess
use simple_commanders_rec,      only: commander_reconstruct3D_distr
use simple_commanders_refine3D, only: commander_refine3D_distr, commander_refine3D_auto
implicit none

type(commander_automask)            :: xautomask
type(commander_postprocess)         :: xpostprocess
type(commander_reconstruct3D_distr) :: xreconstruct3D
type(commander_refine3D_auto)       :: xrefine3D_auto
type(commander_refine3D_distr)      :: xrefine3D_distr

contains

    subroutine exec_refine3D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'automask' )
                call xautomask%execute(cline)
            case( 'postprocess' )
                call xpostprocess%execute(cline)
            case( 'reconstruct3D' )
                call xreconstruct3D%execute(cline)
            case( 'refine3D' )
                if( cline%defined('nrestarts') )then
                call restarted_exec(cline, string('refine3D'), string('simple_exec'))
                else
                    call xrefine3D_distr%execute(cline)
                endif
            case( 'refine3D_auto' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('refine3D_auto'), string('simple_exec'))
                else
                    call xrefine3D_auto%execute(cline)
                endif
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_refine3D_commander

end module simple_exec_refine3D