!@descr: execution of refine2D commanders
module simple_exec_refine2D
use simple_cmdline,                         only: cmdline
use simple_string,                          only: string
use simple_exec_helpers,                    only: restarted_exec
use simple_commanders_project_cls,          only: commander_sample_classes
use simple_commanders_refine2D,             only: commander_ppca_denoise_classes
use simple_commanders_mkcavgs,              only: commander_make_cavgs_distr, commander_bootstrap_cavgs, &
                                                  commander_unbootstrap_cavgs, commander_write_classes
use simple_commanders_solve2D,              only: commander_solve2D, commander_solve2D_chunks
use simple_commanders_cavgs,                only: commander_map_cavgs_selection

implicit none

public :: exec_refine2D_commander
private

type(commander_solve2D)                     :: xsolve2D
type(commander_solve2D_chunks)              :: xsolve2D_chunks
type(commander_make_cavgs_distr)            :: xmake_cavgs_distr
type(commander_bootstrap_cavgs)             :: xbootstrap_cavgs
type(commander_unbootstrap_cavgs)           :: xunbootstrap_cavgs
type(commander_map_cavgs_selection)         :: xmap_cavgs_selection
type(commander_sample_classes)              :: xsample_classes
type(commander_write_classes)               :: xwrite_classes

contains

    subroutine exec_refine2D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'solve2D' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('solve2D'), string('simple_exec'))
                else
                    call xsolve2D%execute(cline)
                endif
            case( 'solve2D_chunks' )
                call xsolve2D_chunks%execute(cline)       
            case( 'make_cavgs' )
                call xmake_cavgs_distr%execute(cline)
            case( 'bootstrap_cavgs' )
                call xbootstrap_cavgs%execute(cline)
            case( 'unbootstrap_cavgs' )
                call xunbootstrap_cavgs%execute(cline)
            case( 'map_cavgs_selection' )
                call xmap_cavgs_selection%execute(cline)
            case( 'sample_classes' )
                call xsample_classes%execute(cline)
            case( 'write_classes' )
                 call xwrite_classes%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_refine2D_commander

end module simple_exec_refine2D
