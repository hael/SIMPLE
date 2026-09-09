!@descr: execution of cluster2D commanders
module simple_exec_cluster2D
use simple_cmdline,                         only: cmdline
use simple_string,                          only: string
use simple_exec_helpers,                    only: restarted_exec
use simple_commanders_project_cls,          only: commander_sample_classes
use simple_commanders_cluster2D,            only: commander_ppca_denoise_classes
use simple_commanders_mkcavgs,              only: commander_make_cavgs_distr, commander_bootstrap_cavgs, &
                                                  commander_unbootstrap_cavgs, commander_write_classes
use simple_commanders_abinitio2D,           only: commander_abinitio2D
use simple_stream_abinitio2D_chunks,        only: stream_abinitio2D_chunks
use simple_commanders_cavgs,                only: commander_map_cavgs_selection

implicit none

public :: exec_cluster2D_commander
private

type(commander_abinitio2D)                  :: xabinitio2D
type(stream_abinitio2D_chunks)              :: xabinitio2D_chunks
type(commander_make_cavgs_distr)            :: xmake_cavgs_distr
type(commander_bootstrap_cavgs)             :: xbootstrap_cavgs
type(commander_unbootstrap_cavgs)           :: xunbootstrap_cavgs
type(commander_map_cavgs_selection)         :: xmap_cavgs_selection
type(commander_sample_classes)              :: xsample_classes
type(commander_write_classes)               :: xwrite_classes

contains

    subroutine exec_cluster2D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'abinitio2D' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('abinitio2D'), string('simple_exec'))
                else
                    call xabinitio2D%execute(cline)
                endif
            case( 'abinitio2D_chunks' )
                call xabinitio2D_chunks%execute(cline)       
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
    end subroutine exec_cluster2D_commander

end module simple_exec_cluster2D
