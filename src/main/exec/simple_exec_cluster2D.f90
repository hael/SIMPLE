!@descr: execution of cluster2D commanders
module simple_exec_cluster2D
use simple_cmdline,                  only: cmdline
use simple_string,                   only: string
use simple_exec_helpers,             only: restarted_exec
use simple_commanders_project_cls,   only: commander_sample_classes
use simple_commanders_cluster2D,     only: commander_cluster2D_autoscale, commander_ppca_denoise_classes
use simple_commanders_mkcavgs,       only: commander_make_cavgs_distr,  commander_write_classes
use simple_commanders_abinitio2D,    only: commander_abinitio2D
use simple_commanders_cleanup2D,     only: commander_cleanup2D
use simple_stream_cluster2D_subsets, only: stream_cluster2D_subsets
use simple_commanders_cavgs,         only: commander_map_cavgs_selection
implicit none

public :: exec_cluster2D_commander
private

type(commander_abinitio2D)                  :: xabinitio2D
type(commander_cleanup2D)                   :: xcleanup2D 
type(commander_cluster2D_autoscale)         :: xcluster2D
type(stream_cluster2D_subsets)              :: xcluster2D_subsets
type(commander_make_cavgs_distr)            :: xmake_cavgs_distr
type(commander_map_cavgs_selection)         :: xmap_cavgs_selection
type(commander_sample_classes)              :: xsample_classes
type(commander_write_classes)               :: xwrite_classes

contains

    subroutine exec_cluster2D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'abinitio2D' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('abinitio2D'), string('simple_exec'))
                else
                    call xabinitio2D%execute(cline)
                endif
            case( 'cleanup2D' )
                call xcleanup2D%execute(cline)
            case( 'cluster2D' )
                call xcluster2D%execute(cline)
            case( 'cluster2D_subsets' )
                call xcluster2D_subsets%execute(cline)
            case( 'make_cavgs' )
                call xmake_cavgs_distr%execute(cline)
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
