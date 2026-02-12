!@descr: "cluster2D" UI api (concrete implementation)
module simple_ui_api_cluster2D
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: abinitio2D
type(ui_program), target :: cleanup2D
type(ui_program), target :: cluster2D
type(ui_program), target :: cluster2D_subsets
type(ui_program), target :: cluster_cavgs
type(ui_program), target :: cluster_cavgs_selection
type(ui_program), target :: cluster_stack
type(ui_program), target :: make_cavgs
type(ui_program), target :: map_cavgs_selection
type(ui_program), target :: match_cavgs
type(ui_program), target :: sample_classes
type(ui_program), target :: select_clusters
type(ui_program), target :: write_classes
type(ui_program), target :: write_mic_filetab

contains

    subroutine register_simple_ui_cluster2D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('abinitio2D',              abinitio2D,              prgtab)
        call add_ui_program('cleanup2D',               cleanup2D,               prgtab)
        call add_ui_program('cluster2D',               cluster2D,               prgtab)
        call add_ui_program('cluster2D_subsets',       cluster2D_subsets,       prgtab)
        call add_ui_program('cluster_cavgs',           cluster_cavgs,           prgtab)
        call add_ui_program('cluster_cavgs_selection', cluster_cavgs_selection, prgtab)
        call add_ui_program('cluster_stack',           cluster_stack,           prgtab)
        call add_ui_program('make_cavgs',              make_cavgs,              prgtab)
        call add_ui_program('map_cavgs_selection',     map_cavgs_selection,     prgtab)
        call add_ui_program('match_cavgs',             match_cavgs,             prgtab)
        call add_ui_program('sample_classes',          sample_classes,          prgtab)
        call add_ui_program('select_clusters',         select_clusters,         prgtab)
        call add_ui_program('write_classes',           write_classes,           prgtab)
        call add_ui_program('write_mic_filetab',       write_mic_filetab,       prgtab)
    end subroutine register_simple_ui_cluster2D

end module simple_ui_api_cluster2D
