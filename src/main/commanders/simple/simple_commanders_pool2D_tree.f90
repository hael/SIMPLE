!@descr: 2D pool analysis
module simple_commanders_pool2D_tree
use simple_commanders_api
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cluster2D,    only: commander_cluster2D
use simple_commanders_project_core, only: commander_selection
use simple_imgarr_utils,            only: dealloc_imgarr, write_imgarr
use simple_block_tree_corr,         only: gen_corr_block_tree_hac
use simple_block_tree_io,           only: write_block_tree
use simple_multi_dendro,            only: multi_dendro
use simple_strategy2D_utils,        only: id_junk_and_prep_cavgs4clust
implicit none

public :: commander_pool2D_tree
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_pool2D_tree
    contains
    procedure :: execute      => exec_pool2D_tree
end type commander_pool2D_tree

contains

    subroutine exec_pool2D_tree( self, cline )
        class(commander_pool2D_tree), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        real, parameter            :: FRAC_SEL = 0.2
        type(parameters)           :: params
        type(sp_project)           :: spproj, spproj2
        type(commander_selection)  :: xselection
        type(commander_abinitio2D) :: xabinitio2D
        type(commander_cluster2D)  :: xcluster2D
        type(multi_dendro)         :: block_tree
        type(string)               :: refs_prep, blocktree_fname
        type(image),   allocatable :: cavg_imgs(:)
        logical,       allocatable :: l_non_junk(:)
        integer,       allocatable :: clsinds(:), clspops(:)
        real,          allocatable :: mm(:,:)
        integer                    :: nptcls_tot, nran, nrefs_non_junk
        ! select FRAC_SEL particles for tree construction
        call params%new(cline)
        call spproj%read(params%projfile)
        call spproj2%copy(spproj)
        nptcls_tot = spproj%os_ptcl2D%get_noris()
        nran = nint(FRAC_SEL * nptcls_tot)
        call cline%set('projfile', params%projfile)
        call cline%set('prune',   'no')
        call cline%set('oritype', 'ptcl2D')
        call cline%set('nran',    nran)
        call xselection%execute(cline)
        ! generate class averages for tree construction
        call cline%set('ncls', params%ncls)
        call cline%set('mkdir', 'no')
        call xabinitio2D%execute(cline)
        ! prep for block-tree generation    
        call params%new(cline)
        call spproj%read(params%projfile)
        call id_junk_and_prep_cavgs4clust(spproj, cavg_imgs, params%mskdiam, clspops, clsinds, l_non_junk, mm )
        nrefs_non_junk = size(cavg_imgs)
        refs_prep = 'pool2D_tree_refs.mrc'
        call write_imgarr(cavg_imgs, refs_prep)
        call cline%set('trs',    5.)
        call cline%set('objfun', 'cc')
        call params%new(cline)
        if( cline%defined('blocktree') )then
            blocktree_fname = params%blocktree
        else
            params%blocktree = 'pool_block_tree.bin'
        endif
        block_tree = gen_corr_block_tree_hac(cavg_imgs, params, params%ncls_sub)
        call write_block_tree(block_tree, params%blocktree)
        call block_tree%kill
        call spproj2%write(params%projfile)
        ! run greedy tree-based search on all particles
        call cline%set('trs', 5.)
        call cline%set('lp', 6.0)
        call cline%set('maxits', 1)
        call cline%set('ncls', nrefs_non_junk)
        call cline%set('refs', refs_prep)
        call cline%set('blocktree', params%blocktree)
        call cline%set('box_crop', params%box_crop)
        call cline%set('refine', 'snhc_ptree')
        call xcluster2D%execute(cline)
        ! cleanup
        call spproj%kill
        call dealloc_imgarr(cavg_imgs)
        call simple_end('**** POOL2D_TREE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_pool2D_tree

end module simple_commanders_pool2D_tree
