!@descr: 2D pool analysis
module simple_commanders_pool2D
use simple_commanders_api
use simple_commanders_abinitio2D, only: commander_abinitio2D
use simple_commanders_cavgs
use simple_commanders_project_core, only: commander_selection
use simple_imgarr_utils, only: read_cavgs_into_imgarr, dealloc_imgarr
use simple_block_tree_corr, only: gen_corr_block_tree_hac
use simple_block_tree_io, only: write_block_tree
use simple_multi_dendro, only: multi_dendro
use simple_abinitio2D_utils_cluster2D 
implicit none

public :: commander_pool2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_pool2D
    contains
        procedure :: execute      => exec_pool2D
end type commander_pool2D

contains

    subroutine exec_pool2D( self, cline )
        class(commander_pool2D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters)                :: params, params_blocktree
        type(sp_project)                :: spproj, spproj_full
        type(cmdline)                   :: cline_selection, cline_abinitio2D, cline_abinitio2D_full, cline_blocktree
        type(commander_selection)       :: xselection
        type(commander_abinitio2D)      :: xabinitio2D
        type(multi_dendro)              :: block_tree
        type(image),        allocatable :: cavg_imgs(:)
        type(string)                    :: refs_fname, blocktree_fname, projfile_orig, projfile_pruned
        integer                         :: nptcls_tot, nran, nrefs
        real                            :: smpd

        if( .not. cline%defined('prune') ) call cline%set('prune', 'no')
        call params%new(cline)
        projfile_orig   = params%projfile
        projfile_pruned = 'pool2D_pruned.simple'
        call simple_copy_file(projfile_orig, projfile_pruned)

        call spproj_full%read(params%projfile)
        nptcls_tot = spproj_full%os_ptcl2D%get_noris()
        nran = nint(0.2 * real(nptcls_tot))
        if( cline%defined('blocktree') )then
            blocktree_fname = cline%get_carg('blocktree')
        else
            blocktree_fname = 'pool_block_tree.bin'
        endif

        cline_selection = cline
        call cline_selection%set('projfile', projfile_pruned)
        call cline_selection%set('prune',   'yes')
        call cline_selection%set('oritype', 'ptcl2D')
        call cline_selection%set('nran',    nran)
        call xselection%execute(cline_selection)

        cline_abinitio2D = cline
        call cline_abinitio2D%set('projfile', projfile_pruned)
        call cline_abinitio2D%set('mkdir', 'no')
        call cline_abinitio2D%set('ncls', params%ncls)
        call xabinitio2D%execute(cline_abinitio2D)
        call params%new(cline_abinitio2D)

        call spproj%read(params%projfile)
        call spproj%get_cavgs_stk(refs_fname, nrefs, smpd, imgkind='cavg')
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        cline_blocktree = cline_abinitio2D
        call cline_blocktree%set('projfile', projfile_pruned)
        call cline_blocktree%set('trs',    5.)
        call cline_blocktree%set('ctf',    'no')
        call cline_blocktree%set('objfun', 'cc')
        call params_blocktree%new(cline_blocktree)
        ! need a cleaning step to remove any empty cavgs if any end up with zero particles after pruning, otherwise the block tree generation will fail
        block_tree = gen_corr_block_tree_hac(cavg_imgs, params_blocktree, params%ncls_sub)

        call write_block_tree(block_tree, blocktree_fname)
        call block_tree%kill
        call dealloc_imgarr(cavg_imgs)
        call spproj%kill
        call spproj_full%kill

        cline_abinitio2D_full = cline
        call cline_abinitio2D_full%set('mkdir',     'yes')
        call cline_abinitio2D_full%set('refine',    'snhc_ptree')
        call cline_abinitio2D_full%set('oritype',   'ptcl2D')
        call cline_abinitio2D_full%set('ncls',      params%ncls)
        call cline_abinitio2D_full%set('projfile',  projfile_orig)
        call cline_abinitio2D_full%set('refs',      refs_fname)
        call cline_abinitio2D_full%set('blocktree', blocktree_fname)
        call xabinitio2D%execute(cline_abinitio2D_full)
    end subroutine exec_pool2D
end module simple_commanders_pool2D
