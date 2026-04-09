!@descr: helper routines shared by tree-guided 2D search strategies
module simple_strategy2D_tree_utils
use simple_core_module_api
use simple_stat, only: power_sampling
use simple_strategy_tree_helpers, only: INVALID_CORR, MAX_NPEAKS, MAX_NTREES, choose_next_child_greedy, choose_next_child_prob, is_invalid_corr
use simple_strategy2D_alloc, only: s2D
use simple_strategy2D_srch, only: strategy2D_srch
implicit none

public :: peak_tree_selection, init_peak_tree_selection, select_peak_trees, &
&descend_tree_prob, descend_tree_greedy, get_tree_for_ref
private
#include "simple_local_flags.inc"

type :: peak_tree_selection
    integer :: ntrees       = 0
    integer :: npeak_target = 0
    integer :: npeak_trees  = 0
    integer :: peak_trees(MAX_NPEAKS)
    real    :: tree_best_corrs(MAX_NTREES)
    real    :: peak_tree_corrs(MAX_NPEAKS)
end type peak_tree_selection

contains

    subroutine init_peak_tree_selection( sel, ntrees, npeak_target )
        type(peak_tree_selection), intent(inout) :: sel
        integer,                   intent(in)    :: ntrees, npeak_target
        if( ntrees < 0 ) THROW_HARD('init_peak_tree_selection: ntrees must be >= 0')
        if( npeak_target < 0 ) THROW_HARD('init_peak_tree_selection: npeak_target must be >= 0')
        if( ntrees > MAX_NTREES ) THROW_HARD('init_peak_tree_selection: ntrees exceeds MAX_NTREES')
        if( npeak_target > MAX_NPEAKS ) THROW_HARD('init_peak_tree_selection: npeak_target exceeds MAX_NPEAKS')
        sel%npeak_trees  = 0
        sel%ntrees       = ntrees
        sel%npeak_target = npeak_target
        if( ntrees > 0 )then
            sel%tree_best_corrs(1:ntrees) = INVALID_CORR
        endif
        if( npeak_target > 0 )then
            sel%peak_trees(1:npeak_target)      = 0
            sel%peak_tree_corrs(1:npeak_target) = INVALID_CORR
        endif
    end subroutine init_peak_tree_selection

    subroutine select_peak_trees( sel )
        type(peak_tree_selection), intent(inout) :: sel
        integer :: loc(1)
        real    :: work_corrs(MAX_NTREES)
        if( sel%ntrees < 0 ) THROW_HARD('select_peak_trees: sel%ntrees must be >= 0')
        if( sel%npeak_target < 0 ) THROW_HARD('select_peak_trees: sel%npeak_target must be >= 0')
        if( sel%ntrees > MAX_NTREES ) THROW_HARD('select_peak_trees: sel%ntrees exceeds MAX_NTREES')
        if( sel%npeak_target > MAX_NPEAKS ) THROW_HARD('select_peak_trees: sel%npeak_target exceeds MAX_NPEAKS')
        sel%npeak_trees                         = 0
        sel%peak_trees(1:sel%npeak_target)      = 0
        sel%peak_tree_corrs(1:sel%npeak_target) = INVALID_CORR
        if( sel%ntrees <= 0 .or. sel%npeak_target <= 0 ) return
        work_corrs(1:sel%ntrees) = sel%tree_best_corrs(1:sel%ntrees)
        do while( sel%npeak_trees < sel%npeak_target )
            loc = maxloc(work_corrs(1:sel%ntrees))
            if( is_invalid_corr(work_corrs(loc(1))) ) exit
            sel%npeak_trees                      = sel%npeak_trees + 1
            sel%peak_trees(sel%npeak_trees)      = loc(1)
            sel%peak_tree_corrs(sel%npeak_trees) = work_corrs(loc(1))
            work_corrs(loc(1)) = INVALID_CORR
        end do
    end subroutine select_peak_trees

    subroutine descend_tree_prob( s, itree, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy2D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref(s, node_root%ref_idx, tree_best_corr, nrefs_tree)
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
            node_cur = s%b_ptr%block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            call eval_child_best(s, itree, node_cur%left_idx,  best_corr_L, nrefs_tree)
            call eval_child_best(s, itree, node_cur%right_idx, best_corr_R, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_L, best_corr_R)
            inode_next = choose_next_child_prob(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_prob

    subroutine descend_tree_greedy( s, itree, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy2D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref(s, node_root%ref_idx, tree_best_corr, nrefs_tree)
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
            node_cur = s%b_ptr%block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            call eval_child_best(s, itree, node_cur%left_idx,  best_corr_L, nrefs_tree)
            call eval_child_best(s, itree, node_cur%right_idx, best_corr_R, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_L, best_corr_R)
            inode_next = choose_next_child_greedy(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_greedy

    integer function get_tree_for_ref( s, iref, ntrees ) result(itree)
        class(strategy2D_srch), intent(in) :: s
        integer,                intent(in) :: iref, ntrees
        itree = 0
        if( iref < 1 ) return
        itree = s%b_ptr%subspace_full2sub_map(iref)
        if( ntrees > 0 )then
            if( itree < 1 .or. itree > ntrees ) itree = 0
        endif
    end function get_tree_for_ref

    subroutine eval_child_best( s, itree, child_idx, best_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy2D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(in)    :: child_idx
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_child
        best_corr = INVALID_CORR
        if( child_idx == 0 ) return
        node_child = s%b_ptr%block_tree%get_node(itree, child_idx)
        if( node_child%ref_idx == 0 ) return
        call eval_tree_ref(s, node_child%ref_idx, best_corr, nrefs_tree)
    end subroutine eval_child_best

    subroutine eval_tree_ref( s, ref_idx, best_corr, nrefs_tree )
        class(strategy2D_srch), intent(inout) :: s
        integer,                intent(in)    :: ref_idx
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        integer :: vec_nrots(s%nrots), inpl_ind, order_ind
        real    :: corr_tmp, inpl_corrs(s%nrots)
        best_corr = INVALID_CORR
        if( ref_idx < 1 ) return
        if( s2D%cls_pops(ref_idx) == 0 ) return
        if( s%l_sh_first )then
            call s%b_ptr%pftc%gen_objfun_vals(ref_idx, s%iptcl, s%xy_first, inpl_corrs)
        else
            call s%b_ptr%pftc%gen_objfun_vals(ref_idx, s%iptcl, [0.,0.], inpl_corrs)
        endif
        call power_sampling( s2D%power, s%nrots, inpl_corrs, vec_nrots, &
                            &s2D%snhc_smpl_ninpl, inpl_ind, order_ind, corr_tmp )
        nrefs_tree = nrefs_tree + 1
        best_corr  = corr_tmp
        call s%store_solution(ref_idx, inpl_ind, corr_tmp)
        if( corr_tmp > s%best_corr )then
            s%best_class = ref_idx
            s%best_rot   = inpl_ind
            s%best_corr  = corr_tmp
        endif
    end subroutine eval_tree_ref

end module simple_strategy2D_tree_utils
