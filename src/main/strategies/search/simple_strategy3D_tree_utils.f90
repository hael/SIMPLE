!@descr: helper routines shared by tree-guided 3D search strategies
! Other multi-state tree strategies to consider:
! *** One-tree, per-state path competition
! Keep the same geometry-selected single tree, but run one descend_tree_prob_fixed_state per state and pick the best final state/path. 
! This costs about the same order as free-state descent, but enforces state consistency along each path.
! *** Top-level multi-state, deep fixed-state
! Use descend_tree_prob for the root and maybe the first 1–2 levels, then lock to the current best state and continue with 
! descend_tree_prob_fixed_state. This is probably the most practical follow-up once the random split starts separating meaningfully.
! *** 2–3 tree beam from geometry
! Still no exhaustive evaluation, but instead of one nearest tree, take the 2 or 3 nearest trees and descend each probabilistically 
! across states. This is much safer if the previous orientation is occasionally off.
! ***  Greedy late-stage version
! A ptree_neigh_states_greedy using descend_tree_greedy instead of descend_tree_prob. Use this after the classes are already fairly stable.
module simple_strategy3D_tree_utils
use simple_core_module_api
use simple_multi_dendro, only: multi_dendro
use simple_strategy3D_alloc
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

real,    parameter :: INVALID_CORR        = -huge(1.0)
real,    parameter :: INVALID_CORR_THRESH = INVALID_CORR / 2.0
integer, parameter :: MAX_NTREES          = 2500
integer, parameter :: MAX_NPEAKS          = 64
integer, parameter :: MAX_TREE_REFS       = 1024

public :: peak_tree_selection, init_peak_tree_selection, select_peak_trees, select_peak_trees_per_state,&
descend_tree_prob, descend_tree_prob_fixed_state, descend_tree_greedy, descend_tree_greedy_fixed_state,&
&get_tree_for_ref, MAX_NTREES, MAX_NPEAKS, MAX_TREE_REFS
private
#include "simple_local_flags.inc"

type :: peak_tree_selection
    integer :: ntrees           = 0
    integer :: npeak_target     = 0
    integer :: nstates          = 1
    integer :: npeaks_per_state = 1
    integer :: npeak_trees      = 0
    integer :: tree_best_states(MAX_NTREES)
    integer :: peak_trees(MAX_NPEAKS)
    integer :: peak_tree_states(MAX_NPEAKS)
    logical :: has_tree_states  = .false.
    logical :: has_peak_states  = .false.
    real    :: tree_best_corrs(MAX_NTREES)
    real    :: peak_tree_corrs(MAX_NPEAKS)
end type peak_tree_selection

contains

    subroutine init_peak_tree_selection( sel, ntrees, npeak_target, need_tree_states, need_peak_states, nstates, npeaks_per_state )
        type(peak_tree_selection), intent(inout) :: sel
        integer,                   intent(in)    :: ntrees, npeak_target
        logical,                   intent(in)    :: need_tree_states, need_peak_states
        integer, optional,         intent(in)    :: nstates, npeaks_per_state
        sel%npeak_trees      = 0
        sel%ntrees           = ntrees
        sel%npeak_target     = npeak_target
        sel%nstates          = 1
        sel%npeaks_per_state = 1
        sel%has_tree_states  = need_tree_states
        sel%has_peak_states  = need_peak_states
        if( present(nstates) )          sel%nstates          = nstates
        if( present(npeaks_per_state) ) sel%npeaks_per_state = npeaks_per_state
        if( ntrees > 0 )then
            sel%tree_best_corrs(1:ntrees) = INVALID_CORR
            if( need_tree_states ) sel%tree_best_states(1:ntrees) = 1
        endif
        if( npeak_target > 0 )then
            sel%peak_trees(1:npeak_target)      = 0
            sel%peak_tree_corrs(1:npeak_target) = INVALID_CORR
            if( need_peak_states ) sel%peak_tree_states(1:npeak_target) = 1
        endif
    end subroutine init_peak_tree_selection

    subroutine select_peak_trees( sel )
        type(peak_tree_selection), intent(inout) :: sel
        integer :: loc(1)
        real    :: work_corrs(MAX_NTREES)
        sel%npeak_trees                         = 0
        sel%peak_trees(1:sel%npeak_target)      = 0
        sel%peak_tree_corrs(1:sel%npeak_target) = INVALID_CORR
        if( sel%has_peak_states ) sel%peak_tree_states(1:sel%npeak_target) = 1
        work_corrs(1:sel%ntrees) = sel%tree_best_corrs(1:sel%ntrees)
        do while( sel%npeak_trees < sel%npeak_target )
            loc = maxloc(work_corrs(1:sel%ntrees))
            if( is_invalid_corr(work_corrs(loc(1))) ) exit
            sel%npeak_trees                      = sel%npeak_trees + 1
            sel%peak_trees(sel%npeak_trees)      = loc(1)
            sel%peak_tree_corrs(sel%npeak_trees) = work_corrs(loc(1))
            if( sel%has_peak_states ) sel%peak_tree_states(sel%npeak_trees) = 1
            work_corrs(loc(1)) = INVALID_CORR
        end do
    end subroutine select_peak_trees

    subroutine select_peak_trees_per_state( sel )
        type(peak_tree_selection), intent(inout) :: sel
        integer :: istate, irank, best_tree, itree
        real    :: best_corr
        real    :: work_corrs(MAX_NTREES)
        sel%npeak_trees = 0
        sel%peak_trees(1:sel%npeak_target)       = 0
        sel%peak_tree_corrs(1:sel%npeak_target)  = INVALID_CORR
        sel%peak_tree_states(1:sel%npeak_target) = 1
        work_corrs(1:sel%ntrees) = sel%tree_best_corrs(1:sel%ntrees)
        state_loop: do istate = 1, sel%nstates
            do irank = 1, sel%npeaks_per_state
                best_corr = INVALID_CORR
                best_tree = 0
                do itree = 1, sel%ntrees
                    if( sel%tree_best_states(itree) /= istate ) cycle
                    if( work_corrs(itree) <= best_corr ) cycle
                    best_corr = work_corrs(itree)
                    best_tree = itree
                end do
                if( best_tree == 0 ) exit
                if( sel%npeak_trees >= sel%npeak_target ) exit state_loop
                sel%npeak_trees = sel%npeak_trees + 1
                sel%peak_trees(sel%npeak_trees)       = best_tree
                sel%peak_tree_corrs(sel%npeak_trees)  = best_corr
                sel%peak_tree_states(sel%npeak_trees) = istate
                work_corrs(best_tree) = INVALID_CORR
            end do
        end do state_loop
    end subroutine select_peak_trees_per_state

    ! Stochastic descent with multi-state evaluation at each node.
    subroutine descend_tree_prob( s, itree, coarse_tree_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        real,                   intent(in)    :: coarse_tree_corr
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref_across_states(s, node_root%ref_idx, tree_best_corr, nrefs_tree)
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

    ! Fixed-state stochastic descent
    subroutine descend_tree_prob_fixed_state( s, itree, coarse_tree_corr, nrefs_tree, istate_fixed )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        real,                   intent(in)    :: coarse_tree_corr
        integer,                intent(inout) :: nrefs_tree
        integer,                intent(in)    :: istate_fixed
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref_fixed_state(s, node_root%ref_idx, istate_fixed, tree_best_corr, nrefs_tree)
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
            node_cur = s%b_ptr%block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            call eval_child_best_fixed_state(s, itree, node_cur%left_idx,  istate_fixed, best_corr_L, nrefs_tree)
            call eval_child_best_fixed_state(s, itree, node_cur%right_idx, istate_fixed, best_corr_R, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_L, best_corr_R)
            inode_next = choose_next_child_prob(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_prob_fixed_state

    ! Greedy descent with multi-state evaluation at each node. This mirrors
    ! srch_eul_bl_tree: always descend via the child with the best local score,
    ! while tracking the best node seen
    subroutine descend_tree_greedy( s, itree, coarse_tree_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        real,                   intent(in)    :: coarse_tree_corr
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref_across_states(s, node_root%ref_idx, tree_best_corr, nrefs_tree)
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

    ! Fixed-state greedy descent.
    subroutine descend_tree_greedy_fixed_state( s, itree, coarse_tree_corr, nrefs_tree, istate_fixed )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        real,                   intent(in)    :: coarse_tree_corr
        integer,                intent(inout) :: nrefs_tree
        integer,                intent(in)    :: istate_fixed
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        node_root = s%b_ptr%block_tree%get_root_node(itree)
        call eval_tree_ref_fixed_state(s, node_root%ref_idx, istate_fixed, tree_best_corr, nrefs_tree)
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
            node_cur = s%b_ptr%block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            call eval_child_best_fixed_state(s, itree, node_cur%left_idx,  istate_fixed, best_corr_L, nrefs_tree)
            call eval_child_best_fixed_state(s, itree, node_cur%right_idx, istate_fixed, best_corr_R, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_L, best_corr_R)
            inode_next = choose_next_child_greedy(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_greedy_fixed_state

    integer function get_tree_for_ref( s, iref, ntrees ) result(itree)
        class(strategy3D_srch), intent(in) :: s
        integer,                intent(in) :: iref, ntrees
        integer :: iproj
        iproj = s3D%proj_space_proj(iref)
        itree = s%b_ptr%subspace_full2sub_map(iproj)
    end function get_tree_for_ref

    subroutine eval_child_best( s, itree, child_idx, best_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(in)    :: child_idx
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_child
        best_corr = INVALID_CORR
        if( child_idx == 0 ) return
        node_child = s%b_ptr%block_tree%get_node(itree, child_idx)
        if( node_child%ref_idx == 0 ) return
        call eval_tree_ref_across_states(s, node_child%ref_idx, best_corr, nrefs_tree)
    end subroutine eval_child_best

    subroutine eval_child_best_fixed_state( s, itree, child_idx, istate_fixed, best_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(in)    :: child_idx
        integer,                intent(in)    :: istate_fixed
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        type(bt_node) :: node_child
        best_corr = INVALID_CORR
        if( child_idx == 0 ) return
        node_child = s%b_ptr%block_tree%get_node(itree, child_idx)
        if( node_child%ref_idx == 0 ) return
        call eval_tree_ref_fixed_state(s, node_child%ref_idx, istate_fixed, best_corr, nrefs_tree)
    end subroutine eval_child_best_fixed_state

    subroutine eval_tree_ref_across_states( s, ref_idx, best_corr, nrefs_tree )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: ref_idx
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        integer :: inds(s%nrots)
        integer :: iref
        integer :: istate
        integer :: loc(1)
        real    :: corr_tmp, inpl_corrs(s%nrots), sorted_corrs(s%nrots)
        best_corr = INVALID_CORR
        do istate = 1, s%nstates
            iref = (istate - 1) * s%p_ptr%nspace + ref_idx
            if( .not. s3D%state_exists(s3D%proj_space_state(iref)) ) cycle
            if( s%p_ptr%l_doshift )then
                call s%b_ptr%pftc%gen_objfun_vals(iref, s%iptcl, s%xy_first, inpl_corrs)
            else
                call s%b_ptr%pftc%gen_objfun_vals(iref, s%iptcl, [0.,0.], inpl_corrs)
            endif
            if( s%p_ptr%l_prob_inpl )then
                loc = angle_sampling( &
                    eulprob_dist_switch(inpl_corrs, s%p_ptr%cc_objfun), &
                    sorted_corrs, inds, &
                    s3D%smpl_inpl_athres(s3D%proj_space_state(iref)), &
                    s%p_ptr%prob_athres )
            else
                loc = maxloc(inpl_corrs)
            endif
            corr_tmp = inpl_corrs(loc(1))
            call s%store_solution(iref, loc(1), corr_tmp)
            nrefs_tree = nrefs_tree + 1
            best_corr  = max(best_corr, corr_tmp)
        end do
    end subroutine eval_tree_ref_across_states

    subroutine eval_tree_ref_fixed_state( s, ref_idx, istate_fixed, best_corr, nrefs_tree )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: ref_idx
        integer,                intent(in)    :: istate_fixed
        real,                   intent(out)   :: best_corr
        integer,                intent(inout) :: nrefs_tree
        integer :: inds(s%nrots)
        integer :: iref
        integer :: loc(1)
        real    :: corr_tmp, inpl_corrs(s%nrots), sorted_corrs(s%nrots)
        iref = (istate_fixed - 1) * s%p_ptr%nspace + ref_idx
        if( .not. s3D%state_exists(s3D%proj_space_state(iref)) )then
            best_corr = INVALID_CORR
            return
        endif
        if( s%p_ptr%l_doshift )then
            call s%b_ptr%pftc%gen_objfun_vals(iref, s%iptcl, s%xy_first, inpl_corrs)
        else
            call s%b_ptr%pftc%gen_objfun_vals(iref, s%iptcl, [0.,0.], inpl_corrs)
        endif
        if( s%p_ptr%l_prob_inpl )then
            loc = angle_sampling( &
                eulprob_dist_switch(inpl_corrs, s%p_ptr%cc_objfun), &
                sorted_corrs, inds, &
                s3D%smpl_inpl_athres(s3D%proj_space_state(iref)), &
                s%p_ptr%prob_athres )
        else
            loc = maxloc(inpl_corrs)
        endif
        corr_tmp = inpl_corrs(loc(1))
        call s%store_solution(iref, loc(1), corr_tmp)
        nrefs_tree = nrefs_tree + 1
        best_corr  = corr_tmp
    end subroutine eval_tree_ref_fixed_state

    subroutine exhaustive_tree_scan( s, itree, tree_best_corr, nrefs_tree )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        real,                   intent(inout) :: tree_best_corr
        integer,                intent(inout) :: nrefs_tree
        integer :: tree_refs(MAX_TREE_REFS)
        integer :: n_tree_refs, iref_tree
        real    :: best_corr_ref
        call s%b_ptr%block_tree%get_tree_refs_static(itree, tree_refs, n_tree_refs)
        do iref_tree = 1, n_tree_refs
            call eval_tree_ref_across_states(s, tree_refs(iref_tree), best_corr_ref, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_ref)
        end do
    end subroutine exhaustive_tree_scan

    subroutine exhaustive_tree_scan_fixed_state( s, itree, istate_fixed, tree_best_corr, nrefs_tree )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(in)    :: istate_fixed
        real,                   intent(inout) :: tree_best_corr
        integer,                intent(inout) :: nrefs_tree
        integer :: tree_refs(MAX_TREE_REFS)
        integer :: n_tree_refs, iref_tree
        real    :: best_corr_ref
        call s%b_ptr%block_tree%get_tree_refs_static(itree, tree_refs, n_tree_refs)
        do iref_tree = 1, n_tree_refs
            call eval_tree_ref_fixed_state(s, tree_refs(iref_tree), istate_fixed, best_corr_ref, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_ref)
        end do
    end subroutine exhaustive_tree_scan_fixed_state

    integer function choose_next_child_prob( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
        integer, intent(in) :: left_idx, right_idx
        real,    intent(in) :: corr_left, corr_right
        real :: cmax
        real :: p_left
        real :: p_right
        inode_next = 0
        if( left_idx == 0 .and. right_idx == 0 ) return
        if( left_idx == 0 )then
            inode_next = right_idx
            return
        endif
        if( right_idx == 0 )then
            inode_next = left_idx
            return
        endif
        if( is_invalid_corr(corr_left) .and. is_invalid_corr(corr_right) )then
            if( sample_two(1.0, 1.0) == 1 )then
                inode_next = left_idx
            else
                inode_next = right_idx
            endif
            return
        endif
        cmax = max(corr_left, corr_right)
        if( is_invalid_corr(corr_left) )then
            p_left = 0.0
        else
            p_left = exp(corr_left - cmax)
        endif
        if( is_invalid_corr(corr_right) )then
            p_right = 0.0
        else
            p_right = exp(corr_right - cmax)
        endif
        if( sample_two(p_left, p_right) == 1 )then
            inode_next = left_idx
        else
            inode_next = right_idx
        endif
    end function choose_next_child_prob

    integer function choose_next_child_greedy( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
        integer, intent(in) :: left_idx, right_idx
        real,    intent(in) :: corr_left, corr_right
        inode_next = 0
        if( left_idx == 0 .and. right_idx == 0 ) return
        if( left_idx == 0 )then
            inode_next = right_idx
            return
        endif
        if( right_idx == 0 )then
            inode_next = left_idx
            return
        endif
        if( is_invalid_corr(corr_left) .and. is_invalid_corr(corr_right) )then
            inode_next = left_idx
            return
        endif
        if( is_invalid_corr(corr_left) )then
            inode_next = right_idx
            return
        endif
        if( is_invalid_corr(corr_right) )then
            inode_next = left_idx
            return
        endif
        inode_next = merge(left_idx, right_idx, corr_left >= corr_right)
    end function choose_next_child_greedy

    logical pure function is_invalid_corr( corr ) result(invalid)
        real, intent(in) :: corr
        invalid = corr <= INVALID_CORR_THRESH
    end function is_invalid_corr

    integer function sample_two( p1, p2 ) result(which)
        real, intent(in) :: p1, p2
        real :: psum
        real :: r
        psum = p1 + p2
        if( psum <= 0.0 )then
            which = merge(1, 2, ran3() < 0.5)
            return
        endif
        r = ran3()
        which = merge(1, 2, r < p1 / psum)
    end function sample_two

end module simple_strategy3D_tree_utils