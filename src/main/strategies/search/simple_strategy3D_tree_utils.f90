!@descr: helper routines shared by tree-guided 3D search strategies
module simple_strategy3D_tree_utils
use simple_core_module_api
use simple_multi_dendro, only: multi_dendro
use simple_strategy3D_alloc
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

real, parameter :: INVALID_CORR        = -huge(1.0)
real, parameter :: INVALID_CORR_THRESH = INVALID_CORR / 2.0

public :: select_peak_trees, descend_tree_prob, descend_tree_prob_fixed_state, &
          &descend_tree_bestfirst, descend_tree_bestfirst_fixed_state, get_tree_for_ref, trace_tree_prob
private
#include "simple_local_flags.inc"

abstract interface
    subroutine tree_ref_eval(ref_idx, iptcl, ithr, best_corr)
        import
        integer, intent(in)  :: ref_idx
        integer, intent(in)  :: iptcl
        integer, intent(in)  :: ithr
        real,    intent(out) :: best_corr
    end subroutine tree_ref_eval
end interface

contains

    subroutine select_peak_trees( tree_best_corrs, peak_trees, peak_tree_corrs, npeak_trees )
        real,    intent(in)  :: tree_best_corrs(:)
        integer, intent(out) :: peak_trees(:)
        real,    intent(out) :: peak_tree_corrs(:)
        integer, intent(out) :: npeak_trees
        integer :: loc(1)
        real    :: work_corrs(size(tree_best_corrs))
        peak_trees      = 0
        peak_tree_corrs = INVALID_CORR
        npeak_trees     = 0
        if( size(tree_best_corrs) == 0 ) return
        if( size(peak_trees)      == 0 ) return
        work_corrs = tree_best_corrs
        do while( npeak_trees < size(peak_trees) )
            loc = maxloc(work_corrs)
            if( is_invalid_corr(work_corrs(loc(1))) ) exit
            npeak_trees                  = npeak_trees + 1
            peak_trees(npeak_trees)      = loc(1)
            peak_tree_corrs(npeak_trees) = work_corrs(loc(1))
            work_corrs(loc(1))           = INVALID_CORR
        end do
    end subroutine select_peak_trees

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

    ! Best-first descent with multi-state evaluation at each node. This mirrors
    ! srch_eul_bl_tree: always descend via the child with the best local score,
    ! while tracking the best node seen
    subroutine descend_tree_bestfirst( s, itree, coarse_tree_corr, nrefs_tree )
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
            inode_next = choose_next_child_bestfirst(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_bestfirst

    ! Fixed-state best-first descent.
    subroutine descend_tree_bestfirst_fixed_state( s, itree, coarse_tree_corr, nrefs_tree, istate_fixed )
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
            inode_next = choose_next_child_bestfirst(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine descend_tree_bestfirst_fixed_state

    integer function get_tree_for_ref( s, iref, ntrees ) result(itree)
        class(strategy3D_srch), intent(in) :: s
        integer,                intent(in) :: iref, ntrees
        integer :: iproj
        iproj = s3D%proj_space_proj(iref)
        itree = s%b_ptr%subspace_full2sub_map(iproj)
    end function get_tree_for_ref

    subroutine trace_tree_prob( block_tree, itree, iptcl, ithr, eval_ref )
        use simple_binary_tree, only: bt_node
        class(multi_dendro), intent(in) :: block_tree
        integer,             intent(in) :: itree
        integer,             intent(in) :: iptcl
        integer,             intent(in) :: ithr
        procedure(tree_ref_eval)        :: eval_ref
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, best_corr_root
        node_root = block_tree%get_root_node(itree)
        if( node_root%node_idx == 0 ) return
        call eval_ref(node_root%ref_idx, iptcl, ithr, best_corr_root)
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( block_tree%is_leaf(itree, inode) ) exit
            node_cur = block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            call eval_child_best_generic(block_tree, itree, node_cur%left_idx,  iptcl, ithr, eval_ref, best_corr_L)
            call eval_child_best_generic(block_tree, itree, node_cur%right_idx, iptcl, ithr, eval_ref, best_corr_R)
            inode_next = choose_next_child_prob(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
    end subroutine trace_tree_prob

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

    subroutine eval_child_best_generic( block_tree, itree, child_idx, iptcl, ithr, eval_ref, best_corr )
        use simple_binary_tree, only: bt_node
        class(multi_dendro), intent(in)    :: block_tree
        integer,             intent(in)    :: itree
        integer,             intent(in)    :: child_idx
        integer,             intent(in)    :: iptcl
        integer,             intent(in)    :: ithr
        procedure(tree_ref_eval)           :: eval_ref
        real,                intent(out)   :: best_corr
        type(bt_node) :: node_child
        best_corr = INVALID_CORR
        if( child_idx == 0 ) return
        node_child = block_tree%get_node(itree, child_idx)
        if( node_child%ref_idx == 0 ) return
        call eval_ref(node_child%ref_idx, iptcl, ithr, best_corr)
    end subroutine eval_child_best_generic

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
            if( s%p_ptr%l_sh_first )then
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
        if( s%p_ptr%l_sh_first )then
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
        integer, allocatable :: tree_refs(:)
        integer              :: iref_tree
        real                 :: best_corr_ref
        tree_refs = s%b_ptr%block_tree%get_tree_refs(itree)
        do iref_tree = 1, size(tree_refs)
            call eval_tree_ref_across_states(s, tree_refs(iref_tree), best_corr_ref, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_ref)
        end do
        if( allocated(tree_refs) ) deallocate(tree_refs)
    end subroutine exhaustive_tree_scan

    subroutine exhaustive_tree_scan_fixed_state( s, itree, istate_fixed, tree_best_corr, nrefs_tree )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: itree
        integer,                intent(in)    :: istate_fixed
        real,                   intent(inout) :: tree_best_corr
        integer,                intent(inout) :: nrefs_tree
        integer, allocatable :: tree_refs(:)
        integer              :: iref_tree
        real                 :: best_corr_ref
        tree_refs = s%b_ptr%block_tree%get_tree_refs(itree)
        do iref_tree = 1, size(tree_refs)
            call eval_tree_ref_fixed_state(s, tree_refs(iref_tree), istate_fixed, best_corr_ref, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_ref)
        end do
        if( allocated(tree_refs) ) deallocate(tree_refs)
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

    integer function choose_next_child_bestfirst( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
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
    end function choose_next_child_bestfirst

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