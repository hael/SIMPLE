!@descr: 3D strategy for coarse subspace search followed by probabilistic block-tree descent
!        Phase 1 mirrors the coarse pass of strategy3D_greedy_sub and ranks unique trees by their
!        best coarse score. Phase 2 evaluates the tree root, descends stochastically using the same
!        in-plane policy as the rest of the strategy, and falls back to exhaustive tree scoring when
!        the stochastic walk fails to improve on the coarse tree score.
module simple_strategy3D_ptree
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: parameters
use simple_oris,            only: oris
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
implicit none

public :: strategy3D_ptree
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_ptree
contains
    procedure :: new         => new_ptree
    procedure :: srch        => srch_ptree
    procedure :: kill        => kill_ptree
    procedure :: oris_assign => oris_assign_ptree
end type strategy3D_ptree

contains

    subroutine new_ptree( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_ptree), intent(inout) :: self
        class(parameters),       intent(in)    :: params
        class(strategy3D_spec),  intent(inout) :: spec
        class(builder),          intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_ptree

    subroutine srch_ptree( self, os, ithr )
        class(strategy3D_ptree), intent(inout) :: self
        class(oris),             intent(inout) :: os
        integer,                 intent(in)    :: ithr
        integer, allocatable :: peak_trees(:)
	    real,    allocatable :: tree_best_corrs(:), peak_tree_corrs(:)
        integer              :: npeak_trees, ntrees, nrefs_coarse, nrefs_tree
        integer              :: isample, iref, itree, npeak_target, loc(1)
        real                 :: inpl_corrs(self%s%nrots), corr_best
        if( os%get_state(self%s%iptcl) <= 0 )then
            call os%reject(self%s%iptcl)
            return
        endif
        self%s%ithr = ithr
        call self%s%prep4srch
        call self%s%inpl_srch_first
        nrefs_coarse = 0
        nrefs_tree   = 0
        ntrees       = self%s%b_ptr%block_tree%get_n_trees()
        if( .not. allocated(self%s%b_ptr%subspace_full2sub_map) )then
            THROW_HARD('Probabilistic tree search requires subspace_full2sub_map. Check builder construction.')
        endif
        if( ntrees > 0 )then
            allocate(tree_best_corrs(ntrees), source=-huge(1.0))
        endif
        ! ----------------------------------------------------------------------
        ! Phase 1: coarse subspace search, matching strategy3D_greedy_sub.
        ! Keep the best coarse score for each tree so phase 2 can search unique
        ! trees without losing coverage when multiple coarse peaks fall in the
        ! same tree.
        ! ----------------------------------------------------------------------
        do isample = 1, self%s%nrefs_sub
            iref = s3D%srch_order_sub(isample, self%s%ithr)
            if( .not. s3D%state_exists(s3D%proj_space_state(iref)) ) cycle
            if( self%s%p_ptr%l_sh_first )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
            endif
            loc       = maxloc(inpl_corrs)
            corr_best = inpl_corrs(loc(1))
            call self%s%store_solution(iref, loc(1), corr_best)
            nrefs_coarse = nrefs_coarse + 1
            if( ntrees > 0 )then
                itree = get_tree_for_ref(self, iref, ntrees)
                tree_best_corrs(itree) = max(tree_best_corrs(itree), corr_best)
            endif
        end do
        ! ----------------------------------------------------------------------
        ! Phase 2: pick the best unique trees from the coarse pass, then do a
        ! corr-guided stochastic descent in each tree. The root is evaluated
        ! explicitly, leaf-only trees are handled naturally, and a failed random
        ! walk falls back to exhaustive tree scoring.
        ! ----------------------------------------------------------------------
        npeak_trees = 0
        if( ntrees > 0 )then
            npeak_target = min(self%s%npeaks, ntrees)
            if( npeak_target > 0 )then
                allocate(peak_trees(npeak_target),      source=0)
                allocate(peak_tree_corrs(npeak_target), source=-huge(1.0))
                call select_peak_trees(tree_best_corrs, peak_trees, peak_tree_corrs, npeak_trees)
                do itree = 1, npeak_trees
                    call descend_tree_prob(self, peak_trees(itree), peak_tree_corrs(itree), nrefs_tree)
                end do
            endif
        endif
        self%s%nrefs_eval = nrefs_coarse + nrefs_tree
        call extract_peak_oris(self%s)
        call self%s%inpl_srch_peaks
        call self%oris_assign
    end subroutine srch_ptree

    subroutine select_peak_trees( tree_best_corrs, peak_trees, peak_tree_corrs, npeak_trees )
        real,    intent(in)  :: tree_best_corrs(:)
        integer, intent(out) :: peak_trees(:)
        real,    intent(out) :: peak_tree_corrs(:)
        integer, intent(out) :: npeak_trees
        integer :: loc(1)
        real    :: work_corrs(size(tree_best_corrs))
        peak_trees      = 0
        peak_tree_corrs = -huge(1.0)
        npeak_trees     = 0
        if( size(tree_best_corrs) == 0 ) return
        if( size(peak_trees)      == 0 ) return
        work_corrs = tree_best_corrs
        do while( npeak_trees < size(peak_trees) )
            loc = maxloc(work_corrs)
            if( work_corrs(loc(1)) <= -huge(1.0) / 2.0 ) exit
            npeak_trees                  = npeak_trees + 1
            peak_trees(npeak_trees)      = loc(1)
            peak_tree_corrs(npeak_trees) = work_corrs(loc(1))
            work_corrs(loc(1))           = -huge(1.0)
        end do
    end subroutine select_peak_trees

    subroutine descend_tree_prob( self, itree, coarse_tree_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_ptree), intent(inout) :: self
        integer,                 intent(in)    :: itree
        real,                    intent(in)    :: coarse_tree_corr
        integer,                 intent(inout) :: nrefs_tree
        type(bt_node) :: node_cur, node_root
        integer       :: inode, inode_next
        real          :: best_corr_L, best_corr_R, tree_best_corr
        logical       :: did_descend
        node_root = self%s%b_ptr%block_tree%get_root_node(itree)
        if( node_root%node_idx == 0 )then
            THROW_HARD('Selected block tree has an invalid root node.')
        endif
        if( node_root%ref_idx == 0 )then
            THROW_HARD('Selected block tree root has no ref_idx.')
        endif
        call eval_tree_ref_across_states(self, node_root%ref_idx, tree_best_corr, nrefs_tree)
        did_descend = .false.
        inode = node_root%node_idx
        do
            if( inode == 0 ) exit
            if( self%s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
            node_cur = self%s%b_ptr%block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            did_descend = .true.
            call eval_child_best(self, itree, node_cur%left_idx,  best_corr_L, nrefs_tree)
            call eval_child_best(self, itree, node_cur%right_idx, best_corr_R, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_L, best_corr_R)
            inode_next = choose_next_child(node_cur%left_idx, node_cur%right_idx, best_corr_L, best_corr_R)
            if( inode_next == 0 ) exit
            inode = inode_next
        end do
        if( did_descend .and. tree_best_corr <= coarse_tree_corr )then
            call exhaustive_tree_scan(self, itree, tree_best_corr, nrefs_tree)
        endif
    end subroutine descend_tree_prob

    subroutine eval_child_best( self, itree, child_idx, best_corr, nrefs_tree )
        use simple_binary_tree, only: bt_node
        class(strategy3D_ptree), intent(inout) :: self
        integer,                 intent(in)    :: itree
        integer,                 intent(in)    :: child_idx
        real,                    intent(out)   :: best_corr
        integer,                 intent(inout) :: nrefs_tree
        type(bt_node) :: node_child
        best_corr = -huge(1.0)
        if( child_idx == 0 ) return
        node_child = self%s%b_ptr%block_tree%get_node(itree, child_idx)
        if( node_child%ref_idx == 0 ) return
        call eval_tree_ref_across_states(self, node_child%ref_idx, best_corr, nrefs_tree)
    end subroutine eval_child_best

    subroutine eval_tree_ref_across_states( self, ref_idx, best_corr, nrefs_tree )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_ptree), intent(inout) :: self
        integer,                 intent(in)    :: ref_idx
        real,                    intent(out)   :: best_corr
        integer,                 intent(inout) :: nrefs_tree
        integer :: inds(self%s%nrots)
        integer :: iref
        integer :: istate
        integer :: loc(1)
        real    :: corr_tmp,  inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots) 
        best_corr = -huge(1.0)
        do istate = 1, self%s%nstates
            iref = (istate - 1) * self%s%p_ptr%nspace + ref_idx
            if( .not. s3D%state_exists(s3D%proj_space_state(iref)) ) cycle
            if( self%s%p_ptr%l_sh_first )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
            endif
            if( self%s%p_ptr%l_prob_inpl )then
                loc = angle_sampling( &
                    eulprob_dist_switch(inpl_corrs, self%s%p_ptr%cc_objfun), &
                    sorted_corrs, inds, &
                    s3D%smpl_inpl_athres(s3D%proj_space_state(iref)), &
                    self%s%p_ptr%prob_athres )
            else
                loc = maxloc(inpl_corrs)
            endif
            corr_tmp = inpl_corrs(loc(1))
            call self%s%store_solution(iref, loc(1), corr_tmp)
            nrefs_tree = nrefs_tree + 1
            best_corr  = max(best_corr, corr_tmp)
        end do
    end subroutine eval_tree_ref_across_states

    subroutine exhaustive_tree_scan( self, itree, tree_best_corr, nrefs_tree )
        class(strategy3D_ptree), intent(inout) :: self
        integer,                 intent(in)    :: itree
        real,                    intent(inout) :: tree_best_corr
        integer,                 intent(inout) :: nrefs_tree
        integer, allocatable :: tree_refs(:)
        integer              :: iref_tree
        real                 :: best_corr_ref
        tree_refs = self%s%b_ptr%block_tree%get_tree_refs(itree)
        do iref_tree = 1, size(tree_refs)
            call eval_tree_ref_across_states(self, tree_refs(iref_tree), best_corr_ref, nrefs_tree)
            tree_best_corr = max(tree_best_corr, best_corr_ref)
        end do
        if( allocated(tree_refs) ) deallocate(tree_refs)
    end subroutine exhaustive_tree_scan

    integer function get_tree_for_ref( self, iref, ntrees ) result(itree)
        class(strategy3D_ptree), intent(in) :: self
        integer,                 intent(in) :: iref, ntrees
        integer :: iproj
        iproj = s3D%proj_space_proj(iref)
        if( iproj < 1 .or. iproj > size(self%s%b_ptr%subspace_full2sub_map) )then
            THROW_HARD('Invalid projection index encountered during block-tree search.')
        endif
        itree = self%s%b_ptr%subspace_full2sub_map(iproj)
        if( itree < 1 .or. itree > ntrees )then
            THROW_HARD('Invalid tree index mapped from reference. Check builder construction.')
        endif
    end function get_tree_for_ref

    integer function choose_next_child( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
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
        if( corr_left <= -huge(1.0) / 2.0 .and. corr_right <= -huge(1.0) / 2.0 )then
            if( sample_two(1.0, 1.0) == 1 )then
                inode_next = left_idx
            else
                inode_next = right_idx
            endif
            return
        endif
        cmax = max(corr_left, corr_right)
        if( corr_left <= -huge(1.0) / 2.0 )then
            p_left = 0.0
        else
            p_left = exp(corr_left - cmax)
        endif
        if( corr_right <= -huge(1.0) / 2.0 )then
            p_right = 0.0
        else
            p_right = exp(corr_right - cmax)
        endif

        if( sample_two(p_left, p_right) == 1 )then
            inode_next = left_idx
        else
            inode_next = right_idx
        endif
    end function choose_next_child

    integer function sample_two( p1, p2 ) result(which)
        real, intent(in) :: p1
        real, intent(in) :: p2
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

    subroutine oris_assign_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_ptree

    subroutine kill_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_ptree

end module simple_strategy3D_ptree