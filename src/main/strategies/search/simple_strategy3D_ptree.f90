!@descr: 3D strategy for coarse subspace search followed by probabilistic block-tree descent
!        Phase 1 mirrors the coarse pass of strategy3D_greedy_sub and ranks unique trees by their
!        best coarse score. Phase 2 evaluates the tree root, descends stochastically using the same
!        in-plane policy as the rest of the strategy, and falls back to exhaustive tree scoring when
!        the stochastic walk fails to improve on the coarse tree score.
module simple_strategy3D_ptree
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_tree_utils,&
                           &only: select_peak_trees, descend_tree_prob_fixed_state, get_tree_for_ref
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
        integer, allocatable :: peak_trees(:), tree_best_states(:), peak_tree_states(:)
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
            allocate(tree_best_states(ntrees), source=1)
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
            if( self%s%p_ptr%l_doshift )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
            endif
            loc       = maxloc(inpl_corrs)
            corr_best = inpl_corrs(loc(1))
            call self%s%store_solution(iref, loc(1), corr_best)
            nrefs_coarse = nrefs_coarse + 1
            if( ntrees > 0 )then
                itree = get_tree_for_ref(self%s, iref, ntrees)
                if( corr_best > tree_best_corrs(itree) )then
                    tree_best_corrs(itree)  = corr_best
                    tree_best_states(itree) = s3D%proj_space_state(iref)
                endif
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
                allocate(peak_tree_states(npeak_target), source=1)
                call select_peak_trees(tree_best_corrs, peak_trees, peak_tree_corrs, npeak_trees)
                peak_tree_states(1:npeak_trees) = tree_best_states(peak_trees(1:npeak_trees))
                do itree = 1, npeak_trees
                    call descend_tree_prob_fixed_state(self%s, peak_trees(itree), peak_tree_corrs(itree), nrefs_tree, peak_tree_states(itree))
                end do
            endif
        endif
        self%s%nrefs_eval = nrefs_coarse + nrefs_tree
        call extract_peak_oris(self%s)
        call self%s%inpl_srch_peaks
        call self%oris_assign
    end subroutine srch_ptree

    subroutine oris_assign_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_ptree

    subroutine kill_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_ptree

end module simple_strategy3D_ptree