!@descr: 3D strategy for coarse subspace search followed by probabilistic block-tree descent
!        Phase 1 mirrors the coarse pass of strategy3D_greedy_sub and ranks unique trees by their
!        best coarse score. Phase 2 evaluates the tree root, descends stochastically using the same
!        in-plane policy as the rest of the strategy.
module simple_strategy3D_ptree
use simple_core_module_api
use simple_strategy3D_alloc,      only: s3D
use simple_strategy3D_utils,      only: extract_peak_ori
use simple_strategy3D_tree_utils, only: init_peak_tree_selection, peak_tree_selection, select_peak_trees, select_peak_trees_per_state, &
    &descend_tree_prob_fixed_state, get_tree_for_ref
use simple_strategy_tree_helpers, only: MAX_NTREES, MAX_NPEAKS
use simple_parameters,            only: parameters
use simple_oris,                  only: oris
use simple_strategy3D,            only: strategy3D
use simple_strategy3D_srch,       only: strategy3D_spec
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
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_ptree), intent(inout) :: self
        class(oris),             intent(inout) :: os
        integer,                 intent(in)    :: ithr
        type(peak_tree_selection) :: peak_sel
        integer :: npeak_trees, ntrees, nrefs_coarse, nrefs_tree
        integer :: isample, iref, itree, npeak_target, loc(1), inds(self%s%nrots)
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots), corr_best
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
        if( ntrees > MAX_NTREES )then
            THROW_HARD('Number of trees exceeds MAX_NTREES; srch_ptree')
        endif
        if( ntrees > 0 )then
            npeak_target = min(self%s%p_ptr%nstates * self%s%npeaks, ntrees, MAX_NPEAKS)
            call init_peak_tree_selection(peak_sel, ntrees, npeak_target, .true., .true., self%s%p_ptr%nstates, self%s%npeaks)
        else
            THROW_HARD('Probabilistic tree search requires at least one block tree.')
        endif
        ! ----------------------------------------------------------------------
        ! Phase 1: coarse subspace search
        ! ----------------------------------------------------------------------
        do isample = 1, self%s%nrefs_sub
            iref = s3D%srch_order_sub(isample, self%s%ithr)
            if( .not. s3D%state_exists(s3D%proj_space_state(iref)) ) cycle
            if( self%s%p_ptr%l_doshift )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
            endif
            loc = angle_sampling( eulprob_dist_switch(inpl_corrs, self%s%p_ptr%cc_objfun), &
                &sorted_corrs, inds, s3D%smpl_inpl_athres(s3D%proj_space_state(iref)), self%s%p_ptr%prob_athres )
            corr_best = inpl_corrs(loc(1))
            call self%s%store_solution(iref, loc(1), corr_best)
            nrefs_coarse = nrefs_coarse + 1
            if( ntrees > 0 )then
                itree = get_tree_for_ref(self%s, iref, ntrees)
                if( corr_best > peak_sel%tree_best_corrs(itree) )then
                    peak_sel%tree_best_corrs(itree)  = corr_best
                    peak_sel%tree_best_states(itree) = s3D%proj_space_state(iref)
                endif
            endif
        end do
        ! ----------------------------------------------------------------------
        ! Phase 2: pick the best trees from the coarse pass, then do a
        ! corr-guided stochastic descent in each tree.
        ! When nstates > 1, select one best peak per state. When nstates == 1,
        ! use the standard selection of top npeaks.
        ! ----------------------------------------------------------------------
        npeak_trees = 0
        if( ntrees > 0 )then
            if( self%s%p_ptr%nstates == 1 )then
                ! Standard selection: top npeaks trees by correlation
                if( npeak_target > 0 )then
                    call select_peak_trees(peak_sel)
                    npeak_trees = peak_sel%npeak_trees
                    if( npeak_trees > 0 )then
                        peak_sel%peak_tree_states(1:npeak_trees) = peak_sel%tree_best_states(peak_sel%peak_trees(1:npeak_trees))
                    endif
                    do itree = 1, npeak_trees
                        call descend_tree_prob_fixed_state(self%s, peak_sel%peak_trees(itree), nrefs_tree, peak_sel%peak_tree_states(itree))
                    end do
                endif
            else
                ! Multi-state selection: npeaks best peaks per state
                if( npeak_target > 0 )then
                    call select_peak_trees_per_state(peak_sel)
                    npeak_trees = peak_sel%npeak_trees
                    do itree = 1, npeak_trees
                        call descend_tree_prob_fixed_state(self%s, peak_sel%peak_trees(itree), nrefs_tree, peak_sel%peak_tree_states(itree))
                    end do
                endif
            endif
        endif
        self%s%nrefs_eval = nrefs_coarse + nrefs_tree
        call self%s%inpl_srch_peaks(min(self%s%npeaks_inpl, self%s%nsolns))
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