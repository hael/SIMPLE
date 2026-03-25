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
        integer, parameter :: MAX_NTREES = 2500, MAX_NPEAKS = 64
        integer :: peak_trees(MAX_NPEAKS), tree_best_states(MAX_NTREES), peak_tree_states(MAX_NPEAKS)
        real    :: tree_best_corrs(MAX_NTREES), peak_tree_corrs(MAX_NPEAKS)
        integer :: npeak_trees, ntrees, nrefs_coarse, nrefs_tree
        integer :: isample, iref, itree, npeak_target, loc(1)
        real    :: inpl_corrs(self%s%nrots), corr_best
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
            tree_best_corrs(1:ntrees) = -huge(1.0)
            tree_best_states(1:ntrees) = 1
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
        ! Phase 2: pick the best trees from the coarse pass, then do a
        ! corr-guided stochastic descent in each tree.
        ! When nstates > 1, select one best peak per state. When nstates == 1,
        ! use the standard selection of top npeaks.
        ! ----------------------------------------------------------------------
        npeak_trees = 0
        if( ntrees > 0 )then
            if( self%s%p_ptr%nstates == 1 )then
                ! Standard selection: top npeaks trees by correlation
                npeak_target = min(self%s%npeaks, ntrees, MAX_NPEAKS)
                if( npeak_target > 0 )then
                    peak_trees(1:npeak_target)      = 0
                    peak_tree_corrs(1:npeak_target) = -huge(1.0)
                    peak_tree_states(1:npeak_target) = 1
                    call select_peak_trees(tree_best_corrs(1:ntrees), peak_trees(1:npeak_target), peak_tree_corrs(1:npeak_target), npeak_trees)
                    if( npeak_trees > 0 )then
                        peak_tree_states(1:npeak_trees) = tree_best_states(peak_trees(1:npeak_trees))
                    endif
                    do itree = 1, npeak_trees
                        call descend_tree_prob_fixed_state(self%s, peak_trees(itree), peak_tree_corrs(itree), nrefs_tree, peak_tree_states(itree))
                    end do
                endif
            else
                ! Multi-state selection: npeaks best peaks per state
                npeak_target = min(self%s%p_ptr%nstates * self%s%npeaks, ntrees, MAX_NPEAKS)
                if( npeak_target > 0 )then
                    peak_trees(1:npeak_target)      = 0
                    peak_tree_corrs(1:npeak_target) = -huge(1.0)
                    peak_tree_states(1:npeak_target) = 1
                    call select_peak_trees_per_state(npeak_trees)
                    do itree = 1, npeak_trees
                        call descend_tree_prob_fixed_state(self%s, peak_trees(itree), peak_tree_corrs(itree), nrefs_tree, peak_tree_states(itree))
                    end do
                endif
            endif
        endif
        self%s%nrefs_eval = nrefs_coarse + nrefs_tree
        call extract_peak_oris(self%s)
        call self%s%inpl_srch_peaks
        call self%oris_assign

    contains

        subroutine select_peak_trees_per_state( npeak_trees_out )
            integer, intent(out) :: npeak_trees_out
            integer :: istate, itree, irank, best_idx
            real    :: state_corrs(MAX_NTREES)
            integer :: state_trees(MAX_NTREES)
            integer :: nstates, npeaks_to_select, npeaks_this_state
            real    :: best_corr
            npeak_trees_out = 0
            nstates = self%s%p_ptr%nstates
            npeaks_to_select = self%s%npeaks
            ! For each state, find the top npeaks trees
            do istate = 1, nstates
                ! Collect all trees for this state
                npeaks_this_state = 0
                do itree = 1, ntrees
                    if( tree_best_states(itree) == istate )then
                        npeaks_this_state = npeaks_this_state + 1
                        state_corrs(npeaks_this_state) = tree_best_corrs(itree)
                        state_trees(npeaks_this_state) = itree
                    endif
                end do
                ! Select top npeaks_to_select for this state
                do irank = 1, min(npeaks_to_select, npeaks_this_state)
                    ! Find the best remaining tree for this state
                    best_corr = -huge(1.0)
                    best_idx = 0
                    do itree = 1, npeaks_this_state
                        if( state_corrs(itree) > best_corr )then
                            best_corr = state_corrs(itree)
                            best_idx = itree
                        endif
                    end do
                    ! Add to results if space available
                    if( best_idx > 0 .and. npeak_trees_out < npeak_target )then
                        npeak_trees_out = npeak_trees_out + 1
                        peak_trees(npeak_trees_out) = state_trees(best_idx)
                        peak_tree_corrs(npeak_trees_out) = state_corrs(best_idx)
                        peak_tree_states(npeak_trees_out) = istate
                        ! Mark as used
                        state_corrs(best_idx) = -huge(1.0)
                    else
                        exit  ! Stop if we run out of space
                    endif
                end do
            end do
        end subroutine select_peak_trees_per_state

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