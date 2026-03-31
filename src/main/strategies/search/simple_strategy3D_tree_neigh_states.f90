!@descr: probabilistic single-tree descent with geometry-based tree selection and
!         multi-state evaluation at each node.
! 3D strategy: select one tree from the previous orientation geometry, then descend
! probabilistically while evaluating all states at every visited node. This is aimed
! at repartitioning particles after an average-volume refinement / random state split:
! the geometry prior comes from the previous assigned orientation, but the state is
! allowed to change freely during the descent.
module simple_strategy3D_tree_neigh_states
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_strategy3D_tree_utils
use simple_strategy_tree_helpers, only: MAX_NTREES
use simple_parameters,            only: parameters
use simple_oris,                  only: oris
use simple_strategy3D,            only: strategy3D
use simple_strategy3D_srch,       only: strategy3D_spec
implicit none

public :: strategy3D_tree_neigh_states
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_tree_neigh_states
contains
    procedure :: new         => new_tree_neigh_states
    procedure :: srch        => srch_tree_neigh_states
    procedure :: kill        => kill_tree_neigh_states
    procedure :: oris_assign => oris_assign_tree_neigh_states
end type strategy3D_tree_neigh_states

contains

    subroutine new_tree_neigh_states( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_tree_neigh_states), intent(inout) :: self
        class(parameters),                    intent(in)    :: params
        class(strategy3D_spec),               intent(inout) :: spec
        class(builder),                       intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_tree_neigh_states

    subroutine srch_tree_neigh_states( self, os, ithr )
        class(strategy3D_tree_neigh_states), intent(inout) :: self
        class(oris),                          intent(inout) :: os
        integer,                              intent(in)    :: ithr
        type(peak_tree_selection) :: peak_sel
        integer :: ntrees, npeak_target, nrefs_tree
        integer :: isub, iproj_full, itree_selected, itree
        real    :: dtmp, inplrotdist
        type(ori) :: o_sub, osym
        if( os%get_state(self%s%iptcl) <= 0 )then
            call os%reject(self%s%iptcl)
            return
        endif
        self%s%ithr = ithr
        call self%s%prep4srch
        call self%s%inpl_srch_first
        nrefs_tree = 0
        ntrees     = self%s%b_ptr%block_tree%get_n_trees()
        if( .not. allocated(self%s%b_ptr%subspace_full2sub_map) )then
            THROW_HARD('tree_neigh_states search requires subspace_full2sub_map. Check builder construction.')
        endif
        if( .not. allocated(self%s%b_ptr%subspace_inds) )then
            THROW_HARD('tree_neigh_states search requires subspace_inds. Check builder construction.')
        endif
        if( ntrees <= 0 )then
            THROW_HARD('tree_neigh_states search requires at least one block tree.')
        endif
        if( ntrees > MAX_NTREES )then
            THROW_HARD('Number of trees exceeds MAX_NTREES; srch_tree_neigh_states')
        endif
        ! ------------------------------------------------------------------
        ! Geometry-only tree choice.
        ! Compute the symmetry-aware angular distance from the previous
        ! orientation to every subspace representative and keep, for each tree,
        ! the representative with the smallest distance. Then select exactly one
        ! tree: the geometry-nearest one.
        ! ------------------------------------------------------------------
        npeak_target = 1
        call init_peak_tree_selection(peak_sel, ntrees, npeak_target, .false., .false.)
        do isub = 1, self%s%p_ptr%nspace_sub
            iproj_full = self%s%b_ptr%subspace_inds(isub)
            if( iproj_full < 1 .or. iproj_full > self%s%p_ptr%nspace )then
                THROW_HARD('subspace index out of bound; tree_neigh_states search')
            endif
            call self%s%b_ptr%eulspace%get_ori(iproj_full, o_sub)
            call self%s%b_ptr%pgrpsyms%sym_dists(self%s%o_prev, o_sub, osym, dtmp, inplrotdist)
            call o_sub%kill
            call osym%kill
            itree = self%s%b_ptr%subspace_full2sub_map(iproj_full)
            if( itree < 1 .or. itree > ntrees )then
                THROW_HARD('tree index out of bound; tree_neigh_states search')
            endif
            ! select_peak_trees maximises, so negate the distance.
            peak_sel%tree_best_corrs(itree) = max(peak_sel%tree_best_corrs(itree), -dtmp)
        enddo
        call select_peak_trees(peak_sel)
        if( peak_sel%npeak_trees <= 0 )then
            THROW_HARD('No geometry-selected tree found; srch_tree_neigh_states')
        endif
        itree_selected = peak_sel%peak_trees(1)

        ! ------------------------------------------------------------------
        ! Single-tree greedy descent with multi-state evaluation.
        ! The selected tree is descended without fixing the state, so the
        ! particle can move between classes while the orientation is refined.
        ! ------------------------------------------------------------------
        call descend_tree_greedy(self%s, itree_selected, nrefs_tree)

        self%s%nrefs_eval = nrefs_tree
        call self%s%inpl_srch_peaks(min(self%s%npeaks_inpl, self%s%nsolns))
        call self%oris_assign
    end subroutine srch_tree_neigh_states

    subroutine oris_assign_tree_neigh_states( self )
        class(strategy3D_tree_neigh_states), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_tree_neigh_states

    subroutine kill_tree_neigh_states( self )
        class(strategy3D_tree_neigh_states), intent(inout) :: self
        call self%s%kill
    end subroutine kill_tree_neigh_states

end module simple_strategy3D_tree_neigh_states
