!@descr: 3D strategy: neighborhood selection by angular distance from previous orientation,
!        followed by probabilistic tree descent within the selected trees.
!        The neighborhood is defined by computing symmetry-aware angular distances from
!        the previous assigned orientation to every subspace (coarse-grid) representative,
!        then selecting the npeaks nearest representatives as tree roots. Each selected
!        tree is then descended probabilistically, with exhaustive fallback when the
!        stochastic walk fails to improve on the previous particle score. States are
!        optimised jointly within each tree via eval_tree_ref_across_states.
!        This mirrors the geometry-based pooling in build_neigh_mask_from_prev_geom,
!        but operates at the tree-search level rather than building a sparse reference
!        graph.
module simple_strategy3D_ptree_neigh
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_tree_utils, &
                           &only: select_peak_trees, descend_tree_prob
use simple_strategy3D_utils
use simple_parameters,      only: parameters
use simple_oris,            only: oris
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
implicit none

public :: strategy3D_ptree_neigh
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_ptree_neigh
contains
    procedure :: new         => new_ptree_neigh
    procedure :: srch        => srch_ptree_neigh
    procedure :: kill        => kill_ptree_neigh
    procedure :: oris_assign => oris_assign_ptree_neigh
end type strategy3D_ptree_neigh

contains

    subroutine new_ptree_neigh( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_ptree_neigh), intent(inout) :: self
        class(parameters),             intent(in)    :: params
        class(strategy3D_spec),        intent(inout) :: spec
        class(builder),                intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_ptree_neigh

    subroutine srch_ptree_neigh( self, os, ithr )
        class(strategy3D_ptree_neigh), intent(inout) :: self
        class(oris),                   intent(inout) :: os
        integer,                       intent(in)    :: ithr
        integer, allocatable :: peak_trees(:)
        real,    allocatable :: neg_dists(:), peak_tree_corrs(:)
        integer :: ntrees, npeak_target, npeak_trees
        integer :: nrefs_tree, isub, iproj_full, itree, i
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
            THROW_HARD('ptree_neigh search requires subspace_full2sub_map. Check builder construction.')
        endif
        if( .not. allocated(self%s%b_ptr%subspace_inds) )then
            THROW_HARD('ptree_neigh search requires subspace_inds. Check builder construction.')
        endif
        if( ntrees <= 0 )then
            THROW_HARD('ptree_neigh search requires at least one block tree.')
        endif
        ! ------------------------------------------------------------------
        ! Geometry-based neighborhood: compute symmetry-aware angular
        ! distance from the previous orientation to every subspace
        ! representative. Store the minimum distance per tree (multiple
        ! subspace reps can map to the same tree in degenerate build cases).
        ! ------------------------------------------------------------------
        npeak_target = min(self%s%npeaks, ntrees)
        allocate(neg_dists(ntrees),             source=-huge(1.0))
        allocate(peak_trees(npeak_target),      source=0)
        allocate(peak_tree_corrs(npeak_target), source=-huge(1.0))
        do isub = 1, self%s%p_ptr%nspace_sub
            iproj_full = self%s%b_ptr%subspace_inds(isub)
            if( iproj_full < 1 .or. iproj_full > self%s%p_ptr%nspace ) THROW_HARD('subspace index out of bound; ptree_neigh search')
            call self%s%b_ptr%eulspace%get_ori(iproj_full, o_sub)
            call self%s%b_ptr%pgrpsyms%sym_dists(self%s%o_prev, o_sub, osym, dtmp, inplrotdist)
            call o_sub%kill
            call osym%kill
            itree = self%s%b_ptr%subspace_full2sub_map(iproj_full)
            if( itree < 1 .or. itree > ntrees ) THROW_HARD('tree index out of bound; ptree_neigh search')
            ! negate: select_peak_trees maximises, we want minimum distance
            neg_dists(itree) = max(neg_dists(itree), -dtmp)
        enddo
        ! ------------------------------------------------------------------
        ! Select the npeak_target trees whose representatives lie closest
        ! (in angular distance) to the previous orientation.
        ! ------------------------------------------------------------------
        call select_peak_trees(neg_dists, peak_trees, peak_tree_corrs, npeak_trees)
        ! ------------------------------------------------------------------
        ! Probabilistic tree descent in each selected tree.
        ! prev_corr serves as the quality bound: if the stochastic walk
        ! fails to improve on the current best assignment, the selected
        ! tree is scanned exhaustively.
        ! ------------------------------------------------------------------
        do i = 1, npeak_trees
            ! this descent consideres all states
            call descend_tree_prob(self%s, peak_trees(i), self%s%prev_corr, nrefs_tree)
        enddo
        self%s%nrefs_eval = nrefs_tree
        call extract_peak_oris(self%s)
        call self%s%inpl_srch_peaks
        call self%oris_assign
    end subroutine srch_ptree_neigh

    subroutine oris_assign_ptree_neigh( self )
        class(strategy3D_ptree_neigh), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_ptree_neigh

    subroutine kill_ptree_neigh( self )
        class(strategy3D_ptree_neigh), intent(inout) :: self
        call self%s%kill
    end subroutine kill_ptree_neigh

end module simple_strategy3D_ptree_neigh
