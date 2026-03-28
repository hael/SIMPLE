!@descr: 2D single-tree strategy: repeated probabilistic descents from the root of
!         the single block tree, followed by shift refinement of top candidates.
!         Assumes block_tree contains exactly one tree (itree = 1). The number of
!         descent trials is controlled by the N_DESCENTS module parameter; each
!         trial independently descends the tree probabilistically, accumulating
!         candidate solutions in the shared class-space arrays. The top candidates
!         are then shift-refined before a final power-sampled class assignment.
module simple_strategy2D_single_ptree
use simple_pftc_srch_api
use simple_strategy2D_alloc
use simple_strategy2D_tree_utils, only: descend_tree_prob
use simple_strategy2D,            only: strategy2D
use simple_strategy2D_srch,       only: strategy2D_spec
use simple_builder,               only: builder
implicit none

public :: strategy2D_single_ptree
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG      = .false.
integer, parameter :: N_DESCENTS = 5 !< number of repeated tree descent trials

type, extends(strategy2D) :: strategy2D_single_ptree
contains
    procedure :: new  => new_single_ptree
    procedure :: srch => srch_single_ptree
    procedure :: kill => kill_single_ptree
end type strategy2D_single_ptree

contains

    subroutine new_single_ptree( self, params, spec, build )
        class(strategy2D_single_ptree), intent(inout) :: self
        class(parameters),              intent(in)    :: params
        class(strategy2D_spec),         intent(inout) :: spec
        class(builder),                 intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_single_ptree

    subroutine srch_single_ptree( self, os )
        class(strategy2D_single_ptree), intent(inout) :: self
        class(oris),                    intent(inout) :: os
        integer :: sorted_cls_inds(self%s%nrefs)
        integer :: itrial, nrefs_tree, class_rank
        if( os%get_state(self%s%iptcl) <= 0 )then
            call os%reject(self%s%iptcl)
            return
        endif
        call self%s%prep4srch(os)
        call self%s%inpl_srch_first
        if( self%s%b_ptr%block_tree%get_n_trees() < 1 )then
            THROW_HARD('single_ptree search requires at least one tree in block_tree.')
        endif
        ! Repeated probabilistic descents from the root of tree 1.
        ! store_solution guards against non-improving overwrites, so each trial
        ! can only improve the per-class candidate corrs in the class-space arrays.
        nrefs_tree = 0
        do itrial = 1, N_DESCENTS
            call descend_tree_prob(self%s, 1, nrefs_tree)
        end do
        self%s%nrefs_eval = nrefs_tree
        ! Shift refinement for top scoring candidates
        call self%s%inpl_srch_peaks(min(s2D%snhc_smpl_ncls, self%s%nsolns))
        ! Class selection via power_sampling over shift-refined scores
        call power_sampling( s2D%power, self%s%nrefs, s2D%class_space_corrs(:, self%s%ithr), &
                            &sorted_cls_inds, s2D%snhc_smpl_ncls, &
                            &self%s%best_class, class_rank, self%s%best_corr )
        self%s%best_rot = 0
        if( self%s%best_class > 0 )then
            self%s%best_rot = s2D%class_space_inplinds(self%s%best_class, self%s%ithr)
        endif
        if( self%s%best_class <= 0 .or. self%s%best_rot <= 0 )then
            self%s%best_class = self%s%prev_class
            self%s%best_rot   = max(1, self%s%prev_rot)
            self%s%best_corr  = self%s%prev_corr
        endif
        ! Final in-plane search on the selected candidate
        call self%s%inpl_srch ! needed because inpl_srch_peaks doesn't store shifts
        call self%s%store_solution(self%s%best_class, self%s%best_rot, self%s%best_corr)
        call self%s%assign_ori(os)
        if( DEBUG ) write(logfhandle,*) '>>> strategy2D_single_ptree::FINISHED SINGLE PTREE SEARCH'
    end subroutine srch_single_ptree

    subroutine kill_single_ptree( self )
        class(strategy2D_single_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_single_ptree

end module simple_strategy2D_single_ptree
