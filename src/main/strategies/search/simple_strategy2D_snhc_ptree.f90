!@descr: 2D hybrid strategy: SNHC coarse class selection followed by probabilistic tree descent
module simple_strategy2D_snhc_ptree
use simple_pftc_srch_api
use simple_strategy2D_alloc
use simple_strategy2D_tree_utils, only: descend_tree_prob, get_tree_for_ref
use simple_strategy2D,            only: strategy2D
use simple_strategy2D_srch,       only: strategy2D_spec
use simple_builder,               only: builder
implicit none

public :: strategy2D_snhc_ptree
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG = .false.

type, extends(strategy2D) :: strategy2D_snhc_ptree
contains
    procedure :: new  => new_snhc_ptree
    procedure :: srch => srch_snhc_ptree
    procedure :: kill => kill_snhc_ptree
end type strategy2D_snhc_ptree

contains

    subroutine new_snhc_ptree( self, params, spec, build )
        class(strategy2D_snhc_ptree), intent(inout) :: self
        class(parameters),             intent(in)    :: params
        class(strategy2D_spec),        intent(inout) :: spec
        class(builder),                intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_snhc_ptree

    subroutine srch_snhc_ptree( self, os )
        class(strategy2D_snhc_ptree), intent(inout) :: self
        class(oris),                  intent(inout) :: os
        real    :: corrs(self%s%nrots), inpl_corr, corr_best_coarse
        real    :: cls_corrs(self%s%nrefs)
        integer :: cls_inpl_inds(self%s%nrefs), sorted_cls_inds(self%s%nrefs)
        integer :: vec_nrots(self%s%nrots)
        integer :: iref, isample, inpl_ind, itree, ntrees, order_ind
        integer :: nrefs_coarse, nrefs_tree, iref_best_coarse
        if( os%get_state(self%s%iptcl) <= 0 )then
            call os%reject(self%s%iptcl)
            return
        endif
        call self%s%prep4srch(os)
        call self%s%inpl_srch_first
        if( .not. allocated(self%s%b_ptr%subspace_full2sub_map) )then
            THROW_HARD('snhc_ptree search requires subspace_full2sub_map. Check builder construction.')
        endif
        ntrees = self%s%b_ptr%block_tree%get_n_trees()
        if( ntrees <= 0 )then
            THROW_HARD('snhc_ptree search requires at least one block tree.')
        endif
        cls_corrs        = -1.
        cls_inpl_inds    = 0
        nrefs_coarse     = 0
        nrefs_tree       = 0
        corr_best_coarse = -huge(1.0)
        iref_best_coarse = self%s%prev_class
        do isample = 1, self%s%nrefs
            iref = s2D%srch_order(self%s%iptcl_batch, isample)
            self%s%nrefs_eval = self%s%nrefs_eval + 1
            if( self%s%nrefs_eval > s2D%snhc_nrefs_bound ) exit
            if( s2D%cls_pops(iref) == 0 ) cycle
            if( self%s%l_sh_first )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         corrs)
            endif
            call power_sampling( s2D%power, self%s%nrots, corrs, vec_nrots, &
                                &s2D%snhc_smpl_ninpl, inpl_ind, order_ind, inpl_corr )
            cls_corrs(iref)     = inpl_corr
            cls_inpl_inds(iref) = inpl_ind
            nrefs_coarse = nrefs_coarse + 1
            if( inpl_corr >= corr_best_coarse )then
                corr_best_coarse  = inpl_corr
                iref_best_coarse  = iref
                self%s%best_class = iref
                self%s%best_corr  = inpl_corr
                self%s%best_rot   = inpl_ind
            endif
        end do
        if( nrefs_coarse > 0 )then
            itree = get_tree_for_ref(self%s, iref_best_coarse, ntrees)
            if( itree > 0 ) call descend_tree_prob(self%s, itree, nrefs_tree, cls_corrs, cls_inpl_inds)
        endif
        ! Shift refinement for top scoring candidates (coarse + tree local extrema)
        call self%s%inpl_srch_peaks(s2D%snhc_smpl_ncls, cls_corrs, cls_inpl_inds)
        ! Class selection via power_sampling over shift-refined scores
        call power_sampling( s2D%power, self%s%nrefs, cls_corrs, sorted_cls_inds, s2D%snhc_smpl_ncls, &
                            &self%s%best_class, self%s%nrefs_eval, self%s%best_corr )
        if( cls_inpl_inds(self%s%best_class) > 0 )then
            self%s%best_rot = cls_inpl_inds(self%s%best_class)
        endif
        if( self%s%best_rot <= 0 )then
            self%s%best_class = self%s%prev_class
            self%s%best_rot   = max(1, self%s%prev_rot)
            self%s%best_corr  = self%s%prev_corr
        endif
        ! Final in-plane search
        call self%s%inpl_srch
        call self%s%store_solution(os, nrefs=min(self%s%nrefs, s2D%snhc_nrefs_bound+1))
        if( DEBUG ) write(logfhandle,*) '>>> strategy2D_snhc_ptree::FINISHED SNHC+PTREE SEARCH'
    end subroutine srch_snhc_ptree

    subroutine kill_snhc_ptree( self )
        class(strategy2D_snhc_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_ptree

end module simple_strategy2D_snhc_ptree
