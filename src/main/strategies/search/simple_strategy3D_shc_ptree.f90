!@descr: 3D hybrid strategy: SHC coarse-node selection followed by probabilistic tree descent
module simple_strategy3D_shc_ptree
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_tree_utils, only: descend_tree_prob_fixed_state, get_tree_for_ref
use simple_strategy3D_utils
use simple_parameters,      only: parameters
use simple_oris,            only: oris
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
implicit none

public :: strategy3D_shc_ptree
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shc_ptree
contains
    procedure :: new         => new_shc_ptree
    procedure :: srch        => srch_shc_ptree
    procedure :: oris_assign => oris_assign_shc_ptree
    procedure :: kill        => kill_shc_ptree
end type strategy3D_shc_ptree

contains

    subroutine new_shc_ptree( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_shc_ptree), intent(inout) :: self
        class(parameters),           intent(in)    :: params
        class(strategy3D_spec),      intent(inout) :: spec
        class(builder),              intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_shc_ptree

    subroutine srch_shc_ptree( self, os, ithr )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_shc_ptree), intent(inout) :: self
        class(oris),                 intent(inout) :: os
        integer,                     intent(in)    :: ithr
        integer :: iref, isample, loc(1), itree, ntrees, inds(self%s%nrots)
        integer :: nrefs_coarse, nrefs_tree, iref_best_coarse
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots), corr_best_coarse
        if( os%get_state(self%s%iptcl) <= 0 )then
            call os%reject(self%s%iptcl)
            return
        endif
        self%s%ithr = ithr
        call self%s%prep4srch
        call self%s%inpl_srch_first
        self%s%nbetter = 0
        nrefs_coarse   = 0
        nrefs_tree     = 0
        corr_best_coarse = -huge(1.0)
        iref_best_coarse = self%s%prev_ref
        ntrees = self%s%b_ptr%block_tree%get_n_trees()
        ! Coarse SHC pass: keep a single selected coarse node/tree.
        do isample = 1, self%s%nrefs_sub
            iref = s3D%srch_order_sub(isample, self%s%ithr)
            if( .not. s3D%state_exists(s3D%proj_space_state(iref)) ) cycle
            if( self%s%p_ptr%l_doshift )then
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
            endif
            loc = angle_sampling(eulprob_dist_switch(inpl_corrs, self%s%p_ptr%cc_objfun), sorted_corrs, inds, s3D%smpl_inpl_athres(s3D%proj_space_state(iref)), self%s%p_ptr%prob_athres)
            call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
            nrefs_coarse = nrefs_coarse + 1
            if( inpl_corrs(loc(1)) > corr_best_coarse )then
                corr_best_coarse = inpl_corrs(loc(1))
                iref_best_coarse = iref
            endif
            if( inpl_corrs(loc(1)) > self%s%prev_corr )then
                self%s%nbetter = self%s%nbetter + 1
                exit
            endif
        end do
        if( nrefs_coarse > 0 )then
            itree = get_tree_for_ref(self%s, iref_best_coarse, ntrees)
            call descend_tree_prob_fixed_state(self%s, itree, corr_best_coarse, nrefs_tree, s3D%proj_space_state(iref_best_coarse))
        endif
        self%s%nrefs_eval = nrefs_coarse ! don't add the tree evaluations, coarse samples only
        call self%s%inpl_srch
        call self%oris_assign
    end subroutine srch_shc_ptree

    subroutine oris_assign_shc_ptree( self )
        class(strategy3D_shc_ptree), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shc_ptree

    subroutine kill_shc_ptree( self )
        class(strategy3D_shc_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_shc_ptree

end module simple_strategy3D_shc_ptree
