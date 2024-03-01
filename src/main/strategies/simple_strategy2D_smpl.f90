! concrete strategy3D: stochastic top sampling
module simple_strategy2D_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_smpl
contains
    procedure :: new  => new_smpl
    procedure :: srch => srch_smpl
    procedure :: kill => kill_smpl
end type strategy2D_smpl

contains

    subroutine new_smpl( self, spec )
        class(strategy2D_smpl), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_smpl

    subroutine srch_smpl( self )
        use simple_regularizer, only: reg_dist_switch
        class(strategy2D_smpl), intent(inout) :: self
        integer :: iref, locs(self%s%nrefs), inds(self%s%nrots), irot
        real    :: inpl_corrs(self%s%nrots), sorted_inpl_corrs(self%s%nrots)
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            ! search
            s2D%cls_corrs(:,self%s%ithr) = TINY
            do iref=1,self%s%nrefs
                if( s2D%cls_pops(iref) == 0 )cycle      
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                irot       = greedy_sampling(reg_dist_switch(inpl_corrs), sorted_inpl_corrs, inds, s2D%smpl_inpl_ns)
                locs(iref) = irot
                s2D%cls_corrs(iref,self%s%ithr) = inpl_corrs(irot)
            enddo
            iref = greedy_sampling(reg_dist_switch(s2D%cls_corrs(:,self%s%ithr)), s2D%smpl_refs_ns)
            if( params_glob%cc_objfun == OBJFUN_CC .and. params_glob%l_kweight_rot )then
                ! back-calculating in-plane angle with k-weighing
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs, kweight=.true.)
                irot = greedy_sampling(reg_dist_switch(inpl_corrs), sorted_inpl_corrs, inds, s2D%smpl_inpl_ns)
                locs(iref) = irot
                s2D%cls_corrs(iref,self%s%ithr) = inpl_corrs(irot)
            endif
            self%s%best_class = iref
            self%s%best_rot   = locs(iref)
            self%s%best_corr  = s2D%cls_corrs(iref,self%s%ithr)
            call self%s%inpl_srch
            self%s%nrefs_eval = self%s%nrefs - s2D%smpl_refs_ns ! for reporting only
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_smpl

    subroutine kill_smpl( self )
        class(strategy2D_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_smpl

end module simple_strategy2D_smpl