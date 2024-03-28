module simple_strategy2D_greedy_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_greedy_smpl
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_greedy_smpl
  contains
    procedure :: new  => new_greedy_smpl
    procedure :: srch => srch_greedy_smpl
    procedure :: kill => kill_greedy_smpl
end type strategy2D_greedy_smpl

contains

    subroutine new_greedy_smpl( self, spec )
        class(strategy2D_greedy_smpl), intent(inout) :: self
        class(strategy2D_spec),   intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_greedy_smpl

    subroutine srch_greedy_smpl( self )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy2D_greedy_smpl), intent(inout) :: self
        integer :: inds(self%s%nrots),loc(1),iref,inpl_ind
        real    :: corrs(self%s%nrots),sorted_corrs(self%s%nrots),inpl_corr,corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            corr = -huge(corr)
            do iref=1,self%s%nrefs
                if( s2D%cls_pops(iref) == 0 )cycle
                ! class best
                call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
                inpl_ind  = angle_sampling(eulprob_dist_switch(corrs), sorted_corrs, inds, s2D%smpl_inpl_athres)
                inpl_corr = corrs(inpl_ind)
                ! updates global best
                if( inpl_corr >= corr )then
                    corr              = inpl_corr
                    self%s%best_class = iref
                    self%s%best_corr  = inpl_corr
                    self%s%best_rot   = inpl_ind
                endif
            end do
            if( params_glob%cc_objfun == OBJFUN_CC .and. params_glob%l_kweight_rot )then
                ! back-calculating in-plane angle with k-weighing
                call pftcc_glob%gencorrs(self%s%best_class, self%s%iptcl, corrs, kweight=.true.)
                self%s%best_rot = angle_sampling(eulprob_dist_switch(corrs), sorted_corrs, inds, s2D%smpl_inpl_athres)
                ! self%s%best_rot  = greedy_sampling(eulprob_dist_switch(corrs), sorted_corrs, inds, s2D%smpl_inpl_ns)
                self%s%best_corr = corrs(inpl_ind)
            endif
            self%s%nrefs_eval = self%s%nrefs
            call self%s%inpl_srch
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_smpl

    subroutine kill_greedy_smpl( self )
        class(strategy2D_greedy_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_smpl

end module simple_strategy2D_greedy_smpl
