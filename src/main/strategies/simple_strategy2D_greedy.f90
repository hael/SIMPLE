module simple_strategy2D_greedy
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_greedy
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_greedy
  contains
    procedure :: new  => new_greedy
    procedure :: srch => srch_greedy
    procedure :: kill => kill_greedy
end type strategy2D_greedy

contains

    subroutine new_greedy( self, spec )
        class(strategy2D_greedy), intent(inout) :: self
        class(strategy2D_spec),   intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_greedy

    subroutine srch_greedy( self )
        class(strategy2D_greedy), intent(inout) :: self
        integer :: iref,inpl_ind
        real    :: corrs(self%s%nrots),inpl_corr,corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! greedy search
            corr = -huge(corr)
            do iref=1,self%s%nrefs
                if( s2D%cls_pops(iref) == 0 )cycle
                ! class best
                if( self%s%l_sh_first )then
                    call pftc_glob%gen_corrs(iref, self%s%iptcl, self%s%xy_first, corrs)
                else
                    call pftc_glob%gen_corrs(iref, self%s%iptcl, corrs)
                endif
                inpl_ind  = maxloc(corrs, dim=1)
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
                if( self%s%l_sh_first )then
                    call pftc_glob%gen_corrs(self%s%best_class, self%s%iptcl, self%s%xy_first, corrs, kweight=.true.)
                else
                    call pftc_glob%gen_corrs(self%s%best_class, self%s%iptcl, corrs, kweight=.true.)
                endif
                self%s%best_rot  = maxloc(corrs, dim=1)
                self%s%best_corr = corrs(inpl_ind)
            endif
            self%s%nrefs_eval = self%s%nrefs
            call self%s%inpl_srch
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy

    subroutine kill_greedy( self )
        class(strategy2D_greedy), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy

end module simple_strategy2D_greedy
