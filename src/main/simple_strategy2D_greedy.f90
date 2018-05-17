module simple_strategy2D_greedy
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! use all in there
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
implicit none

public :: strategy2D_greedy
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_greedy
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
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
        integer :: iref,loc(1),inpl_ind
        real    :: corrs(self%s%nrots),inpl_corr,corr
        if( self%s%a_ptr%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            corr = self%s%prev_corr
            do iref=1,self%s%nrefs
                if( s2D%cls_pops(iref) == 0 )cycle
                call self%s%pftcc_ptr%gencorrs(iref, self%s%iptcl, corrs)
                !call pftcc%gencorrs(iref, self%s%iptcl, corrs)
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)
                if( inpl_corr >= corr )then
                    corr              = inpl_corr
                    self%s%best_class = iref
                    self%s%best_corr  = inpl_corr
                    self%s%best_rot   = inpl_ind
                endif
            end do
            self%s%nrefs_eval = self%s%nrefs
            call self%s%inpl_srch
            call self%s%store_solution
        else
            call self%s%a_ptr%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> STRATEGY2D_GREEDY :: FINISHED STOCHASTIC SEARCH'
    end subroutine srch_greedy

    subroutine kill_greedy( self )
        class(strategy2D_greedy), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy

end module simple_strategy2D_greedy
