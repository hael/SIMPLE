module simple_strategy2D_greedy
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
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
        integer :: iref,inpl_ind
        real    :: corrs(self%s%nrots),inpl_corr,corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            corr = -huge(corr)
            do iref=1,self%s%nrefs
                if( .not.s2D%cls_mask(iref,self%s%ithr) )cycle
                ! if( s2D%cls_pops(iref) == 0 )cycle
                ! class best
                call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
                inpl_ind  = maxloc(corrs, dim=1)
                inpl_corr = corrs(inpl_ind)
                ! updates global best
                if( inpl_corr >= corr )then
                    corr              = inpl_corr
                    self%s%best_class = iref
                    self%s%best_corr  = inpl_corr
                    self%s%best_rot   = inpl_ind
                endif
                ! keep track of visited classes
                if( self%s%l_ptclw )then
                    s2D%cls_searched(iref,self%s%ithr)  = .true.
                    s2D%cls_corrs(iref,self%s%ithr) = inpl_corr
                endif
            end do
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
