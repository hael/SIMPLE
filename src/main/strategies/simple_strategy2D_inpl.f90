module simple_strategy2D_inpl
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_inpl
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_inpl
  contains
    procedure :: new  => new_inpl
    procedure :: srch => srch_inpl
    procedure :: kill => kill_inpl
end type strategy2D_inpl

contains

    subroutine new_inpl( self, spec )
        class(strategy2D_inpl), intent(inout) :: self
        class(strategy2D_spec),   intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_inpl

    subroutine srch_inpl( self )
        class(strategy2D_inpl), intent(inout) :: self
        integer :: inpl_ind
        real    :: corrs(self%s%nrots),inpl_corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! inpl search
            call pftc_glob%gen_objfun_vals(self%s%prev_class, self%s%iptcl, [0.,0.], corrs)
            inpl_ind          = maxloc(corrs, dim=1)
            inpl_corr         = corrs(inpl_ind)
            self%s%best_class = self%s%prev_class
            self%s%best_corr  = inpl_corr
            self%s%best_rot   = inpl_ind
            self%s%nrefs_eval = self%s%nrefs
            call self%s%inpl_srch
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_inpl

    subroutine kill_inpl( self )
        class(strategy2D_inpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_inpl

end module simple_strategy2D_inpl
