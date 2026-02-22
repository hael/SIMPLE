!@descr: 3D strategy for exhaustive projection matching
module simple_strategy3D_greedy
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: parameters
use simple_polarft_calc,     only: pftc_glob
use simple_oris,             only: oris
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
implicit none

public :: strategy3D_greedy
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy
contains
    procedure :: new         => new_greedy
    procedure :: srch        => srch_greedy
    procedure :: kill        => kill_greedy
    procedure :: oris_assign => oris_assign_greedy
end type strategy3D_greedy

contains

    subroutine new_greedy( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_greedy),  intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(strategy3D_spec),    intent(inout) :: spec
        class(builder),    target, intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_greedy

    subroutine srch_greedy( self, os, ithr )
        class(strategy3D_greedy), intent(inout) :: self
        class(oris),              intent(inout) :: os
        integer,                  intent(in)    :: ithr
        integer :: iref, isample, loc(1)
        real    :: inpl_corrs(self%s%nrots)
        if( os%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
             ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! search
            do isample=1,self%s%nrefs
                iref = s3D%srch_order(isample,self%s%ithr) ! set the reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    if( self%s%p_ptr%l_sh_first )then
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                    endif
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nrefs
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy

    subroutine oris_assign_greedy( self )
        class(strategy3D_greedy), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy

    subroutine kill_greedy( self )
        class(strategy3D_greedy), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy

end module simple_strategy3D_greedy
