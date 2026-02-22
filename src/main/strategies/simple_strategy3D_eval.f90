!@descr: 3D strategy for objective function evaluation
module simple_strategy3D_eval
use simple_core_module_api
use simple_strategy3D_utils
use simple_strategy3D_alloc
use simple_parameters,      only: parameters
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
use simple_oris,            only: oris
implicit none

public :: strategy3D_eval
private

type, extends(strategy3D) :: strategy3D_eval
contains
    procedure :: new         => new_eval
    procedure :: srch        => srch_eval
    procedure :: oris_assign => oris_assign_eval
    procedure :: kill        => kill_eval
end type strategy3D_eval

contains

    subroutine new_eval( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_eval),    intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(strategy3D_spec),    intent(inout) :: spec
        class(builder),    target, intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_eval

    subroutine srch_eval( self, os, ithr )
        class(strategy3D_eval), intent(inout) :: self
        class(oris),             intent(inout) :: os
        integer,                intent(in)    :: ithr
        if( os%get_state(self%s%iptcl) > 0 )then
            self%s%ithr = ithr
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            call self%s%store_solution(self%s%prev_ref, self%s%prev_roind, self%s%prev_corr, sh=self%s%prev_shvec)
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_eval

    subroutine oris_assign_eval( self )
        class(strategy3D_eval), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_eval

    subroutine kill_eval( self )
        class(strategy3D_eval), intent(inout) :: self
        call self%s%kill
    end subroutine kill_eval

end module simple_strategy3D_eval
