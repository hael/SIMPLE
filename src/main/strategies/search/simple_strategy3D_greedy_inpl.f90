!@descr: 3D strategy for exhaustive in-plane matching of a single re-projection
module simple_strategy3D_greedy_inpl
use simple_core_module_api
use simple_strategy3D_alloc, only: s3D
use simple_strategy3D_utils, only: extract_peak_ori
use simple_parameters,       only: parameters
use simple_oris,             only: oris
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
implicit none

public :: strategy3D_greedy_inpl
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_inpl
contains
    procedure :: new         => new_greedy_inpl
    procedure :: srch        => srch_greedy_inpl
    procedure :: kill        => kill_greedy_inpl
    procedure :: oris_assign => oris_assign_greedy_inpl
end type strategy3D_greedy_inpl

contains

    subroutine new_greedy_inpl( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_greedy_inpl), intent(inout) :: self
        class(parameters),        intent(in)    :: params
        class(strategy3D_spec),   intent(inout) :: spec
        class(builder),           intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_greedy_inpl

    subroutine srch_greedy_inpl( self, os, ithr )
        class(strategy3D_greedy_inpl), intent(inout) :: self
        class(oris),              intent(inout) :: os
        integer,                  intent(in)    :: ithr
        integer :: iref, isample, loc(1)
        real    :: inpl_corrs(self%s%nrots)
        if( os%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! the single reference is the whole search space
            self%s%nrefs_eval = self%s%nrefs
            ! shift search of previous best reference
            call self%s%store_solution(self%s%prev_ref, self%s%prev_roind, self%s%prev_corr)
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_inpl

    subroutine oris_assign_greedy_inpl( self )
        class(strategy3D_greedy_inpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_inpl

    subroutine kill_greedy_inpl( self )
        class(strategy3D_greedy_inpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_inpl

end module simple_strategy3D_greedy_inpl
