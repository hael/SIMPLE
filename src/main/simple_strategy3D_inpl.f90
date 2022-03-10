! concrete strategy3D: inpl refinement
module simple_strategy3D_inpl
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_inpl
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_inpl
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_inpl
    procedure :: srch        => srch_inpl
    procedure :: kill        => kill_inpl
    procedure :: oris_assign => oris_assign_inpl
end type strategy3D_inpl

contains

    subroutine new_inpl( self, spec )
        class(strategy3D_inpl), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_inpl

    subroutine srch_inpl( self, ithr )
        class(strategy3D_inpl), intent(inout) :: self
        integer,                intent(in)    :: ithr
        real    :: inpl_corrs(self%s%nrots)
        integer :: loc(1)
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! srch
            call pftcc_glob%gencorrs(self%s%prev_proj, self%s%iptcl, inpl_corrs)
            loc = maxloc(inpl_corrs)
            call self%s%store_solution(self%s%prev_proj, loc(1),  inpl_corrs(loc(1)))
            ! in inpl mode, we evaluate one reference
            self%s%nrefs_eval = 1
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign()
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_inpl

    subroutine oris_assign_inpl( self )
        class(strategy3D_inpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_inpl

    subroutine kill_inpl( self )
        class(strategy3D_inpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_inpl

end module simple_strategy3D_inpl
