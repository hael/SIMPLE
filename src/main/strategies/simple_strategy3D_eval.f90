module simple_strategy3D_eval
include 'simple_lib.f08'
use simple_strategy3D_utils
use simple_strategy3D_alloc  ! singleton
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
use simple_builder,          only: build_glob
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

    subroutine new_eval( self, spec )
        class(strategy3D_eval), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_eval

    subroutine srch_eval( self, ithr )
        class(strategy3D_eval), intent(inout) :: self
        integer,                intent(in)    :: ithr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            self%s%ithr = ithr
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            call self%s%store_solution(self%s%prev_ref, self%s%prev_roind, self%s%prev_corr, sh=self%s%prev_shvec)
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
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
