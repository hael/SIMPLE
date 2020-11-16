module simple_strategy2D_eval
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
implicit none

public :: strategy2D_eval
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_eval
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
contains
    procedure :: new  => new_eval
    procedure :: srch => srch_eval
    procedure :: kill => kill_eval
end type strategy2D_eval

contains

    subroutine new_eval( self, spec )
        class(strategy2D_eval), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_eval

    subroutine srch_eval( self )
        class(strategy2D_eval), intent(inout) :: self
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        if( DEBUG ) write(logfhandle,*) '>>> strategy2D_srch::FINISHED EVALUATION'
    end subroutine srch_eval

    subroutine kill_eval( self )
        class(strategy2D_eval), intent(inout) :: self
        call self%s%kill
    end subroutine kill_eval

end module simple_strategy2D_eval
