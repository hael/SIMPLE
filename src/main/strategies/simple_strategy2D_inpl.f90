module simple_strategy2D_inpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
implicit none

public :: strategy2D_inpl
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_inpl
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
contains
    procedure :: new  => new_inpl
    procedure :: srch => srch_inpl
    procedure :: kill => kill_inpl
end type strategy2D_inpl

contains

    subroutine new_inpl( self, spec )
        class(strategy2D_inpl), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_inpl

    subroutine srch_inpl( self )
        class(strategy2D_inpl), intent(inout) :: self
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            ! greedy in-plane update
            call self%s%inpl_srch
            self%s%nrefs_eval = self%s%nrefs
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        if( DEBUG ) write(logfhandle,*) '>>> strategy2D_srch::FINISHED GREEDY IN-PLANE SEARCH'
    end subroutine srch_inpl

    subroutine kill_inpl( self )
        class(strategy2D_inpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_inpl

end module simple_strategy2D_inpl
