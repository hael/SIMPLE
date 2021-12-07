! concrete strategy3D: continuous single-state refinement
module simple_strategy3D_cont
include 'simple_lib.f08'
use simple_strategy3D_alloc    ! singleton
use simple_strategy3D,         only: strategy3D
use simple_strategy3D_srch,    only: strategy3D_srch, strategy3D_spec
use simple_pftcc_orisrch_grad, only: pftcc_orisrch_grad
use simple_ori,                only: ori
use simple_parameters,         only: params_glob
use simple_builder,            only: build_glob
implicit none

public :: strategy3D_cont
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cont
    type(strategy3D_srch)    :: s
    type(strategy3D_spec)    :: spec
    type(pftcc_orisrch_grad) :: cont_srch
    type(ori) :: o
    integer   :: irot
    real      :: corr
contains
    procedure :: new         => new_cont
    procedure :: srch        => srch_cont
    procedure :: oris_assign => oris_assign_cont
    procedure :: kill        => kill_cont
end type strategy3D_cont

contains

    subroutine new_cont( self, spec )
        class(strategy3D_cont), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
        call self%cont_srch%new
    end subroutine new_cont

    subroutine srch_cont( self, ithr )
        class(strategy3D_cont), intent(inout) :: self
        integer,                intent(in)    :: ithr
        real, allocatable :: cxy(:)
        logical :: found_better
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! initialize
            call self%s%prep4srch
            call self%cont_srch%set_particle(self%s%iptcl)
            call build_glob%spproj_field%get_ori(self%s%iptcl, self%o)
            cxy       = self%cont_srch%minimize(self%o, params_glob%athres/2., params_glob%trs, found_better)
            self%corr = cxy(1)
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_cont

    subroutine oris_assign_cont( self )
        use simple_ori,  only: ori
        class(strategy3D_cont), intent(inout) :: self
        ! no update of convergence stats, only ori info
        call build_glob%spproj_field%set_euler(self%s%iptcl, self%o%get_euler())
        call build_glob%spproj_field%set_shift(self%s%iptcl, self%o%get_2Dshift()) ! vector addition in pftcc_orisrch_grad class
    end subroutine oris_assign_cont

    subroutine kill_cont( self )
        class(strategy3D_cont),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_cont

end module simple_strategy3D_cont
