! concrete strategy3D: continuous stochastic neighborhood hill climbing
module simple_strategy3D_snhcc
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: params_glob
use simple_builder,         only: build_glob
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
use simple_cartft_corrcalc, only: cartftcc_glob
use simple_ori,             only: ori
implicit none

public :: strategy3D_snhcc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_snhcc
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_snhcc
    procedure          :: srch        => srch_snhcc
    procedure          :: oris_assign => oris_assign_snhcc
    procedure          :: kill        => kill_snhcc
end type strategy3D_snhcc

contains

    subroutine new_snhcc( self, spec )
        class(strategy3D_snhcc), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_snhcc

    subroutine srch_snhcc( self, ithr )
        class(strategy3D_snhcc), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        integer   :: isample
        type(ori) :: o, o_best
        real      :: corr, corr_best
        o = self%s%o_prev
        o_best = o
        ! continuous sochastic search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4_cont_srch
            self%s%nbetter = 0
            corr_best      = -1.
            do isample=1,self%spec%szsn
                ! make a random rotation matrix within the assymetric unit
                call build_glob%pgrpsyms%rnd_euler(self%s%o_prev, self%s%athres, o)
                ! calculate Cartesian corr
                call cartftcc_glob%project_and_correlate(self%s%iptcl, o, corr)
                ! update condition
                if( corr >= corr_best )then
                    corr_best = corr
                    o_best    = o
                endif
            end do
            call build_glob%spproj_field%set_ori(self%s%iptcl, o_best)
            self%s%nrefs_eval = self%spec%szsn
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_snhcc

    subroutine oris_assign_snhcc( self )
        class(strategy3D_snhcc), intent(inout) :: self
    end subroutine oris_assign_snhcc

    subroutine kill_snhcc( self )
        class(strategy3D_snhcc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhcc

end module simple_strategy3D_snhcc
