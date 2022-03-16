! concrete strategy3D: refinement
module simple_strategy3D_proj
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_proj
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_proj
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_proj
    procedure          :: srch        => srch_proj
    procedure          :: oris_assign => oris_assign_proj
    procedure          :: kill        => kill_proj
end type strategy3D_proj

contains

    subroutine new_proj( self, spec )
        class(strategy3D_proj), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_proj

    subroutine srch_proj( self, ithr )
        class(strategy3D_proj), intent(inout) :: self
        integer,                intent(in)    :: ithr
        integer :: iref, isample
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            do isample=1,self%s%nrefs
                iref = s3D%srch_order(self%s%ithr,isample) ! set the stochastic reference index
                call per_ref_srch                          ! actual search
                ! exit condition
                if( self%s%nbetter > 0 ) exit
            end do
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif

    contains

        subroutine per_ref_srch
            real :: corr
            if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                corr = pftcc_glob%gencorr_for_rot(iref, self%s%iptcl, self%s%prev_shvec, self%s%prev_roind)
                call self%s%store_solution(iref, self%s%prev_roind, corr)
                ! update nbetter to keep track of how many improving solutions we have identified
                if( corr > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
            end if
        end subroutine per_ref_srch

    end subroutine srch_proj

    subroutine oris_assign_proj( self )
        class(strategy3D_proj), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_proj

    subroutine kill_proj( self )
        class(strategy3D_proj),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_proj

end module simple_strategy3D_proj
