! concrete strategy3D: probabilistic multi-state refinement
module simple_strategy3D_multi
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_multi
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_multi
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_multi
    procedure          :: srch        => srch_multi
    procedure          :: oris_assign => oris_assign_multi
    procedure          :: kill        => kill_multi
end type strategy3D_multi

contains

    subroutine new_multi( self, spec )
        class(strategy3D_multi), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_multi

    subroutine srch_multi( self, ithr )
        class(strategy3D_multi), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        integer :: iref,isample,nrefs
        real    :: inpl_corrs(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            if( self%s%neigh )then
                nrefs = self%s%nnnrefs
            else
                nrefs = self%s%nrefs
            endif
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            do isample=1,nrefs
                iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                call per_ref_srch                           ! actual search
                ! exit condition
                if( self%s%nbetter > 0 ) exit
            end do
            call sort_corrs(self%s)  ! sort in correlation projection direction space
            call self%s%inpl_srch    ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif

    contains

        subroutine per_ref_srch
            integer :: loc(1)
            if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                ! identify the top scoring in-plane angle
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                loc            = maxloc(inpl_corrs)
                call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)), .true.)
                ! update nbetter to keep track of how many improving solutions we have identified
                if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
            end if
        end subroutine per_ref_srch

    end subroutine srch_multi

    subroutine oris_assign_multi( self )
        class(strategy3D_multi), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_multi

    subroutine kill_multi( self )
        class(strategy3D_multi),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_multi

end module simple_strategy3D_multi
