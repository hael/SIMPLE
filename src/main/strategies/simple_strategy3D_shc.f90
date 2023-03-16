! concrete strategy3D: refinement
module simple_strategy3D_shc
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_shc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shc
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_shc
    procedure          :: srch        => srch_shc
    procedure          :: oris_assign => oris_assign_shc
    procedure          :: kill        => kill_shc
end type strategy3D_shc

contains

    subroutine new_shc( self, spec )
        class(strategy3D_shc), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shc

    subroutine srch_shc( self, ithr )
        class(strategy3D_shc), intent(inout) :: self
        integer,               intent(in)    :: ithr
        integer :: iref, isample, loc(1)
        real    :: inpl_corrs(self%s%nrots)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            if( params_glob%l_obj_reg )then
                do isample=1,self%s%nrefs
                    iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                    if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                        ! identify the top scoring in-plane angle
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs, s3D%srch_order)
                        loc = maxloc(inpl_corrs)
                        call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                        ! update nbetter to keep track of how many improving solutions we have identified
                        if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                        ! keep track of how many references we are evaluating
                        self%s%nrefs_eval = self%s%nrefs_eval + 1
                    end if
                    ! exit condition
                    if( self%s%nbetter > 0 ) exit
                end do
            else
                do isample=1,self%s%nrefs
                    iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                    if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                        ! identify the top scoring in-plane angle
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                        loc = maxloc(inpl_corrs)
                        call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                        ! update nbetter to keep track of how many improving solutions we have identified
                        if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                        ! keep track of how many references we are evaluating
                        self%s%nrefs_eval = self%s%nrefs_eval + 1
                    end if
                    ! exit condition
                    if( self%s%nbetter > 0 ) exit
                end do
            endif
            call self%s%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shc

    subroutine oris_assign_shc( self )
        class(strategy3D_shc), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shc

    subroutine kill_shc( self )
        class(strategy3D_shc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_shc

end module simple_strategy3D_shc
