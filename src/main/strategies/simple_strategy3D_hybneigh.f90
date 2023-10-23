! concrete strategy3D: neighbourhood refinement
module simple_strategy3D_hybneigh
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_hybneigh
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_hybneigh
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_hybneigh
    procedure          :: srch        => srch_hybneigh
    procedure          :: oris_assign => oris_assign_hybneigh
    procedure          :: kill        => kill_hybneigh
end type strategy3D_hybneigh

contains

    subroutine new_hybneigh( self, spec )
        class(strategy3D_hybneigh), intent(inout) :: self
        class(strategy3D_spec),     intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_hybneigh

    subroutine srch_hybneigh( self, ithr )
        class(strategy3D_hybneigh), intent(inout) :: self
        integer,                    intent(in)    :: ithr
        type(ori) :: o
        integer   :: iref, isample, iproj, minnrefs, loc(1)
        real      :: inpl_corrs(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            call build_glob%spproj_field%get_ori(self%s%iptcl, o)
            lnns = .false.
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, params_glob%athres, lnns)
            self%s%nnn = count(lnns)
            minnrefs   = ceiling(real(self%s%nnn) * NEIGH_MINFRAC)
            ! search
            do isample=1,self%s%nrefs
                iref  = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                iproj = iref - (self%s%prev_state - 1) * params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                    ! update nbetter to keep track of how many improving solutions we have identified
                    if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                endif
                ! exit condition
                if( self%s%nbetter > 0 .and. self%s%nrefs_eval >= minnrefs ) exit
            end do
            ! doing shc when there is no candidate in the neighbor
            if( self%s%nbetter == 0 )then
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
                self%s%nrefs_eval = self%s%nrefs_eval / 2
            endif
            call self%s%inpl_srch   ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_hybneigh

    subroutine oris_assign_hybneigh( self )
        class(strategy3D_hybneigh), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_hybneigh

    subroutine kill_hybneigh( self )
        class(strategy3D_hybneigh), intent(inout) :: self
        call self%s%kill
    end subroutine kill_hybneigh

end module simple_strategy3D_hybneigh
