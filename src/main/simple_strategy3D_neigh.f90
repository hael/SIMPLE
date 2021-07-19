! concrete strategy3D: neighbourhood refinement
module simple_strategy3D_neigh
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_neigh
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_neigh
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_neigh
    procedure          :: srch        => srch_neigh
    procedure          :: oris_assign => oris_assign_neigh
    procedure          :: kill        => kill_neigh
end type strategy3D_neigh

contains

    subroutine new_neigh( self, spec )
        class(strategy3D_neigh), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_neigh

    subroutine srch_neigh( self, ithr )
        use simple_ori, only: ori
        class(strategy3D_neigh), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        type(ori) :: o
        integer   :: iref, iproj, minnrefs
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
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, params_glob%athres, lnns)
            self%s%nnn = count(lnns)
            minnrefs   = ceiling(real(self%s%nnn) * NEIGH_MINFRAC)
            ! search
            do iproj=1,params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                iref = (self%s%prev_state - 1) * params_glob%nspace + iproj
                call per_ref_srch
                ! exit condition
                if( self%s%nbetter > 0 .and. self%s%nrefs_eval >= minnrefs ) exit
            end do
            call sort_corrs(self%s) ! sort in correlation projection direction space
            call self%s%inpl_srch   ! search shifts
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif

    contains

        subroutine per_ref_srch
            integer :: loc(1)
            if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                ! identify the top scoring in-plane angle
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                loc = maxloc(inpl_corrs)
                call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)), .true.)
                ! update nbetter to keep track of how many improving solutions we have identified
                if( inpl_corrs(loc(1)) > self%s%prev_corr ) self%s%nbetter = self%s%nbetter + 1
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
            endif
        end subroutine per_ref_srch

    end subroutine srch_neigh

    subroutine oris_assign_neigh( self )
        class(strategy3D_neigh), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_neigh

    subroutine kill_neigh( self )
        class(strategy3D_neigh),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_neigh

end module simple_strategy3D_neigh
