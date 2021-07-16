! concrete strategy3D: probabilistic multi-state refinement
module simple_strategy3D_neigh_multi
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_neigh_multi
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_neigh_multi
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_neigh_multi
    procedure          :: srch        => srch_neigh_multi
    procedure          :: oris_assign => oris_assign_neigh_multi
    procedure          :: kill        => kill_neigh_multi
end type strategy3D_neigh_multi

contains

    subroutine new_neigh_multi( self, spec )
        class(strategy3D_neigh_multi), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_neigh_multi

    subroutine srch_neigh_multi( self, ithr )
        use simple_ori, only: ori
        class(strategy3D_neigh_multi), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        type(ori) :: o
        integer   :: iref,nrefs,iproj
        real      :: inpl_corrs(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            nrefs = self%s%nrefs
            call build_glob%spproj_field%get_ori(self%s%iptcl, o)
            call build_glob%eulspace%nearest_proj_neighbors(o, params_glob%nnn, lnns)
            ! search
            do iproj=1,params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                iref = (self%s%prev_state - 1)*params_glob%nspace + iproj
                call per_ref_srch
            end do
            self%s%nrefs_eval = nrefs
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
                ! calculate in-plane correlations
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                ! identify the top scoring in-plane angle
                loc = maxloc(inpl_corrs)
                ! stash
                call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)), .true.)
            endif
        end subroutine per_ref_srch

    end subroutine srch_neigh_multi

    subroutine oris_assign_neigh_multi( self )
        class(strategy3D_neigh_multi), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_neigh_multi

    subroutine kill_neigh_multi( self )
        class(strategy3D_neigh_multi),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_neigh_multi

end module simple_strategy3D_neigh_multi
