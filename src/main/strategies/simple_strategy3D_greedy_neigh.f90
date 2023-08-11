! concrete strategy3D: neighbourhood refinement
module simple_strategy3D_greedy_neigh
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_greedy_neigh
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_neigh
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_greedy_neigh
    procedure          :: srch        => srch_greedy_neigh
    procedure          :: oris_assign => oris_assign_greedy_neigh
    procedure          :: kill        => kill_greedy_neigh
end type strategy3D_greedy_neigh

contains

    subroutine new_greedy_neigh( self, spec )
        class(strategy3D_greedy_neigh), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_greedy_neigh

    subroutine srch_greedy_neigh( self, ithr )
        class(strategy3D_greedy_neigh), intent(inout) :: self
        integer,                        intent(in)    :: ithr
        type(ori) :: o
        integer   :: iref, iproj, loc(1)
        real      :: inpl_corrs(self%s%nrots)
        logical   :: lnns(params_glob%nspace)
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! initialize, ctd
            call build_glob%spproj_field%get_ori(self%s%iptcl, o)
            lnns = .false.
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, params_glob%athres, lnns)
            self%s%nnn = count(lnns)
            ! search
            do iproj=1,params_glob%nspace
                if( .not. lnns(iproj) ) cycle
                iref = (self%s%prev_state - 1) * params_glob%nspace + iproj
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nnn
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_neigh

    subroutine oris_assign_greedy_neigh( self )
        class(strategy3D_greedy_neigh), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_neigh

    subroutine kill_greedy_neigh( self )
        class(strategy3D_greedy_neigh),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_neigh

end module simple_strategy3D_greedy_neigh
