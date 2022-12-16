! concrete strategy3D: greedy refinement
module simple_strategy3D_greedy_sub
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_greedy_sub
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_sub
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_greedy_sub
    procedure :: srch        => srch_greedy_sub
    procedure :: kill        => kill_greedy_sub
    procedure :: oris_assign => oris_assign_greedy_sub
end type strategy3D_greedy_sub

contains

    subroutine new_greedy_sub( self, spec )
        class(strategy3D_greedy_sub), intent(inout) :: self
        class(strategy3D_spec),   intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_greedy_sub

    subroutine srch_greedy_sub( self, ithr )
        class(strategy3D_greedy_sub), intent(inout) :: self
        integer,                      intent(in)    :: ithr
        integer   :: iref, isample, loc(1)
        real      :: inpl_corrs(self%s%nrots)
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! search
            do isample=1,self%s%nrefs_sub
                iref = s3D%srch_order_sub(self%s%ithr,isample) ! set the reference index
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    ! identify the top scoring in-plane angle
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nrefs_sub
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign()
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_sub

    subroutine oris_assign_greedy_sub( self )
        class(strategy3D_greedy_sub), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_sub

    subroutine kill_greedy_sub( self )
        class(strategy3D_greedy_sub), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_sub

end module simple_strategy3D_greedy_sub
