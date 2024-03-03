! concrete strategy3D: greedy refinement
module simple_strategy3D_shift
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_shift
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shift
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_shift
    procedure :: srch        => srch_shift
    procedure :: kill        => kill_shift
    procedure :: oris_assign => oris_assign_shift
end type strategy3D_shift

contains

    subroutine new_shift( self, spec )
        class(strategy3D_shift), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shift

    subroutine srch_shift( self, ithr )
        class(strategy3D_shift), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        integer   :: loc(1), irot_best
        real      :: inpl_corrs(self%s%nrots), score_best, xsh, ysh
        logical   :: score_set
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! srch
            score_set = .false.
            xsh = -params_glob%trs
            do while(xsh <= params_glob%trs)
                ysh = -params_glob%trs
                do while( ysh <= params_glob%trs)
                    call pftcc_glob%gencorrs(self%s%prev_ref, self%s%iptcl, [xsh,ysh], inpl_corrs)
                    loc = maxloc(inpl_corrs)
                    if( score_set )then
                        if( inpl_corrs(loc(1)) > score_best )then
                            score_best = inpl_corrs(loc(1))
                            irot_best  = loc(1)
                        endif
                    else
                        score_best = inpl_corrs(loc(1))
                        irot_best  = loc(1)
                    endif
                    ysh = ysh + 0.5
                end do
                xsh = xsh + 0.5
            end do
            call self%s%store_solution(self%s%prev_ref, irot_best, score_best)
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nrefs
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare orientation
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shift

    subroutine oris_assign_shift( self )
        class(strategy3D_shift), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_shift

    subroutine kill_shift( self )
        class(strategy3D_shift), intent(inout) :: self
        call self%s%kill
    end subroutine kill_shift

end module simple_strategy3D_shift
