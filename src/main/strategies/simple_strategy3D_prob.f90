! concrete strategy3D: greedy refinement
module simple_strategy3D_prob
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_prob
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_prob
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_prob
    procedure :: srch        => srch_prob
    procedure :: kill        => kill_prob
    procedure :: oris_assign => oris_assign_prob
end type strategy3D_prob

contains

    subroutine new_prob( self, spec )
        class(strategy3D_prob), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_prob

    subroutine srch_prob( self, ithr )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: ithr
        integer :: iref, iptcl, ir, irot
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            iptcl = self%s%iptcl
            iref  = self%spec%reg_obj%ptcl_ref_map(iptcl)
            irot  = self%spec%reg_obj%ptcl_loc_map(iptcl)
            do ir = 1, self%s%nrefs
                call self%s%store_solution(ir, 1, 0.)
            enddo
            call self%s%store_solution(iref, self%spec%reg_obj%ref_ptcl_tab(iptcl, iref, irot)%loc,&
                                            &self%spec%reg_obj%ref_ptcl_tab(iptcl, iref, irot)%prob,&
                                            &self%spec%reg_obj%ref_ptcl_tab(iptcl, iref, irot)%sh,&
                                            &self%spec%reg_obj%ref_ptcl_tab(iptcl, iref, irot)%prob)
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = self%s%nrefs
            ! prepare orientation
            call self%oris_assign()
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_prob

    subroutine oris_assign_prob( self )
        class(strategy3D_prob), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_prob

    subroutine kill_prob( self )
        class(strategy3D_prob), intent(inout) :: self
        call self%s%kill
    end subroutine kill_prob

end module simple_strategy3D_prob
