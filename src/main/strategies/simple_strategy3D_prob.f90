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
        integer :: iref, iptcl, irot
        real    :: corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            iptcl = self%s%iptcl
            iref  = self%spec%reg_obj%ptcl_ref_map(iptcl)
            corr  =     self%spec%reg_obj%dist_loc_tab(iref, iptcl, 1)
            irot  = int(self%spec%reg_obj%dist_loc_tab(iref, iptcl, 2))
            call self%s%store_solution(iref, irot, corr)
            if( self%s%doshift )then
                call self%s%inpl_srch
                irot = s3D%proj_space_inplinds(self%s%ithr, iref)
                corr = s3D%proj_space_corrs(self%s%ithr, iref)
                call assign_ori(self%s, iref, irot, corr, sh_in=s3D%proj_space_shift(:,iref,self%s%ithr))
            else
                call assign_ori(self%s, iref, irot, corr)
            endif
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
