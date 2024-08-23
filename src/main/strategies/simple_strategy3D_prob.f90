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
        use simple_eul_prob_tab, only: eulprob_corr_switch
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: ithr
        integer   :: iproj, iptcl, iptcl_map, irot, istate, iref
        real      :: corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            iptcl     = self%s%iptcl
            iptcl_map = self%s%iptcl_map
            istate    =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%istate
            iproj     =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%iproj
            corr      = eulprob_corr_switch(self%spec%eulprob_obj_part%assgn_map(iptcl_map)%dist)
            irot      =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%inpl
            iref      = (istate-1)*params_glob%nspace + iproj
            if( self%s%doshift )then
                if( self%spec%eulprob_obj_part%assgn_map(iptcl_map)%has_sh )then
                    call assign_ori(self%s, iref, irot, corr,&
                    &[self%spec%eulprob_obj_part%assgn_map(iptcl_map)%x,&
                    & self%spec%eulprob_obj_part%assgn_map(iptcl_map)%y], corr)
                else
                    call self%s%inpl_srch(ref=iref, irot_in=irot)
                    ! checking if shift search is good
                    if( s3D%proj_space_inplinds(iref, ithr) < 1 )then
                        call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
                    else
                        irot = s3D%proj_space_inplinds(iref, ithr)
                        corr = s3D%proj_space_corrs(   iref, ithr)
                        call assign_ori(self%s, iref, irot, corr, s3D%proj_space_shift(:,iref,ithr), corr)
                    endif
                endif
            else
                call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
            endif
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_prob

    subroutine oris_assign_prob( self )
        class(strategy3D_prob), intent(inout) :: self
    end subroutine oris_assign_prob

    subroutine kill_prob( self )
        class(strategy3D_prob), intent(inout) :: self
        call self%s%kill
    end subroutine kill_prob

end module simple_strategy3D_prob
