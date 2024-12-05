! concrete strategy3D: probabilistic in states refinement
module simple_strategy3D_prob_state
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_prob_state
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_prob_state
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_prob_state
    procedure :: srch        => srch_prob_state
    procedure :: kill        => kill_prob_state
    procedure :: oris_assign => oris_assign_prob_state
end type strategy3D_prob_state

contains

    subroutine new_prob_state( self, spec )
        class(strategy3D_prob_state), intent(inout) :: self
        class(strategy3D_spec),       intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_prob_state

    subroutine srch_prob_state( self, ithr )
        use simple_eul_prob_tab, only: eulprob_corr_switch
        class(strategy3D_prob_state), intent(inout) :: self
        integer,                      intent(in)    :: ithr
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
                    call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
                endif
            else
                call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
            endif
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_prob_state

    subroutine oris_assign_prob_state( self )
        class(strategy3D_prob_state), intent(inout) :: self
    end subroutine oris_assign_prob_state

    subroutine kill_prob_state( self )
        class(strategy3D_prob_state), intent(inout) :: self
        call self%s%kill
    end subroutine kill_prob_state

end module simple_strategy3D_prob_state
