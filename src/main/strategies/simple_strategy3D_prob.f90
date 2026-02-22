!@descr: 3D strategy for probabilistic projection matching
module simple_strategy3D_prob
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: parameters
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
use simple_oris,             only: oris
implicit none

public :: strategy3D_prob
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_prob
contains
    procedure :: new         => new_prob
    procedure :: srch        => srch_prob
    procedure :: kill        => kill_prob
    procedure :: oris_assign => oris_assign_prob
end type strategy3D_prob

contains

    subroutine new_prob( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_prob), intent(inout) :: self
        class(parameters),      intent(in)    :: params
        class(strategy3D_spec), intent(inout) :: spec
        class(builder),         intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_prob

    subroutine srch_prob( self, os, ithr )
        use simple_eul_prob_tab, only: eulprob_corr_switch
        class(strategy3D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer,                intent(in)    :: ithr
        integer :: iproj, iptcl_map, irot, istate, iref, rot
        real    :: corr
        if( os%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            self%s%nrefs_eval = self%s%nrefs
            iptcl_map = self%s%iptcl_map
            istate    =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%istate
            iproj     =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%iproj
            corr      = eulprob_corr_switch(self%spec%eulprob_obj_part%assgn_map(iptcl_map)%dist, self%s%p_ptr%cc_objfun)
            irot      =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%inpl
            iref      = (istate-1)*self%s%p_ptr%nspace + iproj
            if( self%s%doshift )then
                if( self%s%p_ptr%l_sh_first .or. self%s%p_ptr%l_prob_sh )then
                    if( self%spec%eulprob_obj_part%assgn_map(iptcl_map)%has_sh )then
                        call assign_ori(self%s, iref, irot, corr,&
                        &[self%spec%eulprob_obj_part%assgn_map(iptcl_map)%x,&
                        & self%spec%eulprob_obj_part%assgn_map(iptcl_map)%y], corr)
                    else
                        call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
                    endif
                else
                    ! take care of the shift with current assignment
                    rot = irot
                    call self%s%inpl_srch(ref=iref, irot_in=rot)
                    if( rot > 0 )then
                        call assign_ori(self%s, iref, rot, corr, s3D%proj_space_shift(:,iref,self%s%ithr), corr)
                    else
                        call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
                    endif
                endif
            else
                call assign_ori(self%s, iref, irot, corr, [0.,0.], corr)
            endif
        else
            call os%reject(self%s%iptcl)
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
