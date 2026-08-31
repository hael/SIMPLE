!@descr: 3D strategy for probabilistic projection matching
module simple_strategy3D_prob
use simple_core_module_api
use simple_strategy3D_utils, only: assign_ori
use simple_parameters,      only: parameters
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
use simple_oris,            only: oris
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
        use simple_eul_prob_tab_utils, only: eulprob_corr_switch
        class(strategy3D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer,                intent(in)    :: ithr
        integer :: iproj, iptcl_map, irot, istate, iref
        real    :: corr, frac, sh(2), inpl_coord, fixed_euler(3), assigned_euler(3)
        logical :: inpl_valid, l_fixed_projection
        if( os%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4prob
            self%s%nrefs_eval = self%s%nrefs
            iptcl_map = self%s%iptcl_map
            istate    =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%istate
            if( istate < 1 )then
                ! particle received no valid assignment (e.g. a state emptied mid-run);
                ! keep its current orientation unchanged and move on
                return
            endif
            l_fixed_projection = trim(self%s%p_ptr%multivol_mode) == 'input_oris_fixed' .and. &
                &trim(self%s%p_ptr%refine) == 'prob_state'
            if( l_fixed_projection ) fixed_euler = self%s%b_ptr%spproj_field%get_euler(self%s%iptcl)
            iproj     =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%iproj
            corr      = eulprob_corr_switch(self%spec%eulprob_obj_part%assgn_map(iptcl_map)%dist, self%s%p_ptr%cc_objfun)
            irot      =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%inpl
            frac      =                     self%spec%eulprob_obj_part%assgn_map(iptcl_map)%frac
            if( frac <= 0. ) frac = 100.
            iref      = (istate-1)*self%s%p_ptr%nspace + iproj
            sh = 0.
            if( self%s%doshift .and. &
                &self%spec%eulprob_obj_part%assgn_map(iptcl_map)%has_sh )then
                sh = [self%spec%eulprob_obj_part%assgn_map(iptcl_map)%x,&
                    & self%spec%eulprob_obj_part%assgn_map(iptcl_map)%y]
            endif
            if( self%s%uses_continuous_refinement() )then
                call self%s%refine_assignment_continuously(iref, irot, corr, sh, &
                    &inpl_coord, inpl_valid)
                call assign_ori(self%s, iref, irot, corr, sh, inpl_coord, inpl_valid)
            else
                call assign_ori(self%s, iref, irot, corr, sh)
            endif
            if( l_fixed_projection )then
                assigned_euler      = self%s%b_ptr%spproj_field%get_euler(self%s%iptcl)
                assigned_euler(1:2) = fixed_euler(1:2)
                call self%s%b_ptr%spproj_field%set_euler(self%s%iptcl, assigned_euler)
                call self%s%b_ptr%spproj_field%set(self%s%iptcl, 'mi_proj', 1.)
            endif
            call self%s%b_ptr%spproj_field%set(self%s%iptcl, 'frac', frac)
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
