!@descr: 3D strategy for probabilistic projection matching
module simple_strategy3D_prob
use simple_core_module_api
use simple_strategy3D_utils, only: assign_ori, candidate_pose, resolve_inplane_e3
use simple_eul_prob_tab_utils, only: prob_insert_pose, record_shift_state
use simple_parameters,      only: parameters
use simple_strategy3D_pftc,  only: strategy3D_pftc
use simple_strategy3D_srch, only: strategy3D_spec
use simple_oris,            only: oris
implicit none

public :: strategy3D_prob
private

#include "simple_local_flags.inc"

type, extends(strategy3D_pftc) :: strategy3D_prob
    logical :: assignment_valid = .false.
contains
    procedure :: new         => new_prob
    procedure :: srch        => srch_prob
    procedure :: kill        => kill_prob
    procedure :: oris_assign => oris_assign_prob
    procedure :: has_valid_assignment
    procedure, private :: resolve_extra_poses
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
        self%assignment_valid = .false.
    end subroutine new_prob

    subroutine srch_prob( self, os, ithr )
        use simple_eul_prob_tab_utils, only: eulprob_corr_switch
        class(strategy3D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer,                intent(in)    :: ithr
        integer :: iproj, iptcl_map, irot, istate, iref
        real    :: corr, frac, sh(2), inpl_coord
        logical :: inpl_valid
        self%assignment_valid = .false.
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
            ! the weighted insertion's extra poses of this particle (topk > 1), now, while prev_shvec is the shift
            ! the table's native increments refer to; they take the continuous in-plane route when the committed pose does
            if( associated(self%spec%extra_poses) )then
                if( self%spec%extra_poses%has(self%s%iptcl) ) call self%resolve_extra_poses(iptcl_map,&
                    &self%s%uses_continuous_refinement())
            endif
            if( self%s%uses_continuous_refinement() )then
                call self%s%refine_assignment_continuously(iref, irot, corr, sh, &
                    &inpl_coord, inpl_valid)
                call assign_ori(self%s, iref, irot, corr, sh, inpl_coord, inpl_valid)
            else
                call assign_ori(self%s, iref, irot, corr, sh)
            endif
            ! the dynamic centring's state (shift_sigma > 0): this iteration's move and the particle's prior width for
            ! its next search
            if( self%s%p_ptr%l_centring ) call record_shift_state(self%s%b_ptr%spproj_field, self%s%iptcl, norm2(sh),&
                &self%spec%eulprob_obj_part%assgn_map(iptcl_map)%sdev, self%s%p_ptr%shift_sigma, self%s%p_ptr%smpd,&
                &self%s%p_ptr%box, self%s%p_ptr%trs, self%s%b_ptr%pftc%get_kfromto())
            call self%s%b_ptr%spproj_field%set(self%s%iptcl, 'frac', frac)
            self%assignment_valid = .true.
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_prob

    !> The particle's extra poses for the weighted insertion: the part's table entries (reference, in-plane index,
    !! native shift increment, weight) as Euler triplets and absolute stored shifts through candidate_pose, the
    !! committed pose's route, into the particle's slots of the part's structure (disjoint slots per particle: safe
    !! across the search's threads). With l_polish each extra first takes the continuous in-plane route from its own
    !! cell and shift, as the committed pose does after it; its weight stays the table's
    subroutine resolve_extra_poses( self, iptcl_map, l_polish )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iptcl_map
        logical,                intent(in)    :: l_polish
        type(prob_insert_pose) :: pose
        integer :: iloc, k, ipose, inpl
        real    :: euler(3), shvec(2), sh(2), corr, inpl_coord
        logical :: inpl_valid
        iloc = self%spec%extra_poses%loc(self%s%iptcl)
        do k = 1,self%spec%extra_poses%nextra(iloc)
            pose       = self%spec%eulprob_obj_part%insert_tab%get_extra(iptcl_map, k)
            ipose      = self%spec%extra_poses%first(iloc) + k - 1
            inpl       = pose%inpl
            sh         = [pose%x, pose%y]
            inpl_coord = real(inpl)
            inpl_valid = .false.
            if( l_polish )then
                corr = 0. ! re-scored by the solve, only its finiteness is checked
                call self%s%refine_extra_pose_continuously(pose%iref, inpl, corr, sh, inpl_coord, inpl_valid)
            endif
            call candidate_pose(self%s, pose%iref, inpl, sh, euler, shvec)
            if( l_polish ) euler(3) = resolve_inplane_e3(inpl, inpl_coord, inpl_valid,&
                &self%s%uses_continuous_refinement(), self%s%b_ptr%pftc%get_dang())
            self%spec%extra_poses%euler(:,ipose) = euler
            self%spec%extra_poses%shift(:,ipose) = shvec
            self%spec%extra_poses%w(ipose)       = pose%w
        enddo
    end subroutine resolve_extra_poses

    !> True only when this search committed a valid probabilistic assignment.
    !! A missing assignment must not seed a later local pose refinement from
    !! the particle's stale incoming orientation.
    pure logical function has_valid_assignment(self) result(valid)
        class(strategy3D_prob), intent(in) :: self
        valid = self%assignment_valid
    end function has_valid_assignment

    subroutine oris_assign_prob( self )
        class(strategy3D_prob), intent(inout) :: self
    end subroutine oris_assign_prob

    subroutine kill_prob( self )
        class(strategy3D_prob), intent(inout) :: self
        call self%s%kill
        self%assignment_valid = .false.
    end subroutine kill_prob

end module simple_strategy3D_prob
