!@descr: 3D strategy for probabilistic projection matching
module simple_strategy3D_prob
use simple_core_module_api
use simple_strategy3D_utils,   only: assign_ori, candidate_pose, resolve_inplane_e3
use simple_eul_prob_tab_utils, only: prob_insert_pose, record_shift_state, merge_view_basins, temperature_or_argmin
use simple_parameters,         only: parameters
use simple_strategy3D_pftc,    only: strategy3D_pftc
use simple_strategy3D_srch,    only: strategy3D_spec, view_solve
use simple_ori_utils,          only: euler2m
use simple_ori,                only: ori
use simple_oris,               only: oris
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
    procedure, private :: resolve_extra_poses_view
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
        real    :: corr, frac, sh(2), inpl_coord, sh_rec(2), sdev_rec
        logical :: inpl_valid, l_view_extras
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
            ! the table's native increments refer to; they take the continuous in-plane route when the committed pose
            ! does. Under cont_view they follow the committed pose, with the view, as basins (resolve_extra_poses_view)
            l_view_extras = .false.
            if( associated(self%spec%extra_poses) )then
                if( self%spec%extra_poses%has(self%s%iptcl) )then
                    if( self%s%uses_view_refinement() )then
                        l_view_extras = .true.
                    else
                        call self%resolve_extra_poses(iptcl_map, self%s%uses_continuous_refinement())
                    endif
                endif
            endif
            if( self%s%uses_continuous_refinement() )then
                call self%s%refine_assignment_continuously(iref, irot, corr, sh, &
                    &inpl_coord, inpl_valid)
                call assign_ori(self%s, iref, irot, corr, sh, inpl_coord, inpl_valid)
                ! an improving view solve's pose replaces the drawn cell's (cont_view)
                call self%s%commit_view
            else
                call assign_ori(self%s, iref, irot, corr, sh)
            endif
            ! the committed pose's shift increment and posterior width; an extra that replaces the pose replaces them
            sh_rec   = sh
            sdev_rec = self%spec%eulprob_obj_part%assgn_map(iptcl_map)%sdev
            if( l_view_extras ) call self%resolve_extra_poses_view(iptcl_map, sh_rec, sdev_rec)
            ! the dynamic centring's state (shift_sigma > 0): this iteration's move and the particle's prior width for
            ! its next search, from the pose the particle keeps
            if( self%s%p_ptr%l_centring ) call record_shift_state(self%s%b_ptr%spproj_field, self%s%iptcl, norm2(sh_rec),&
                &sdev_rec, self%s%p_ptr%shift_sigma, self%s%p_ptr%smpd,&
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
    !! across the search's threads). With l_cont_refine each extra first takes the continuous in-plane route from its own
    !! cell and shift, as the committed pose does after it; its weight stays the table's
    subroutine resolve_extra_poses( self, iptcl_map, l_cont_refine )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iptcl_map
        logical,                intent(in)    :: l_cont_refine
        type(prob_insert_pose) :: pose
        type(view_solve)       :: vs
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
            if( l_cont_refine )then
                corr = 0. ! re-scored by the solve, only its finiteness is checked
                call self%s%refine_extra_pose_continuously(pose%iref, inpl, corr, sh, inpl_coord, inpl_valid, vs)
            endif
            call candidate_pose(self%s, pose%iref, inpl, sh, euler, shvec)
            if( l_cont_refine ) euler(3) = resolve_inplane_e3(inpl, inpl_coord, inpl_valid,&
                &self%s%uses_continuous_refinement(), self%s%b_ptr%pftc%get_dang())
            self%spec%extra_poses%euler(:,ipose) = euler
            self%spec%extra_poses%shift(:,ipose) = shvec
            self%spec%extra_poses%w(ipose)       = pose%w
        enddo
    end subroutine resolve_extra_poses

    !> The particle's extra poses under cont_view (section 6 of polar_continuous_view.md): each of the table's top
    !! candidates takes the view solve from its own cell and shift; the committed pose (the project's, after commit_view)
    !! and the valid refined extras are merged into basins within the insertion cone of one another (merge_view_basins);
    !! each basin enters at its best pose with its members' posterior weights added, and the committed pose becomes
    !! its basin's best when a refined extra beats it (the same hypothesis at its better optimum); the dynamic
    !! centring's move sh_rec and width sdev_rec then become that extra's (its refined shift increment, its table
    !! entry's posterior width). A committed pose whose view solve was invalid enters alone
    subroutine resolve_extra_poses_view( self, iptcl_map, sh_rec, sdev_rec )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iptcl_map
        real,                   intent(inout) :: sh_rec(2), sdev_rec
        type(prob_insert_pose)        :: pose
        type(view_solve), allocatable :: vs(:)
        type(view_solve)              :: vs_c
        type(ori)                     :: o
        real,    allocatable :: rmats(:,:,:), eulers(:,:), shifts(:,:), scores(:), w(:), shincs(:,:), sdevs(:)
        integer, allocatable :: basin(:)
        logical, allocatable :: valid(:), best(:)
        real    :: euler(3), shvec(2), sh(2), corr, inpl_coord, wc
        integer :: iloc, nextra, k, ipose, inpl, istate, nbasins, kbest, ninserted, nmerged
        logical :: inpl_valid
        iloc   = self%spec%extra_poses%loc(self%s%iptcl)
        nextra = self%spec%extra_poses%nextra(iloc)
        allocate(vs(nextra), rmats(3,3,0:nextra), eulers(3,0:nextra), shifts(2,0:nextra), scores(0:nextra), w(0:nextra), &
            &basin(0:nextra), valid(0:nextra), best(0:nextra), shincs(2,nextra), sdevs(nextra))
        do k = 1,nextra
            self%spec%extra_poses%w(self%spec%extra_poses%first(iloc) + k - 1) = 0.
        enddo
        ! the committed pose as the project holds it after assign_ori and commit_view
        call self%s%get_view_solution(vs_c)
        call self%s%b_ptr%spproj_field%get_ori(self%s%iptcl, o)
        istate        = o%get_state()
        eulers(:,0)   = o%get_euler()
        rmats(:,:,0)  = euler2m(eulers(:,0))
        shifts(:,0)   = o%get_2Dshift()
        scores(0)     = vs_c%score
        valid(0)      = vs_c%evaluated
        do k = 1,nextra
            pose = self%spec%eulprob_obj_part%insert_tab%get_extra(iptcl_map, k)
            inpl = pose%inpl
            sh   = [pose%x, pose%y]
            corr = 0. ! re-scored by the solve, only its finiteness is checked
            call self%s%refine_extra_pose_continuously(pose%iref, inpl, corr, sh, inpl_coord, inpl_valid, vs(k))
            call candidate_pose(self%s, pose%iref, inpl, sh, euler, shvec)
            euler(3) = resolve_inplane_e3(inpl, inpl_coord, inpl_valid, .true., self%s%b_ptr%pftc%get_dang())
            if( vs(k)%committed ) euler = self%s%view_pose(euler, vs(k))
            eulers(:,k)  = euler
            rmats(:,:,k) = euler2m(euler)
            shifts(:,k)  = shvec
            ! the increment the stored shift took (candidate_pose: none without a shift search)
            shincs(:,k)  = 0.
            if( self%s%doshift ) shincs(:,k) = sh
            sdevs(k)     = pose%sdev
            scores(k)    = vs(k)%score
            valid(k)     = vs(k)%evaluated
        enddo
        call merge_view_basins(self%spec%eulprob_obj_part%insert_cone, istate, &
            &temperature_or_argmin(self%spec%eulprob_obj_part%temps(iptcl_map)), nextra, rmats, shifts, scores, valid, &
            &PROB_INSERT_WMIN, basin, nbasins, w, best)
        wc        = 1.
        ninserted = 0
        nmerged   = 0
        if( nbasins > 0 )then
            ! the committed basin: its weight, and its best member's pose when an extra beats the committed one
            kbest = 0
            do k = 1,nextra
                if( basin(k) == basin(0) .and. best(k) ) kbest = k
            enddo
            if( kbest > 0 )then
                wc = w(kbest)
                call self%s%b_ptr%spproj_field%set_euler(self%s%iptcl, eulers(:,kbest))
                call self%s%b_ptr%spproj_field%set_shift(self%s%iptcl, shifts(:,kbest))
                call self%s%b_ptr%spproj_field%set(self%s%iptcl, 'corr', scores(kbest))
                call o%set_euler(eulers(:,kbest))
                call self%s%b_ptr%spproj_field%set(self%s%iptcl, 'proj', real(self%s%b_ptr%eulspace%find_closest_proj(o)))
                call self%s%set_view(vs(kbest))
                sh_rec   = shincs(:,kbest)
                sdev_rec = sdevs(kbest)
            else
                wc = w(0)
            endif
            ! the other basins at their best poses, the extras merged into a basin counted
            do k = 1,nextra
                if( basin(k) == 0 ) cycle
                if( basin(k) == basin(0) .or. .not. best(k) )then
                    nmerged = nmerged + 1
                    cycle
                endif
                ninserted = ninserted + 1
                ipose = self%spec%extra_poses%first(iloc) + ninserted - 1
                self%spec%extra_poses%euler(:,ipose) = eulers(:,k)
                self%spec%extra_poses%shift(:,ipose) = shifts(:,k)
                self%spec%extra_poses%w(ipose)       = w(k)
            enddo
        endif
        self%spec%extra_poses%wc(iloc) = wc
        call self%s%set_view_extras_counts(count(valid(1:nextra)), ninserted, nmerged, wc)
        call o%kill
    end subroutine resolve_extra_poses_view

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
