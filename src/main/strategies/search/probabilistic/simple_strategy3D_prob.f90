!@descr: 3D strategy for probabilistic projection matching
module simple_strategy3D_prob
use simple_core_module_api
use simple_strategy3D_utils, only: assign_ori, candidate_pose, resolve_inplane_e3
use simple_eul_prob_tab_utils, only: prob_insert_pose, record_shift_state, prob_inpl_counts, posterior_temperature, &
    &temperature_or_argmin, prob_rnd_scope, inpl_quadrature, inpl_view_weights, inpl_prior_penalty, inpl_entry_pack, &
    &inpl_particle_counts
use simple_rnd,              only: rnd_keyed, rnd_draw_key, RND_FAM_INPLQ
use simple_polarft_calc,     only: SHERRSQ
use simple_strategy3D_alloc, only: s3D
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
    type(prob_inpl_counts) :: inpl_counts  !< the in-plane posterior insertion's counts of this particle
contains
    procedure :: new         => new_prob
    procedure :: srch        => srch_prob
    procedure :: kill        => kill_prob
    procedure :: oris_assign => oris_assign_prob
    procedure :: has_valid_assignment
    procedure :: add_inpl_counts
    procedure, private :: resolve_extra_poses
    procedure, private :: inpl_active
    procedure, private :: inpl_pose_row
    procedure, private :: insert_inpl_posterior
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
        self%inpl_counts      = prob_inpl_counts()
    end subroutine new_prob

    subroutine srch_prob( self, os, ithr )
        use simple_eul_prob_tab_utils, only: eulprob_corr_switch
        class(strategy3D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer,                intent(in)    :: ithr
        integer :: iproj, iptcl_map, irot, istate, iref
        real    :: corr, frac, sh(2), inpl_coord, euler(3)
        logical :: inpl_valid, l_inpl
        self%assignment_valid = .false.
        ! the in-plane posterior insertion: one particle of the table, its committed view weight 1 unless reweighted
        l_inpl = self%inpl_active()
        if( l_inpl ) self%inpl_counts%view_weight_sum = 1.d0
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
                    &self%s%uses_continuous_refinement(), l_inpl)
            endif
            if( self%s%uses_continuous_refinement() )then
                call self%s%refine_assignment_continuously(iref, irot, corr, sh, &
                    &inpl_coord, inpl_valid)
                call assign_ori(self%s, iref, irot, corr, sh, inpl_coord, inpl_valid, euler_out=euler)
            else
                call assign_ori(self%s, iref, irot, corr, sh, euler_out=euler)
            endif
            if( l_inpl ) call self%insert_inpl_posterior(iptcl_map, iref, irot, sh, euler)
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
    !! across the search's threads). With l_cont_refine each extra first takes the continuous in-plane route from its own
    !! cell and shift, as the committed pose does after it; its weight stays the table's (l_inpl: reweighed after)
    subroutine resolve_extra_poses( self, iptcl_map, l_cont_refine, l_inpl )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iptcl_map
        logical,                intent(in)    :: l_cont_refine
        logical,                intent(in)    :: l_inpl !< the in-plane posterior insertion: each extra's row at its solve
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
            if( l_cont_refine )then
                corr = 0. ! re-scored by the solve, only its finiteness is checked
                call self%s%refine_extra_pose_continuously(pose%iref, inpl, corr, sh, inpl_coord, inpl_valid)
            endif
            call candidate_pose(self%s, pose%iref, inpl, sh, euler, shvec)
            if( l_cont_refine ) euler(3) = resolve_inplane_e3(inpl, inpl_coord, inpl_valid,&
                &self%s%uses_continuous_refinement(), self%s%b_ptr%pftc%get_dang())
            self%spec%extra_poses%euler(:,ipose) = euler
            self%spec%extra_poses%shift(:,ipose) = shvec
            self%spec%extra_poses%w(ipose)       = pose%w
            if( l_inpl ) call self%inpl_pose_row(pose%iref, inpl, sh, [pose%cx, pose%cy], iloc, k)
        enddo
    end subroutine resolve_extra_poses

    !> whether the in-plane posterior insertion runs for this particle: the key, and the part's blocks reserved
    logical function inpl_active( self ) result( l_inpl )
        class(strategy3D_prob), intent(in) :: self
        l_inpl = .false.
        if( .not. self%s%p_ptr%l_inpl_posterior ) return
        if( .not. associated(self%spec%extra_poses) ) return
        if( .not. allocated(self%spec%extra_poses%blk_nent) ) return
        l_inpl = self%spec%extra_poses%has(self%s%iptcl)
    end function inpl_active

    !> Pose k (0 the committed one) of the particle at local index iloc: its in-plane row at the solved increment sh (the
    !! stored frame) and bin, at the physical temperature; the quadrature, keyed on the particle and the reference, into
    !! the pose's block, and the row's marginal terms and the shift prior's term into the thread's work
    subroutine inpl_pose_row( self, iref, bin, sh, centre, iloc, k )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iref, bin, iloc, k
        real,                   intent(in)    :: sh(2), centre(2)
        real     :: dists(self%s%nrots), w(INPL_POSTERIOR_BUDGET), sh_row(2), temp, sigma2_window
        real(dp) :: v
        integer  :: inds(INPL_POSTERIOR_BUDGET), nent, j, iblk, ithr, iptcl
        iptcl  = self%s%iptcl
        ithr   = self%s%ithr
        sh_row = 0.
        if( self%s%doshift ) sh_row = sh
        if( sum(sh_row**2) > SHERRSQ )then
            call self%s%b_ptr%pftc%memoize_shifted_ptcl(iptcl, sh_row, s3D%inpl_bank(:,:,ithr))
            call self%s%b_ptr%pftc%gen_euclid_dists_bank(iref, iptcl, s3D%inpl_bank(:,:,ithr), dists)
        else
            call self%s%b_ptr%pftc%gen_euclid_dists(iref, iptcl, sh_row, dists)
        endif
        temp = temperature_or_argmin(posterior_temperature(self%s%b_ptr%pftc%wsqsum_ptcl_now(iptcl), &
            &self%s%b_ptr%pftc%get_pftsz()))
        v    = rnd_keyed(rnd_draw_key(prob_rnd_scope(self%s%p_ptr%which_iter, .false.), iptcl, RND_FAM_INPLQ, &
            &reference=iref))
        call inpl_quadrature(self%s%nrots, dists, temp, INPL_POSTERIOR_BUDGET, v, nent, inds, w, &
            &s3D%inpl_neff(k+1,ithr), s3D%inpl_lnz(k+1,ithr), s3D%inpl_dmin(k+1,ithr))
        sigma2_window = self%spec%eulprob_obj_part%get_shift_prior_sigma2()
        s3D%inpl_pen(k+1,ithr)  = inpl_prior_penalty(sh_row, centre, self%s%p_ptr%l_centring, &
            &self%s%p_ptr%shift_sigma, self%s%p_ptr%smpd, sigma2_window)
        s3D%inpl_nent(k+1,ithr) = nent
        iblk = self%spec%extra_poses%block(iloc, k)
        self%spec%extra_poses%blk_nent(iblk) = nent
        self%spec%extra_poses%blk_bin(iblk)  = bin
        do j = 1,nent
            self%spec%extra_poses%blk_ent(j,iblk) = inpl_entry_pack(inds(j), nint(w(j) * real(INPL_POSTERIOR_BUDGET)))
        enddo
    end subroutine inpl_pose_row

    !> The committed pose's row after its solve, then the view weights of all the particle's poses (inpl_view_weights) into
    !! the slots the insertion reads; euler is the committed pose about its direction (assign_ori). A pose without a
    !! valid row leaves the whole particle to the insertion without the key (its blocks emptied, the table's weights kept)
    subroutine insert_inpl_posterior( self, iptcl_map, iref, bin, sh, euler )
        class(strategy3D_prob), intent(inout) :: self
        integer,                intent(in)    :: iptcl_map, iref, bin
        real,                   intent(in)    :: sh(2), euler(3)
        integer :: iloc, npose, k, ithr
        real    :: temp
        iloc  = self%spec%extra_poses%loc(self%s%iptcl)
        npose = 1 + self%spec%extra_poses%nextra(iloc)
        ithr  = self%s%ithr
        self%spec%extra_poses%committed_euler(:,iloc) = euler
        call self%inpl_pose_row(iref, bin, sh, self%spec%eulprob_obj_part%get_prior_centre(iptcl_map), iloc, 0)
        if( any(s3D%inpl_nent(1:npose,ithr) < 1) )then
            do k = 0,npose - 1
                self%spec%extra_poses%blk_nent(self%spec%extra_poses%block(iloc, k)) = 0
            enddo
            self%inpl_counts%view_weight_sum = real(self%spec%extra_poses%wc(iloc), dp)
            return
        endif
        temp = temperature_or_argmin(posterior_temperature(self%s%b_ptr%pftc%wsqsum_ptcl_now(self%s%iptcl), &
            &self%s%b_ptr%pftc%get_pftsz()))
        call inpl_view_weights(npose, s3D%inpl_lnz(1:npose,ithr), s3D%inpl_dmin(1:npose,ithr), temp, &
            &s3D%inpl_pen(1:npose,ithr), PROB_INSERT_WMIN, s3D%inpl_m(1:npose,ithr))
        self%spec%extra_poses%wc(iloc) = s3D%inpl_m(1,ithr)
        do k = 1,npose - 1
            self%spec%extra_poses%w(self%spec%extra_poses%first(iloc) + k - 1) = s3D%inpl_m(k+1,ithr)
        enddo
        self%inpl_counts = inpl_particle_counts(self%spec%extra_poses, iloc, s3D%inpl_neff(1:npose,ithr))
    end subroutine insert_inpl_posterior

    !> this particle's in-plane posterior counts into its slot of the pass (the matcher sums them in particle order)
    subroutine add_inpl_counts( self, counts )
        class(strategy3D_prob), intent(in)    :: self
        type(prob_inpl_counts), intent(inout) :: counts
        call counts%add(self%inpl_counts)
    end subroutine add_inpl_counts

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
