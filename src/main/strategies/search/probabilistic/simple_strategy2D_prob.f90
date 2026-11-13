!@descr: 2D strategy for probabilistic class assignment (the class drawn from the particle's row by the matcher's pipeline)
module simple_strategy2D_prob
use simple_pftc_srch_api
use simple_polarft_calc,     only: SHERRSQ
use simple_eul_prob_tab_utils, only: record_shift_state, eulprob_corr_switch, posterior_temperature, &
    &temperature_or_argmin, prob_rnd_scope, inpl_quadrature
use simple_rnd,              only: rnd_keyed, rnd_draw_key, RND_FAM_INPLQ
use simple_strategy2D_alloc, only: s2D
use simple_parameters,       only: parameters
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_oris,             only: oris
use simple_builder,          only: builder
implicit none

public :: strategy2D_prob
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_prob
contains
    procedure :: new  => new_prob
    procedure :: srch => srch_prob
    procedure :: kill => kill_prob
end type strategy2D_prob

contains

    subroutine new_prob( self, params, spec, build )
        class(strategy2D_prob), intent(inout) :: self
        class(parameters),      intent(in)    :: params
        class(strategy2D_spec), intent(inout) :: spec
        class(builder),         intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_prob

    subroutine srch_prob( self, os )
        class(strategy2D_prob), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer :: iptcl_map, icls, inpl
        real    :: corr, frac
        if( os%get_state(self%s%iptcl) > 0 )then
            if( .not. associated(self%spec%eulprob_obj_part2D) ) THROW_HARD('strategy2D_prob requires eulprob_obj_part2D')
            call self%s%prep4srch(os)
            self%s%nrefs_eval = self%s%nrefs
            iptcl_map = self%s%iptcl_map
            icls      = self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%icls
            inpl      = self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%inpl
            frac      = self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%frac
            if( frac <= 0. ) frac = 100.
            self%s%nrefs_eval = max(1, min(self%s%nrefs, nint(real(self%s%nrefs) * frac / 100.)))
            corr      = eulprob_corr_switch(self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%dist, self%s%p_ptr%cc_objfun)
            self%s%best_class = icls
            self%s%best_rot   = inpl
            self%s%best_corr  = corr
            self%s%best_shvec = [0.,0.]
            if( self%s%p_ptr%l_doshift .and. self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%has_sh )then
                self%s%best_shvec = [self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%x, &
                                  self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%y]
            endif
            call self%s%refine_selected_continuously
            ! the dynamic centring's state: this iteration's move and the particle's prior width for its next search
            if( self%s%p_ptr%l_centring ) call record_shift_state(self%s%b_ptr%spproj_field, self%s%iptcl,&
                &norm2(self%s%best_shvec), self%spec%eulprob_obj_part2D%assgn_map(iptcl_map)%sdev, self%s%p_ptr%shift_sigma,&
                &self%s%p_ptr%smpd, self%s%p_ptr%box, self%s%p_ptr%trs, self%s%b_ptr%pftc%get_kfromto())
            call self%s%store_solution(self%s%best_class, self%s%best_rot, self%s%best_corr)
            call self%s%assign_ori(os)
            call os%set(self%s%iptcl, 'frac', frac)
            if( associated(self%spec%inpl_ins) ) call insert_inplane(self, os)
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_prob

    !> The in-plane insertion of the committed particle (topk_inpl = K > 1): the row of its committed class over every
    !! rotation at its committed shift (the particle shifted once in its own frame), the K-point quadrature of that
    !! row's posterior at the particle's own temperature with the particle's keyed phase, the entries into the
    !! batch's side channel. The committed rotation's entries sit at the committed (polished) angle; a quadrature in
    !! one rotation leaves the particle at the committed pose, weight 1. The pose in the field is untouched
    subroutine insert_inplane( self, os )
        class(strategy2D_prob), intent(inout) :: self
        class(oris),            intent(in)    :: os
        real    :: dists(self%s%nrots), e3(self%s%p_ptr%topk_inpl), w(self%s%p_ptr%topk_inpl)
        real    :: shvec(2), e3c, temp, neff
        real(dp) :: v
        integer :: inds(self%s%p_ptr%topk_inpl), nent, icls, irot_c, iptcl, j, ithr
        iptcl  = self%s%iptcl
        icls   = nint(os%get(iptcl, 'class'))
        irot_c = nint(os%get(iptcl, 'inpl'))
        ! the committed pose is the polished pose of the class the search committed; another class or rotation in
        ! the field (a fallback to the previous class) leaves the particle at its pose alone
        if( icls /= self%s%best_class .or. irot_c /= self%s%best_rot ) return
        e3c   = os%e3get(iptcl)
        shvec = self%s%best_shvec
        ithr  = self%s%ithr
        if( sum(shvec**2) > SHERRSQ )then
            call self%s%b_ptr%pftc%memoize_shifted_ptcl(iptcl, shvec, s2D%inpl_bank(:,:,ithr))
            call self%s%b_ptr%pftc%gen_euclid_dists_bank(icls, iptcl, s2D%inpl_bank(:,:,ithr), dists)
        else
            call self%s%b_ptr%pftc%gen_euclid_dists(icls, iptcl, shvec, dists)
        endif
        temp = temperature_or_argmin(posterior_temperature(self%s%b_ptr%pftc%wsqsum_ptcl_now(iptcl), &
            &self%s%b_ptr%pftc%get_pftsz()))
        v    = rnd_keyed(rnd_draw_key(prob_rnd_scope(self%s%p_ptr%which_iter, .true.), iptcl, RND_FAM_INPLQ))
        call inpl_quadrature(self%s%nrots, dists, temp, self%s%p_ptr%topk_inpl, v, nent, inds, w, neff)
        if( nent < 1 ) return
        if( nent == 1 )then
            call self%spec%inpl_ins%set_single(self%spec%iptcl_batch, neff)
            return
        endif
        do j = 1,nent
            if( inds(j) == irot_c )then
                e3(j) = e3c
            else
                e3(j) = 360. - self%s%b_ptr%pftc%get_rot(inds(j))
            endif
        enddo
        call self%spec%inpl_ins%set_spread(self%spec%iptcl_batch, nent, e3(1:nent), w(1:nent), neff)
    end subroutine insert_inplane

    subroutine kill_prob( self )
        class(strategy2D_prob), intent(inout) :: self
        call self%s%kill
    end subroutine kill_prob

end module simple_strategy2D_prob
