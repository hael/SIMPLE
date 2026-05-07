!@descr: 2D strategy for in-plane refinement with probabilistic sampling
module simple_strategy2D_inpl_smpl
use simple_pftc_srch_api
use simple_strategy2D_alloc, only: s2D
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: builder
use simple_type_defs,        only: OBJFUN_EUCLID
implicit none

public :: strategy2D_inpl_smpl
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_inpl_smpl
    contains
    procedure :: new  => new_inpl_smpl
    procedure :: srch => srch_inpl_smpl
    procedure :: kill => kill_inpl_smpl
end type strategy2D_inpl_smpl

contains

    subroutine new_inpl_smpl( self, params, spec, build )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        class(parameters),           intent(in)    :: params
        class(strategy2D_spec),      intent(inout) :: spec
        class(builder),              intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_inpl_smpl

    subroutine srch_inpl_smpl( self, os )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        class(oris),                 intent(inout) :: os
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots), offsets(2,self%s%nrots), cxy(3)
        real    :: inpl_dist
        integer :: sorted_inds(self%s%nrots), isample, inpl_ind, order_ind
        logical :: l_prob_objfun
        if( os%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch(os)
            l_prob_objfun = (self%s%p_ptr%cc_objfun == OBJFUN_EUCLID)
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! In-plane sampling
            if( .not. (l_prob_objfun .and. .not. s2D%do_inplsrch(self%s%iptcl_batch)) )then
                if( self%s%l_sh_first )then
                    call self%s%b_ptr%pftc%gen_objfun_vals(self%s%best_class, self%s%iptcl, self%s%xy_first, inpl_corrs)
                else
                    call self%s%b_ptr%pftc%gen_objfun_vals(self%s%best_class, self%s%iptcl, [0.,0.],         inpl_corrs)
                endif
            endif
            ! Shift search
            if( s2D%do_inplsrch(self%s%iptcl_batch) )then
                ! Shift search for top-ranking subset
                call self%s%grad_shsrch_obj2%set_indices(self%s%best_class, self%s%iptcl)
                sorted_corrs = inpl_corrs
                sorted_inds  = (/(inpl_ind,inpl_ind=1,self%s%nrots)/)
                call hpsort(sorted_corrs, sorted_inds)
                inpl_corrs = -1.
                offsets   = 0.
                do isample = self%s%nrots-s2D%snhc_smpl_ninpl+1,self%s%nrots
                    inpl_ind = sorted_inds(isample)
                    if( self%s%l_sh_first )then
                        cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%s%xy_first)
                        if( inpl_ind == 0 )then
                            inpl_ind = sorted_inds(isample)
                            cxy(1)   = real(self%s%b_ptr%pftc%gen_corr_for_rot_8(self%s%best_class, self%s%iptcl, real(self%s%xy_first,dp), inpl_ind))
                            cxy(2:3) = self%s%xy_first_rot
                        endif
                    else
                        cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                        if( inpl_ind == 0 )then
                            inpl_ind = sorted_inds(isample)
                            cxy      = [real(self%s%b_ptr%pftc%gen_corr_for_rot_8(self%s%best_class, self%s%iptcl, inpl_ind)), 0.,0.]
                        endif
                    endif
                    inpl_corrs(inpl_ind) = cxy(1)
                    offsets(:,inpl_ind) = cxy(2:3)
                enddo
                ! In-plane selection
                call power_sampling( s2D%power, self%s%nrots, inpl_corrs, sorted_inds, s2D%snhc_smpl_ninpl,&
                                    &self%s%best_rot, order_ind, self%s%best_corr )
                self%s%best_shvec = offsets(:,self%s%best_rot)
            else
                ! In-plane selection
                if( l_prob_objfun )then
                    if( self%s%l_sh_first )then
                        call self%s%b_ptr%pftc%gen_prob_power_objfun_val(self%s%best_class, self%s%iptcl,&
                            &self%s%xy_first, s2D%power, s2D%snhc_smpl_ninpl, inpl_dist, self%s%best_corr,&
                            &self%s%best_rot, inpl_corrs, sorted_inds)
                    else
                        call self%s%b_ptr%pftc%gen_prob_power_objfun_val(self%s%best_class, self%s%iptcl,&
                            &[0.,0.], s2D%power, s2D%snhc_smpl_ninpl, inpl_dist, self%s%best_corr,&
                            &self%s%best_rot, inpl_corrs, sorted_inds)
                    endif
                else
                    call power_sampling( s2D%power, self%s%nrots, inpl_corrs, sorted_inds, s2D%snhc_smpl_ninpl,&
                                        &self%s%best_rot, order_ind, self%s%best_corr )
                endif
                self%s%best_shvec = 0.
            endif
            self%s%nrefs_eval = self%s%nrefs
            ! Updates solution
            call self%s%store_solution(self%s%best_class, self%s%best_rot, self%s%best_corr)
            call self%s%assign_ori(os)
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_inpl_smpl

    subroutine kill_inpl_smpl( self )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_inpl_smpl

end module simple_strategy2D_inpl_smpl
