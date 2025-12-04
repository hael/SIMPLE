module simple_strategy2D_greedy_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_parameters,       only: params_glob
use simple_eul_prob_tab2D,   only: squared_sampling

implicit none

public :: strategy2D_greedy_smpl
private

#include "simple_local_flags.inc"

type, extends(strategy2D) :: strategy2D_greedy_smpl
  contains
    procedure :: new  => new_greedy_smpl
    procedure :: srch => srch_greedy_smpl
    procedure :: kill => kill_greedy_smpl
end type strategy2D_greedy_smpl

contains

    subroutine new_greedy_smpl( self, spec )
        class(strategy2D_greedy_smpl), intent(inout) :: self
        class(strategy2D_spec),   intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_greedy_smpl

    subroutine srch_greedy_smpl( self )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy2D_greedy_smpl), intent(inout) :: self
        integer :: refs_inds(self%s%nrefs), refs_inplinds(self%s%nrefs), inds(self%s%nrots)
        integer :: iref, inpl_ind, isample, order_ind
        real    :: refs_corrs(self%s%nrefs), inpl_corrs(self%s%nrots), cxy(3)
        real    :: inpl_corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! Shift search on previous best reference
            call self%s%inpl_srch_first
            ! Class search
            do iref = 1,self%s%nrefs
                refs_inds(iref) = iref
                if( s2D%cls_pops(iref) == 0 )then
                    refs_corrs(iref)    = -1.
                    refs_inplinds(iref) = 0
                else
                    ! in-plane sampling
                    if( self%s%l_sh_first )then
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                    endif
                    call squared_sampling(self%s%nrots, inpl_corrs, inds,&
                                        &s2D%smpl_ninpl, inpl_ind, order_ind, inpl_corr)
                    refs_inplinds(iref) = inpl_ind
                    refs_corrs(iref)    = inpl_corr
                endif
            enddo
            self%s%best_class = maxloc(refs_corrs,dim=1)
            self%s%best_corr  = refs_corrs(self%s%best_class)
            self%s%best_rot   = refs_inplinds(self%s%best_class)
            ! Performs shift search for top scoring subset
            if( s2D%do_inplsrch(self%s%iptcl_batch) )then
                call hpsort(refs_corrs, refs_inds)
                self%s%best_corr = -1.
                do isample = self%s%nrefs-s2D%smpl_ncls+1,self%s%nrefs
                    iref = refs_inds(isample)
                    if( s2D%cls_pops(iref) == 0 ) cycle
                    call self%s%grad_shsrch_obj2%set_indices(iref, self%s%iptcl)
                    if( self%s%grad_shsrch_obj2%does_opt_angle() )then
                        if( self%s%l_sh_first )then
                            cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%s%xy_first)
                            if( inpl_ind == 0 )then
                                inpl_ind = refs_inplinds(iref)
                                cxy      = [refs_corrs(isample), self%s%xy_first_rot(1), self%s%xy_first_rot(2)]
                            endif
                        else
                            cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                            if( inpl_ind == 0 )then
                                inpl_ind = refs_inplinds(iref)
                                cxy      = [refs_corrs(isample), 0., 0.]
                            endif
                        endif
                    else
                        inpl_ind = refs_inplinds(iref)
                        if( self%s%l_sh_first )then
                            cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%s%xy_first)
                            if( inpl_ind == 0 )then
                                inpl_ind = refs_inplinds(iref)
                                cxy(1)   = real(pftc_glob%gen_corr_for_rot_8(iref, self%s%iptcl, real(self%s%xy_first,dp), inpl_ind))
                                cxy(2:3) = self%s%xy_first_rot
                            endif
                        else
                            cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                            if( inpl_ind == 0 )then
                                inpl_ind = refs_inplinds(iref)
                                cxy      = [real(pftc_glob%gen_corr_for_rot_8(iref, self%s%iptcl, inpl_ind)), 0.,0.]
                            endif
                        endif
                    endif
                    ! Select best score
                    if( cxy(1) > self%s%best_corr )then
                        self%s%best_class = iref
                        self%s%best_corr  = cxy(1)
                        self%s%best_rot   = inpl_ind
                        self%s%best_shvec = cxy(2:3)
                    endif
                enddo
            endif
            call self%s%inpl_srch
            self%s%nrefs_eval = self%s%nrefs
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedy_smpl

    subroutine kill_greedy_smpl( self )
        class(strategy2D_greedy_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_smpl

end module simple_strategy2D_greedy_smpl
