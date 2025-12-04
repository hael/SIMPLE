module simple_strategy2D_inpl_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_parameters,       only: params_glob
use simple_eul_prob_tab2D,   only: power_sampling
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

    subroutine new_inpl_smpl( self, spec )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_inpl_smpl

    subroutine srch_inpl_smpl( self )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        real    :: inpl_corrs(self%s%nrots), sorted_corrs(self%s%nrots), offsets(2,self%s%nrots), cxy(3)
        integer :: sorted_inds(self%s%nrots), isample, inpl_ind, order_ind
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! In-plane sampling
            if( self%s%l_sh_first )then
                call pftc_glob%gen_objfun_vals(self%s%best_class, self%s%iptcl, self%s%xy_first, inpl_corrs)
            else
                call pftc_glob%gen_objfun_vals(self%s%best_class, self%s%iptcl, [0.,0.],         inpl_corrs)
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
                            cxy(1)   = real(pftc_glob%gen_corr_for_rot_8(self%s%best_class, self%s%iptcl, real(self%s%xy_first,dp), inpl_ind))
                            cxy(2:3) = self%s%xy_first_rot
                        endif
                    else
                        cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                        if( inpl_ind == 0 )then
                            inpl_ind = sorted_inds(isample)
                            cxy      = [real(pftc_glob%gen_corr_for_rot_8(self%s%best_class, self%s%iptcl, inpl_ind)), 0.,0.]
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
                call power_sampling( s2D%power, self%s%nrots, inpl_corrs, sorted_inds, s2D%snhc_smpl_ninpl,&
                                    &self%s%best_rot, order_ind, self%s%best_corr )
                self%s%best_shvec = 0.
            endif
            self%s%nrefs_eval = self%s%nrefs
            ! Updates solution
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_inpl_smpl

    subroutine kill_inpl_smpl( self )
        class(strategy2D_inpl_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_inpl_smpl

end module simple_strategy2D_inpl_smpl
