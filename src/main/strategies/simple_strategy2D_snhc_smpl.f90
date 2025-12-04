module simple_strategy2D_snhc_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_eul_prob_tab2D,   only: power_sampling
implicit none

public :: strategy2D_snhc_smpl
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_snhc_smpl
    contains
    procedure :: new  => new_snhc_smpl
    procedure :: srch => srch_snhc_smpl
    procedure :: kill => kill_snhc_smpl
end type strategy2D_snhc_smpl

contains

    subroutine new_snhc_smpl( self, spec )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_snhc_smpl

    subroutine srch_snhc_smpl( self )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        real    :: inpl_corrs(self%s%nrots), sorted_cls_corrs(self%s%nrefs), cls_corrs(self%s%nrefs)
        real    :: cxy(3), inpl_corr
        integer :: cls_inpl_inds(self%s%nrefs), vec_nrots(self%s%nrots), sorted_cls_inds(self%s%nrefs)
        integer :: iref, isample, inpl_ind, order_ind
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! Class search
            cls_corrs     = -1.
            cls_inpl_inds = 0
            do isample = 1,self%s%nrefs
                ! stochastic reference index
                iref = s2D%srch_order(self%s%iptcl_batch, isample)
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! neighbourhood size
                if(self%s%nrefs_eval > s2D%snhc_nrefs_bound) exit
                if( s2D%cls_pops(iref) == 0 )cycle
                ! In-plane sampling
                if( self%s%l_sh_first )then
                    call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                else
                    call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                endif
                call power_sampling( s2D%power, self%s%nrots, inpl_corrs, vec_nrots,&
                                    &s2D%snhc_smpl_ninpl, inpl_ind, order_ind, inpl_corr )
                cls_corrs(iref)     = inpl_corr
                cls_inpl_inds(iref) = inpl_ind
            end do
            ! Performs shift search for top scoring subset
            if( s2D%do_inplsrch(self%s%iptcl_batch) )then
                sorted_cls_corrs = cls_corrs
                sorted_cls_inds  = (/(iref,iref=1,self%s%nrefs)/)
                call hpsort(sorted_cls_corrs, sorted_cls_inds)
                cls_corrs = -1.
                do isample = self%s%nrefs-s2D%snhc_smpl_ncls+1,self%s%nrefs
                    iref = sorted_cls_inds(isample)
                    if( s2D%cls_pops(iref) == 0 ) cycle
                    call self%s%grad_shsrch_obj2%set_indices(iref, self%s%iptcl)
                    inpl_ind = cls_inpl_inds(iref)
                    if( self%s%l_sh_first )then
                        cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%s%xy_first)
                        if( inpl_ind == 0 )then
                            inpl_ind = cls_inpl_inds(iref)
                            cxy(1)   = real(pftc_glob%gen_corr_for_rot_8(iref, self%s%iptcl, real(self%s%xy_first,dp), inpl_ind))
                            cxy(2:3) = self%s%xy_first_rot
                        endif
                    else
                        cxy = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                        if( inpl_ind == 0 )then
                            inpl_ind = cls_inpl_inds(iref)
                            cxy      = [real(pftc_glob%gen_corr_for_rot_8(iref, self%s%iptcl, inpl_ind)), 0.,0.]
                        endif
                    endif
                    cls_corrs(iref) = cxy(1)
                enddo
            endif
            ! Class selection
            call power_sampling( s2D%power, self%s%nrefs, cls_corrs, sorted_cls_inds, s2D%snhc_smpl_ncls,&
                                &self%s%best_class, self%s%nrefs_eval, self%s%best_corr )
            ! In-plane angle
            self%s%best_rot = cls_inpl_inds(self%s%best_class)
            ! In-plane search
            call self%s%inpl_srch
            ! Updates solution
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_snhc_smpl

    subroutine kill_snhc_smpl( self )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_smpl

end module simple_strategy2D_snhc_smpl
    