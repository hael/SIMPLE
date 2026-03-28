!@descr: 2D strategy for stochastic neighborhood hill climbing with probabilistic in-plane search
module simple_strategy2D_snhc_smpl
use simple_pftc_srch_api
use simple_strategy2D_alloc
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: builder
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

    subroutine new_snhc_smpl( self, params, spec, build )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        class(parameters),           intent(in)    :: params
        class(strategy2D_spec),      intent(inout) :: spec
        class(builder),              intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_snhc_smpl

    subroutine srch_snhc_smpl( self, os )   
        class(strategy2D_snhc_smpl), intent(inout) :: self
        class(oris),                 intent(inout) :: os
        real    :: inpl_corrs(self%s%nrots)
        real    :: inpl_corr
        integer :: sorted_cls_inds(self%s%nrefs), vec_nrots(self%s%nrots)
        integer :: iref, isample, inpl_ind, order_ind, class_rank
        integer :: nrefs_coarse_eval
        if( os%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch(os)
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! Class search
            nrefs_coarse_eval = 0
            do isample = 1,self%s%nrefs
                ! stochastic reference index
                iref = s2D%srch_order(self%s%iptcl_batch, isample)
                ! keep track of how many references we are evaluating
                nrefs_coarse_eval = nrefs_coarse_eval + 1
                ! neighbourhood size
                if(nrefs_coarse_eval > s2D%snhc_nrefs_bound) exit
                if( s2D%cls_pops(iref) == 0 )cycle
                ! In-plane sampling
                if( self%s%l_sh_first )then
                    call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                else
                    call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                endif
                call power_sampling( s2D%power, self%s%nrots, inpl_corrs, vec_nrots,&
                                    &s2D%snhc_smpl_ninpl, inpl_ind, order_ind, inpl_corr )
                call self%s%store_solution(iref, inpl_ind, inpl_corr)
            end do
            ! Performs shift search for top scoring subset
            call self%s%inpl_srch_peaks(min(s2D%snhc_smpl_ncls, self%s%nsolns))
            ! Class selection
            call power_sampling( s2D%power, self%s%nrefs, s2D%class_space_corrs(:, self%s%ithr), &
                                &sorted_cls_inds, s2D%snhc_smpl_ncls, &
                                &self%s%best_class, class_rank, self%s%best_corr )
            self%s%nrefs_eval = nrefs_coarse_eval
            ! In-plane angle
            self%s%best_rot = s2D%class_space_inplinds(self%s%best_class, self%s%ithr)
            ! In-plane search
            call self%s%inpl_srch ! needed because inpl_srch_peaks doesn't store shifts
            ! Updates solution
            call self%s%store_solution(self%s%best_class, self%s%best_rot, self%s%best_corr)
            call self%s%assign_ori(os)
        else
            call os%reject(self%s%iptcl)
        endif
    end subroutine srch_snhc_smpl

    subroutine kill_snhc_smpl( self )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_smpl

end module simple_strategy2D_snhc_smpl
    