! concrete strategy3D: refinement
module simple_strategy3D_snhc_smpl
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
use simple_eul_prob_tab2D,   only: neighfrac2nsmpl, power_sampling
use simple_polarft_calc, only: pftc_glob
use simple_decay_funs,       only: extremal_decay
implicit none

public :: strategy3D_snhc_smpl
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_snhc_smpl
contains
    procedure :: new         => new_snhc_smpl
    procedure :: srch        => srch_snhc_smpl
    procedure :: oris_assign => oris_assign_snhc_smpl
    procedure :: kill        => kill_snhc_smpl
end type strategy3D_snhc_smpl

contains

    subroutine new_snhc_smpl( self, spec )
        class(strategy3D_snhc_smpl), intent(inout) :: self
        class(strategy3D_spec),      intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_snhc_smpl

    subroutine srch_snhc_smpl( self, ithr )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy3D_snhc_smpl), intent(inout) :: self
        integer,                     intent(in)    :: ithr
        real    :: sorted_corrs(self%s%nrefs), inpl_corrs(self%s%nrots)
        real    :: cxy(3), inpl_corr, neigh_frac, power
        integer :: vec_nrots(self%s%nrots), sorted_inds(self%s%nrefs)
        integer :: iref, isample, nrefs_bound, inpl_ind, order_ind, smpl_nrefs, smpl_ninpl
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! # of references to evaluate (extremal optimization)
            neigh_frac  = extremal_decay( params_glob%extr_iter, params_glob%extr_lim )
            nrefs_bound = max(2,min(self%s%nrefs, nint(real(self%s%nrefs)*(1.-neigh_frac))))
            ! # for out-of-plane sampling
            smpl_nrefs = neighfrac2nsmpl(neigh_frac, self%s%nrefs)
            ! # for in-plane sampling
            smpl_ninpl = neighfrac2nsmpl(neigh_frac, self%s%nrots)
            ! greediness of sampling
            power = merge(EXTR_POWER, POST_EXTR_POWER, params_glob%extr_iter<=params_glob%extr_lim)
            ! initialize
            self%s%nrefs_eval =  0
            ! search
            do isample=1,self%s%nrefs
                ! set the stochastic reference index
                iref = s3D%srch_order(isample,self%s%ithr)
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! neighbourhood size
                if(self%s%nrefs_eval > nrefs_bound) exit
                ! empty space
                if( .not.s3D%state_exists(s3D%proj_space_state(iref)) )cycle
                ! In-plane sampling
                call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.], inpl_corrs)
                call power_sampling( power, self%s%nrots, inpl_corrs, vec_nrots,&
                                    &smpl_ninpl, inpl_ind, order_ind, inpl_corr )
                call self%s%store_solution(iref, inpl_ind, inpl_corr)
            enddo
            ! Performs shift search for top scoring subset
            if( self%s%doshift )then
                ! Sort
                sorted_corrs = s3D%proj_space_corrs(:, self%s%ithr)
                sorted_inds  = (/(iref, iref=1,self%s%nrefs)/)
                call hpsort(sorted_corrs, sorted_inds)
                ! Subset offset search
                s3D%proj_space_corrs(:,self%s%ithr) = -1.0
                do isample = self%s%nrefs-smpl_nrefs+1, self%s%nrefs
                    iref = sorted_inds(isample)
                    call self%s%grad_shsrch_obj2%set_indices(iref, self%s%iptcl)
                    inpl_ind = s3D%proj_space_inplinds(iref,self%s%ithr)
                    cxy      = self%s%grad_shsrch_obj2%minimize(irot=inpl_ind)
                    if( inpl_ind == 0 )then
                        inpl_ind = s3D%proj_space_inplinds(iref,self%s%ithr)
                        cxy      = [real(pftc_glob%gen_corr_for_rot_8(iref, self%s%iptcl, inpl_ind)), 0.,0.]
                    endif
                    call self%s%store_solution(iref, inpl_ind, cxy(1), sh=cxy(2:3))
                enddo
            endif
            ! Projection direction samling
            sorted_corrs = s3D%proj_space_corrs(:,self%s%ithr)
            call power_sampling( power, self%s%nrefs, sorted_corrs, sorted_inds, smpl_nrefs,&
                                &iref, self%s%nrefs_eval, inpl_corr )
            ! In-plane search
            call self%s%inpl_srch(ref=iref)
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_snhc_smpl

    subroutine oris_assign_snhc_smpl( self )
        class(strategy3D_snhc_smpl), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_snhc_smpl

    subroutine kill_snhc_smpl( self )
        class(strategy3D_snhc_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_smpl

end module simple_strategy3D_snhc_smpl
