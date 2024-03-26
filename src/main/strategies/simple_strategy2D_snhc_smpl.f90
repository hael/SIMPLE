module simple_strategy2D_snhc_smpl
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
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
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        class(strategy2D_snhc_smpl), intent(inout) :: self
        integer :: inds(self%s%nrots), iref, isample, inpl_ind, class_glob, inpl_glob, nrefs_bound, nrefs
        real    :: inpl_corrs(self%s%nrots), sorted_inpl_corrs(self%s%nrots), inpl_corr, cc_glob
        logical :: found_better
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            cc_glob      = -huge(cc_glob)
            found_better = .false.
            nrefs        = self%s%nrefs
            nrefs_bound  = min(nrefs, nint(real(nrefs)*(1.-self%spec%stoch_bound)))
            nrefs_bound  = max(2, nrefs_bound)
            do isample=1,self%s%nrefs
                ! stochastic reference index
                iref = s2D%srch_order(self%s%iptcl_map, isample)
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! neighbourhood size
                if(self%s%nrefs_eval > nrefs_bound) exit
                ! passes empty classes
                if( s2D%cls_pops(iref) == 0 )cycle
                ! multinomal in-plane update
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                inpl_ind  = angle_sampling(eulprob_dist_switch(inpl_corrs), sorted_inpl_corrs, inds, s2D%smpl_inpl_athres)
                inpl_corr = inpl_corrs(inpl_ind)
                ! keep track of global best
                if( inpl_corr > cc_glob )then
                    cc_glob       = inpl_corr
                    class_glob    = iref
                    inpl_glob     = inpl_ind
                endif
                ! improvement found (hill climbing over classes)
                if( inpl_corr > self%s%prev_corr ) found_better = .true.
                ! keep track of visited classes
                if( self%s%l_ptclw )then
                    s2D%cls_searched(iref,self%s%ithr) = .true.
                    s2D%cls_corrs(iref,self%s%ithr)    = inpl_corr
                endif
                ! first improvement heuristic
                if( found_better ) exit
            end do
            ! updates solution
            self%s%best_class = class_glob
            self%s%best_corr  = cc_glob
            self%s%best_rot   = inpl_glob
            if( params_glob%cc_objfun == OBJFUN_CC .and. params_glob%l_kweight_rot )then
                ! back-calculating in-plane angle with k-weighing
                if( found_better )then
                    call pftcc_glob%gencorrs(self%s%prev_class, self%s%iptcl, inpl_corrs, kweight=.true.)
                    self%s%prev_corr = inpl_corrs(self%s%prev_rot) ! updated threshold
                    call pftcc_glob%gencorrs(self%s%best_class, self%s%iptcl, inpl_corrs, kweight=.true.)
                    inpl_ind  = angle_sampling(eulprob_dist_switch(inpl_corrs), sorted_inpl_corrs, inds, s2D%smpl_inpl_athres)
                    inpl_corr = inpl_corrs(inpl_ind)
                    if( inpl_corr > self%s%prev_corr )then
                        ! improvement found
                    else
                        ! defaults to best
                        inpl_ind = maxloc(inpl_corrs, dim=1)
                    endif
                    self%s%best_rot  = inpl_ind
                    self%s%best_corr = inpl_corr
                else
                    call pftcc_glob%gencorrs(self%s%best_class, self%s%iptcl, inpl_corrs, kweight=.true.)
                    self%s%best_rot  = maxloc(inpl_corrs, dim=1)
                    self%s%best_corr = inpl_corrs(self%s%best_rot)
                endif
            endif
            call self%s%inpl_srch
            nrefs_bound = min(nrefs, nrefs_bound+1)
            call self%s%store_solution(nrefs=nrefs_bound)
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_snhc_smpl

    subroutine kill_snhc_smpl( self )
        class(strategy2D_snhc_smpl), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_smpl

end module simple_strategy2D_snhc_smpl
    