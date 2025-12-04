module simple_strategy2D_snhc
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_calc, only: pftc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy2D_snhc
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_snhc
  contains
    procedure :: new  => new_snhc
    procedure :: srch => srch_snhc
    procedure :: kill => kill_snhc
end type strategy2D_snhc

contains

    subroutine new_snhc( self, spec )
        class(strategy2D_snhc), intent(inout) :: self
        class(strategy2D_spec), intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_snhc

    subroutine srch_snhc( self )
        class(strategy2D_snhc), intent(inout) :: self
        integer :: iref, isample, inpl_ind, class_glob, inpl_glob
        real    :: corrs(self%s%nrots), inpl_corr, cc_glob
        logical :: found_better
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! Prep
            call self%s%prep4srch
            ! Shift search on previous best reference
            call self%s%inpl_srch_first
            ! Class search
            cc_glob      = -huge(cc_glob)
            found_better = .false.
            do isample=1,self%s%nrefs
                ! stochastic reference index
                iref = s2D%srch_order(self%s%iptcl_batch, isample)
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! neighbourhood size
                if(self%s%nrefs_eval > s2D%snhc_nrefs_bound) exit
                ! passes empty classes
                if( s2D%cls_pops(iref) == 0 )cycle
                ! shc update
                if( self%s%l_sh_first )then
                    call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, corrs)
                else
                    call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         corrs)
                endif
                inpl_ind = shcloc(self%s%nrots, corrs, self%s%prev_corr)
                if( inpl_ind == 0 )then
                    ! update inpl_ind & inpl_corr to greedy best
                    inpl_ind  = maxloc(corrs, dim=1)
                    inpl_corr = corrs(inpl_ind)
                else
                    ! use the parameters selected by SHC condition
                    inpl_corr         = corrs(inpl_ind)
                    self%s%best_class = iref
                    self%s%best_corr  = inpl_corr
                    self%s%best_rot   = inpl_ind
                    found_better      = .true.
                endif
                ! keep track of global best
                if( inpl_corr > cc_glob )then
                    cc_glob       = inpl_corr
                    class_glob    = iref
                    inpl_glob     = inpl_ind
                endif
                ! first improvement heuristic
                if( found_better ) exit
            end do
            if( found_better )then
                ! best ref has already been updated
            else
                ! use the globally best parameters
                self%s%best_class = class_glob
                self%s%best_corr  = cc_glob
                self%s%best_rot   = inpl_glob
            endif
            call self%s%inpl_srch
            call self%s%store_solution(nrefs=min(self%s%nrefs, s2D%snhc_nrefs_bound+1))
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        if( DEBUG ) write(logfhandle,*) '>>> strategy2D_srch::FINISHED STOCHASTIC NEIGH SEARCH'
    end subroutine srch_snhc

    subroutine kill_snhc( self )
        class(strategy2D_snhc), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc

end module simple_strategy2D_snhc
