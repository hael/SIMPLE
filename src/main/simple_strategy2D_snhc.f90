module simple_strategy2D_snhc
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy2D_snhc
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_snhc
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
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
        integer :: iref, isample, inpl_ind, class_glob, inpl_glob, nrefs_bound, nrefs
        real    :: corrs(self%s%nrots), inpl_corr, cc_glob
        logical :: found_better
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            call self%s%prep4srch
            cc_glob       = -1.
            found_better  = .false.
            nrefs         = count(s2D%cls_chunk==self%s%chunk_id)
            nrefs_bound   = min(nrefs, nint(real(nrefs)*(1.-self%spec%stoch_bound)))
            nrefs_bound   = max(2, nrefs_bound)
            do isample=1,self%s%nrefs
                ! stochastic reference index
                iref = s2D%srch_order(self%s%iptcl_map, isample)
                ! check whether we are in the correct chunk
                if( s2D%cls_chunk(iref) /= self%s%chunk_id )cycle
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! neighbourhood size
                if(self%s%nrefs_eval > nrefs_bound) exit
                ! passes empty classes
                if( s2D%cls_pops(iref) == 0 )cycle
                ! shc update
                call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
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
            nrefs_bound = min(nrefs, nrefs_bound+1)
            call self%s%store_solution(nrefs=nrefs_bound)
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        if( DEBUG ) print *, '>>> strategy2D_srch::FINISHED STOCHASTIC NEIGH SEARCH'
    end subroutine srch_snhc

    subroutine kill_snhc( self )
        class(strategy2D_snhc), intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc

end module simple_strategy2D_snhc
