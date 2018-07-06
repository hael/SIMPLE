module simple_strategy2D_stochasticneigh
include 'simple_lib.f08'
use simple_strategy2D_alloc  ! singleton
use simple_strategy2D,       only: strategy2D
use simple_strategy2D_srch,  only: strategy2D_srch, strategy2D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy2D_stochasticneigh
private

logical, parameter :: DEBUG   = .false.

type, extends(strategy2D) :: strategy2D_stochasticneigh
    type(strategy2D_srch) :: s
    type(strategy2D_spec) :: spec
contains
    procedure :: new  => new_stochasticneigh
    procedure :: srch => srch_stochasticneigh
    procedure :: kill => kill_stochasticneigh
end type strategy2D_stochasticneigh

contains

    subroutine new_stochasticneigh( self, spec )
        class(strategy2D_stochasticneigh), intent(inout) :: self
        class(strategy2D_spec),            intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_stochasticneigh

    subroutine srch_stochasticneigh( self )
        class(strategy2D_stochasticneigh), intent(inout) :: self
        integer :: iref, isample, inpl_ind, class_glob, inpl_glob, nrefs_bound
        real    :: corrs(self%s%nrots), inpl_corr, cc_glob
        logical :: found_better, do_inplsrch, glob_best_set, do_shc
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            do_inplsrch       = .true.
            cc_glob           = -1.
            glob_best_set     = .false.
            found_better      = .false.
            self%s%nrefs_eval = 0
            call self%s%prep4srch
            ! objective function based logics for performing extremal search
            if( pftcc_glob%get_objfun() == 2 )then
                ! based on b-factor for objfun=ccres (the lower the better)
                do_shc = (self%spec%extr_bound < 0.) .or. (self%s%prev_bfac < self%spec%extr_bound)
            else
                ! based on correlation for objfun=cc (the higher the better)
                do_shc = (self%spec%extr_bound < 0.) .or. (self%s%prev_corr > self%spec%extr_bound)
            endif
            if( do_shc )then
                ! SHC move
                do isample=1,self%s%nrefs
                    iref = s2D%srch_order(self%s%iptcl_map, isample)
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                    ! passes empty classes
                    if( s2D%cls_pops(iref) == 0 )cycle
                    ! shc update
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
                    inpl_ind  = shcloc(self%s%nrots, corrs, self%s%prev_corr)
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
                        glob_best_set = .true.
                    endif
                    ! first improvement heuristic
                    if( found_better ) exit
                end do
                if( found_better )then
                    ! best ref has already been updated
                else
                    if( glob_best_set )then
                        ! use the globally best parameters
                        self%s%best_class = class_glob
                        self%s%best_corr  = cc_glob
                        self%s%best_rot   = inpl_glob
                    else
                        ! keep the old parameters
                        self%s%best_class = self%s%prev_class
                        self%s%best_corr  = self%s%prev_corr
                        self%s%best_rot   = self%s%prev_rot
                    endif
                endif
            else
                ! quasi-random move: best of 10% of the search space
                nrefs_bound = max(2, min(nint(0.1*real(self%s%nrefs)), self%s%nrefs))
                do_inplsrch = .false.
                do isample=1,self%s%nrefs
                    ! stochatsic reference index
                    iref = s2D%srch_order(self%s%iptcl_map, isample)
                    ! passes empty classes
                    if( s2D%cls_pops(iref) == 0 )cycle
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                    ! restricted neighbourhood
                    if( self%s%nrefs_eval > nrefs_bound )exit
                    ! shc update
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs)
                    inpl_ind  = maxloc(corrs, dim=1)
                    inpl_corr = corrs(inpl_ind)
                    ! keep track of global best
                    if( inpl_corr > cc_glob )then
                        self%s%best_class = iref
                        self%s%best_corr  = inpl_corr
                        self%s%best_rot   = inpl_ind
                        cc_glob = inpl_corr
                    endif
                end do
                self%s%nrefs_eval = min(isample,nrefs_bound)
            endif
            if( do_inplsrch )call self%s%inpl_srch
            call self%s%store_solution
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        if( DEBUG ) print *, '>>> strategy2D_srch::FINISHED STOCHASTIC NEIGH SEARCH'
    end subroutine srch_stochasticneigh

    subroutine kill_stochasticneigh( self )
        class(strategy2D_stochasticneigh), intent(inout) :: self
        call self%s%kill
    end subroutine kill_stochasticneigh

end module simple_strategy2D_stochasticneigh
