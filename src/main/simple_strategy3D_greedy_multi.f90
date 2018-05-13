! concrete strategy3D: greedy multi-state refinement
module simple_strategy3D_greedy_multi
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
implicit none

public :: strategy3D_greedy_multi
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedy_multi
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_greedy_multi
    procedure :: srch        => srch_greedy_multi
    procedure :: kill        => kill_greedy_multi
    procedure :: oris_assign => oris_assign_greedy_multi
end type strategy3D_greedy_multi

contains

    subroutine new_greedy_multi( self, spec, npeaks )
        class(strategy3D_greedy_multi), intent(inout) :: self
        class(strategy3D_spec),         intent(inout) :: spec
        integer,                   intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_greedy_multi

    subroutine srch_greedy_multi( self )
        class(strategy3D_greedy_multi), intent(inout) :: self
        integer :: iref,isample,nrefs,target_projs(self%s%npeaks_grid), inpl_ind(1), state
        real    :: corrs(self%s%nrefs), inpl_corrs(self%s%nrots), inpl_corr
        if( self%s%a_ptr%get_state(self%s%iptcl) > 0 )then
            if( self%s%neigh )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. associated(self%spec%grid_projs) )&
                &stop 'need optional grid_projs 4 subspace srch; strategy3D_greedy_multi :: srch_greedy_multi'
                call self%s%greedy_subspace_srch(self%spec%grid_projs, target_projs)
                ! initialize
                call self%s%prep4srch(self%spec%nnmat, target_projs)
                nrefs = self%s%nnnrefs
            else
                ! initialize
                call self%s%prep4srch()
                nrefs = self%s%nrefs
            endif
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            proj_space_corrs(self%s%iptcl_map,:) = -1.
            ! search
            do isample=1,nrefs
                iref = srch_order(self%s%iptcl_map,isample) ! set the reference index
                call per_ref_srch                         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%s%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%s%iptcl_map,:))
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare weights & orientation
            call self%oris_assign()
        else
            call self%s%a_ptr%reject(self%s%iptcl)
        endif
        DebugPrint   '>>> STRATEGY3D_GREEDY_MULTI :: FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch
                state = proj_space_state(self%s%iptcl_map,iref)
                if( state_exists(state) )then    
                    call self%s%pftcc_ptr%gencorrs(iref, self%s%iptcl, inpl_corrs) ! in-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)                                ! greedy in-plane index
                    inpl_corr  = inpl_corrs(inpl_ind(1))                           ! max in plane correlation
                    call self%s%store_solution(iref, iref, inpl_ind(1), inpl_corr)
                endif
            end subroutine per_ref_srch

    end subroutine srch_greedy_multi

    subroutine oris_assign_greedy_multi( self )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_greedy_multi), intent(inout) :: self
        type(ori) :: osym
        real      :: corrs(self%s%npeaks), ws(self%s%npeaks), state_ws(self%s%nstates)
        real      :: wcorr, frac, ang_sdev, dist, dist_inpl, euldist
        integer   :: best_loc(1), neff_states, state
        logical   :: included(self%s%npeaks)
        ! extract peak info
        call extract_peaks( self%s, corrs, state )
        ! stochastic weights
        call corrs2softmax_weights( self%s, corrs, self%spec%pp%tau, ws, included, best_loc, wcorr )
        ! state reweighting
        call states_reweight( self%s, corrs, ws, state_ws, included, state, best_loc, wcorr )
        ! B factors
        call fit_bfactors( self%s, ws )
        ! angular standard deviation
        ang_sdev = estimate_ang_sdev( self%s, best_loc )
        ! angular distances
        call self%s%se_ptr%sym_dists( self%s%a_ptr%get_ori(self%s%iptcl),&
            &o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl )
        ! generate convergence stats
        call convergence_stats_multi( self%s, best_loc, euldist )
        ! fraction of search space scanned
        neff_states = count(state_exists)
        if( self%s%neigh )then
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nnn * neff_states)
        else
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs * neff_states)
        endif
        ! set the distances before we update the orientation
        if( self%s%a_ptr%isthere(self%s%iptcl,'dist') )then
            call self%s%a_ptr%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*self%s%a_ptr%get(self%s%iptcl,'dist'))
        else
            call self%s%a_ptr%set(self%s%iptcl, 'dist', euldist)
        endif
        call self%s%a_ptr%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call self%s%a_ptr%set_euler(self%s%iptcl, o_peaks(self%s%iptcl)%get_euler(best_loc(1)))
        call self%s%a_ptr%set_shift(self%s%iptcl, o_peaks(self%s%iptcl)%get_2Dshift(best_loc(1)))
        call self%s%a_ptr%set(self%s%iptcl, 'state',     real(state))
        call self%s%a_ptr%set(self%s%iptcl, 'frac',      frac)
        call self%s%a_ptr%set(self%s%iptcl, 'corr',      wcorr)
        call self%s%a_ptr%set(self%s%iptcl, 'specscore', self%s%specscore)
        call self%s%a_ptr%set(self%s%iptcl, 'ow',        o_peaks(self%s%iptcl)%get(best_loc(1),'ow')  )
        call self%s%a_ptr%set(self%s%iptcl, 'proj',      o_peaks(self%s%iptcl)%get(best_loc(1),'proj'))
        call self%s%a_ptr%set(self%s%iptcl, 'sdev',      ang_sdev)
        call self%s%a_ptr%set(self%s%iptcl, 'npeaks',    real(self%s%npeaks_eff))
        DebugPrint  '>>> STRATEGY3D_GREEDY_MULTI :: EXECUTED ORIS_ASSIGN_GREEDY_MULTI'
    end subroutine oris_assign_greedy_multi

    subroutine kill_greedy_multi( self )
        class(strategy3D_greedy_multi), intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_multi

end module simple_strategy3D_greedy_multi
