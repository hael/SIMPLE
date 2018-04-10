! concrete strategy3D: probabilistic single-state refinement
module simple_strategy3D_single
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D_multi, only: strategy3D_multi
implicit none

public :: strategy3D_single
private

logical, parameter :: DEBUG = .false.

type, extends(strategy3D_multi) :: strategy3D_single
contains
    procedure :: oris_assign => oris_assign_single
end type strategy3D_single

contains

    subroutine oris_assign_single( self )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_single), intent(inout) :: self
        type(ori) :: osym
        real      :: corrs(self%s%npeaks), ws(self%s%npeaks)
        real      :: wcorr, frac,  ang_sdev, dist, dist_inpl, euldist
        integer   :: best_loc(1)
        logical   :: included(self%s%npeaks)
        ! extract peak info
        call extract_peaks( self%s, corrs )
        ! stochastic weights
        call corrs2softmax_weights( self%s, corrs, self%spec%pp%tau, ws, included, best_loc, wcorr )
        ! B factors
        call fit_bfactors( self%s, ws )
        ! angular standard deviation
        ang_sdev = estimate_ang_sdev( self%s, best_loc )
        ! angular distances
        call self%s%se_ptr%sym_dists( self%s%a_ptr%get_ori(self%s%iptcl),&
            &o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl )
        ! generate convergence stats
        call convergence_stats_single( self%s, best_loc, euldist )
        ! fraction of search space scanned
        if( self%s%neigh )then
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nnn)
        else
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs)
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
        call self%s%a_ptr%set(self%s%iptcl, 'state',     1.)
        call self%s%a_ptr%set(self%s%iptcl, 'frac',      frac)
        call self%s%a_ptr%set(self%s%iptcl, 'corr',      wcorr)
        call self%s%a_ptr%set(self%s%iptcl, 'specscore', self%s%specscore)
        call self%s%a_ptr%set(self%s%iptcl, 'ow',        o_peaks(self%s%iptcl)%get(best_loc(1),'ow'))
        call self%s%a_ptr%set(self%s%iptcl, 'proj',      o_peaks(self%s%iptcl)%get(best_loc(1),'proj'))
        call self%s%a_ptr%set(self%s%iptcl, 'sdev',      ang_sdev)
        call self%s%a_ptr%set(self%s%iptcl, 'npeaks',    real(self%s%npeaks_eff))
        if( DEBUG ) print *,  '>>> STRATEGY3D_SINGLE :: EXECUTED ORIS_ASSIGN_SINGLE'
    end subroutine oris_assign_single

end module simple_strategy3D_single
