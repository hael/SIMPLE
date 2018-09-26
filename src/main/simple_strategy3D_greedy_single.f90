! concrete strategy3D: greedy single-state refinement
module simple_strategy3D_greedy_single
use simple_strategy3D_alloc         ! use all in there
use simple_strategy3D_utils         ! use all in there
use simple_strategy3D_greedy_multi, only: strategy3D_greedy_multi
use simple_parameters,              only: params_glob
use simple_builder,                 only: build_glob
implicit none

public :: strategy3D_greedy_single
private
#include "simple_local_flags.inc"

type, extends(strategy3D_greedy_multi) :: strategy3D_greedy_single

contains
    procedure :: oris_assign => oris_assign_greedy_single
end type strategy3D_greedy_single

contains

    subroutine oris_assign_greedy_single( self )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_greedy_single), intent(inout) :: self
        type(ori) :: osym
        real      :: corrs(self%s%npeaks), ws(self%s%npeaks)
        real      :: wcorr, frac, ang_spread, dist_inpl, euldist
        real      :: shwmean, shwstdev        
        integer   :: best_loc(1)
        logical   :: included(self%s%npeaks)
        ! extract peak info
        call extract_peaks(self%s, corrs)
        ! stochastic weights
        call corrs2softmax_weights(self%s, self%s%npeaks, corrs, ws, included, best_loc, wcorr)
        ! angular standard deviation
        ang_spread = estimate_ang_spread(self%s)
        call estimate_shift_increment(self%s, shwmean, shwstdev)        
        ! angular distances
        call build_glob%pgrpsyms%sym_dists( build_glob%spproj_field%get_ori(self%s%iptcl),&
            & s3D%o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl )
        ! fraction of search space scanned
        if( self%s%neigh )then
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nnn)
        else
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs)
        endif
        ! set the distances before we update the orientation
        if( build_glob%spproj_field%isthere(self%s%iptcl,'dist') )then
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(self%s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call build_glob%spproj_field%set_euler(self%s%iptcl,  s3D%o_peaks(self%s%iptcl)%get_euler(best_loc(1)))
        call build_glob%spproj_field%set_shift(self%s%iptcl,  s3D%o_peaks(self%s%iptcl)%get_2Dshift(best_loc(1)))
        call build_glob%spproj_field%set(self%s%iptcl, 'state',     1.)
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      frac)
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',      wcorr)
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call build_glob%spproj_field%set(self%s%iptcl, 'ow',        s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'ow'))
        call build_glob%spproj_field%set(self%s%iptcl, 'proj',      s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'proj'))
        call build_glob%spproj_field%set(self%s%iptcl, 'inpl',      s3D%o_peaks(self%s%iptcl)%get(best_loc(1),'inpl'))
        call build_glob%spproj_field%set(self%s%iptcl, 'spread',    ang_spread)
        call build_glob%spproj_field%set(self%s%iptcl, 'shwmean',   shwmean)
        call build_glob%spproj_field%set(self%s%iptcl, 'shwstdev',  shwstdev)
        call build_glob%spproj_field%set(self%s%iptcl, 'npeaks',    real(self%s%npeaks_eff))
    end subroutine oris_assign_greedy_single

end module simple_strategy3D_greedy_single
