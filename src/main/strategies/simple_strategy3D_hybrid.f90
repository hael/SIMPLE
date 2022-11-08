! concrete strategy3D: hybrid greedy,discrete,polar/stochastic,continuous,Cartesian refinement
module simple_strategy3D_hybrid
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_parameters,       only: params_glob
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_cartft_corrcalc,  only: cftcc_glob
implicit none

public :: strategy3D_hybrid
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_hybrid
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_hybrid
    procedure :: srch        => srch_hybrid
    procedure :: kill        => kill_hybrid
    procedure :: oris_assign => oris_assign_hybrid
end type strategy3D_hybrid

contains

    subroutine new_hybrid( self, spec )
        class(strategy3D_hybrid), intent(inout) :: self
        class(strategy3D_spec),   intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_hybrid

    subroutine srch_hybrid( self, ithr )
        class(strategy3D_hybrid), intent(inout) :: self
        integer,                  intent(in)    :: ithr
        integer, allocatable :: peak_inds(:)
        type(ori) :: osym, oi, oj, o, obest
        integer   :: iref, ipeak, isample, loc, cnt, npeaks, i, j, nevals(2), inds(self%s%nrefs)
        real      :: inpl_corrs(self%s%nrots), corrs(self%s%nrefs), angdist
        real      :: euldist, dist_inpl, sdev, var, cc_peak_avg, cc_nonpeak_avg
        real      :: corr, corr_best, cxy(3), shvec(2), shvec_incr(2), shvec_best(2)
        logical   :: peaks(self%s%nrefs), got_better, l_within_lims
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! DISCRETE GREEDY SEARCH
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! search
            corrs = -1.
            do iref= 1,self%s%nrefs
                ! store the reference index 
                inds(iref) = iref
                ! identify the top scoring in-plane angle
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                loc = maxloc(inpl_corrs,dim=1)
                ! store for peak detection
                corrs(iref) = inpl_corrs(loc)
                ! store for continuous refinement
                call self%s%eulspace%e3set(iref, 360. - pftcc_glob%get_rot(loc))
            end do
            ! detect peaks
            call self%s%eulspace%detect_peaks(s3D%proj_space_nnmat, corrs, peaks)
            npeaks = count(peaks)
            call build_glob%spproj_field%set(self%s%iptcl, 'npeaks', real(npeaks))
            angdist = 0.
            if( npeaks > 0 )then
                ! report average correlation for peak and non-peak distributions
                call build_glob%spproj_field%set(self%s%iptcl, 'cc_peak',    sum(corrs, mask=     peaks) / real(npeaks))
                call build_glob%spproj_field%set(self%s%iptcl, 'cc_nonpeak', sum(corrs, mask=.not.peaks) / real(count(.not.peaks)))
                ! calculate average angular distance between peaks
                angdist = 0.
                cnt     = 0
                do i = 1,self%s%nrefs
                    if( .not. peaks(i) ) cycle
                    call self%s%eulspace%get_ori(i, oi)
                    do j = 1,self%s%nrefs
                        if( .not. peaks(j) .or. i == j ) cycle
                        call self%s%eulspace%get_ori(j, oj)
                        call build_glob%pgrpsyms%sym_dists(oi, oj, osym, euldist, dist_inpl)
                        angdist = angdist + euldist
                        cnt = cnt + 1
                    end do
                end do
                if( cnt > 0 ) angdist = angdist / real(cnt)
            endif
            call build_glob%spproj_field%set(self%s%iptcl, 'dist_peaks', angdist)

            ! CONTINUOUS STOCHASTIC REFINEMENT
            ! prep
            call self%s%prep4_cont_srch
            ! transfer critical per-particle params
            o          = self%s%o_prev
            obest      = self%s%o_prev
            ! zero shifts because particle is shifted to its previous origin
            shvec_best = 0.
            ! currently the best correlation is the previous one
            corr_best  = self%s%prev_corr
            got_better = .false.
            ! point loc to the best discrete orientation
            loc        = maxloc(corrs,dim=1)
            ! extract the peak indices
            if( npeaks > 1 ) peak_inds  = pack(inds, mask=peaks) 
            do isample=1,self%s%nsample_neigh
                if( npeaks > 1 )then
                    ipeak = peak_inds(irnd_uni(npeaks))
                else
                    ipeak = loc
                endif
                ! make a random rotation matrix neighboring ipeak within the assymetric unit
                call obest%set_euler(self%s%eulspace%get_euler(ipeak))
                call build_glob%pgrpsyms%rnd_euler(obest, self%s%athres, o)
                if( self%s%nsample_trs > 0 )then
                    ! calculate Cartesian corr and simultaneously stochastically search shifts (Gaussian sampling)
                    shvec = shvec_best
                    call cftcc_glob%project_and_srch_shifts(self%s%iptcl, o, self%s%nsample_trs, params_glob%trs, shvec, corr)
                    if( corr > corr_best )then
                        corr_best  = corr
                        shvec_best = shvec
                        obest      = o
                        got_better = .true.
                    endif
                else
                    ! calculate Cartesian corr
                    corr = cftcc_glob%project_and_correlate(self%s%iptcl, o)
                    if( corr > corr_best )then
                        corr_best  = corr
                        obest      = o
                        got_better = .true.
                    endif
                endif
            end do
            if( allocated(peak_inds) ) deallocate(peak_inds)
            if( got_better )then
                call build_glob%pgrpsyms%sym_dists(self%s%o_prev, obest, osym, euldist, dist_inpl)
                call obest%set('dist',      euldist)
                call obest%set('dist_inpl', dist_inpl)
                call obest%set('corr',      corr_best)
                call obest%set('better',      1.0)
                call obest%set('frac',      100.0)
                call build_glob%spproj_field%set_ori(self%s%iptcl, obest)
                if( self%s%nsample_trs > 0 )then ! we did stochastic shift search
                    ! since particle image is shifted in the Cartesian formulation and we apply 
                    ! with negative sign in rec3D the sign of the increment found needs to be negative
                    shvec_incr = - shvec_best
                    shvec      =   self%s%prev_shvec + shvec_incr
                    where( abs(shvec) < 1e-6 ) shvec = 0.
                    call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
                    call build_glob%spproj_field%set(self%s%iptcl, 'shincarg', arg(shvec_incr))
                endif
            else
                call build_glob%spproj_field%set(self%s%iptcl, 'better', 0.)
            endif
            if( self%s%doshift .and. got_better ) then
                ! Cartesian shift search using L-BFGS-B with analytical derivatives
                call cftcc_glob%prep4shift_srch(self%s%iptcl, obest)
                l_within_lims = .true.
                if( self%s%nsample_trs > 0 )then ! we have a prior shift vector
                    ! check that it is within the limits
                    if( any(shvec_best < - params_glob%trs) ) l_within_lims = .false. 
                    if( any(shvec_best >   params_glob%trs) ) l_within_lims = .false.
                endif
                if( self%s%nsample_trs > 0 .and. l_within_lims )then ! refine the stochastic solution
                    cxy = self%s%shift_srch_cart(nevals, shvec_best)
                    ! report # cost function and gradient evaluations
                    call build_glob%spproj_field%set(self%s%iptcl, 'nevals',  real(nevals(1)))
                    call build_glob%spproj_field%set(self%s%iptcl, 'ngevals', real(nevals(2)))
                else if( self%s%nsample_trs == 0 )then ! shifts were not searched before
                    cxy = self%s%shift_srch_cart(nevals)
                    ! report # cost function and gradient evaluations
                    call build_glob%spproj_field%set(self%s%iptcl, 'nevals',  real(nevals(1)))
                    call build_glob%spproj_field%set(self%s%iptcl, 'ngevals', real(nevals(2)))
                else
                    ! wait with L-BFGS-B refinement until the particle is within the bounds
                    return
                endif
                shvec      = 0.
                shvec_incr = 0.
                if( cxy(1) >= corr_best )then
                    shvec      = self%s%prev_shvec
                    ! since particle image is shifted in the Cartesian formulation and we apply 
                    ! with negative sign in rec3D the sign of the increment found needs to be negative
                    shvec_incr = - cxy(2:3) 
                    shvec      = shvec + shvec_incr
                end if
                where( abs(shvec) < 1e-6 ) shvec = 0.
                call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
                call build_glob%spproj_field%set(self%s%iptcl, 'shincarg', arg(shvec_incr))
            endif
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_hybrid

    subroutine oris_assign_hybrid( self )
        class(strategy3D_hybrid), intent(inout) :: self
    end subroutine oris_assign_hybrid

    subroutine kill_hybrid( self )
        class(strategy3D_hybrid),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_hybrid

end module simple_strategy3D_hybrid
