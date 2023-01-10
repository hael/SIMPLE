module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: extract_peak_ori, extract_peak_oris
private
#include "simple_local_flags.inc"

contains

    subroutine extract_peak_ori( s )
        class(strategy3D_srch), intent(inout) :: s
        type(ori) :: osym, o_prev, o_new
        integer   :: ref, inpl, state, neff_states, loc(1), nrefs_eval, nrefs_tot
        real      :: shvec(2), shvec_incr(2), mi_state, euldist, dist_inpl, corr, mi_proj, frac, pw
        logical   :: l_multistates
        ! stash previous ori
        call build_glob%spproj_field%get_ori(s%iptcl, o_prev)
        ! reference (proj)
        loc = maxloc(s3D%proj_space_corrs(s%ithr,:))
        ref = loc(1)
        if( ref < 1 .or. ref > s%nrefs ) THROW_HARD('ref index: '//int2str(ref)//' out of bound; extract_peak_ori')
        call build_glob%spproj_field%set(s%iptcl, 'proj', real(s3D%proj_space_proj(ref)))
        ! in-plane (inpl)
        inpl = s3D%proj_space_inplinds(s%ithr, ref)
        call build_glob%spproj_field%set(s%iptcl, 'inpl', real(inpl))
        ! Euler angle
        call build_glob%spproj_field%set_euler(s%iptcl, s3D%proj_space_euls(:,ref,s%ithr))
        ! shift
        shvec      = s%prev_shvec
        shvec_incr = 0.
        if( s%doshift ) then
            shvec_incr = s3D%proj_space_shift(:,ref,s%ithr)
            shvec      = shvec + shvec_incr
        end if
        where( abs(shvec) < 1e-6 ) shvec = 0.
        call build_glob%spproj_field%set_shift(s%iptcl, shvec)
        call build_glob%spproj_field%set(s%iptcl, 'shincarg', arg(shvec_incr))
        ! state
        state = 1
        l_multistates = s%nstates > 1
        if( l_multistates )then
            state = s3D%proj_space_state(ref)
            if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; extract_peak_ori')
        endif
        mi_state = 0.
        if( s%prev_state == state ) mi_state = 1.
        if( l_multistates )then
            call build_glob%spproj_field%set(s%iptcl, 'state',  real(state))
            call build_glob%spproj_field%set(s%iptcl, 'mi_state', mi_state)
        else
            call build_glob%spproj_field%set(s%iptcl, 'state',    1.)
            call build_glob%spproj_field%set(s%iptcl, 'mi_state', 1.)
        endif
        ! correlation
        corr = s3D%proj_space_corrs(s%ithr,ref)
        if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
            if( corr < 0. ) corr = 0.
        end if
        call build_glob%spproj_field%set(s%iptcl, 'corr', corr)
        ! angular distances
        call build_glob%spproj_field%get_ori(s%iptcl, o_new)
        call build_glob%pgrpsyms%sym_dists(o_prev, o_new, osym, euldist, dist_inpl)
        if( build_glob%spproj_field%isthere(s%iptcl,'dist') )then
            call build_glob%spproj_field%set(s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(s%iptcl, 'dist_inpl', dist_inpl)
        ! CONVERGENCE STATS
        ! projection direction overlap
        mi_proj  = 0.
        if( euldist < 0.5 ) mi_proj  = 1.
        call build_glob%spproj_field%set(s%iptcl, 'mi_proj', mi_proj)
        ! fraction of search space scanned
        neff_states = 1
        if( l_multistates ) neff_states = count(s3D%state_exists)
        if( s%l_neigh )then
            nrefs_tot  = s%nnn * neff_states
            if( s%nnn > 1 )then
                nrefs_eval = s%nrefs_eval
            else
                nrefs_eval = nrefs_tot  ! the case of global srch
            endif
        else if( s%l_greedy )then
            nrefs_tot  = s%nprojs * neff_states
            nrefs_eval = nrefs_tot
        else
            nrefs_eval = s%nrefs_eval
            nrefs_tot  = s%nprojs * neff_states
        endif
        frac = 100.0 * real(nrefs_eval) / real(nrefs_tot)
        call build_glob%spproj_field%set(s%iptcl, 'frac', frac)
        ! weight
        pw = 1.0
        if( s%l_ptclw ) call calc_ori_weight(s, ref, pw)
        call build_glob%spproj_field%set(s%iptcl, 'w', pw)
        ! destruct
        call osym%kill
        call o_prev%kill
        call o_new%kill
    end subroutine extract_peak_ori

    subroutine extract_peak_oris( s )
        class(strategy3D_srch), intent(inout) :: s
        integer   :: refs(s%npeaks), state, ipeak
        real      :: shvec(2), corr
        logical   :: l_multistates
        refs = maxnloc(s3D%proj_space_corrs(s%ithr,:), s%npeaks)
        if( any(refs < 1) .or. any(refs > s%nrefs) ) THROW_HARD('at least one refs index out of bound; extract_peak_oris')
        l_multistates = s%nstates > 1
        do ipeak = 1, s%npeaks
            call s%opeaks%set_euler(ipeak, s3D%proj_space_euls(:,refs(ipeak),s%ithr))
            shvec = s%prev_shvec
            if( s%doshift ) shvec = shvec + s3D%proj_space_shift(:,refs(ipeak),s%ithr)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            call s%opeaks%set_shift(ipeak, shvec)
            state = 1
            if( l_multistates )then
                state = s3D%proj_space_state(refs(ipeak))
                if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; extract_peak_oris')
            endif
            call s%opeaks%set_state(ipeak, state)
            corr = s3D%proj_space_corrs(s%ithr,refs(ipeak))
            if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
                if( corr < 0. ) corr = 0.
            end if
            call s%opeaks%set(ipeak, 'corr', corr)
        end do
    end subroutine extract_peak_oris

    subroutine calc_ori_weight( s, ref, pw )
        class(strategy3D_srch), intent(in)  :: s
        integer,                intent(in)  :: ref
        real,                   intent(out) :: pw
        ! real(dp) :: sumw, diff2, max_diff2, best_score, sum_diff, score_t, score
        real, allocatable :: scores(:), peak(:), nonpeak(:)
        real    :: best_score, score_t, ksstat, prob
        integer :: sz_peak, sz_nonpeak
        ! integer :: iref, npix
        pw = 1.0
        ! if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
        !     npix      = pftcc_glob%get_npix()
        !     max_diff2 = corr2distweight(s3D%proj_space_corrs(s%ithr,ref), npix, params_glob%tau)
        !     sumw      = 0.
        !     do iref = 1,s%nrefs
        !         if( s3D%proj_space_mask(iref,s%ithr) )then
        !             diff2 = corr2distweight(s3D%proj_space_corrs(s%ithr,iref), npix, params_glob%tau) - max_diff2
        !             sumw = sumw + exp(-diff2)
        !         endif
        !     enddo
        !     pw = max(0.,min(1.,real(1.d0 / sumw)))
        ! else
        !     ! objective function is exp(-euclid/denom) in [0,1] where 1 is best 
        !     best_score = real(s3D%proj_space_corrs(s%ithr,ref),kind=dp) ! best score as identified by stochastic search
        !     score_t    = real(calc_score_thres(s%nrefs, s3D%proj_space_corrs(s%ithr,:), NPEAKS_ORIW),kind=dp)
        !     sumw       = 0.d0
        !     do iref = 1,s%nrefs
        !         if( s3D%proj_space_mask(iref,s%ithr) )then ! only summing over references that have been evaluated
        !             score = real(s3D%proj_space_corrs(s%ithr,iref),kind=dp)
        !             if( score >= score_t )then
        !                 ! the argument to the exponential function is always negative
        !                 diff2 = score - best_score
        !                 ! hence, for the best score the exponent will be 1 and < 1 for all others
        !                 if( ref == iref )then
        !                     sumw  = sumw + 1.d0
        !                 else if( diff2 < 0.d0 )then
        !                     sumw  = sumw + exp(diff2)
        !                 endif
        !             endif
        !         endif
        !     enddo
        !     ! this normalization ensures that particles that do not show a distinct peak are down-weighted
        !     pw = max(0.,min(1.,real(1.d0 / sumw)))
        ! endif
        ! objective functions are exp(-euclid/denom) in [0,1] where 1 is best or corr [-1,1] where 1 is best
        ! best_score = s3D%proj_space_corrs(s%ithr,ref) ! best score as identified by stochastic search
        scores     = pack(s3D%proj_space_corrs(s%ithr,:), s3D%proj_space_mask(:,s%ithr))
        score_t    = calc_score_thres(s%nrefs, scores, NPEAKS_ORIW)
        peak       = pack(scores, mask=scores >= score_t)
        nonpeak    = pack(scores, mask=scores <  score_t)
        sz_peak    = 0
        sz_nonpeak = 0
        if( allocated(peak)    ) sz_peak    = size(peak)
        if( allocated(nonpeak) ) sz_nonpeak = size(nonpeak)
        if( sz_peak > 0 .and. sz_nonpeak > 0 )then
            call kstwo(peak, sz_peak, nonpeak, sz_nonpeak, ksstat, prob) ! prob is probability that distributions are the same
            ! pw = best_score * (1. - prob)
            pw = 1. - prob
        endif
        if( allocated(scores)  ) deallocate(scores)
        if( allocated(peak)    ) deallocate(peak)
        if( allocated(nonpeak) ) deallocate(nonpeak)
    end subroutine calc_ori_weight

end module simple_strategy3D_utils
