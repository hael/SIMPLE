module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: extract_peak_ori
private
#include "simple_local_flags.inc"

contains

    subroutine extract_peak_ori( s )
        use simple_ori, only: ori
        class(strategy3D_srch), intent(inout) :: s
        type(ori) :: osym, o_prev, o_new
        integer   :: ref, inpl, state, neff_states, loc(1)
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
        ! specscore
        call build_glob%spproj_field%set(s%iptcl, 'specscore', s%specscore)
        ! weight
        pw = 1.0
        if( s%l_ptclw ) call calc_ori_weight(s, ref, pw)
        call build_glob%spproj_field%set(s%iptcl, 'w', pw)
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
            if( s%nnn > 1 )then
                frac = 100. * real(s%nrefs_eval) / real(s%nnn * neff_states)
            else
                ! the case of global srch
                frac = 100.
            endif
        else if( s%l_greedy .or. s%l_cont )then
            frac = 100.
        else
            frac = 100. * real(s%nrefs_eval) / real(s%nprojs * neff_states)
        endif
        call build_glob%spproj_field%set(s%iptcl, 'frac', frac)
        ! destruct
        call osym%kill
        call o_prev%kill
        call o_new%kill
    end subroutine extract_peak_ori

    subroutine calc_ori_weight( s, ref, pw )
        class(strategy3D_srch), intent(in)  :: s
        integer,                intent(in)  :: ref
        real,                   intent(out) :: pw
        real(dp) :: sumw, diff2, max_diff2
        integer  :: iref
        pw = 1.0
        if( params_glob%cc_objfun /= OBJFUN_EUCLID ) return
        max_diff2 = s3D%proj_space_corrs(s%ithr,ref)
        sumw      = 0.d0
        do iref = 1,s%nrefs
            if( s3D%proj_space_mask(iref,s%ithr) )then
                diff2 = real(max_diff2 - s3D%proj_space_corrs(s%ithr,iref),dp)
                if( diff2 < 700.d0 ) sumw = sumw + exp(-diff2)
            endif
        enddo
        pw = min(1.0,real(1.d0 / sumw))
    end subroutine calc_ori_weight

end module simple_strategy3D_utils
