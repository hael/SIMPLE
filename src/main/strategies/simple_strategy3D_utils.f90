module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_calc, only: pftc_glob
implicit none

public :: extract_peak_ori, extract_peak_oris, assign_ori
private
#include "simple_local_flags.inc"

contains

    subroutine assign_ori( s, ref, inpl, corr, sh, w )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: ref, inpl  ! ref here is multi-state ref index
        real,                   intent(in)    :: corr
        real,                   intent(in)    :: sh(2)
        real,         optional, intent(in)    :: w
        type(ori) :: osym, o_prev, o_new
        integer   :: state, neff_states, nrefs_eval, nrefs_tot
        real      :: shvec(2), shvec_incr(2), mi_state, euldist, dist_inpl, mi_proj, frac, pw
        logical   :: l_multistates
        s3D%proj_space_euls(3,ref,s%ithr) = 360. - pftc_glob%get_rot(inpl)
        ! stash previous ori
        call build_glob%spproj_field%get_ori(s%iptcl, o_prev)
        ! reference (proj)
        if( ref < 1 .or. ref > s%nrefs ) THROW_HARD('ref index: '//int2str(ref)//' out of bound; assign_ori')
        call build_glob%spproj_field%set(s%iptcl, 'proj', real(s3D%proj_space_proj(ref)))
        ! in-plane (inpl)
        call build_glob%spproj_field%set(s%iptcl, 'inpl', real(inpl))
        ! Euler angle
        call build_glob%spproj_field%set_euler(s%iptcl, s3D%proj_space_euls(:,ref,s%ithr))
        ! shift
        shvec      = s%prev_shvec
        shvec_incr = 0.
        if( s%doshift ) then
            shvec_incr = sh
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
            if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; assign_ori')
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
            select case(trim(s%refine))
                case('shc_neigh')
                    nrefs_tot  = s%nprojs_sub * neff_states
                    nrefs_eval = s%nrefs_eval
                case DEFAULT
                    nrefs_tot  = s%nnn * neff_states
                    if( s%nnn > 1 )then
                        nrefs_eval = s%nrefs_eval
                    else
                        nrefs_eval = nrefs_tot  ! the case of global srch
                    endif
            end select
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
        pw = s3D%proj_space_w(ref, s%ithr)
        if( (trim(params_glob%ptclw) .eq. 'yes') .and. present(w) ) pw = w
        call build_glob%spproj_field%set(s%iptcl, 'w', pw)
        ! destruct
        call osym%kill
        call o_prev%kill
        call o_new%kill
    end subroutine assign_ori

    subroutine extract_peak_ori( s )
        class(strategy3D_srch), intent(inout) :: s
        integer :: ref, inpl, loc(1)
        real    :: corr, sh(2)
        loc  = maxloc(s3D%proj_space_corrs(:,s%ithr))
        ref  = loc(1)
        corr = s3D%proj_space_corrs(   ref,s%ithr)
        sh   = s3D%proj_space_shift(:, ref,s%ithr)
        inpl = s3D%proj_space_inplinds(ref,s%ithr)
        if( params_glob%cc_objfun == OBJFUN_CC .and. corr < 0. ) corr = 0.
        call assign_ori( s, ref, inpl, corr, sh )
    end subroutine extract_peak_ori

    subroutine extract_peak_oris( s )
        class(strategy3D_srch), intent(inout) :: s
        integer   :: refs(s%npeaks), state, ipeak
        real      :: shvec(2), corr
        logical   :: l_multistates
        refs = maxnloc(s3D%proj_space_corrs(:,s%ithr), s%npeaks)
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
            corr = s3D%proj_space_corrs(refs(ipeak),s%ithr)
            if( params_glob%cc_objfun == OBJFUN_CC )then
                if( corr < 0. ) corr = 0.
            end if
            call s%opeaks%set(ipeak, 'corr', corr)
        end do
    end subroutine extract_peak_oris

end module simple_strategy3D_utils
