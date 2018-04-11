module simple_strategy3D_utils
use simple_defs              ! use all in there
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_srch,  only: strategy3D_srch
implicit none

contains

    subroutine extract_peaks( s, corrs, state )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: corrs(s%npeaks)
        integer, optional,      intent(out)   :: state
        integer :: ipeak, cnt, ref
        real    :: shvec(2)
        logical :: state_present
        state_present = present(state)
        do ipeak = 1, s%npeaks
            cnt = s%nrefs - s%npeaks + ipeak
            ref = proj_space_inds(s%iptcl_map, cnt)
            if( ref < 1 .or. ref > s%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; strategy3D_utils :: extract_peaks'
            endif
            if( state_present )then
                state = proj_space_state(s%iptcl_map,ref)
                if( .not. state_exists(state) )then
                    print *, 'empty state: ', state
                    stop 'strategy3D_utils :: extract_peak'
                endif
            endif
            ! add shift
            shvec = s%prev_shvec
            if( s%doshift )shvec = shvec + proj_space_shift(s%iptcl_map,ref,1:2)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = proj_space_corrs(s%iptcl_map,ref)
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            if( state_present )then
                call o_peaks(s%iptcl)%set(ipeak, 'state', real(state))
            else
                call o_peaks(s%iptcl)%set(ipeak, 'state', 1.)
            endif
            call o_peaks(s%iptcl)%set(ipeak, 'proj',  real(proj_space_proj(s%iptcl_map,ref)))
            call o_peaks(s%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call o_peaks(s%iptcl)%set_euler(ipeak, proj_space_euls(s%iptcl_map,ref,1:3))
            call o_peaks(s%iptcl)%set_shift(ipeak, shvec)
        enddo
    end subroutine extract_peaks

    subroutine corrs2softmax_weights( s, corrs, tau, ws, included, best_loc, wcorr )
        use simple_ori, only: ori
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: corrs(s%npeaks), tau
        real,                   intent(out)   :: ws(s%npeaks), wcorr
        logical,                intent(out)   :: included(s%npeaks)
        integer,                intent(out)   :: best_loc(1)
        type(ori) :: o, o_best
        real      :: dists(s%npeaks), arg4softmax(s%npeaks), ang_dist, inpl_dist
        real      :: ang_weights(s%npeaks), ang_lim
        integer   :: ipeak
        if( s%npeaks == 1 )then
            best_loc(1)  = 1
            ws(1)        = 1.
            included(1)  = .true.
            s%npeaks_eff = 1
            wcorr        = corrs(1)
        else
            best_loc = maxloc(corrs)
            ! reweighting according to angular distance to best peak
            ang_lim = deg2rad( s%se_ptr%srchrange_theta()/2. )
            o_best  = o_peaks(s%iptcl)%get_ori(best_loc(1))
            do ipeak = 1,s%npeaks
                if( ipeak == best_loc(1) )then
                    ang_weights(ipeak) = 1.0
                    cycle
                endif
                if( trim(s%se_ptr%get_pgrp()).eq.'c1' )then
                    ang_dist = o_peaks(s%iptcl)%get_ori(ipeak).euldist.o_best
                else
                    call s%se_ptr%sym_dists( o_best, o_peaks(s%iptcl)%get_ori(ipeak), o, ang_dist, inpl_dist)
                    ang_dist = deg2rad( ang_dist )
                endif
                if( ang_dist > ang_lim )then
                    ang_weights(ipeak) = 0.
                else
                    ang_weights(ipeak) = max(0., cos(ang_dist/ang_lim*PIO2))**2.
                endif
            enddo
            ! convert correlations to distances
            dists = 1.0 - corrs
            ! scale distances with TAU
            dists = dists / tau
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate softmax weights
            ws = exp(arg4softmax) * ang_weights
            ! normalise
            ws = ws / sum(ws)
            ! threshold weights
            included   = (ws >= SOFTMAXW_THRESH)
            s%npeaks_eff = count(included)
            where( .not. included ) ws = 0.
            ! weighted corr
            wcorr = sum(ws*corrs,mask=included)
        endif
        ! update npeaks individual weights
        call o_peaks(s%iptcl)%set_all('ow', ws)
    end subroutine corrs2softmax_weights

    subroutine states_reweight( s, corrs, ws, state_ws, included, state, best_loc, wcorr )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: corrs(s%npeaks)
        real,                   intent(inout) :: ws(s%npeaks)
        real,                   intent(out)   :: state_ws(s%nstates), wcorr
        logical,                intent(inout) :: included(s%npeaks)
        integer,                intent(out)   :: state, best_loc(1)
        integer :: istate, loc(1), ipeak, states(s%npeaks)
        if( s%npeaks > 1 .and. s%nstates > 1 )then
            ! states weights
            do ipeak = 1, s%npeaks
                corrs(ipeak)  = o_peaks(s%iptcl)%get(ipeak,'corr')
                states(ipeak) = nint(o_peaks(s%iptcl)%get(ipeak,'state'))
            enddo
            ! greedy state assignment
            do istate = 1, s%nstates
                state_ws(istate) = sum(ws,mask=(states==istate))
            enddo
            loc   = maxloc(state_ws)
            state = loc(1)
            ! in-state re-weighing
            included = included .and. (states==state)
            where( .not. included ) ws = 0.
            ws = ws / sum(ws, mask=included)
            best_loc = maxloc(ws)
            ! weighted corr
            wcorr    = sum(ws*corrs, mask=included)
            ! update npeaks individual weights
            call o_peaks(s%iptcl)%set_all('ow', ws)
        endif
    end subroutine states_reweight

    subroutine fit_bfactors( s, ws )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: ws(s%npeaks)
        integer :: ipeak, cnt, ref, roind
        real    :: shvec(2), bfacs(s%npeaks), bfac
        if( s%pftcc_ptr%objfun_is_ccres() )then
            bfacs = 0.
            do ipeak = 1, s%npeaks
                if( ws(ipeak) > TINY .or. s%npeaks == 1 )then
                    cnt   = s%nrefs - s%npeaks + ipeak
                    ref   = proj_space_inds(s%iptcl_map, cnt)
                    shvec = 0.
                    if( s%doshift )shvec = proj_space_shift(s%iptcl_map, ref, 1:2)
                    roind        = s%pftcc_ptr%get_roind(360. - proj_space_euls(s%iptcl_map, ref, 3))
                    bfacs(ipeak) = s%pftcc_ptr%fit_bfac(ref, s%iptcl, roind, shvec)
                else
                    bfacs(ipeak) = 0.
                endif
                call o_peaks(s%iptcl)%set(ipeak, 'bfac', bfacs(ipeak))
            enddo
            bfac = sum(ws * bfacs, mask=(ws>TINY))
            call s%a_ptr%set(s%iptcl, 'bfac',  bfac )
        endif
    end subroutine fit_bfactors

    function estimate_ang_sdev( s, best_loc ) result( ang_sdev )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: best_loc(1)
        real       :: ang_sdev, dist, inpl_dist
        integer    :: ipeak
        type(ori)  :: osym
        type(oris) :: sym_os
        ang_sdev = 0.
        if( trim(s%se_ptr%get_pgrp()).eq.'c1' )then
            ang_sdev = o_peaks(s%iptcl)%ang_sdev(s%nstates, s%npeaks)
        else
            if( s%npeaks > 2 )then
                sym_os = o_peaks(s%iptcl)
                do ipeak = 1, s%npeaks
                    if( ipeak == best_loc(1) )cycle
                    call s%se_ptr%sym_dists( o_peaks(s%iptcl)%get_ori(best_loc(1)),&
                        &o_peaks(s%iptcl)%get_ori(ipeak), osym, dist, inpl_dist)
                    call sym_os%set_ori(ipeak, osym)
                enddo
                ang_sdev = sym_os%ang_sdev(s%nstates, s%npeaks)
            endif
        endif
    end function estimate_ang_sdev

    subroutine convergence_stats_single( s, best_loc, euldist )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: best_loc(1)
        real,                   intent(in)    :: euldist
        integer :: roind
        real    :: mi_proj, mi_inpl, mi_joint
        roind    = s%pftcc_ptr%get_roind(360. - o_peaks(s%iptcl)%e3get(best_loc(1)))
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj = 1.
            mi_joint = mi_joint + 1.
        endif
        if( s%prev_roind == roind )then
            mi_inpl  = 1.
            mi_joint = mi_joint + 1.
        endif
        mi_joint = mi_joint/2.
        call s%a_ptr%set(s%iptcl, 'mi_proj',   mi_proj)
        call s%a_ptr%set(s%iptcl, 'mi_inpl',   mi_inpl)
        call s%a_ptr%set(s%iptcl, 'mi_state',  1.)
        call s%a_ptr%set(s%iptcl, 'mi_joint',  mi_joint)
    end subroutine convergence_stats_single

    subroutine convergence_stats_multi( s, best_loc, euldist )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: best_loc(1)
        real,                   intent(in)    :: euldist
        integer :: roind, state
        real    :: mi_proj, mi_inpl, mi_joint, mi_state
        roind    = s%pftcc_ptr%get_roind(360. - o_peaks(s%iptcl)%e3get(best_loc(1)))
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj = 1.
            mi_joint = mi_joint + 1.
        endif
        if( s%prev_roind == roind )then
            mi_inpl  = 1.
            mi_joint = mi_joint + 1.
        endif
        mi_state = 0.
        state = nint( o_peaks(s%iptcl)%get(best_loc(1), 'state') )
        if( .not. state_exists(state) )then
            print *, 'empty state: ', state
            stop 'strategy3D_utils; convergence_stats_multi'
        endif
        if( s%prev_state == state )then
            mi_state = 1.
            mi_joint = mi_joint + mi_state
        endif
        mi_joint = mi_joint/3.
        ! set the overlaps
        call s%a_ptr%set(s%iptcl, 'mi_proj',   mi_proj)
        call s%a_ptr%set(s%iptcl, 'mi_inpl',   mi_inpl)
        call s%a_ptr%set(s%iptcl, 'mi_state',  mi_state)
        call s%a_ptr%set(s%iptcl, 'mi_joint',  mi_joint)
    end subroutine convergence_stats_multi

end module simple_strategy3D_utils
