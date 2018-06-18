module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
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
            ref = s3D%proj_space_inds(s%iptcl_map, cnt)
            if( ref < 1 .or. ref > s%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; strategy3D_utils :: extract_peaks'
            endif
            if( state_present )then
                state = s3D%proj_space_state(s%iptcl_map,ref)
                if( .not. s3D%state_exists(state) )then
                    print *, 'empty state: ', state
                    stop 'strategy3D_utils :: extract_peak'
                endif
            endif
            ! add shift
            shvec = s%prev_shvec
            if( s%doshift )shvec = shvec + s3D%proj_space_shift(s%iptcl_map,ref,1:2)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = s3D%proj_space_corrs(s%iptcl_map,ref)
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            if( state_present )then
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', real(state))
            else
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', 1.)
            endif
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'proj',  real(s3D%proj_space_proj(s%iptcl_map,ref)))
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call s3D%o_peaks(s%iptcl)%set_euler(ipeak, s3D%proj_space_euls(s%iptcl_map,ref,1:3))
            call s3D%o_peaks(s%iptcl)%set_shift(ipeak, shvec)
        enddo
    end subroutine extract_peaks

    subroutine prob_select_peak( s, tau, updatecnt )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: tau, updatecnt
        real    :: corrs(s%npeaks), shvecs(s%npeaks,2), euls(s%npeaks,3), pvec(s%npeaks)
        integer :: cnt, inds(s%npeaks), refs(s%npeaks), states(s%npeaks), projs(s%npeaks)
        integer :: prob_peak, ipeak, rank, loc(1)
        real    :: bound, rnd, dists(s%npeaks), arg4softmax(s%npeaks)
        do ipeak = 1, s%npeaks
            ! stash peak index
            inds(ipeak) = ipeak
            ! stash reference index
            cnt = s%nrefs - s%npeaks + ipeak
            refs(ipeak) = s3D%proj_space_inds(s%iptcl_map, cnt)
            if( refs(ipeak) < 1 .or. refs(ipeak) > s%nrefs )then
                print *, 'refs(ipeak): ', refs(ipeak)
                stop 'refs(ipeak) index out of bound; strategy3D_utils :: prob_select_peak'
            endif
            ! stash state index
            if( params_glob%nstates > 1 )then
                states(ipeak) = s3D%proj_space_state(s%iptcl_map,refs(ipeak))
                if( .not. s3D%state_exists(states(ipeak)) )then
                    print *, 'empty states(ipeak): ', states(ipeak)
                    stop 'strategy3D_utils :: prob_select_peak'
                endif
            else
                states(ipeak) = 1
            endif
            ! stash shift (obtained with vector addition)
            shvecs(ipeak,:) = s%prev_shvec
            if( s%doshift ) shvecs(ipeak,:) = shvecs(ipeak,:) + s3D%proj_space_shift(s%iptcl_map,refs(ipeak),1:2)
            where( abs(shvecs(ipeak,:)) < 1e-6 ) shvecs(ipeak,:) = 0.
            ! stash corr
            corrs(ipeak) = s3D%proj_space_corrs(s%iptcl_map,refs(ipeak))
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            ! stash proj index
            projs(ipeak) = s3D%proj_space_proj(s%iptcl_map,refs(ipeak))
            ! stash Euler
            euls(ipeak,:) = s3D%proj_space_euls(s%iptcl_map,refs(ipeak),1:3)
        end do
        if( updatecnt < 5.0 )then
            ! multinomal peak selection
            ! convert correlations to distances
            dists = 1.0 - corrs
            ! scale distances with TAU
            dists = dists / tau
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate probabilities
            pvec = exp(arg4softmax)
            ! normalise
            pvec = pvec / sum(pvec)
            ! sample
            call hpsort(pvec, inds)
            rank = 0
            rnd  = ran3()
            do prob_peak=s%npeaks,1,-1 ! we want to start in the high end
                bound = sum(pvec(prob_peak:s%npeaks))
                rank = rank + 1
                if( rnd <= bound ) exit
            enddo
            prob_peak = inds(prob_peak) ! translate to unsorted
        else
            ! greedy peak selection
            loc = maxloc(corrs)
            prob_peak = loc(1)
        endif
        ! update o_peaks
        call s3D%o_peaks(s%iptcl)%set(1, 'state', real(states(prob_peak)))
        call s3D%o_peaks(s%iptcl)%set(1, 'proj',  real(projs(prob_peak)))
        call s3D%o_peaks(s%iptcl)%set(1, 'corr',  corrs(prob_peak))
        call s3D%o_peaks(s%iptcl)%set(1, 'ow',    1.0)
        call s3D%o_peaks(s%iptcl)%set_euler(1,    euls(prob_peak,:))
        call s3D%o_peaks(s%iptcl)%set_shift(1,    shvecs(prob_peak,:))
    end subroutine prob_select_peak

    subroutine corrsweights( s, corrs, tau, ws, included, best_loc, wcorr )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: corrs(s%npeaks), tau
        real,                   intent(out)   :: ws(s%npeaks), wcorr
        logical,                intent(out)   :: included(s%npeaks)
        integer,                intent(out)   :: best_loc(1)
        integer :: i, ind
        if( s%npeaks == 1 )then
            best_loc(1)  = 1
            ws(1)        = 1.
            included(1)  = .true.
            s%npeaks_eff = 1
            wcorr        = corrs(1)
        else
            best_loc = maxloc(corrs)
            included = .true.
            where( corrs<TINY )included = .false.
            s%npeaks_eff = count(included)
            where( included )
                ws = corrs / sum(corrs, mask=included)
            else where
                ws = 0.
            end where
            wcorr = sum(ws*corrs,mask=included)
        endif
        ! update npeaks individual weights
        call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
    end subroutine corrsweights

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
            ang_lim = deg2rad( build_glob%pgrpsyms%srchrange_theta()/2. )
            o_best  = s3D%o_peaks(s%iptcl)%get_ori(best_loc(1))
            do ipeak = 1,s%npeaks
                if( ipeak == best_loc(1) )then
                    ang_weights(ipeak) = 1.0
                    cycle
                endif
                if( trim(build_glob%pgrpsyms%get_pgrp()).eq.'c1' )then
                    ang_dist = s3D%o_peaks(s%iptcl)%get_ori(ipeak).euldist.o_best
                else
                    call build_glob%pgrpsyms%sym_dists( o_best, s3D%o_peaks(s%iptcl)%get_ori(ipeak), o, ang_dist, inpl_dist)
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
        call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
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
                corrs(ipeak)  = s3D%o_peaks(s%iptcl)%get(ipeak,'corr')
                states(ipeak) = nint(s3D%o_peaks(s%iptcl)%get(ipeak,'state'))
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
            call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
        endif
    end subroutine states_reweight

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
        if( trim(build_glob%pgrpsyms%get_pgrp()).eq.'c1' )then
            ang_sdev = s3D%o_peaks(s%iptcl)%ang_sdev(s%nstates, s%npeaks)
        else
            if( s%npeaks > 2 )then
                sym_os = s3D%o_peaks(s%iptcl)
                do ipeak = 1, s%npeaks
                    if( ipeak == best_loc(1) )cycle
                    call build_glob%pgrpsyms%sym_dists( s3D%o_peaks(s%iptcl)%get_ori(best_loc(1)),&
                        &s3D%o_peaks(s%iptcl)%get_ori(ipeak), osym, dist, inpl_dist)
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
        roind    = pftcc_glob%get_roind(360. - s3D%o_peaks(s%iptcl)%e3get(best_loc(1)))
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
        call build_glob%spproj_field%set(s%iptcl, 'mi_proj',   mi_proj)
        call build_glob%spproj_field%set(s%iptcl, 'mi_inpl',   mi_inpl)
        call build_glob%spproj_field%set(s%iptcl, 'mi_state',  1.)
        call build_glob%spproj_field%set(s%iptcl, 'mi_joint',  mi_joint)
    end subroutine convergence_stats_single

    subroutine convergence_stats_multi( s, best_loc, euldist )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: best_loc(1)
        real,                   intent(in)    :: euldist
        integer :: roind, state
        real    :: mi_proj, mi_inpl, mi_joint, mi_state
        roind    = pftcc_glob%get_roind(360. - s3D%o_peaks(s%iptcl)%e3get(best_loc(1)))
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
        state = nint( s3D%o_peaks(s%iptcl)%get(best_loc(1), 'state') )
        if( .not. s3D%state_exists(state) )then
            print *, 'empty state: ', state
            stop 'strategy3D_utils; convergence_stats_multi'
        endif
        if( s%prev_state == state )then
            mi_state = 1.
            mi_joint = mi_joint + mi_state
        endif
        mi_joint = mi_joint/3.
        ! set the overlaps
        call build_glob%spproj_field%set(s%iptcl, 'mi_proj',   mi_proj)
        call build_glob%spproj_field%set(s%iptcl, 'mi_inpl',   mi_inpl)
        call build_glob%spproj_field%set(s%iptcl, 'mi_state',  mi_state)
        call build_glob%spproj_field%set(s%iptcl, 'mi_joint',  mi_joint)
    end subroutine convergence_stats_multi

end module simple_strategy3D_utils
