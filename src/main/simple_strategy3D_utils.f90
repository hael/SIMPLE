module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

contains

    subroutine extract_peaks( s, corrs, multistates )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: corrs(s%npeaks * MAXNINPLPEAKS)
        logical, optional,      intent(in)    :: multistates
        integer :: ipeak, cnt, ref, inpl, state
        real    :: shvec(2)
        logical :: l_multistates
        if( present(multistates) )then
            l_multistates = multistates
        else
            l_multistates = .false.
        endif
        cnt = 0
        do ipeak = 1, s%npeaks
            ref = s3D%proj_space_refinds(s%ithr, s%nrefs - s%npeaks + ipeak)
            if( ref < 1 .or. ref > s%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; strategy3D_utils :: extract_peaks'
            endif
            if( l_multistates )then
                state = s3D%proj_space_state(ref)
                if( .not. s3D%state_exists(state) )then
                    print *, 'empty state: ', state
                    stop 'strategy3D_utils :: extract_peak'
                endif
            endif
            do inpl=1,MAXNINPLPEAKS
                cnt = cnt + 1
                ! add shift
                shvec = s%prev_shvec
                if( s%doshift ) shvec = shvec + s3D%proj_space_shift(s%ithr,ref,inpl,1:2)
                where( abs(shvec) < 1e-6 ) shvec = 0.
                ! transfer to solution set
                corrs(cnt) = s3D%proj_space_corrs(s%ithr,ref,inpl)
                if (.not. pftcc_glob%is_euclid(s%iptcl)) then
                    if( corrs(cnt) < 0. ) corrs(cnt) = 0.
                end if
                if( l_multistates )then
                    call s3D%o_peaks(s%iptcl)%set(cnt, 'state', real(state))
                else
                    call s3D%o_peaks(s%iptcl)%set(cnt, 'state', 1.)
                endif
                call s3D%o_peaks(s%iptcl)%set(cnt, 'proj',  real(s3D%proj_space_proj(ref)))
                call s3D%o_peaks(s%iptcl)%set(cnt, 'corr',  corrs(cnt))
                call s3D%o_peaks(s%iptcl)%set_euler(cnt, s3D%proj_space_euls(s%ithr,ref,inpl,1:3))
                call s3D%o_peaks(s%iptcl)%set_shift(cnt, shvec)
            end do
        enddo
    end subroutine extract_peaks

    ! assumes hard orientation assignment. Hence inpl=1
    subroutine prob_select_peak( s, tau, updatecnt )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: tau, updatecnt
        real    :: corrs(s%npeaks), shvecs(s%npeaks,2), euls(s%npeaks,3), pvec(s%npeaks)
        integer :: inds(s%npeaks), refs(s%npeaks), states(s%npeaks), projs(s%npeaks)
        integer :: prob_peak, ipeak, rank, loc(1)
        real    :: bound, rnd, dists(s%npeaks), arg4softmax(s%npeaks)
        do ipeak = 1, s%npeaks
            ! stash peak index
            inds(ipeak) = ipeak
            ! stash reference index
            refs(ipeak) = s3D%proj_space_refinds(s%ithr, s%nrefs - s%npeaks + ipeak)
            if( refs(ipeak) < 1 .or. refs(ipeak) > s%nrefs )then
                print *, 'refs(ipeak): ', refs(ipeak)
                stop 'refs(ipeak) index out of bound; strategy3D_utils :: prob_select_peak'
            endif
            ! stash state index
            if( params_glob%nstates > 1 )then
                states(ipeak) = s3D%proj_space_state(refs(ipeak))
                if( .not. s3D%state_exists(states(ipeak)) )then
                    print *, 'empty states(ipeak): ', states(ipeak)
                    stop 'strategy3D_utils :: prob_select_peak'
                endif
            else
                states(ipeak) = 1
            endif
            ! stash shift (obtained with vector addition)
            shvecs(ipeak,:) = s%prev_shvec
            if( s%doshift ) shvecs(ipeak,:) = shvecs(ipeak,:) + s3D%proj_space_shift(s%ithr,refs(ipeak),1,1:2)
            where( abs(shvecs(ipeak,:)) < 1e-6 ) shvecs(ipeak,:) = 0.
            ! stash corr
            corrs(ipeak) = s3D%proj_space_corrs(s%ithr,refs(ipeak),1)
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            ! stash proj index
            projs(ipeak) = s3D%proj_space_proj(refs(ipeak))
            ! stash Euler
            euls(ipeak,:) = s3D%proj_space_euls(s%ithr,refs(ipeak),1,1:3)
        end do
        if( updatecnt < 5.0 )then
            if (.not. pftcc_glob%is_euclid(s%iptcl)) then
                ! multinomal peak selection
                ! convert correlations to distances
                dists = 1.0 - corrs
                ! scale distances with TAU
                dists = dists / tau
            else
                dists = - corrs / params_glob%sigma2_fudge
            end if
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

    subroutine corrs2softmax_weights( s, n, corrs, tau, ws, included, best_loc, wcorr )
        use simple_ori, only: ori
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: n
        real,                   intent(in)    :: corrs(n), tau
        real,                   intent(out)   :: ws(n), wcorr
        logical,                intent(out)   :: included(n)
        integer,                intent(out)   :: best_loc(1)
        type(ori) :: o, o_best
        real      :: dists(n), arg4softmax(n), ang_dist, inpl_dist
        integer   :: ipeak
        if( n == 1 )then
            best_loc(1)  = 1
            ws(1)        = 1.
            included(1)  = .true.
            s%npeaks_eff = 1
            wcorr        = corrs(1)
        else
            best_loc = maxloc(corrs)
            if (.not. pftcc_glob%is_euclid(s%iptcl)) then
                ! convert correlations to distances
                dists = 1.0 - corrs
                ! scale distances with TAU
                dists = dists / tau
            else
                dists = - corrs / params_glob%sigma2_fudge
            end if
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate softmax weights
            ws = exp(arg4softmax)
            if (.not. pftcc_glob%is_euclid(s%iptcl)) then
                where( corrs <= TINY ) ws = 0.
            end if
            ! normalise
            ws = ws / sum(ws)
            ! threshold weights
            included = (ws >= SOFTMAXW_THRESH)
            s%npeaks_eff = count(included)
            if( s%npeaks_eff > 0 )then
                ! exclusion
                where( .not. included ) ws = 0.
                ! weighted corr
                wcorr = sum(ws*corrs,mask=included)
            else
                ! defaults to best when all fall below threshold
                s%npeaks_eff          = 1
                wcorr                 = corrs(best_loc(1))
                ws                    = 0.
                ws(best_loc(1))       = 1.
                included(best_loc(1)) = .true.
            endif
        endif
        ! update npeaks individual weights
        call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
    end subroutine corrs2softmax_weights

    subroutine states_reweight( s, npeaks, ws, included, state, best_loc, wcorr )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: npeaks
        real,                   intent(inout) :: ws(npeaks)
        real,                   intent(out)   :: wcorr
        logical,                intent(inout) :: included(npeaks)
        integer,                intent(out)   :: state, best_loc(1)
        integer :: istate, ipeak, states(npeaks)
        real    :: corrs(npeaks), state_ws(s%nstates)
        if( npeaks > 1 .and. s%nstates > 1 )then
            ! states weights
            do ipeak = 1, npeaks
                corrs(ipeak)  = s3D%o_peaks(s%iptcl)%get(ipeak,'corr')
                states(ipeak) = s3D%o_peaks(s%iptcl)%get_state(ipeak)
            enddo
            ! greedy state assignment
            state_ws = 0.
            do istate = 1, s%nstates
                state_ws(istate) = sum(ws,mask=(states==istate))
            enddo
            state = maxloc(state_ws, dim=1)
            ! in-state re-weighing
            included     = included .and. (states==state)
            s%npeaks_eff = count(included)
            where( .not. included ) ws = 0.
            ws       = ws / sum(ws, mask=included)
            best_loc = maxloc(ws)
            ! weighted corr
            wcorr = sum(ws*corrs, mask=included)
            ! update individual weights
            call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
        else if( npeaks == 1 .and. s%nstates > 1 )then
            best_loc = 1
            state    = s3D%o_peaks(s%iptcl)%get_state(1)
            wcorr    = s3D%o_peaks(s%iptcl)%get(1,'corr')
        endif
    end subroutine states_reweight

    function estimate_ang_spread( s ) result( ang_spread )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_srch), intent(inout) :: s
        integer    :: ipeak, jpeak, states(s3D%o_peaks(s%iptcl)%get_noris())
        integer    :: best_state, loc(1), npeaks, cnt, i
        type(ori)  :: osym
        real       :: ang_spread, dist, inpl_dist, ws(s3D%o_peaks(s%iptcl)%get_noris()), ave, dev, var, ep
        real       :: dists(s3D%o_peaks(s%iptcl)%get_noris() * (s3D%o_peaks(s%iptcl)%get_noris() - 1) / 2)
        logical    :: multi_states
        ang_spread = 0.
        npeaks     = s3D%o_peaks(s%iptcl)%get_noris()
        if( npeaks < 3 ) return ! need at least 3
        ! gather weights & states
        do ipeak=1,npeaks
            ws(ipeak)     = s3D%o_peaks(s%iptcl)%get(ipeak, 'ow')
            states(ipeak) = s3D%o_peaks(s%iptcl)%get_state(ipeak)
        enddo
        if( count(ws > TINY) < 3 ) return ! need at least 3
        ! multi-states or not?
        multi_states = .false.
        if( s%nstates > 1 ) multi_states = .true.
        ! best state is?
        loc = maxloc(ws)
        best_state = states(loc(1))
        if( multi_states )then
            if( count(states == best_state) < 3 ) return ! need at least 3
        endif
        ! loop over valid peaks and find maximum angular distance as metric for angular spread
        cnt = 0
        do ipeak = 1, npeaks - 1
            ! take care of ow > 0. requirement
            if( ws(ipeak) <= TINY ) cycle
            ! take care of multi-state case
            if( multi_states )then
                if( states(ipeak) /= best_state ) cycle
            endif
            do jpeak = ipeak + 1, npeaks
                ! take care of ow > 0. requirement
                if( ws(jpeak) <= TINY ) cycle
                ! take care of multi-state case
                if( multi_states )then
                    if( states(jpeak) /= best_state ) cycle
                endif
                ! incr dist counter
                cnt = cnt + 1
                ! to minimize the angular distance over the symmetry operations
                ! this routine takes care of the c1 case
                call build_glob%pgrpsyms%sym_dists( s3D%o_peaks(s%iptcl)%get_ori(ipeak),&
                    &s3D%o_peaks(s%iptcl)%get_ori(jpeak), osym, dist, inpl_dist)
                dists(cnt) = dist
            enddo
        enddo
        if( cnt < 3 ) return ! need at least 3
        ! calculate standard deviation of distances as measure of angular spread
        ave = sum(dists(:cnt)) / real(cnt)
        ep  = 0.
        var = 0.
        do i=1,cnt
            dev = dists(i) - ave
            ep  = ep + dev
            var = var + dev * dev
        end do
        var = (var - ep * ep / real(cnt))/(real(cnt) - 1.) ! corrected two-pass formula
        ang_spread = 0.
        if( var > 0. ) ang_spread = sqrt(var)
    end function estimate_ang_spread

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

    subroutine sort_corrs( s )
        class(strategy3D_srch), intent(inout) :: s
        real                                  :: corrs(s%nrefs)
        real                                  :: areal
        corrs = s3D%proj_space_corrs(s%ithr,:,1)
        where(.not. s3D%proj_space_corrs_srchd(s%ithr,:)) corrs = -HUGE(areal)
        call hpsort(corrs, s3D%proj_space_refinds(s%ithr,:))
    end subroutine sort_corrs

end module simple_strategy3D_utils
