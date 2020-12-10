module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: extract_peaks, extract_peaks_dev, calc_ori_weights, update_softmax_weights
public :: states_reweight, estimate_ang_spread, estimate_shift_increment
public :: set_state_overlap, sort_corrs
private
#include "simple_local_flags.inc"

contains

    subroutine extract_peaks( s, corrs, multistates )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: corrs(s%npeaks)
        logical, optional,      intent(in)    :: multistates
        integer :: ipeak, ref, inpl, state
        real    :: shvec(2), shvec_incr(2)
        logical :: l_multistates
        if( present(multistates) )then
            l_multistates = multistates
        else
            l_multistates = .false.
        endif
        do ipeak = 1, s%npeaks
            ref = s3D%proj_space_refinds_sorted(s%ithr, s%nrefsmaxinpl - s%npeaks + ipeak)
            if( ref < 1 .or. ref > s%nrefs ) THROW_HARD('ref index: '//int2str(ref)//' out of bound; extract_peaks')
            if( l_multistates )then
                state = s3D%proj_space_state(ref)
                if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; extract_peaks')
            endif
            inpl = s3D%proj_space_inplinds_sorted(s%ithr, s%nrefsmaxinpl - s%npeaks + ipeak)
            ! add shift
            shvec      = s%prev_shvec
            shvec_incr = 0.
            if( s%doshift ) then
                shvec_incr = s3D%proj_space_shift(s%ithr,ref,inpl,1:2)
                shvec      = shvec + shvec_incr
            end if
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = s3D%proj_space_corrs(s%ithr,ref,inpl)
            if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
                if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            end if
            if( l_multistates )then
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', real(state))
            else
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', 1.)
            endif
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'proj',  real(s3D%proj_space_proj(ref)))
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'inpl',  real(inpl))
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call s3D%o_peaks(s%iptcl)%set_euler(ipeak, s3D%proj_space_euls(s%ithr,ref,inpl,1:3))
            call s3D%o_peaks(s%iptcl)%set_shift(ipeak, shvec)
            call s3D%o_peaks(s%iptcl)%set_shift_incr(ipeak, shvec_incr)
        enddo
    end subroutine extract_peaks

    subroutine extract_peaks_dev( s, corrs, multistates )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: corrs(s%npeaks)
        logical, optional,      intent(in)    :: multistates
        integer :: order(s3D%inplpeaks%n), ipeak, iref, inpl, state, ind, i,j
        real    :: euls(3), shvec(2), shvec_incr(2)
        logical :: l_multistates
        if( present(multistates) )then
            l_multistates = multistates
        else
            l_multistates = .false.
        endif
        if( l_multistates ) THROW_HARD('extract_peaks_dev: unsupported multi-states')
        do i = 1,s3D%inplpeaks%n
            order(i) = i
        enddo
        call hpsort(s3D%inplpeaks%ccs(:,s%ithr), order)
        ! stash best s%npeaks
        ipeak = 0
        do i = s3D%inplpeaks%n,s3D%inplpeaks%n-s%npeaks+1,-1
            ind   = order(i)
            ipeak = ipeak + 1
            ! euler angles
            iref = s3D%inplpeaks%eulinds(1,ind,s%ithr)
            if( iref < 1 .or. iref > s%nrefs )then
                print *,'iptcl i ind ithr: ',s%iptcl, i, ind, s%ithr
                print *,'eulinds: ',s3D%inplpeaks%eulinds(:,ind,s%ithr)
                print *,'shifts:  ',s3D%inplpeaks%shifts(:,ind,s%ithr)
                do j = 1,s3D%inplpeaks%n
                    ind = order(j)
                    print *,j,s3D%inplpeaks%ccs(j,s%ithr),order(j),s3D%inplpeaks%eulinds(:,ind,s%ithr),s3D%inplpeaks%shifts(:,ind,s%ithr),s3D%inplpeaks%ccs(j,s%ithr)
                enddo
                ! do j = params_glob%kfromto(1),params_glob%kstop
                !     print *,j,csq_fast(pftcc_glob%heap_vars(s%ithr)%pft_ref(:,j)),pftcc_glob%heap_vars(s%ithr)%pft_ref(:,j)
                ! enddo
                ! ind = pftcc_glob%pinds(s%iptcl)
                do j = params_glob%kfromto(1),params_glob%kstop
                    ! print *,j,sum(csq_fast(pftcc_glob%pfts_ptcls(:,j,ind))),pftcc_glob%pfts_ptcls(:,j,ind)
                enddo
                THROW_HARD('ref index: '//int2str(iref)//' out of bound; extract_peaks_dev')
            endif
            if( l_multistates )then
                state = s3D%proj_space_state(iref)
                if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; extract_peaks')
            endif
            euls    = build_glob%eulspace%get_euler(s3D%proj_space_proj(iref))
            inpl    = s3D%inplpeaks%eulinds(2,ind,s%ithr)
            euls(3) = 360. - pftcc_glob%get_rot(inpl)
            ! add shift
            shvec = s%prev_shvec
            shvec_incr = 0.
            if( s%doshift ) then
                shvec_incr = s3D%inplpeaks%shifts(:,ind,s%ithr)
                shvec      = shvec + shvec_incr
            end if
            ! transfer to solution set
            corrs(ipeak) = s3D%inplpeaks%ccs(i,s%ithr)
            if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
                if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            end if
            if( l_multistates )then
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', real(state))
            else
                call s3D%o_peaks(s%iptcl)%set(ipeak, 'state', 1.)
            endif
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'proj',  real(s3D%proj_space_proj(iref)))
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'inpl',  real(inpl))
            call s3D%o_peaks(s%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call s3D%o_peaks(s%iptcl)%set_euler(ipeak, euls)
            call s3D%o_peaks(s%iptcl)%set_shift(ipeak, shvec)
            call s3D%o_peaks(s%iptcl)%set_shift_incr(ipeak, shvec_incr)
        enddo
    end subroutine extract_peaks_dev

    ! (1) Calculate normalised softmax weights per particle
    ! (2) Use Otsu's algorithm to detect correlation peaks
    subroutine calc_ori_weights( s, corrs, ws, best_loc, best_corr )
        use simple_ori, only: ori
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(in)    :: corrs(s%npeaks)
        real,                   intent(out)   :: ws(s%npeaks), best_corr
        integer,                intent(out)   :: best_loc(1)
        real,    allocatable :: ws_nonzero(:)
        real    :: wsum, thres
        s%npeaks_eff = 1
        if( s%npeaks == 1 )then
            best_loc(1)  = 1
            ws(1)        = 1.
            best_corr    = corrs(best_loc(1))
        else
            ! find highest corr pos
            best_loc  = maxloc(corrs)
            best_corr = corrs(best_loc(1))
            ! calculate weights
            call calc_ori_weights_here(s%npeaks, corrs, ws)
            if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
                ! define a threshold using Otsu's algorithm
                ws_nonzero = pack(ws, mask=ws > TINY)
                call otsu(ws_nonzero, thres)
                s%npeaks_eff = max(1,count(ws_nonzero > thres))
                ! limit the number of peaks
                if( s%npeaks_eff > MAXNPEAKS )then
                    call hpsort(ws_nonzero)
                    call reverse(ws_nonzero) ! largest first
                    thres = ws_nonzero(MAXNPEAKS + 1)
                    s%npeaks_eff = MAXNPEAKS
                endif
                ! zero weights below the threshold and re-normalize
                if( s%npeaks_eff == 1 )then
                    where(ws < maxval(ws) ) ws = 0. ! always one nonzero weight
                else
                    where(ws <= thres) ws = 0.
                endif
                wsum = sum(ws)
                if( wsum > TINY )then
                    ws = ws / wsum
                else
                    ws = 0.
                endif
                select case(params_glob%wcrit_enum)
                    case(RANK_SUM_CRIT,RANK_CEN_CRIT,RANK_EXP_CRIT,RANK_INV_CRIT)
                        if( params_glob%wcrit_enum == RANK_EXP_CRIT )then
                            call conv2rank_weights(s%npeaks, ws, params_glob%wcrit_enum, RANKW_EXP)
                        else
                            call conv2rank_weights(s%npeaks, ws, params_glob%wcrit_enum)
                        endif
                end select
            else
                s%npeaks_eff = max(1,count(ws > TINY))
                ! limit the number of peaks
                if( s%npeaks_eff > MAXNPEAKS )then
                    ws_nonzero   = pack(ws, mask=ws > TINY)
                    call hpsort(ws_nonzero)
                    call reverse(ws_nonzero) ! largest first
                    thres = ws_nonzero(MAXNPEAKS + 1)
                    where(ws <= thres) ws = 0.
                else if( s%npeaks_eff == 1 )then
                    ws = 0.
                    ws((best_loc(1))) = 1.
                endif
                ! wsum = sum(ws)
                ! if( wsum > TINY )then
                !     ws = ws / wsum
                ! else
                !     ws = 0.
                ! endif
            endif
        endif
        ! update npeaks individual weights
        call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
    end subroutine calc_ori_weights

    subroutine update_softmax_weights( iptcl, npeaks )
        integer, intent(in) :: iptcl, npeaks
        real, allocatable   :: corrs(:)
        real :: ws(npeaks)
        corrs = s3D%o_peaks(iptcl)%get_all('corr')
        call calc_ori_weights_here(npeaks, corrs, ws )
        call s3D%o_peaks(iptcl)%set_all('ow', ws)
    end subroutine update_softmax_weights

    subroutine calc_ori_weights_here( npeaks, corrs, ws )
        integer, intent(in)    :: npeaks
        real,    intent(in)    :: corrs(npeaks)
        real,    intent(inout) :: ws(npeaks)
        real    :: dists(npeaks), arg4softmax(npeaks)
        real    :: wsum
        if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            ! subtracts minimum distance
            dists = corrs / params_glob%sigma2_fudge
            dists = dists - maxval(dists)
            ! exponential weights
            ws = exp(dists)
        else
            ! convert correlations to distances
            dists = 1.0 - corrs
            ! scale distances with TAU
            dists = dists / params_glob%tau
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate softmax weights
            ws = exp(arg4softmax)
            where( corrs <= TINY ) ws = 0.
        endif
        ! critical for performance to normalize here as well
        wsum = sum(ws)
        if( wsum > TINY )then
            ws = ws / wsum
        else
            ws = 0.
        endif
    end subroutine calc_ori_weights_here

    subroutine states_reweight( s, ws, state, best_loc )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(inout) :: ws(s%npeaks)
        integer,                intent(out)   :: state, best_loc(1)
        integer :: istate, ipeak, states(s%npeaks)
        real    :: corrs(s%npeaks), state_ws(s%nstates), wsum
        if( s%npeaks > 1 .and. s%nstates > 1 )then
            ! states weights
            do ipeak = 1, s%npeaks
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
            where( .not. states==state ) ws = 0.
            ! critical for performance to normalize here as well
            wsum = sum(ws)
            if( wsum > TINY )then
                ws = ws / wsum
            else
                ws = 0.
            endif
            best_loc = maxloc(ws)
            ! update individual weights
            call s3D%o_peaks(s%iptcl)%set_all('ow', ws)
        else if( s%npeaks == 1 .and. s%nstates > 1 )then
            best_loc = 1
            state    = s3D%o_peaks(s%iptcl)%get_state(1)
        endif
    end subroutine states_reweight

    function estimate_ang_spread( s ) result( ang_spread )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_srch), intent(inout) :: s
        integer    :: ipeak, jpeak, states(s3D%o_peaks(s%iptcl)%get_noris())
        integer    :: best_state, loc(1), npeaks, cnt, i
        type(ori)  :: osym, o1, o2
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
                call s3D%o_peaks(s%iptcl)%get_ori(ipeak, o1)
                call s3D%o_peaks(s%iptcl)%get_ori(jpeak, o2)
                call build_glob%pgrpsyms%sym_dists( o1, o2, osym, dist, inpl_dist)
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
        call osym%kill
        call o1%kill
        call o2%kill
    end function estimate_ang_spread

    subroutine estimate_shift_increment( s, shwmean, shwstdev )
        class(strategy3D_srch), intent(inout) :: s
        real,                   intent(out)   :: shwmean, shwstdev
        integer    :: ipeak, states(s3D%o_peaks(s%iptcl)%get_noris())
        integer    :: best_state, npeaks, cnt, i
        real       :: ws(s3D%o_peaks(s%iptcl)%get_noris()), dev_w, var_w
        real       :: shift_incrs(s3D%o_peaks(s%iptcl)%get_noris()), ws_here(s3D%o_peaks(s%iptcl)%get_noris())
        logical    :: multi_states
        shwmean  = 0.
        shwstdev = 0.
        npeaks   = s3D%o_peaks(s%iptcl)%get_noris()
        if( npeaks < 1 ) return ! need at least 1
        ! gather weights & states
        do ipeak=1,npeaks
            ws(ipeak)     = s3D%o_peaks(s%iptcl)%get(ipeak, 'ow')
            states(ipeak) = s3D%o_peaks(s%iptcl)%get_state(ipeak)
        enddo
        if( count(ws > TINY) < 1 ) return ! need at least 1
        ! multi-states or not?
        multi_states = s%nstates > 1
        ! best state is?
        best_state = states(maxloc(ws,dim=1))
        if( multi_states )then
            if( count(states == best_state) < 1 ) return ! need at least 1
        endif
        ! loop over valid peaks and find shift increment
        cnt = 0
        do ipeak = 1, max(1, npeaks - 1)
            ! take care of ow > 0. requirement
            if( ws(ipeak) <= TINY ) cycle
            ! take care of multi-state case
            if( multi_states )then
                if( states(ipeak) /= best_state ) cycle
            endif
            ! incr shift increment counter
            cnt = cnt + 1
            ws_here(cnt)     = ws(ipeak)
            shift_incrs(cnt) = arg(s3D%o_peaks(s%iptcl)%get_2Dshift_incr(ipeak))
        enddo
        if( cnt < 1 ) return ! need at least 1
        ! calculate standard deviation of distances as measure of angular spread
        if( cnt == 1)then
            shwmean  = shift_incrs(1)
            shwstdev = 0.
        else
            shwmean = sum(shift_incrs(:cnt)*ws_here(:cnt)) / sum(ws_here(:cnt))
            var_w   = 0.
            do i=1,cnt
                dev_w = shift_incrs(i) - shwmean
                var_w = var_w + ws_here(i) * dev_w * dev_w
            end do
            var_w    = var_w/((real(cnt) - 1.)/real(cnt)*sum(ws_here(:cnt))) ! corrected two-pass formula
            if( var_w > 0. ) shwstdev = sqrt(var_w)
        endif
    end subroutine estimate_shift_increment

    subroutine set_state_overlap( s, best_loc )
        class(strategy3D_srch), intent(inout) :: s
        integer,                intent(in)    :: best_loc(1)
        integer :: state
        real    :: mi_state
        mi_state = 0.
        state = nint( s3D%o_peaks(s%iptcl)%get(best_loc(1), 'state') )
        if( .not. s3D%state_exists(state) )then
            THROW_HARD('empty state: '//int2str(state)//' set_state_overlap')
        endif
        if( s%prev_state == state ) mi_state = 1.
        call build_glob%spproj_field%set(s%iptcl, 'mi_state',  mi_state)
    end subroutine set_state_overlap

    subroutine sort_corrs( s )
        class(strategy3D_srch), intent(inout) :: s
        real    :: corrs(s%nrefs*params_glob%ninplpeaks), corrs_highest(s%nrefs)
        integer :: proj_space_tmp(s%nrefs*params_glob%ninplpeaks)
        integer :: i, j, arange(2), idx_array(s%nrefs*params_glob%ninplpeaks)
        integer :: asequence(params_glob%ninplpeaks)
        real    :: areal
        asequence = (/(j, j=1,params_glob%ninplpeaks)/)
        do i = 1,s%nrefs
            arange(1) = (i-1)*params_glob%ninplpeaks+1
            arange(2) = i*params_glob%ninplpeaks
            s3D%proj_space_refinds_sorted        (s%ithr,arange(1):arange(2)) = i
            s3D%proj_space_inplinds_sorted       (s%ithr,arange(1):arange(2)) = asequence(:)
            s3D%proj_space_refinds_sorted_highest(s%ithr,                  i) = i
            if (.not. s3D%proj_space_corrs_srchd(s%ithr,i)) then
                corrs(arange(1):arange(2)) = -HUGE(areal)
                corrs_highest(i)           = -HUGE(areal)
            else
                corrs(arange(1):arange(2)) = s3D%proj_space_corrs(s%ithr,i,:)
                corrs_highest(i)           = s3D%proj_space_corrs(s%ithr,i,1)
            end if
        end do
        do j = 1,s%nrefs*params_glob%ninplpeaks
            idx_array(j) = j
        end do
        call hpsort(corrs, idx_array)
        proj_space_tmp(:) = s3D%proj_space_refinds_sorted(s%ithr,:)
        do j = 1,s%nrefs*params_glob%ninplpeaks
            s3D%proj_space_refinds_sorted(s%ithr,j) = proj_space_tmp(idx_array(j))
        end do
        proj_space_tmp(:) = s3D%proj_space_inplinds_sorted(s%ithr,:)
        do j = 1,s%nrefs*params_glob%ninplpeaks
            s3D%proj_space_inplinds_sorted(s%ithr,j) = proj_space_tmp(idx_array(j))
        end do
        call hpsort(corrs_highest, s3D%proj_space_refinds_sorted_highest(s%ithr, :))
    end subroutine sort_corrs

end module simple_strategy3D_utils
