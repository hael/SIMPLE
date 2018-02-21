module simple_strategy3D_greedy_multi
use simple_strategy3D_alloc
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
implicit none

public :: strategy3D_greedy_multi
private

logical, parameter :: DEBUG = .false.

type, extends(strategy3D) :: strategy3D_greedy_multi
    type(strategy3D_srch) :: s
contains
    procedure          :: new         => new_greedy_multi
    procedure          :: srch        => srch_greedy_multi
    procedure          :: kill        => kill_greedy_multi
    procedure          :: oris_assign => oris_assign_greedy_multi
end type strategy3D_greedy_multi

contains

    subroutine new_greedy_multi( self, spec )
        class(strategy3D_greedy_multi), intent(inout) :: self
        class(strategy3D_spec),          intent(inout) :: spec
        call self%s%new( spec )
    end subroutine new_greedy_multi

    !>  \brief  greedy hill-climbing
    subroutine srch_greedy_multi( self, spec )
        class(strategy3D_greedy_multi), intent(inout) :: self
        class(strategy3D_spec),          intent(inout) :: spec
        integer :: iref,isample,nrefs,target_projs(self%s%npeaks_grid), inpl_ind(1), state
        real    :: corrs(self%s%nrefs), inpl_corrs(self%s%nrots), inpl_corr
        if( self%s%a_ptr%get_state(self%s%iptcl) > 0 )then
            if( self%s%neigh )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. associated(spec%grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: greedy_srch'
                call self%s%greedy_subspace_srch(spec%grid_projs, target_projs)
                ! initialize
                call self%s%prep4srch(spec%nnmat, target_projs)
                nrefs = self%s%nnnrefs
            else
                ! initialize
                call self%s%prep4srch()
                nrefs = self%s%nrefs
            endif
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            proj_space_corrs(self%s%iptcl_map,:) = -1.
            ! search
            do isample=1,nrefs
                iref = srch_order(self%s%iptcl_map,isample) ! set the reference index
                call per_ref_srch                         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%s%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%s%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%s%iptcl_map,:))
            ! take care of the in-planes
            call self%s%inpl_srch
            ! prepare weights & orientation
            call self%oris_assign()
        else
            call self%s%a_ptr%reject(self%s%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch
                state = proj_space_state(self%s%iptcl_map,iref)
                if( state_exists(state) )then
                    call self%s%pftcc_ptr%gencorrs(iref, self%s%iptcl, inpl_corrs) ! In-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)                            ! greedy in-plane index
                    inpl_corr  = inpl_corrs(inpl_ind(1))                       ! max in plane correlation
                    call self%s%store_solution(iref, iref, inpl_ind(1), inpl_corr)
                endif
            end subroutine per_ref_srch

    end subroutine srch_greedy_multi

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine oris_assign_greedy_multi( self )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(strategy3D_greedy_multi), intent(inout) :: self
        type(ori)  :: osym
        type(oris) :: sym_os
        real       :: shvec(2), corrs(self%s%npeaks), ws(self%s%npeaks), dists(self%s%npeaks)
        real       :: arg4softmax(self%s%npeaks), state_ws(self%s%nstates), bfacs(self%s%npeaks)
        real       :: mi_proj, mi_inpl, mi_state, dist_inpl, wcorr, frac
        real       :: ang_sdev, dist, inpl_dist, euldist, mi_joint, bfac
        integer    :: best_loc(1), loc(1), states(self%s%npeaks)
        integer    :: s, ipeak, cnt, ref, state, roind, neff_states
        logical    :: included(self%s%npeaks)
        ! empty states
        neff_states = count(state_exists)
        ! init npeaks
        do ipeak = 1, self%s%npeaks
            cnt = self%s%nrefs - self%s%npeaks + ipeak
            ref = proj_space_inds(self%s%iptcl_map, cnt)
            if( ref < 1 .or. ref > self%s%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            endif
            state = proj_space_state(self%s%iptcl_map,ref)
            if( .not. state_exists(state) )then
                print *, 'empty state: ', state
                stop 'simple_prime3D_srch::prep_npeaks_oris'
            endif
            ! add shift
            shvec = self%s%prev_shvec
            if( self%s%doshift )shvec = shvec + proj_space_shift(self%s%iptcl_map,ref,1:2)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = proj_space_corrs(self%s%iptcl_map,ref)
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            call o_peaks(self%s%iptcl)%set(ipeak, 'state', real(state))
            call o_peaks(self%s%iptcl)%set(ipeak, 'proj',  real(proj_space_proj(self%s%iptcl_map,ref)))
            call o_peaks(self%s%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call o_peaks(self%s%iptcl)%set_euler(ipeak, proj_space_euls(self%s%iptcl_map,ref,1:3))
            call o_peaks(self%s%iptcl)%set_shift(ipeak, shvec)
        enddo
        best_loc = maxloc(corrs)
        ! stochastic weights
        if( self%s%npeaks == 1 )then
            ws(1) = 1.
            call o_peaks(self%s%iptcl)%set(1,'ow',1.0)
            wcorr = o_peaks(self%s%iptcl)%get(1,'corr')
        else
            ! convert correlations to distances
            dists = 1.0 - corrs
            ! scale distances with TAU
            dists = dists / TAU
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate softmax weights
            ws = exp(arg4softmax)
            ws = ws / sum(ws)
            ! threshold weights
            included = (ws >= SOFTMAXW_THRESH)
            self%s%npeaks_eff = count(included)
            where( .not. included ) ws = 0.
            ! weighted corr
            wcorr = sum(ws*corrs,mask=included)
            ! update npeaks individual weights
            call o_peaks(self%s%iptcl)%set_all('ow', ws)
        endif
        if( self%s%npeaks > 1 .and. self%s%nstates > 1 )then
            ! states weights
            do ipeak = 1, self%s%npeaks
                corrs(ipeak)  = o_peaks(self%s%iptcl)%get(ipeak,'corr')
                states(ipeak) = nint(o_peaks(self%s%iptcl)%get(ipeak,'state'))
            enddo
            ! greedy state assignment
            do s = 1, self%s%nstates
                state_ws(s) = sum(ws,mask=(states==s))
            enddo
            loc      = maxloc(state_ws)
            state    = loc(1)
            ! in-state re-weighing
            included = included .and. (states==state)
            where( .not.included )ws = 0.
            ws       = ws / sum(ws, mask=included)
            best_loc = maxloc(ws)
            ! weighted corr
            wcorr    = sum(ws*corrs, mask=included)
            ! update npeaks individual weights
            call o_peaks(self%s%iptcl)%set_all('ow', ws)
        endif
        ! B factors
        if( self%s%pftcc_ptr%objfun_is_ccres() )then
            bfacs = 0.
            do ipeak = 1, self%s%npeaks
                if( ws(ipeak) > TINY .or. self%s%npeaks==1 )then
                    cnt   = self%s%nrefs - self%s%npeaks + ipeak
                    ref   = proj_space_inds(self%s%iptcl_map, cnt)
                    shvec = 0.
                    if( self%s%doshift )shvec = proj_space_shift(self%s%iptcl_map, ref, 1:2)
                    roind        = self%s%pftcc_ptr%get_roind(360. - proj_space_euls(self%s%iptcl_map, ref, 3))
                    bfacs(ipeak) = self%s%pftcc_ptr%fit_bfac(ref, self%s%iptcl, roind, shvec)
                else
                    bfacs(ipeak) = 0.
                endif
                call o_peaks(self%s%iptcl)%set(ipeak, 'bfac', bfacs(ipeak))
            enddo
            bfac = sum(ws * bfacs, mask=(ws>TINY))
            call self%s%a_ptr%set(self%s%iptcl, 'bfac',  bfac )
        endif
        ! angular standard deviation
        ang_sdev = 0.
        if( trim(self%s%se_ptr%get_pgrp()).eq.'c1' )then
            ang_sdev = o_peaks(self%s%iptcl)%ang_sdev(self%s%nstates, self%s%npeaks)
        else
            if( self%s%npeaks > 2 )then
                sym_os = o_peaks(self%s%iptcl)
                do ipeak = 1, self%s%npeaks
                    if( ipeak == best_loc(1) )cycle
                    call self%s%se_ptr%sym_dists( o_peaks(self%s%iptcl)%get_ori(best_loc(1)),&
                        &o_peaks(self%s%iptcl)%get_ori(ipeak), osym, dist, inpl_dist)
                    call sym_os%set_ori(ipeak, osym)
                enddo
                ang_sdev = sym_os%ang_sdev(self%s%nstates, self%s%npeaks)
            endif
        endif
        ! angular distances
        call self%s%se_ptr%sym_dists( self%s%a_ptr%get_ori(self%s%iptcl),&
            &o_peaks(self%s%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl)
        ! convergence parameters
        roind    = self%s%pftcc_ptr%get_roind(360. - o_peaks(self%s%iptcl)%e3get(best_loc(1)))
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj = 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%s%prev_roind == roind )then
            mi_inpl  = 1.
            mi_joint = mi_joint + 1.
        endif
        ! states convergence
        mi_state = 0.
        state = nint( o_peaks(self%s%iptcl)%get(best_loc(1), 'state') )
        if( .not. state_exists(state) )then
            print *, 'empty state: ', state
            stop 'simple_prime3d_srch; update_best'
        endif
        if(self%s%nstates > 1)then
            if( self%s%prev_state == state )then
                mi_state = 1.
                mi_joint = mi_joint + mi_state
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        ! fraction search space
        if( self%s%neigh )then
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nnn * neff_states)
        else
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs * neff_states)
        endif
        ! set the overlaps
        call self%s%a_ptr%set(self%s%iptcl, 'mi_proj',  mi_proj )
        call self%s%a_ptr%set(self%s%iptcl, 'mi_inpl',  mi_inpl )
        call self%s%a_ptr%set(self%s%iptcl, 'mi_state', mi_state)
        call self%s%a_ptr%set(self%s%iptcl, 'mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call self%s%a_ptr%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*self%s%a_ptr%get(self%s%iptcl,'dist'))
        call self%s%a_ptr%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call self%s%a_ptr%set_euler(self%s%iptcl, o_peaks(self%s%iptcl)%get_euler(best_loc(1)))
        call self%s%a_ptr%set_shift(self%s%iptcl, o_peaks(self%s%iptcl)%get_2Dshift(best_loc(1)))
        call self%s%a_ptr%set(self%s%iptcl, 'state', real(state))
        call self%s%a_ptr%set(self%s%iptcl, 'frac', frac )
        call self%s%a_ptr%set(self%s%iptcl, 'corr', wcorr )
        call self%s%a_ptr%set(self%s%iptcl, 'specscore', self%s%specscore)
        call self%s%a_ptr%set(self%s%iptcl, 'ow',    o_peaks(self%s%iptcl)%get(best_loc(1),'ow')   )
        call self%s%a_ptr%set(self%s%iptcl, 'proj',  o_peaks(self%s%iptcl)%get(best_loc(1),'proj') )
        call self%s%a_ptr%set(self%s%iptcl, 'sdev',  ang_sdev )
        call self%s%a_ptr%set(self%s%iptcl, 'npeaks', real(self%s%npeaks_eff) )
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine oris_assign_greedy_multi

    subroutine kill_greedy_multi( self )
        class(strategy3D_greedy_multi),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedy_multi

end module simple_strategy3D_greedy_multi
