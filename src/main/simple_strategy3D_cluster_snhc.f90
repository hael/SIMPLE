! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster_snhc
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy3D_cluster_snhc
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cluster_snhc
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_cluster3D_snhc
    procedure :: srch        => srch_cluster3D_snhc
    procedure :: oris_assign => oris_assign_cluster3D_snhc
    procedure :: kill        => kill_cluster3D_snhc
end type strategy3D_cluster_snhc

contains

    subroutine new_cluster3D_snhc( self, spec, npeaks )
        class(strategy3D_cluster_snhc), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        integer,                   intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_cluster3D_snhc

    ! subroutine srch_cluster3D_snhc( self )
    !     use simple_rnd, only: shcloc, irnd_uni
    !     class(strategy3D_cluster_snhc),   intent(inout) :: self
    !     integer :: iref, state, nstates_bound, nstates
    !     real    :: corrs(self%s%nstates), corrs_inpl(self%s%nrots), shvec(2)
    !     real    :: corr, mi_state, frac, mi_inpl, bfac
    !     logical :: avail(self%s%nstates)
    !     self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
    !     if( self%s%prev_state > 0 )then
    !         ! prep4srch
    !         self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
    !         self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
    !         ! B-factor memoization
    !         if( params_glob%l_bfac_static )then
    !             bfac = params_glob%bfac_static
    !         else
    !             bfac = pftcc_glob%fit_bfac(self%s%prev_ref, self%s%iptcl, self%s%prev_roind, [0.,0.])
    !         endif
    !         if( params_glob%cc_objfun == OBJFUN_RES )call pftcc_glob%memoize_bfac(self%s%iptcl, bfac)
    !         ! specscore
    !         self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
    !         ! neighbourhood
    !         avail = s3D%state_exists
    !         if( self%spec%extr_score_thresh < 0.02 )then
    !             ! full search space
    !             nstates_bound = self%s%nstates
    !         else
    !             ! stochastic search space
    !             if( ran3() > self%spec%extr_score_thresh )then
    !                 nstates_bound = self%s%nstates
    !             else
    !                 ! nstates_bound = irnd_uni(self%s%nstates-1)
    !                 nstates_bound = 1
    !             endif
    !         endif
    !         if( nstates_bound == 1 )then
    !             ! pick a random one that is not the previous
    !             state = irnd_uni(self%s%nstates)
    !             do while( (state==self%s%prev_state) .and. .not.avail(state) )
    !                 state = irnd_uni(self%s%nstates)
    !             enddo
    !             self%s%nrefs_eval = 1
    !         else
    !             if( nstates_bound < self%s%nstates )then
    !                 ! reject previous state
    !                 avail(self%s%prev_state) = .false.
    !                 ! randomly rejects states
    !                 do while( count(avail) > nstates_bound )
    !                     state = irnd_uni(self%s%nstates)
    !                     if( avail(state) ) avail(state) = .false.
    !                 enddo
    !             endif
    !             ! evaluate correlations
    !             do state=1,self%s%nstates
    !                 if( avail(state) .or. state==self%s%prev_state )then
    !                     iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
    !                     call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs_inpl)
    !                     ! corrs(state) = corrs_inpl(self%s%prev_roind)
    !                     corrs(state) = maxval(corrs_inpl)
    !                 endif
    !             enddo
    !             self%s%prev_corr = corrs(self%s%prev_state)
    !             where( .not.avail )corrs = -1.
    !             ! state shc move
    !             state = shcloc(self%s%nstates, corrs, self%s%prev_corr)
    !             if( state == 0 )then
    !                 ! no improvement found, pick best.
    !                 state = maxloc(corrs, dim=1)
    !                 self%s%nrefs_eval = nstates_bound
    !             else
    !                 self%s%nrefs_eval = count((corrs<=self%s%prev_corr).and.avail)
    !             endif
    !         endif
    !         ! parameters & correlation update
    !         iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
    !         corr = pftcc_glob%gencorr_for_rot_8(iref, self%s%iptcl, [0.d0,0.d0], self%s%prev_roind)
    !         if( self%s%nstates == nstates_bound )then
    !             ! greedy in-plane after stochastic neighbourhood phase
    !             s3D%proj_space_corrs(self%s%iptcl_map,iref,1)         = corr
    !             s3D%proj_space_refinds(self%s%iptcl_map,self%s%nrefs) = iref
    !             call self%s%inpl_srch
    !             if( s3D%proj_space_corrs(self%s%iptcl_map,iref,1) > corr )then
    !                 corr  = s3D%proj_space_corrs(self%s%iptcl_map,iref,1)
    !                 shvec = build_glob%spproj_field%get_2Dshift(self%s%iptcl) + s3D%proj_space_shift(self%s%iptcl_map,iref,1,:)
    !                 call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
    !                 call build_glob%spproj_field%e3set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map,iref,1,3))
    !             endif
    !         endif
    !         ! updates peak and orientation
    !         mi_inpl  = 0.
    !         if( self%s%prev_roind==pftcc_glob%get_roind(360.-s3D%proj_space_euls(self%s%iptcl_map,iref,1,3)) )mi_inpl = 0.
    !         mi_state = 0.
    !         if( state == self%s%prev_state )mi_state = 1.
    !         frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
    !         call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
    !         call build_glob%spproj_field%set(self%s%iptcl,'state',    real(state))
    !         call build_glob%spproj_field%set(self%s%iptcl,'corr',     corr)
    !         call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  1.)
    !         call build_glob%spproj_field%set(self%s%iptcl,'mi_inpl',  mi_inpl)
    !         call build_glob%spproj_field%set(self%s%iptcl,'mi_state', mi_state)
    !         call build_glob%spproj_field%set(self%s%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
    !         call build_glob%spproj_field%set(self%s%iptcl,'w',        1.)
    !         call build_glob%spproj_field%set(self%s%iptcl,'bfac',     bfac)
    !         call build_glob%spproj_field%set(self%s%iptcl,'specscore',self%s%specscore)
    !         call s3D%o_peaks(self%s%iptcl)%set(1,'state', real(state))
    !         call s3D%o_peaks(self%s%iptcl)%set(1,'corr',  corr)
    !         call s3D%o_peaks(self%s%iptcl)%set(1,'w',     1.)
    !         call s3D%o_peaks(self%s%iptcl)%set(1,'ow',    1.)
    !         call s3D%o_peaks(self%s%iptcl)%set_euler(1, build_glob%spproj_field%get_euler(self%s%iptcl))
    !         call s3D%o_peaks(self%s%iptcl)%set_shift(1, build_glob%spproj_field%get_2Dshift(self%s%iptcl))
    !     else
    !         call build_glob%spproj_field%reject(self%s%iptcl)
    !     endif
    !     DebugPrint   '>>> STRATEGY3D_CLUSTER :: FINISHED SRCH_CLUSTER3D'
    ! end subroutine srch_cluster3D_snhc

    subroutine srch_cluster3D_snhc( self )
        use simple_rnd, only: shcloc, irnd_uni
        class(strategy3D_cluster_snhc),   intent(inout) :: self
        integer :: iproj, iref, state, ran_state
        real    :: corrs(self%s%nstates), corrs_inpl(self%s%nrots)
        real    :: shvec(2), corr, mi_state, frac, mi_inpl, mi_proj, bfac
        logical :: do_snhc, do_greedy_inpl, avail(self%s%nstates)
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! local prep
            self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
            ! specscore & B-factor
            if( params_glob%l_bfac_static )then
                bfac = params_glob%bfac_static
            else
                bfac = pftcc_glob%fit_bfac(self%s%prev_ref, self%s%iptcl, self%s%prev_roind, [0.,0.])
            endif
            if( params_glob%cc_objfun == OBJFUN_RES )call pftcc_glob%memoize_bfac(self%s%iptcl, bfac)
            self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
            ! search decisions
            self%s%prev_corr = build_glob%spproj_field%get(self%s%iptcl, 'corr')
            do_snhc          = self%spec%do_extr .and. self%s%prev_corr < self%spec%extr_score_thresh
            !do_greedy_inpl   = .not.self%spec%do_extr
            !do_snhc        = ran3() < self%spec%extr_score_thresh
            !write(*,*)self%s%iptcl, self%spec%extr_score_thresh, do_snhc
            do_greedy_inpl = .false.
            ! do_greedy_inpl = self%spec%extr_score_thresh < 0.02
            ! evaluate all correlations
            corrs = -1.
            do state = 1, self%s%nstates
                if( .not. s3D%state_exists(state) ) cycle
                iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
                call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs_inpl)
                ! to somewhat take account for alignment errors due to state mixing
                ! replaced corrs(state) = corrs_inpl(self%s%prev_roind) with:
                corrs(state) = maxval(corrs_inpl)
            enddo
            self%s%prev_corr = corrs(self%s%prev_state)
            ! make moves
            mi_state = 0.
            mi_inpl  = 1.
            mi_proj  = 1.
            if( do_snhc )then
                ! state randomization
                avail = s3D%state_exists
                avail(self%s%prev_state) = .false.
                ran_state = irnd_uni(self%s%nstates)
                do while(.not.avail(ran_state))
                    ran_state = irnd_uni(self%s%nstates)
                enddo
                state = ran_state
                corr  = corrs(state)
                if( corr > self%s%prev_corr )then
                    ! accept
                    self%s%nrefs_eval = 1
                else
                    if( count(avail) == 0 )then
                        self%s%nrefs_eval = 1
                    else
                        avail(state) = .false.
                        ran_state    = irnd_uni(self%s%nstates)
                        do while(.not.avail(ran_state))
                            ran_state = irnd_uni(self%s%nstates)
                        enddo
                        if(corrs(ran_state) > corr)then
                            state = ran_state
                            corr  = corrs(state)
                        endif
                        self%s%nrefs_eval = 2
                    endif
                endif
                ! state randomization
                ! state = irnd_uni(self%s%nstates)
                ! do while(state == self%s%prev_state .or. .not.s3D%state_exists(state))
                !     state = irnd_uni(self%s%nstates)
                ! enddo
                ! corr = corrs(state)
                ! self%s%nrefs_eval = 1
            else
                ! SHC state optimization
                if( .not.self%spec%do_extr .or. ran3() > 0.5 )then
                    state = maxloc(corrs, dim=1)
                    self%s%nrefs_eval = self%s%nstates
                else
                    state = shcloc(self%s%nstates, corrs, self%s%prev_corr)
                    self%s%nrefs_eval = count(corrs <= self%s%prev_corr)
                endif
                corr = corrs(state)
                if( self%s%prev_state .eq. state ) mi_state = 1.
                if( do_greedy_inpl )then
                    ! greedy moves after extremal optimization complete
                    if( mi_state > 0.5 )then
                        ! greedy in plane, need to recalculate bound for in-plane search
                        iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
                        corr = pftcc_glob%gencorr_for_rot_8(iref, self%s%iptcl, [0.d0,0.d0], self%s%prev_roind)
                        s3D%proj_space_corrs(self%s%iptcl_map,iref,1)         = corr ! inpl = 1
                        s3D%proj_space_refinds(self%s%iptcl_map,self%s%nrefs) = iref
                        call self%s%inpl_srch
                        if( s3D%proj_space_corrs(self%s%iptcl_map,iref,1) > corr )then
                            corr  = s3D%proj_space_corrs(self%s%iptcl_map,iref,1)
                            shvec = build_glob%spproj_field%get_2Dshift(self%s%iptcl) + s3D%proj_space_shift(self%s%iptcl_map,iref,1,:) ! inpl = 1
                            call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
                            call build_glob%spproj_field%e3set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref,1,3))           ! inpl = 1
                        endif
                        if( self%s%prev_roind .ne. pftcc_glob%get_roind(360.-s3D%proj_space_euls(self%s%iptcl_map, iref,1,3)) ) mi_inpl = 0.
                    endif
                    call build_glob%spproj_field%set(self%s%iptcl,'proj', real(s3D%proj_space_proj(self%s%iptcl_map,iref)))
                endif
            endif
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
            call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
            call build_glob%spproj_field%set(self%s%iptcl,'state',    real(state))
            call build_glob%spproj_field%set(self%s%iptcl,'corr',     corr)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  mi_proj)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_inpl',  mi_inpl)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_state', mi_state)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
            call build_glob%spproj_field%set(self%s%iptcl,'w',        1.)
            call build_glob%spproj_field%set(self%s%iptcl,'bfac',     bfac)
            call build_glob%spproj_field%set(self%s%iptcl,'specscore',self%s%specscore)
            call s3D%o_peaks(self%s%iptcl)%set(1,'state', real(state))
            call s3D%o_peaks(self%s%iptcl)%set(1,'corr',  corr)
            call s3D%o_peaks(self%s%iptcl)%set(1,'w',     1.)
            call s3D%o_peaks(self%s%iptcl)%set(1,'ow',    1.)
            call s3D%o_peaks(self%s%iptcl)%set_euler(1, build_glob%spproj_field%get_euler(self%s%iptcl))
            call s3D%o_peaks(self%s%iptcl)%set_shift(1, build_glob%spproj_field%get_2Dshift(self%s%iptcl))
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint   '>>> STRATEGY3D_CLUSTER :: FINISHED SRCH_CLUSTER3D'
    end subroutine srch_cluster3D_snhc

    !>  placeholder
    subroutine oris_assign_cluster3D_snhc( self )
        class(strategy3D_cluster_snhc), intent(inout) :: self
        ! nothing to do
    end subroutine oris_assign_cluster3D_snhc

    subroutine kill_cluster3D_snhc( self )
        class(strategy3D_cluster_snhc), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D_snhc

end module simple_strategy3D_cluster_snhc
