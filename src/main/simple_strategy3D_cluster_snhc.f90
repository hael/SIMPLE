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

    subroutine srch_cluster3D_snhc( self, ithr )
        use simple_rnd, only: shcloc, irnd_uni
        class(strategy3D_cluster_snhc),   intent(inout) :: self
        integer,                          intent(in)    :: ithr
        integer :: iref, state, ran_state
        real    :: corrs(self%s%nstates)
        real    :: shvec(2), corr, mi_state, frac, mi_inpl, bfac
        logical :: do_snhc, do_greedy_inpl, avail(self%s%nstates)
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! local prep
            call prep_strategy3D_thread(self%s%ithr)
            self%s%prev_proj  = build_glob%eulspace%find_closest_proj(build_glob%spproj_field%get_ori(self%s%iptcl))
            self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + self%s%prev_proj
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
            do_greedy_inpl   = .not.self%spec%do_extr
            ! evaluate all correlations
            corrs = -1.
            do state = 1, self%s%nstates
                if( .not. s3D%state_exists(state) ) cycle
                iref = (state-1) * self%s%nprojs + self%s%prev_proj
                ! to somewhat take account alignment errors due to state mixing
                ! call pftcc_glob%gencorrs(iref, self%s%iptcl, corrs_inpl)
                ! corrs(state) = maxval(corrs_inpl)
                ! to not take account alignment errors due to state mixing
                corrs(state) = real(pftcc_glob%gencorr_for_rot_8(iref, self%s%iptcl, [0.d0,0.d0], self%s%prev_roind))
            enddo
            self%s%prev_corr = corrs(self%s%prev_state)
            ! make moves
            mi_state = 0.
            mi_inpl  = 1.
            if( do_snhc )then
                ! state randomization, best of two
                avail = s3D%state_exists
                !avail(self%s%prev_state) = .false. ! without replacement
                ran_state = irnd_uni(self%s%nstates)
                do while(.not.avail(ran_state))
                    ran_state = irnd_uni(self%s%nstates)
                enddo
                state        = ran_state
                corr         = corrs(state)
                avail(state) = .false.
                if( count(avail) == 0 )then
                    self%s%nrefs_eval = 1
                else
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
                if( state.eq.self%s%prev_state )mi_state = 1.
            else
                ! shc in state
                ! state = shcloc(self%s%nstates, corrs, self%s%prev_corr)
                ! self%s%nrefs_eval = count(corrs <= self%s%prev_corr)
                ! greedy in state
                state = maxloc(corrs, dim=1)
                self%s%nrefs_eval = self%s%nstates
                corr = corrs(state)
                if( state.eq.self%s%prev_state )mi_state = 1.
                if( do_greedy_inpl )then
                    ! greedy in-plane moves after extremal optimization complete
                    if( mi_state > 0.5 )then
                        iref = (state-1) * self%s%nprojs + self%s%prev_proj
                        corr = corrs(state)
                        s3D%proj_space_corrs(self%s%ithr,iref,1)         = corr
                        s3D%proj_space_refinds(self%s%ithr,self%s%nrefs) = iref
                        call self%s%inpl_srch
                        if( s3D%proj_space_corrs(self%s%ithr,iref,1) > corr )then
                            corr  = s3D%proj_space_corrs(self%s%ithr,iref,1)
                            shvec = build_glob%spproj_field%get_2Dshift(self%s%iptcl) + s3D%proj_space_shift(self%s%ithr,iref,1,:) ! inpl = 1
                            call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
                            call build_glob%spproj_field%e3set(self%s%iptcl, s3D%proj_space_euls(self%s%ithr, iref,1,3))           ! inpl = 1
                        endif
                        if( self%s%prev_roind .ne. pftcc_glob%get_roind(360.-s3D%proj_space_euls(self%s%ithr, iref,1,3)) ) mi_inpl = 0.
                    endif
                    call build_glob%spproj_field%set(self%s%iptcl,'proj', real(self%s%prev_proj))
                endif
            endif
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
            call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
            call build_glob%spproj_field%set(self%s%iptcl,'state',    real(state))
            call build_glob%spproj_field%set(self%s%iptcl,'corr',     corr)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  1.)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_inpl',  mi_inpl)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_state', mi_state)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
            call build_glob%spproj_field%set(self%s%iptcl,'w',        1.)
            call build_glob%spproj_field%set(self%s%iptcl,'npeaks',   1.)
            call build_glob%spproj_field%set(self%s%iptcl,'bfac',     bfac)
            call build_glob%spproj_field%set(self%s%iptcl,'specscore',self%s%specscore)
            call s3D%o_peaks(self%s%iptcl)%set(1,'state', real(state))
            call s3D%o_peaks(self%s%iptcl)%set(1,'corr',  corr)
            call s3D%o_peaks(self%s%iptcl)%set(1,'w',     1.)
            call s3D%o_peaks(self%s%iptcl)%set(1,'ow',    1.)
            call s3D%o_peaks(self%s%iptcl)%set(1,'npeaks',1.)
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
