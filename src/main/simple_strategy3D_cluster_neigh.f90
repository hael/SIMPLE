! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster_neigh
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy3D_cluster_neigh
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cluster_neigh
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_cluster3D_neigh
    procedure :: srch        => srch_cluster3D_neigh
    procedure :: oris_assign => oris_assign_cluster3D_neigh
    procedure :: kill        => kill_cluster3D_neigh
end type strategy3D_cluster_neigh

contains

    subroutine new_cluster3D_neigh( self, spec, npeaks )
        class(strategy3D_cluster_neigh), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        integer,                   intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_cluster3D_neigh

    subroutine srch_cluster3D_neigh( self )
        use simple_rnd, only: shcloc, irnd_uni
        class(strategy3D_cluster_neigh),   intent(inout) :: self
        integer :: iref, state, nstates_bound, nstates
        real    :: corrs(self%s%nstates), corrs_inpl(self%s%nrots), shvec(2)
        real    :: corr, mi_state, frac, mi_inpl, bfac
        logical :: avail(self%s%nstates)
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! prep4srch
            self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
            ! B-factor memoization
            if( params_glob%l_bfac_static )then
                bfac = params_glob%bfac_static
            else
                bfac = pftcc_glob%fit_bfac(self%s%prev_ref, self%s%iptcl, self%s%prev_roind, [0.,0.])
            endif
            if( pftcc_glob%get_objfun() == 2 )call pftcc_glob%memoize_bfac(self%s%iptcl, bfac)
            ! specscore
            self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
            ! stochastic states neighbourhood
            avail = s3D%state_exists
            if( self%spec%extr_score_thresh < 0.05 )then
                ! full search space
                nstates_bound = self%s%nstates
            else
                ! stochastic search space bound
                nstates_bound = min(self%s%nstates-1, nint(real(count(s3D%state_exists))*(1.-self%spec%extr_score_thresh)) )
                nstates_bound = max(2, nstates_bound)
                ! removes previous state
                ! avail(self%s%prev_state) = .false.
                ! randomly rejects states
                do while( count(avail) > nstates_bound )
                    state = irnd_uni(self%s%nstates)
                    if( avail(state) ) avail(state) = .false.
                enddo
                write(*,*)self%s%iptcl,self%spec%extr_score_thresh,count(avail),nstates_bound,self%s%prev_state,avail
            endif
            ! evaluate correlations
            corrs = -1.
            do state=1,self%s%nstates
                if( .not.avail(state) ) cycle
                ! dual resolution limit scheme:
                ! state assignement employs correlations up to kstop_grid
                iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
                call pftcc_glob%gencorrs(iref, self%s%iptcl, self%s%kstop_grid, corrs_inpl)
                corrs(state) = corrs_inpl(self%s%prev_roind)
            enddo
            self%s%prev_corr = corrs(self%s%prev_state)
            ! state shc move
            state = shcloc(self%s%nstates, corrs, self%s%prev_corr)
            if( state == 0 )then
                ! no improvement found, pick best.
                state = maxloc(corrs, dim=1)
                self%s%nrefs_eval = self%s%nstates
            else
                self%s%nrefs_eval = count(corrs <= self%s%prev_corr)
            endif
            ! parameters & correlation update
            iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
            corr = pftcc_glob%gencorr_for_rot_8(iref, self%s%iptcl, [0.d0,0.d0], self%s%prev_roind)
            if( self%s%nstates == nstates_bound )then
                ! greedy in-plane after stochastic neighbourhood phase
                s3D%proj_space_corrs(self%s%iptcl_map,iref,1)         = corr
                s3D%proj_space_refinds(self%s%iptcl_map,self%s%nrefs) = iref
                call self%s%inpl_srch
                if( s3D%proj_space_corrs(self%s%iptcl_map,iref,1) > corr )then
                    corr  = s3D%proj_space_corrs(self%s%iptcl_map,iref,1)
                    shvec = build_glob%spproj_field%get_2Dshift(self%s%iptcl) + s3D%proj_space_shift(self%s%iptcl_map,iref,1,:)
                    call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
                    call build_glob%spproj_field%e3set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map,iref,1,3))
                endif
            endif
            ! updates peak and orientation
            mi_inpl  = 0.
            if( self%s%prev_roind==pftcc_glob%get_roind(360.-s3D%proj_space_euls(self%s%iptcl_map,iref,1,3)) )mi_inpl = 0.
            mi_state = 0.
            if( state == self%s%prev_state )mi_state = 1.
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
            call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
            call build_glob%spproj_field%set(self%s%iptcl,'state',    real(state))
            call build_glob%spproj_field%set(self%s%iptcl,'corr',     corr)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  1.)
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
    end subroutine srch_cluster3D_neigh

    !>  placeholder
    subroutine oris_assign_cluster3D_neigh( self )
        class(strategy3D_cluster_neigh), intent(inout) :: self
        ! nothing to do
    end subroutine oris_assign_cluster3D_neigh

    subroutine kill_cluster3D_neigh( self )
        class(strategy3D_cluster_neigh), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D_neigh

end module simple_strategy3D_cluster_neigh
