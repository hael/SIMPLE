! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: strategy3D_cluster
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cluster
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure :: new         => new_cluster3D
    procedure :: srch        => srch_cluster3D
    procedure :: oris_assign => oris_assign_cluster3D
    procedure :: kill        => kill_cluster3D
end type strategy3D_cluster

contains

    subroutine new_cluster3D( self, spec, npeaks )
        class(strategy3D_cluster), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        integer,                   intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_cluster3D

    subroutine srch_cluster3D( self )
        use simple_rnd, only: shcloc, irnd_uni
        class(strategy3D_cluster),   intent(inout) :: self
        integer :: sym_projs(self%s%nstates), loc(1), iproj, iref, isym, state
        real    :: corrs(self%s%nstates), corrs_sym(self%s%nsym), corrs_inpl(self%s%nrots)
        real    :: shvec(2), corr, mi_state, frac, mi_inpl, mi_proj, bfac, score4extr
        logical :: hetsym
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! local version of prep4srch
            hetsym            = associated(self%spec%symmat)
            self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
            self%s%prev_corr  = build_glob%spproj_field%get(self%s%iptcl, 'corr')
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
            ! extremal optimization score
            score4extr = self%s%prev_corr
            ! B-factor memoization
            if( params_glob%l_bfac_static )then
                bfac = params_glob%bfac_static
            else
                bfac = pftcc_glob%fit_bfac(self%s%prev_ref, self%s%iptcl, self%s%prev_roind, [0.,0.])
            endif
            if( pftcc_glob%get_objfun() == 2 ) call pftcc_glob%memoize_bfac(self%s%iptcl, bfac)
            ! specscore
            self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
            ! evaluate all correlations
            corrs = -1.
            do state = 1, self%s%nstates
                if( .not. s3D%state_exists(state) ) cycle
                if( hetsym )then
                    ! greedy in symmetric unit
                    do isym = 1, self%s%nsym
                        iproj = self%spec%symmat(s3D%prev_proj(self%s%iptcl_map), isym)
                        iref  = (state-1) * self%s%nprojs + iproj
                        ! dual resolution limit scheme:
                        ! state assignement employs correlations up to kstop_grid
                        ! image alignment to state employs correlations up to kfromto(2)
                        call pftcc_glob%gencorrs(iref, self%s%iptcl, self%s%kstop_grid, corrs_inpl)
                        corrs_sym(isym) = corrs_inpl(self%s%prev_roind)
                    enddo
                    loc              = maxloc(corrs_sym)
                    isym             = loc(1)
                    corrs(state)     = corrs_sym(isym)
                    sym_projs(state) = self%spec%symmat(s3D%prev_proj(self%s%iptcl_map), isym)
                else
                    iref = (state-1) * self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
                    ! dual resolution limit scheme:
                    ! state assignement employs correlations up to kstop_grid
                    ! image alignment to state employs correlations up to kfromto(2)
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, self%s%kstop_grid, corrs_inpl)
                    ! replaced:
                    ! corrs(state) = corrs_inpl(self%s%prev_roind)
                    ! with:
                    corrs(state) = maxval(corrs_inpl)
                    ! to somewhat take account for alignment errors due to state mixing
                endif
            enddo
            ! make moves
            mi_state = 0.
            mi_inpl  = 1.
            mi_proj  = 1.
            if( score4extr < self%spec%extr_score_thresh )then
                ! state randomization
                state = irnd_uni(self%s%nstates)
                do while(state == self%s%prev_state .or. .not.s3D%state_exists(state))
                    state = irnd_uni(self%s%nstates)
                enddo
                corr = corrs(state)
                self%s%nrefs_eval = 1
            else
                ! SHC state optimization
                self%s%prev_corr  = corrs(self%s%prev_state)
                state             = shcloc(self%s%nstates, corrs, self%s%prev_corr)
                corr              = corrs(state)
                self%s%nrefs_eval = count(corrs <= self%s%prev_corr)
                if( self%s%prev_state .eq. state ) mi_state = 1.
                if( self%spec%do_extr )then
                    ! extremal optimization
                    if( hetsym )then
                        ! greedy in symmetric unit
                        if( sym_projs(state) .ne. s3D%prev_proj(self%s%iptcl_map) )then
                            mi_proj = 0.
                            iref    = (state-1)*self%s%nprojs + sym_projs(state)
                            call build_glob%spproj_field%e1set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map,iref,1,1)) ! inpl = 1
                            call build_glob%spproj_field%e2set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map,iref,1,2)) ! inpl = 1
                        endif
                    endif
                else
                    ! greedy moves after extremal optimization complete
                    if( mi_state > 0.5 )then
                        if( hetsym )then
                            ! greedy in symmetric units
                            iproj = sym_projs(state)
                            if( iproj .ne. s3D%prev_proj(self%s%iptcl_map) )then
                                mi_proj = 0.
                                iref    = (state-1) * self%s%nprojs + iproj
                                call build_glob%spproj_field%e1set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 1, 1)) ! inpl = 1
                                call build_glob%spproj_field%e2set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 1, 2)) ! inpl = 1
                            endif
                        else
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
    end subroutine srch_cluster3D

    subroutine oris_assign_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        ! nothing to do
    end subroutine oris_assign_cluster3D

    subroutine kill_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D

end module simple_strategy3D_cluster
