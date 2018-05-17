! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster
use simple_strategy3D_alloc ! use all in there
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
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
        real    :: shvec(2), corr, mi_state, frac, mi_inpl, mi_proj, w, bfac
        logical :: hetsym
        self%s%prev_state = self%s%a_ptr%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            hetsym            = associated(self%spec%symmat)
            self%s%prev_roind = self%s%pftcc_ptr%get_roind(360.-self%s%a_ptr%e3get(self%s%iptcl))
            self%s%prev_corr  = self%s%a_ptr%get(self%s%iptcl, 'corr')
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
            w = 1.
            ! B-factor memoization
            if( self%s%pftcc_ptr%objfun_is_ccres() )then
                bfac = self%s%pftcc_ptr%fit_bfac(self%s%prev_ref, self%s%iptcl, self%s%prev_roind, [0.,0.])
                call self%s%pftcc_ptr%memoize_bfac(self%s%iptcl, bfac)
                call self%s%a_ptr%set(self%s%iptcl,'bfac', bfac)
            endif
            ! evaluate all correlations
            corrs = -1.
            do state = 1, self%s%nstates
                if( .not.s3D%state_exists(state) ) cycle
                if( hetsym )then
                    ! greedy in symmetric unit
                    do isym = 1, self%s%nsym
                        iproj = self%spec%symmat(s3D%prev_proj(self%s%iptcl_map), isym)
                        iref  = (state-1) * self%s%nprojs + iproj
                        call self%s%pftcc_ptr%gencorrs(iref, self%s%iptcl, corrs_inpl)
                        corrs_sym(isym) = corrs_inpl(self%s%prev_roind)
                    enddo
                    loc              = maxloc(corrs_sym)
                    isym             = loc(1)
                    corrs(state)     = corrs_sym(isym)
                    sym_projs(state) = self%spec%symmat(s3D%prev_proj(self%s%iptcl_map), isym)
                else
                    iref = (state-1) * self%s%nprojs +s3D%prev_proj(self%s%iptcl_map)
                    call self%s%pftcc_ptr%gencorrs(iref, self%s%iptcl, corrs_inpl)
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
            if( self%s%prev_corr < self%spec%corr_thresh )then
                ! state randomization
                state = irnd_uni(self%s%nstates)
                do while(state == self%s%prev_state .or. .not.s3D%state_exists(state))
                    state = irnd_uni(self%s%nstates)
                enddo
                corr              = corrs(state)
                self%s%nrefs_eval = 1
            else
                ! SHC state optimization
                self%s%prev_corr = corrs(self%s%prev_state)
                state            = shcloc(self%s%nstates, corrs, self%s%prev_corr)
                corr             = corrs(state)
                self%s%nrefs_eval = count(corrs <= self%s%prev_corr)
                if( self%s%prev_state .eq. state ) mi_state = 1.
                if( self%spec%do_extr )then
                    ! extremal optimization
                    if( hetsym )then
                        ! greedy in symmetric unit
                        if( sym_projs(state) .ne. s3D%prev_proj(self%s%iptcl_map) )then
                            mi_proj = 0.
                            iref    = (state-1)*self%s%nprojs + sym_projs(state)
                            call self%s%a_ptr%e1set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 1))
                            call self%s%a_ptr%e2set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 2))
                        endif
                    endif
                else
                    ! greedy moves after extremal optimization complete
                    if( mi_state > 0.5 )then
                        if( hetsym )then
                            ! greedy in symmetric units
                            iproj   = sym_projs(state)
                            if( iproj .ne. s3D%prev_proj(self%s%iptcl_map) )then
                                mi_proj = 0.
                                iref    = (state-1) * self%s%nprojs + iproj
                                call self%s%a_ptr%e1set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 1))
                                call self%s%a_ptr%e2set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 2))
                            endif
                        else
                            ! greedy in plane
                            iref = (state-1)*self%s%nprojs + s3D%prev_proj(self%s%iptcl_map)
                            s3D%proj_space_corrs(self%s%iptcl_map,iref)        = corr
                            s3D%proj_space_inds(self%s%iptcl_map,self%s%nrefs) = iref
                            call self%s%inpl_srch
                            if( s3D%proj_space_corrs(self%s%iptcl_map,iref) > corr )then
                                corr  = s3D%proj_space_corrs(self%s%iptcl_map,iref)
                                shvec = self%s%a_ptr%get_2Dshift(self%s%iptcl) + s3D%proj_space_shift(self%s%iptcl_map,iref,:)
                                call self%s%a_ptr%set_shift(self%s%iptcl, shvec)
                                call self%s%a_ptr%e3set(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map, iref, 3))
                            endif
                            if( self%s%prev_roind .ne. self%s%pftcc_ptr%get_roind(360.-s3D%proj_space_euls(self%s%iptcl_map, iref, 3)) )then
                                mi_inpl = 0.
                            endif
                        endif
                    endif
                    call self%s%a_ptr%set(self%s%iptcl,'proj', real(s3D%proj_space_proj(self%s%iptcl_map,iref)))
                endif
            endif
            ! updates peaks and orientation orientation
            frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
            call self%s%a_ptr%set(self%s%iptcl,'frac',     frac)
            call self%s%a_ptr%set(self%s%iptcl,'state',    real(state))
            call self%s%a_ptr%set(self%s%iptcl,'corr',     corr)
            call self%s%a_ptr%set(self%s%iptcl,'mi_proj',  mi_proj)
            call self%s%a_ptr%set(self%s%iptcl,'mi_inpl',  mi_inpl)
            call self%s%a_ptr%set(self%s%iptcl,'mi_state', mi_state)
            call self%s%a_ptr%set(self%s%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
            call self%s%a_ptr%set(self%s%iptcl,'w',        w)
            call s3D%o_peaks(self%s%iptcl)%set(1,'state', real(state))
            call s3D%o_peaks(self%s%iptcl)%set(1,'corr',  corr)
            call s3D%o_peaks(self%s%iptcl)%set(1,'w',     w)
            call s3D%o_peaks(self%s%iptcl)%set(1,'ow',    1.)
            call s3D%o_peaks(self%s%iptcl)%set_euler(1, self%s%a_ptr%get_euler(self%s%iptcl))
            call s3D%o_peaks(self%s%iptcl)%set_shift(1, self%s%a_ptr%get_2Dshift(self%s%iptcl))
        else
            call self%s%a_ptr%reject(self%s%iptcl)
        endif
        DebugPrint   '>>> STRATEGY3D_CLUSTER :: FINISHED SRCH_CLUSTER3D'
    end subroutine srch_cluster3D

    !>  placeholder
    subroutine oris_assign_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        ! nothing to do
    end subroutine oris_assign_cluster3D

    subroutine kill_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D

end module simple_strategy3D_cluster
