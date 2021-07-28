! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
use simple_rnd,              only: shcloc, irnd_uni
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

    subroutine new_cluster3D( self, spec )
        class(strategy3D_cluster), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        call self%s%new( spec )
        self%spec = spec
    end subroutine new_cluster3D

    !>  \brief search driver
    subroutine srch_cluster3D( self, ithr )
        use simple_ori, only: ori
        class(strategy3D_cluster), intent(inout) :: self
        integer,                   intent(in)    :: ithr
        type(ori) :: o
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! prep
            self%s%ithr = ithr
            call prep_strategy3D_thread(self%s%ithr)
            self%s%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%s%iptcl))
            self%s%prev_corr  = build_glob%spproj_field%get(self%s%iptcl, 'corr') ! score for EO
            call build_glob%spproj_field%get_ori(self%s%iptcl, o)
            self%s%prev_proj  = build_glob%eulspace%find_closest_proj(o)
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + self%s%prev_proj
            self%s%prev_shvec = build_glob%spproj_field%get_2Dshift(self%s%iptcl)
            ! specscore
            self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
            call build_glob%spproj_field%set(self%s%iptcl,'specscore',self%s%specscore)
            ! fork
            if( associated(self%spec%symmat) )then
                call symsrch_cluster3D(self%s, self%spec)
            else
                call statesrch_cluster3D(self%s, self%spec)
            endif
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        call o%kill
    end subroutine srch_cluster3D

    subroutine statesrch_cluster3D(s, spec)
        type(strategy3D_srch), intent(inout) :: s
        type(strategy3D_spec), intent(in)    :: spec
        integer :: ind,iproj,iref,state,inpl,i ! neigh_projs(s%nstates), ineigh
        real    :: corrs(s%nstates),corrs_inpl(s%nrots), corr,mi_state,mi_proj
        corrs       = -1.
        ! neigh_projs = 0
        ! evaluate all correlations
        do state=1,s%nstates
            if( .not. s3D%state_exists(state) ) cycle
            ! if( s%neigh )then
            !     ! greedy neighbourhood search
            !     do ineigh=1,s%nnn
            !         iproj = build_glob%nnmat(s%prev_proj,ineigh)
            !         iref  = (state-1)*s%nprojs + iproj
            !         call pftcc_glob%gencorrs(iref, s%iptcl, corrs_inpl)
            !         inpl = maxloc(corrs_inpl,dim=1)
            !         if( corrs_inpl(inpl) > corrs(state) )then
            !             corrs(state)       = corrs_inpl(inpl)
            !             neigh_projs(state) = iproj
            !         endif
            !     enddo
            ! else
                iref = (state-1)*s%nprojs + s%prev_proj
                call pftcc_glob%gencorrs(iref, s%iptcl, corrs_inpl)
                corrs(state) = maxval(corrs_inpl)
            ! endif
        enddo
        ! make moves
        if(s%prev_corr < spec%extr_score_thresh)then
            ! state randomization
            state = irnd_uni(s%nstates)
            if( s%nstates > 2 )then
                do while( state==s%prev_state .or. .not.s3D%state_exists(state) )
                    state = irnd_uni(s%nstates)
                enddo
            else
                do while( .not.s3D%state_exists(state) )
                    state = irnd_uni(s%nstates)
                enddo
            endif
            corr         = corrs(state)
            iproj        = s%prev_proj
            s%nrefs_eval = 1
        else
            ! SHC state optimization
            s%prev_corr  = corrs(s%prev_state)
            state        = shcloc(s%nstates, corrs, s%prev_corr)
            corr         = corrs(state)
            s%nrefs_eval = count(corrs <= s%prev_corr)
            ! if( s%neigh )then
            !     iproj = neigh_projs(state)
            !     iref  = (state-1)*s%nprojs + iproj
            !     call build_glob%spproj_field%e1set(s%iptcl,s3D%proj_space_euls(s%ithr,iref,1))
            !     call build_glob%spproj_field%e2set(s%iptcl,s3D%proj_space_euls(s%ithr,iref,2))
            ! else
                iproj = s%prev_proj
            ! endif
            if( state == s%prev_state ) call greedy_inplsrch(s, corr, state, iproj)
        endif
        ! reporting & convergence
        mi_state = merge(1.,0., state==s%prev_state)
        mi_proj  = merge(1.,0., iproj==s%prev_proj)
        call build_glob%spproj_field%set(s%iptcl,'state',    real(state))
        call build_glob%spproj_field%set(s%iptcl,'corr',     corr)
        call build_glob%spproj_field%set(s%iptcl,'proj',     real(iproj))
        call build_glob%spproj_field%set(s%iptcl,'mi_proj',  mi_proj)
        call build_glob%spproj_field%set(s%iptcl,'mi_state', mi_state)
    end subroutine statesrch_cluster3D

    subroutine symsrch_cluster3D(s, spec)
        type(strategy3D_srch), intent(inout) :: s
        type(strategy3D_spec), intent(in)    :: spec
        integer :: projs(s%nstates),iproj,iref,isym,state,iproj_sym ! ineigh
        real    :: corrs(s%nstates),sym_corrs(s%nsym),corrs_inpl(s%nrots)
        real    :: corr,mi_state,mi_proj
        corrs = -1.
        projs = 0
        ! evaluate all correlations
        do state = 1, s%nstates
            if( .not. s3D%state_exists(state) ) cycle
            ! if( s%neigh )then
            !     ! greedy in symmetric unit & neighbourhood
            !     do ineigh=1,s%nnn
            !         iproj = build_glob%nnmat(s%prev_proj,ineigh)
            !         do isym=1,s%nsym
            !             iproj_sym = spec%symmat(iproj,isym)
            !             iref      = (state-1)*s%nprojs+iproj_sym
            !             call pftcc_glob%gencorrs(iref,s%iptcl,corrs_inpl)
            !             corr = maxval(corrs_inpl)
            !             if( corr > corrs(state) )then
            !                 corrs(state) = corr
            !                 projs(state) = iproj_sym
            !             endif
            !         enddo
            !     enddo
            ! else
                ! greedy in symmetric unit
                do isym=1,s%nsym
                    iproj = spec%symmat(s%prev_proj, isym)
                    iref  = (state-1)*s%nprojs+iproj
                    call pftcc_glob%gencorrs(iref,s%iptcl,corrs_inpl)
                    sym_corrs(isym) = maxval(corrs_inpl)
                enddo
                isym = maxloc(sym_corrs, dim=1)
                corrs(state) = sym_corrs(isym)
                projs(state) = spec%symmat(s%prev_proj, isym)
            ! endif
        enddo
        ! makes move
        if( s%prev_corr < spec%extr_score_thresh )then
            ! state randomization
            state = irnd_uni(s%nstates)
            if( s%nstates > 2 )then
                do while( state==s%prev_state .or. .not.s3D%state_exists(state))
                    state = irnd_uni(s%nstates)
                enddo
            else
                do while(.not.s3D%state_exists(state))
                    state = irnd_uni(s%nstates)
                enddo
            endif
            s%nrefs_eval = 1
            iproj = s%prev_proj
            corr  = corrs(state)
        else
            ! SHC state optimization
            s%prev_corr  = corrs(s%prev_state)
            state        = shcloc(s%nstates, corrs, s%prev_corr)
            s%nrefs_eval = count(corrs <= s%prev_corr)
            iproj        = projs(state)
            iref         = (state-1)*s%nprojs+iproj
            call build_glob%spproj_field%e1set(s%iptcl,s3D%proj_space_euls(s%ithr,iref,1))
            call build_glob%spproj_field%e2set(s%iptcl,s3D%proj_space_euls(s%ithr,iref,2))
            corr = corrs(state)
            if(state == s%prev_state) call greedy_inplsrch(s, corr, state, iproj)
        endif
        mi_state = merge(1.,0., state==s%prev_state)
        mi_proj  = merge(1.,0., iproj==s%prev_proj)
        call build_glob%spproj_field%set(s%iptcl,'state',    real(state))
        call build_glob%spproj_field%set(s%iptcl,'corr',     corr)
        call build_glob%spproj_field%set(s%iptcl,'proj',     real(iproj))
        call build_glob%spproj_field%set(s%iptcl,'mi_proj',  mi_proj)
        call build_glob%spproj_field%set(s%iptcl,'mi_state', mi_state)
    end subroutine symsrch_cluster3D

    subroutine greedy_inplsrch(s, corr, istate, iproj)
        type(strategy3D_srch), intent(inout) :: s
        real,                  intent(inout) :: corr
        integer,               intent(in)    :: istate, iproj
        integer :: iref
        iref = (istate-1)*s%nprojs+iproj
        s3D%proj_space_refinds_sorted_highest(s%ithr,s%nrefs) = iref ! inpl angle + shift
        s3D%proj_space_corrs_srchd(s%ithr,iref) = .true.
        call s%inpl_srch
        if( s3D%proj_space_corrs(s%ithr,iref) > corr )then
            corr = s3D%proj_space_corrs(s%ithr,iref)
            call build_glob%spproj_field%set_shift(s%iptcl,s%prev_shvec+s3D%proj_space_shift(s%ithr,iref,:))
            call build_glob%spproj_field%e3set(s%iptcl,s3D%proj_space_euls(s%ithr,iref,3))
        endif
    end subroutine greedy_inplsrch

    subroutine oris_assign_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        real :: frac
        frac = 100.*real(self%s%nrefs_eval) / real(self%s%nstates)
        call build_glob%spproj_field%set(self%s%iptcl,'frac',frac)
        call build_glob%spproj_field%set(self%s%iptcl,'w',   1.)
    end subroutine oris_assign_cluster3D

    subroutine kill_cluster3D( self )
        class(strategy3D_cluster), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D

end module simple_strategy3D_cluster
