! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_clustersoft
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils  ! use all in there
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_builder,          only: build_glob
use simple_polarft_corrcalc, only: pftcc_glob
!use simple_parameters,       only: params_glob
use simple_rnd,              only: shcloc, irnd_uni
use simple_ori,              only: ori
implicit none

public :: strategy3D_clustersoft
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_clustersoft
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_clustersoft3D
    procedure          :: srch        => srch_clustersoft3D
    procedure          :: oris_assign => oris_assign_clustersoft3D
    procedure          :: kill        => kill_clustersoft3D
end type strategy3D_clustersoft

contains

    subroutine new_clustersoft3D( self, spec, npeaks )
        class(strategy3D_clustersoft), intent(inout) :: self
        class(strategy3D_spec),    intent(inout) :: spec
        integer,                   intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_clustersoft3D

    subroutine srch_clustersoft3D( self, ithr )
        class(strategy3D_clustersoft), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        type(ori) :: o
        real      :: wcorrs(self%s%nstates), mi_state, frac
        !real      :: corrs(self%s%nstates,s3D%o_peaks(self%s%iptcl)%get_noris())
        integer   :: s, npeaks
        ! execute search
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! set thread index
            ! self%s%ithr = ithr
            ! prep
            o = build_glob%spproj_field%get_ori(self%s%iptcl)
            self%s%prev_proj  = build_glob%eulspace%find_closest_proj(o)
            self%s%prev_ref   = (self%s%prev_state-1)*self%s%nprojs + self%s%prev_proj
            self%s%prev_roind = pftcc_glob%get_roind(360.-o%e3get())
            self%s%prev_shvec = o%get_2Dshift()
            self%s%prev_corr  = build_glob%spproj_field%get(self%s%iptcl, 'corr')
            ! particle prep
            call pftcc_glob%prep_matchfilt(self%s%iptcl,self%s%prev_ref,self%s%prev_roind)
            ! calc specscore
            self%s%specscore = pftcc_glob%specscore(self%s%prev_ref, self%s%iptcl, self%s%prev_roind)
            ! other init
            npeaks = s3D%o_peaks(self%s%iptcl)%get_noris()
            ! extremal optimisation fork
            if( self%s%prev_corr < self%spec%extr_score_thresh )then
                ! state randomization
                s = irnd_uni(self%s%nstates)
                if( self%s%nstates > 2 )then
                    do while( s==self%s%prev_state .or. .not.s3D%state_exists(s) )
                        s = irnd_uni(self%s%nstates)
                    enddo
                else
                    do while( .not.s3D%state_exists(s) )
                        s = irnd_uni(self%s%nstates)
                    enddo
                endif
                self%s%nrefs_eval = 1
                ! weighted correlation update
                call update_corr(s)
            else
                ! exhaustive neighbourhood evaluation & shc condition
                do s = 1,self%s%nstates
                    call update_corr(s)
                enddo
                self%s%prev_corr = wcorrs(self%s%prev_state)
                s = shcloc(self%s%nstates, wcorrs, self%s%prev_corr)
                self%s%nrefs_eval = count(wcorrs <= self%s%prev_corr)
            endif
            ! assignments
            frac     = 100. * real(self%s%nrefs_eval) / real(self%s%nstates)
            mi_state = merge(1.,0., s==self%s%prev_state)
            call build_glob%spproj_field%set(self%s%iptcl,'w',        1.) !?
            call build_glob%spproj_field%set(self%s%iptcl,'state',    real(s))
            call build_glob%spproj_field%set(self%s%iptcl,'corr',     wcorrs(s))
            call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  1.)
            call build_glob%spproj_field%set(self%s%iptcl,'mi_state', mi_state)
            call build_glob%spproj_field%set(self%s%iptcl,'specscore', self%s%specscore)
            call s3D%o_peaks(self%s%iptcl)%set_all2single('state',real(s))
            call o%kill
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_clustersoft :: FINISHED SEARCH'
        contains

            subroutine update_corr( istate )
                integer, intent(in) :: istate
                real      :: shvec(2), corr, ow
                integer   :: ipeak, roind, iref, iproj
                wcorrs(istate) = 0.
                if( .not.s3D%state_exists(s) )then
                    ! corrs(istate)  = 0.
                else
                    do ipeak = 1,npeaks
                        ow = s3D%o_peaks(self%s%iptcl)%get(ipeak,'ow')
                        if( ow < TINY )cycle
                        o     = s3D%o_peaks(self%s%iptcl)%get_ori(ipeak)
                        roind = pftcc_glob%get_roind(360.-o%e3get())
                        iproj = build_glob%eulspace%find_closest_proj(o)
                        shvec = s3D%o_peaks(self%s%iptcl)%get_2Dshift(ipeak) - self%s%prev_shvec
                        iref  = (istate-1)*self%s%nprojs + iproj
                        corr  = real(pftcc_glob%gencorr_for_rot_8(iref, self%s%iptcl, real(shvec,kind=dp), roind), kind=sp)
                        wcorrs(istate) = wcorrs(istate) + ow*corr
                        !corrs(istate)  = corr
                    enddo
                endif
            end subroutine update_corr

    end subroutine srch_clustersoft3D

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine oris_assign_clustersoft3D( self )
        class(strategy3D_clustersoft), intent(inout) :: self
        DebugPrint   '>>> STRATEGY3D_clustersoft :: EXECUTED ORIS_ASSIGN_MULTI'
    end subroutine oris_assign_clustersoft3D

    subroutine kill_clustersoft3D( self )
        class(strategy3D_clustersoft), intent(inout) :: self
        call self%s%kill
    end subroutine kill_clustersoft3D

    subroutine store_solution( self, iref, inpl_ind, corr )
        class(strategy3D_clustersoft), intent(in) :: self
        integer,                       intent(in) :: iref, inpl_ind
        real,                          intent(in) :: corr
        s3D%proj_space_inplinds(self%s%ithr,iref,1) = inpl_ind
        s3D%proj_space_euls(self%s%ithr,iref,1,1)   = 360. - pftcc_glob%get_rot(inpl_ind)
        s3D%proj_space_corrs(self%s%ithr,iref,1)    = corr
    end subroutine store_solution

end module simple_strategy3D_clustersoft
