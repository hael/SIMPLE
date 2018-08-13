! concrete strategy3D: stochastic state quantization (3D clustering)
module simple_strategy3D_cluster_snhc
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! use all in there
use simple_strategy3D_utils
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
        class(strategy3D_cluster_snhc), intent(inout) :: self
        integer,                        intent(in)    :: ithr
        real    :: inpl_corrs(self%s%nrots)
        integer :: iref, isample, nrefs
        self%s%prev_state = build_glob%spproj_field%get_state(self%s%iptcl)
        if( self%s%prev_state > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch(build_glob%nnmat)
            self%s%nrefs_eval = 0
            nrefs = min(self%s%nnnrefs,nint(self%s%nnnrefs*(1.-self%spec%extr_score_thresh)))
            ! correlations
            do isample=1,nrefs
                ! set the stochastic reference index
                iref = s3D%srch_order(self%s%ithr,isample)
                ! search
                call per_ref_srch
            enddo
            ! sort in correlation projection direction space
            call sort_corrs(self%s)
            ! take care of the in-planes
            if(nrefs == self%s%nnnrefs) call self%s%inpl_srch
            ! prepare weights & orientation
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint   '>>> strategy3D_cluster_snhc :: FINISHED SRCH_cluster3D_snhc'
        contains

            subroutine per_ref_srch
                integer :: loc(3)
                if( s3D%state_exists( s3D%proj_space_state(iref) ) )then
                    ! identify the 3 top scoring in-planes
                    call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                    loc = max3loc(inpl_corrs)
                    call self%s%store_solution(iref, loc, [inpl_corrs(loc(1)), inpl_corrs(loc(2)), inpl_corrs(loc(3))], .true.)
                    ! keep track of how many references we are evaluating
                    self%s%nrefs_eval = self%s%nrefs_eval + 1
                endif
            end subroutine per_ref_srch

    end subroutine srch_cluster3D_snhc

    subroutine oris_assign_cluster3D_snhc( self )
        class(strategy3D_cluster_snhc), intent(inout) :: self
        integer :: loc(1), iproj, iref, istate, roind
        real    :: shvec(2), euls(3), corr, mi_state, frac, mi_inpl, mi_proj
        iref     = s3D%proj_space_refinds_sorted_highest(self%s%ithr,self%s%nrefs)
        iproj    = s3D%proj_space_proj(iref)
        corr     = s3D%proj_space_corrs(self%s%ithr,iref,1)
        istate   = s3D%proj_space_state(iref)
        roind    = s3D%proj_space_inplinds(self%s%ithr,iref,1)
        euls     = s3D%proj_space_euls(self%s%ithr,iref,1,1:3)
        shvec    = s3D%proj_space_shift(self%s%ithr,iref,1,1:2)
        frac     = 100.*real(self%s%nrefs_eval) / real(self%s%nnnrefs)
        mi_state = 0.
        mi_inpl  = 0.
        mi_proj  = 0.
        if(self%s%prev_roind==roind)  mi_inpl  = 1.
        if(self%s%prev_state==istate) mi_state = 1.
        if(self%s%prev_proj ==iproj)   mi_proj  = 1.
        call build_glob%spproj_field%set(self%s%iptcl,'frac',     frac)
        call build_glob%spproj_field%set(self%s%iptcl,'state',    real(istate))
        call build_glob%spproj_field%set(self%s%iptcl,'corr',     corr)
        call build_glob%spproj_field%set(self%s%iptcl,'mi_proj',  mi_proj)
        call build_glob%spproj_field%set(self%s%iptcl,'mi_inpl',  mi_inpl)
        call build_glob%spproj_field%set(self%s%iptcl,'mi_state', mi_state)
        call build_glob%spproj_field%set(self%s%iptcl,'mi_joint', (mi_state+mi_inpl+mi_inpl)/3.)
        call build_glob%spproj_field%set(self%s%iptcl,'npeaks',   1.)
        call build_glob%spproj_field%set(self%s%iptcl,'ow',       1.)
        call build_glob%spproj_field%set(self%s%iptcl,'w',        1.)
        call build_glob%spproj_field%set(self%s%iptcl,'specscore',self%s%specscore)
        call build_glob%spproj_field%set_euler(self%s%iptcl, euls)
        call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
        call s3D%o_peaks(self%s%iptcl)%set(1,'state', real(istate))
        call s3D%o_peaks(self%s%iptcl)%set(1,'corr',  corr)
        call s3D%o_peaks(self%s%iptcl)%set(1,'w',     1.)
        call s3D%o_peaks(self%s%iptcl)%set(1,'ow',    1.)
        call s3D%o_peaks(self%s%iptcl)%set_euler(1, euls)
        call s3D%o_peaks(self%s%iptcl)%set_shift(1, shvec)
    end subroutine oris_assign_cluster3D_snhc

    subroutine kill_cluster3D_snhc( self )
        class(strategy3D_cluster_snhc), intent(inout) :: self
        call self%s%kill
    end subroutine kill_cluster3D_snhc

end module simple_strategy3D_cluster_snhc
