! concrete strategy3D: stochastic neighbourhood hill-climbing
module simple_strategy3D_snhc_single
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_builder,          only: build_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_snhc_single
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_snhc_single
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_snhc_single
    procedure          :: srch        => srch_snhc_single
    procedure          :: oris_assign => oris_assign_snhc_single
    procedure          :: kill        => kill_snhc_single
end type strategy3D_snhc_single

contains

    subroutine new_snhc_single( self, spec, npeaks )
        class(strategy3D_snhc_single), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        integer,                       intent(in)    :: npeaks
        call self%s%new( spec, npeaks )
        self%spec = spec
    end subroutine new_snhc_single

    subroutine srch_snhc_single( self )
        class(strategy3D_snhc_single), intent(inout) :: self
        integer :: iref, isample, loc(1), inpl_ind
        real    :: corrs(self%s%nrefs), inpl_corrs(self%s%nrots), inpl_corr
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! initialize
            call self%s%prep4srch
            self%s%nbetter    = 0
            self%s%nrefs_eval = 0
            s3D%proj_space_corrs(self%s%iptcl_map,:) = -1.
            ! search
            do isample=1,self%spec%szsn
                iref = s3D%srch_order(self%s%iptcl_map,isample) ! set the stochastic reference index
                call per_ref_srch                           ! actual search
            end do
            self%s%nrefs_eval = self%spec%szsn
            ! sort in correlation projection direction space
            corrs = s3D%proj_space_corrs(self%s%iptcl_map,:)
            call hpsort(corrs, s3D%proj_space_inds(self%s%iptcl_map,:))
            ! output
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_SNHC_SINGLE :: FINISHED SEARCH'

        contains

            subroutine per_ref_srch
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs) ! In-plane correlations
                loc        = maxloc(inpl_corrs)               ! greedy in-plane
                inpl_ind   = loc(1)                           ! in-plane angle index
                inpl_corr  = inpl_corrs(inpl_ind)             ! max in plane correlation
                call self%s%store_solution(iref, iref, inpl_ind, inpl_corr)
            end subroutine per_ref_srch

    end subroutine srch_snhc_single

    subroutine oris_assign_snhc_single( self )
        use simple_ori,  only: ori
        class(strategy3D_snhc_single), intent(inout) :: self
        type(ori)  :: osym
        real       :: dist_inpl, corr, frac, euldist, bfac
        integer    :: ref, roind
        ! orientation parameters
        ref = s3D%proj_space_inds(self%s%iptcl_map, self%s%nrefs)
        if( ref < 1 .or. ref > self%s%nrefs )then
            print *, 'ref: ', ref
            write(*,*)'ref index out of bound; strategy3d_snhc_single :: oris_assign_snhc_single'
            stop 'ref index out of bound; strategy3d_snhc_single :: oris_assign_snhc_single'
        endif
        roind = pftcc_glob%get_roind(360. - s3D%proj_space_euls(self%s%iptcl_map, ref, 3))
        ! transfer to solution set
        corr = max(0., s3D%proj_space_corrs(self%s%iptcl_map,ref))
        call s3D%o_peaks(self%s%iptcl)%set(1, 'state', 1.)
        call s3D%o_peaks(self%s%iptcl)%set(1, 'proj',  real(s3D%proj_space_proj(self%s%iptcl_map,ref)))
        call s3D%o_peaks(self%s%iptcl)%set(1, 'corr',  corr)
        call s3D%o_peaks(self%s%iptcl)%set_euler(1, s3D%proj_space_euls(self%s%iptcl_map,ref,1:3))
        call s3D%o_peaks(self%s%iptcl)%set_shift(1, [0.,0.]) ! no shift search in snhc
        call s3D%o_peaks(self%s%iptcl)%set(1, 'ow', 1.0)
        ! B factor
        if( pftcc_glob%objfun_is_ccres() )then
            bfac  = pftcc_glob%fit_bfac(ref, self%s%iptcl, roind, [0.,0.])
            call build_glob%spproj_field%set(self%s%iptcl, 'bfac',  bfac )
        endif
        ! angular distances
        call build_glob%pgrpsyms%sym_dists( build_glob%spproj_field%get_ori(self%s%iptcl),&
            &s3D%o_peaks(self%s%iptcl)%get_ori(1), osym, euldist, dist_inpl)
        ! fraction search space
        frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs)
        ! set the overlaps
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_proj',   0.)
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_inpl',   0.)
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_state',  1.)
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_joint',  0.)
        if( build_glob%spproj_field%isthere(self%s%iptcl,'dist') )then
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(self%s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call build_glob%spproj_field%set(self%s%iptcl, 'state',     1.)
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      frac)
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',      corr)
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call build_glob%spproj_field%set(self%s%iptcl, 'proj',      s3D%o_peaks(self%s%iptcl)%get(1,'proj'))
        call build_glob%spproj_field%set(self%s%iptcl, 'sdev',      0.)
        call build_glob%spproj_field%set(self%s%iptcl, 'npeaks',    1.)
        call build_glob%spproj_field%set_euler(self%s%iptcl, s3D%proj_space_euls(self%s%iptcl_map,ref,1:3))
        call build_glob%spproj_field%set_shift(self%s%iptcl, [0.,0.]) ! no shift search in snhc
        DebugPrint  '>>> STRATEGY3D_SNHC_SINGLE :: EXECUTED ORIS_ASSIGN_SNHC_SINGLE'
    end subroutine oris_assign_snhc_single

    subroutine kill_snhc_single( self )
        class(strategy3D_snhc_single),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc_single

end module simple_strategy3D_snhc_single
