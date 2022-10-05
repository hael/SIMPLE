! concrete strategy3D: stochastic neighbourhood hill-climbing
module simple_strategy3D_snhc
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_srch, strategy3D_spec
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: strategy3D_snhc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_snhc
    private
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_snhc
    procedure          :: srch        => srch_snhc
    procedure          :: oris_assign => oris_assign_snhc
    procedure          :: kill        => kill_snhc
end type strategy3D_snhc

contains

    subroutine new_snhc( self, spec )
        class(strategy3D_snhc), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_snhc

    subroutine srch_snhc( self, ithr )
        class(strategy3D_snhc), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        integer :: iref, isample
        real    ::  inpl_corrs(self%s%nrots)
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! initialize
            call self%s%prep4srch
            self%s%nbetter = 0
            ! search
            do isample=1,self%spec%szsn
                iref = s3D%srch_order(self%s%ithr,isample)  ! set the stochastic reference index
                call per_ref_srch                           ! actual search
            end do
            self%s%nrefs_eval = self%spec%szsn
            ! output
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif

        contains

            subroutine per_ref_srch
                integer :: loc(1)
                ! calculate in-plane correlations
                call pftcc_glob%gencorrs(iref, self%s%iptcl, inpl_corrs)
                ! identify the top scoring in-plane angle
                loc = maxloc(inpl_corrs)
                ! stash
                call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
            end subroutine per_ref_srch

    end subroutine srch_snhc

    subroutine oris_assign_snhc( self )
        use simple_ori, only: ori
        class(strategy3D_snhc), intent(inout) :: self
        type(ori)  :: osym, o1, o2
        real       :: dist_inpl, corr, frac, euldist
        integer    :: ref, roind, loc(1)
        ! stash prev ori
        call build_glob%spproj_field%get_ori(self%s%iptcl, o1)
        ! orientation parameters
        loc = maxloc(s3D%proj_space_corrs(self%s%ithr,:))
        ref = loc(1)
        if( ref < 1 .or. ref > self%s%nrefs )then
            THROW_HARD('ref index: '//int2str(ref)//' out of bound; oris_assign_snhc')
        endif
        roind = pftcc_glob%get_roind(360. - s3D%proj_space_euls(3,ref,self%s%ithr))
        ! transfer to spproj_field
        corr = max(0., s3D%proj_space_corrs(self%s%ithr,ref))
        call build_glob%spproj_field%set(self%s%iptcl, 'state', 1.)
        call build_glob%spproj_field%set(self%s%iptcl, 'proj',  real(s3D%proj_space_proj(ref)))
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',  corr)
        call build_glob%spproj_field%set_euler(self%s%iptcl, s3D%proj_space_euls(:,ref,self%s%ithr))
        call build_glob%spproj_field%set_shift(self%s%iptcl, [0.,0.]) ! no shift search in snhc
        ! angular distances
        call build_glob%spproj_field%get_ori(self%s%iptcl, o2)
        call build_glob%pgrpsyms%sym_dists( o1, o2, osym, euldist, dist_inpl)
        ! fraction search space
        frac = 100.*real(self%s%nrefs_eval) / real(self%s%nprojs)
        ! set the overlaps
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_proj',   0.)
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_state',  1.)
        if( build_glob%spproj_field%isthere(self%s%iptcl,'dist') )then
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(self%s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      frac)
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call osym%kill
        call o1%kill
        call o2%kill
    end subroutine oris_assign_snhc

    subroutine kill_snhc( self )
        class(strategy3D_snhc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_snhc

end module simple_strategy3D_snhc
