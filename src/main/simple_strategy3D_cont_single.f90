! concrete strategy3D: continuous single-state refinement
module simple_strategy3D_cont_single
include 'simple_lib.f08'
use simple_strategy3D_alloc    ! singleton
use simple_strategy3D,         only: strategy3D
use simple_strategy3D_srch,    only: strategy3D_srch, strategy3D_spec
use simple_pftcc_orisrch_grad, only: pftcc_orisrch_grad
use simple_ori,                only: ori
use simple_parameters,         only: params_glob
use simple_builder,            only: build_glob
implicit none

public :: strategy3D_cont_single
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cont_single
    type(strategy3D_srch)    :: s
    type(strategy3D_spec)    :: spec
    type(pftcc_orisrch_grad) :: cont_srch
    type(ori) :: o
    integer   :: irot
    real      :: corr
contains
    procedure :: new         => new_cont_single
    procedure :: srch        => srch_cont_single
    procedure :: oris_assign => oris_assign_cont_single
    procedure :: kill        => kill_cont_single
end type strategy3D_cont_single

contains

    subroutine new_cont_single( self, spec, npeaks )
        class(strategy3D_cont_single), intent(inout) :: self
        class(strategy3D_spec),        intent(inout) :: spec
        integer,                       intent(in)    :: npeaks
        call self%s%new(spec, npeaks)
        self%spec = spec
        call self%cont_srch%new
    end subroutine new_cont_single

    subroutine srch_cont_single( self, ithr )
        class(strategy3D_cont_single), intent(inout) :: self
        integer,                       intent(in)    :: ithr
        real, allocatable :: cxy(:)
        logical :: found_better
        ! execute search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! initialize
            call self%s%prep4srch
            call self%cont_srch%set_particle(self%s%iptcl)
            call build_glob%spproj_field%get_ori(self%s%iptcl, self%o)
            cxy       = self%cont_srch%minimize(self%o, NPEAKSATHRES/2.0, params_glob%trs, found_better)
            self%corr = cxy(1)
            if( .not. found_better )then
                ! put back the original one
                call build_glob%spproj_field%get_ori(self%s%iptcl, self%o)
            endif
            ! prepare weights and orientations
            call self%oris_assign
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_cont_single

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine oris_assign_cont_single( self )
        use simple_ori,  only: ori
        class(strategy3D_cont_single), intent(inout) :: self
        type(ori) :: osym, o_tmp
        real      :: dist_inpl, euldist, mi_proj, frac
        ! angular distances
        call build_glob%spproj_field%get_ori(self%s%iptcl, o_tmp)
        call build_glob%pgrpsyms%sym_dists(o_tmp, self%o, osym, euldist, dist_inpl)
        ! generate convergence stats
        mi_proj  = 0.
        if( euldist < 0.5 ) mi_proj  = 1.
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_proj',   mi_proj)
        call build_glob%spproj_field%set(self%s%iptcl, 'mi_state',  1.)
        ! fraction of search space scanned
        frac = 100.
        ! set the distances before we update the orientation
        if( build_glob%spproj_field%isthere(self%s%iptcl,'dist') )then
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(self%s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(self%s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        call build_glob%spproj_field%set_euler(self%s%iptcl, self%o%get_euler())
        call build_glob%spproj_field%set_shift(self%s%iptcl, self%o%get_2Dshift())
        call build_glob%spproj_field%set(self%s%iptcl, 'frac',      frac)
        call build_glob%spproj_field%set(self%s%iptcl, 'state',     1.)
        call build_glob%spproj_field%set(self%s%iptcl, 'corr',      self%corr)
        call build_glob%spproj_field%set(self%s%iptcl, 'specscore', self%s%specscore)
        call build_glob%spproj_field%set(self%s%iptcl, 'ow',        1.0)
        call build_glob%spproj_field%set(self%s%iptcl, 'spread',    0.)
        call build_glob%spproj_field%set(self%s%iptcl, 'npeaks',    1.)
        ! transfer data to o_peaks
        call build_glob%spproj_field%get_ori(self%s%iptcl, o_tmp)
        call s3D%o_peaks(self%s%iptcl)%set_ori(1, o_tmp)
        call osym%kill
        call o_tmp%kill
    end subroutine oris_assign_cont_single

    subroutine kill_cont_single( self )
        class(strategy3D_cont_single),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_cont_single

end module simple_strategy3D_cont_single
