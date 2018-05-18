! concrete strategy3D: continuous single-state refinement
module simple_strategy3D_cont_single
include 'simple_lib.f08'
use simple_strategy3D_alloc    ! singleton
use simple_strategy3D_utils,   only: fit_bfactors
use simple_strategy3D,         only: strategy3D
use simple_strategy3D_srch,    only: strategy3D_srch, strategy3D_spec
use simple_pftcc_orisrch_grad, only: pftcc_orisrch_grad
use simple_ori,              only: ori
use simple_singletons
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
        call self%cont_srch%new(spec%ppftcc)
    end subroutine new_cont_single

    subroutine srch_cont_single( self )
        class(strategy3D_cont_single), intent(inout) :: self
        real, allocatable :: cxy(:)
        ! execute search
        if( b%a%get_state(self%s%iptcl) > 0 )then
            ! initialize
            call self%s%prep4srch()
            call self%cont_srch%set_particle(self%s%iptcl)
            self%o    = b%a%get_ori(self%s%iptcl)
            cxy       = self%cont_srch%minimize(self%o, NPEAKSATHRES/2.0, p%trs)
            self%corr = cxy(1)
            ! prepare weights and orientations
            call self%oris_assign
        else
            call b%a%reject(self%s%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_CONT_SINGLE :: FINISHED CONTINUOUS SEARCH'
    end subroutine srch_cont_single

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine oris_assign_cont_single( self )
        use simple_ori,  only: ori
        class(strategy3D_cont_single), intent(inout) :: self
        type(ori) :: osym
        real      :: ws(1), dist_inpl, euldist, mi_proj, mi_inpl, mi_joint, frac, ang_sdev
        ! B factors
        ws = 1.0
        call fit_bfactors(self%s, ws)
        ! angular standard deviation
        ang_sdev = 0.
        ! angular distances
        call b%se%sym_dists(b%a%get_ori(self%s%iptcl), self%o, osym, euldist, dist_inpl)
        ! generate convergence stats
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj  = 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%irot == 0 .or. self%s%prev_roind == self%irot )then
            mi_inpl  = 1.
            mi_joint = mi_joint + 1.
        endif
        mi_joint = mi_joint / 2.
        call b%a%set(self%s%iptcl, 'mi_proj',   mi_proj)
        call b%a%set(self%s%iptcl, 'mi_inpl',   mi_inpl)
        call b%a%set(self%s%iptcl, 'mi_state',  1.)
        call b%a%set(self%s%iptcl, 'mi_joint',  mi_joint)
        ! fraction of search space scanned
        frac = 100.
        ! set the distances before we update the orientation
        if( b%a%isthere(self%s%iptcl,'dist') )then
            call b%a%set(self%s%iptcl, 'dist', 0.5*euldist + 0.5*b%a%get(self%s%iptcl,'dist'))
        else
            call b%a%set(self%s%iptcl, 'dist', euldist)
        endif
        call b%a%set(self%s%iptcl, 'dist_inpl', dist_inpl)
        call b%a%set_euler(self%s%iptcl, self%o%get_euler())
        call b%a%set_shift(self%s%iptcl, self%o%get_2Dshift())
        call b%a%set(self%s%iptcl, 'frac',      frac)
        call b%a%set(self%s%iptcl, 'state',     1.)
        call b%a%set(self%s%iptcl, 'corr',      self%corr)
        call b%a%set(self%s%iptcl, 'specscore', self%s%specscore)
        call b%a%set(self%s%iptcl, 'ow',        1.0)
        call b%a%set(self%s%iptcl, 'sdev',      0.)
        call b%a%set(self%s%iptcl, 'npeaks',    1.)
        ! transfer data to o_peaks
        call s3D%o_peaks(self%s%iptcl)%set_ori(1,b%a%get_ori(self%s%iptcl))
    end subroutine oris_assign_cont_single

    subroutine kill_cont_single( self )
        class(strategy3D_cont_single),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_cont_single

end module simple_strategy3D_cont_single
