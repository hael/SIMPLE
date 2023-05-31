! rotational origin shift alignment used in the initial model regularization, gradient based minimizer
module simple_pftcc_reg_shsrch
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: pftcc_reg_shsrch
private
#include "simple_local_flags.inc"

type :: pftcc_reg_shsrch
    private
    integer  :: reference    = 0       !< reference pft
    integer  :: particle     = 0       !< particle pft
    integer  :: nrots        = 0       !< # rotations
    integer  :: maxits       = 100     !< max # iterations
    integer  :: cur_inpl_idx = 0       !< index of inplane angle for shift search
contains
    procedure :: new         => grad_shsrch_new
    procedure :: set_indices => grad_shsrch_set_indices
    procedure :: minimize    => grad_shsrch_minimize
    procedure :: kill        => grad_shsrch_kill
end type pftcc_reg_shsrch

contains

    subroutine grad_shsrch_new( self, maxits )
        use simple_opt_factory, only: opt_factory
        class(pftcc_reg_shsrch),   intent(inout) :: self           !< instance
        integer,          optional, intent(in)    :: maxits         !< maximum iterations
        call self%kill
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%nrots = pftcc_glob%get_nrots()
    end subroutine grad_shsrch_new

    !> set indicies for shift search
    subroutine grad_shsrch_set_indices( self, ref, ptcl )
        class(pftcc_reg_shsrch), intent(inout) :: self
        integer,                  intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine grad_shsrch_set_indices

    !> minimisation routine
    function grad_shsrch_minimize( self, irot ) result( cxy )
        class(pftcc_reg_shsrch), intent(inout) :: self
        integer,                  intent(inout) :: irot
        real     :: lowest_shift(2), cxy(3), lowest_cost, rotmat(2,2)
        real(dp) :: lowest_cost_overall
        self%cur_inpl_idx = irot
        lowest_cost_overall = -pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, [0.d0,0.d0], self%cur_inpl_idx)
        ! [TODO] shift search with gradient and get the updated lowest_cost, lowest_shift
        if( lowest_cost < lowest_cost_overall )then
            cxy(1)  = - real(lowest_cost)  ! correlation
            cxy(2:) =   lowest_shift               ! shift
            ! rotate the shift vector to the frame of reference
            call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
        else
            irot = 0 ! to communicate that a better solution was not found
        endif
    end function grad_shsrch_minimize

    subroutine grad_shsrch_kill( self )
        class(pftcc_reg_shsrch), intent(inout) :: self
    end subroutine grad_shsrch_kill
end module simple_pftcc_reg_shsrch