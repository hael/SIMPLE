module simple_motion_anisocor
include 'simple_lib.f08'
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_image,       only: image
implicit none
private
public :: motion_anisocor, POLY_DIM, TOL
#include "simple_local_flags.inc"

real,    parameter                  :: TOL      = 1e-7      !< tolerance parameter
integer, parameter                  :: POLY_DIM = 2    !12  !< dimensionality of polynomial dome model

type, abstract, public :: motion_anisocor
    type(opt_spec)                  :: ospec                !< optimizer specification object
    class(optimizer),   pointer     :: nlopt      => null() !< pointer to nonlinear optimizer
    real                            :: maxHWshift = 0.      !< maximum half-width of shift
    class(image),       pointer     :: reference  => null() !< reference image ptr
    class(image),       pointer     :: frame      => null() !< particle image ptr
    integer                         :: ldim(2)              !< dimensions of reference, particle
    integer                         :: ldim_out(2)          !< dimensions of output image
    real                            :: motion_correctftol
    real                            :: motion_correctgtol
contains
    procedure(motion_anisocor_calc_T_out), deferred :: calc_T_out
    procedure(motion_anisocor_new       ), deferred :: new
    procedure(motion_anisocor_kill      ), deferred :: kill
    procedure(motion_anisocor_minimize  ), deferred :: minimize
end type motion_anisocor

interface
    subroutine motion_anisocor_calc_T_out( self, a, frame_in, frame_out )
        import motion_anisocor, POLY_DIM, image
        class(motion_anisocor), intent(inout) :: self
        real(kind=8),           intent(in)    :: a(POLY_DIM)
        class(image), target,   intent(in)    :: frame_in
        class(image), target,   intent(inout) :: frame_out
    end subroutine motion_anisocor_calc_T_out

    subroutine motion_anisocor_new( self, motion_correct_ftol, motion_correct_gtol )
        import motion_anisocor
        class(motion_anisocor), intent(inout) :: self
        real        , optional, intent(in)    :: motion_correct_ftol, motion_correct_gtol
    end subroutine motion_anisocor_new

    subroutine motion_anisocor_kill( self )
        import motion_anisocor
        class(motion_anisocor), intent(inout) :: self
    end subroutine motion_anisocor_kill

    function motion_anisocor_minimize( self, ref, frame, corr, regu ) result(cxy)
        import motion_anisocor, POLY_DIM, image
        class(motion_anisocor), intent(inout) :: self
        class(image), target,   intent(in)    :: ref, frame
        real(kind=8),           intent(out)   :: corr, regu
        real(kind=8) :: cxy(POLY_DIM+1)
    end function motion_anisocor_minimize

end interface

end module simple_motion_anisocor
