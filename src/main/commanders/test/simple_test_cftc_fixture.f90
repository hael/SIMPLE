!@descr: cartft test fixture shared by testers: one particle's objective and five-parameter pose gradient
module simple_test_cftc_fixture
use simple_core_module_api, only: dp
use simple_cartft_calc,     only: cartft_calc
implicit none
private
public :: objective_gradient

contains

    !> pose_normal_terms without its Gauss-Newton matrix: the objective (1 - cc or the normalized loss L)
    !! and its derivatives in the three right-increment rotation components (radians) and two shifts (pixels)
    subroutine objective_gradient( calc, state, iseven, iptcl, rotmat, shift, objective, gradient )
        class(cartft_calc), intent(in)  :: calc
        integer,            intent(in)  :: state, iptcl
        logical,            intent(in)  :: iseven
        real(dp),           intent(in)  :: rotmat(3,3), shift(2)
        real(dp),           intent(out) :: objective, gradient(5)
        real(dp) :: hessian(5,5)
        call calc%pose_normal_terms(state, iseven, iptcl, rotmat, shift, objective, gradient, hessian)
    end subroutine objective_gradient

end module simple_test_cftc_fixture
