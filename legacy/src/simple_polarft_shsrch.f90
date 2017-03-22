module simple_polarft_shsrch
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_polarft,     only: polarft
implicit none


public :: polarft_shsrch_init, polarft_shsrch_set_rot, polarft_shsrch_minimize
private

type(opt_factory)         :: ofac              !< optimizer factory
type(opt_spec)            :: ospec             !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()     !< pointer to nonlinear optimizer
class(polarft), pointer   :: reference=>null() !< reference pft
class(polarft), pointer   :: particle=>null()  !< particle pft
integer                   :: rot               !< in-plane rotation

contains
    
    subroutine polarft_shsrch_init( ref, ptcl, lims )
        class(polarft), intent(in), target :: ref, ptcl
        real, intent(in) :: lims(2,2)
        ! set pointers
        reference => ref 
        particle  => ptcl
        ! make optimizer spec
        call ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=3)
        ! set optimizer cost function
        call ospec%set_costfun(polarft_shsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
    end subroutine
    
    subroutine polarft_shsrch_set_rot( r )
        integer, intent(in) :: r
        rot = r
    end subroutine

    ! function polarft_shsrch_grad_cost( vec, D ) result( grad_cost )
    !   integer, intent(in) :: D
    !   real, intent(in)    :: vec(D)
    !   real :: grad_cost
    !   !TODO: later or never
    ! end function polarft_shsrch_grad_cost
    
    function polarft_shsrch_cost( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real :: cost
        cost = -reference%corr_shifted(particle, rot, vec)
    end function
    
    function polarft_shsrch_minimize( ) result( cxy )
        use simple_math, only: rotmat2d
        real :: cxy(3), m(2,2)
        ospec%x = 0.
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation 
        ! rotate the shift vector to the frame of reference
        cxy(2:) = ospec%x ! shift
        m = rotmat2d(reference%get_rot(rot)) 
        cxy(2:) = matmul(cxy(2:),m)
    end function
    
end module simple_polarft_shsrch
