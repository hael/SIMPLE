! shift search with L-BFGS-B using expanded Fourier transforms (used in motion_correct)
module simple_ftexp_shsrch
include 'simple_lib.f08'
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded, ft_exp_reset_tmp_pointers
use simple_image,       only: image
implicit none

public :: ftexp_shsrch_init, ftexp_shsrch_set_ptrs, ftexp_shsrch_minimize, test_ftexp_shsrch
private

type(opt_factory)             :: ofac              !< optimizer factory
type(opt_spec)                :: ospec             !< optimizer specification object
class(optimizer),   pointer   :: nlopt=>null()     !< pointer to nonlinear optimizer
class(ft_expanded), pointer   :: reference=>null() !< reference pft
class(ft_expanded), pointer   :: particle=>null()  !< particle pft
real,    parameter            :: TOL=1e-4          !< tolerance parameter
integer, parameter            :: MAXITS=30         !< maximum number of iterations
real                          :: maxHWshift = 0.   !< maximum half-width of shift
real                          :: motion_correctftol = 1e-4
real                          :: motion_correctgtol = 1e-4

contains

    !> Initialise  ftexp_shsrch
    subroutine ftexp_shsrch_init( ref, ptcl, trs, motion_correct_ftol, motion_correct_gtol )
        class(ft_expanded), target, intent(in) :: ref, ptcl
        real,                       intent(in) :: trs
        real,             optional, intent(in) :: motion_correct_ftol, motion_correct_gtol
        if( present(motion_correct_ftol) )then
            motion_correctftol = motion_correct_ftol
        else
            motion_correctftol = TOL
        end if
        if( present(motion_correct_gtol) )then
            motion_correctgtol = motion_correct_gtol
        else
            motion_correctgtol = TOL
        end if
        reference  => ref
        particle   => ptcl
        maxHWshift =  trs
    end subroutine ftexp_shsrch_init

    subroutine ftexp_shsrch_set_ptrs( ref, ptcl )
        class(ft_expanded), intent(in), target :: ref, ptcl
        reference => ref
        particle  => ptcl
    end subroutine ftexp_shsrch_set_ptrs

    !> Cost function, double precision
    function ftexp_shsrch_cost_8( fun_self, vec, D ) result(cost)
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
        cost = -reference%corr_shifted_8(particle, -vec)
    end function ftexp_shsrch_cost_8

    !> Gradient function, double precision
    subroutine ftexp_shsrch_gcost_8( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(dp), intent(inout) :: vec( D )
        real(dp), intent(out)   :: grad( D )
        call reference%corr_gshifted_8(particle, -vec, grad)
    end subroutine ftexp_shsrch_gcost_8

    !> Gradient & cost function, double precision
    subroutine ftexp_shsrch_fdfcost_8( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        call reference%corr_fdfshifted_8(particle, -vec, f, grad)
        f    = f    * (-1.0_dp)
    end subroutine ftexp_shsrch_fdfcost_8

    !> Main search routine
    function ftexp_shsrch_minimize( prev_corr, prev_shift ) result( cxy )
        real, optional, intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3), lims(2,2) ! max_shift, maxlim
        class(*), pointer :: fun_self => null()
        lims(1,1) = prev_shift(1) - maxHWshift
        lims(1,2) = prev_shift(1) + maxHWshift
        lims(2,1) = prev_shift(2) - maxHWshift
        lims(2,2) = prev_shift(2) + maxHWshift
        call ospec%specify('lbfgsb', 2, ftol=motion_correctftol, gtol=motion_correctgtol, limits=lims) ! nrestarts=nrestarts
        call ospec%set_costfun_8(ftexp_shsrch_cost_8)
        call ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
        call ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(nlopt) )then
            call nlopt%kill
            deallocate(nlopt)
        end if
        call ofac%new(ospec, nlopt)
        ! set initial solution to previous shift
        ospec%x = prev_shift
        call nlopt%minimize(ospec, fun_self, cxy(1))
        call reference%corr_normalize(particle, cxy(1))
        call ft_exp_reset_tmp_pointers
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= TOL )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
    end function ftexp_shsrch_minimize

    subroutine test_ftexp_shsrch
        type(image)       :: img_ref, img_ptcl
        type(ft_expanded) :: ftexp_ref, ftexp_ptcl
        real, parameter   :: TRS=5.0, lp=6., hp=100.
        real              :: cxy(3), x, y, lims(2,2)
        integer           :: i
        lims(:,1) = -TRS
        lims(:,2) = TRS
        call img_ref%new([100,100,1], 2.)
        call img_ptcl%new([100,100,1], 2.)
        call img_ref%square(20)
        call img_ref%fft()
        call ftexp_ref%new(img_ref, hp, lp)
        call ftexp_ptcl%new(img_ptcl, hp, lp)
        call ftexp_shsrch_init(ftexp_ref, ftexp_ptcl, trs)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            img_ptcl = img_ref
            call img_ptcl%shift([x,y,0.])
            call ftexp_ptcl%new(img_ptcl, hp, lp)
            cxy = ftexp_shsrch_minimize()
            if( cxy(1) < 0.995 )then
                stop 'shift alignment failed; test_ftexp_shsrch; simple_ftexp_shsrch'
            endif
        end do
        write(*,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftexp_shsrch

end module simple_ftexp_shsrch
