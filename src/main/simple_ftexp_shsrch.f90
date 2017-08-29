! shift search using expanded Fourier transforms (used in unblur)
module simple_ftexp_shsrch
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded
use simple_image,       only: image
implicit none

public :: ftexp_shsrch_init, ftexp_shsrch_reset_ptrs, ftexp_shsrch_minimize, test_ftexp_shsrch
private

type(opt_factory)             :: ofac              !< optimizer factory
type(opt_spec)                :: ospec             !< optimizer specification object
class(optimizer),   pointer   :: nlopt=>null()     !< pointer to nonlinear optimizer
class(ft_expanded), pointer   :: reference=>null() !< reference pft
class(ft_expanded), pointer   :: particle=>null()  !< particle pft
character(len=:), allocatable :: opt_str           !< optimiser string descriptor
integer                       :: nrestarts=3       !< number of optimisation restarts
real,    parameter            :: TOL=1e-4          !< tolerance parameter
integer, parameter            :: MAXITS=30         !< maximum number of iterations

contains

    !> Initialise  ftexp_shsrch
    subroutine ftexp_shsrch_init( ref, ptcl, lims, opt, nrestarts_in )
        class(ft_expanded), target, intent(in) :: ref, ptcl
        real,                       intent(in) :: lims(2,2)
        character(len=*), optional, intent(in) :: opt
        integer,          optional, intent(in) :: nrestarts_in
        ! set pointers
        reference => ref
        particle  => ptcl
        ! set opt_str & nrestarts
        if( allocated(opt_str) ) deallocate(opt_str)
        if( present(opt) )then
            allocate(opt_str, source=opt)
        else
            allocate(opt_str, source='simplex')
        endif
        nrestarts = 1
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer spec
        call ospec%specify(opt_str, 2, ftol=TOL, gtol=TOL,&
        limits=lims, nrestarts=nrestarts )
        ! set optimizer cost function
        call ospec%set_costfun(ftexp_shsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
    end subroutine ftexp_shsrch_init

    subroutine ftexp_shsrch_reset_ptrs( ref, ptcl )
        class(ft_expanded), intent(in), target :: ref, ptcl
        ! re-set pointers
        reference => ref
        particle  => ptcl
    end subroutine ftexp_shsrch_reset_ptrs

    !> Cost function
    function ftexp_shsrch_cost( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real :: cost
        cost = -reference%corr_shifted(particle, -vec)
    end function ftexp_shsrch_cost
    
    !> Main search routine
    function ftexp_shsrch_minimize( prev_corr, prev_shift ) result( cxy )
        real, optional, intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3), maxshift, maxlim, lims(2,2)
        if( present(prev_shift) )then
            maxshift = maxval(abs(prev_shift))
            maxlim   = maxval(abs(ospec%limits))
            if( maxshift > maxlim )then
                ! re-specify the limits
                lims(:,1) = -maxshift
                lims(:,2) = maxshift
                call ospec%specify(opt_str, 2, ftol=TOL, gtol=TOL, limits=lims, nrestarts=nrestarts )
                ! set optimizer cost function
                call ospec%set_costfun(ftexp_shsrch_cost)
                ! generate optimizer object with the factory
                call ofac%new(ospec, nlopt)
            endif
            ospec%x = prev_shift
        else
            ospec%x = 0.
        endif
        call nlopt%minimize(ospec, cxy(1))
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
        use simple_rnd, only: ran3
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
        call img_ref%fwd_ft
        call ftexp_ref%new(img_ref, hp, lp)
        call ftexp_ptcl%new(img_ptcl, hp, lp)
        call ftexp_shsrch_init(ftexp_ref, ftexp_ptcl, lims, 'simplex', 3)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            call img_ref%shift([x,y,0.], imgout=img_ptcl)
            call ftexp_ptcl%new(img_ptcl, hp, lp)
            cxy = ftexp_shsrch_minimize()
            if( cxy(1) < 0.995 )then
                stop 'shift alignment failed; test_ftexp_shsrch; simple_ftexp_shsrch'
            endif
        end do
        write(*,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftexp_shsrch

end module simple_ftexp_shsrch
