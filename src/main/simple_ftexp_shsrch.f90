! shift search with L-BFGS-B using expanded Fourier transforms (used in motion_correct)
module simple_ftexp_shsrch
include 'simple_lib.f08'
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded

implicit none

public :: ftexp_shsrch, test_ftexp_shsrch
private
#include "simple_local_flags.inc"

real,    parameter :: TOL    = 1e-4 !< tolerance parameter
integer, parameter :: MAXITS = 30   !< maximum number of iterations

type :: ftexp_shsrch
    type(opt_spec)              :: ospec                     !< optimizer specification object
    class(optimizer),   pointer :: nlopt     => null()       !< pointer to nonlinear optimizer
    class(ft_expanded), pointer :: reference => null()       !< reference ft_exp
    class(ft_expanded), pointer :: particle  => null()       !< particle ft_exp
    real                        :: maxHWshift         = 0.   !< maximum half-width of shift
    real                        :: motion_correctftol = 1e-4 !< function error tolerance
    real                        :: motion_correctgtol = 1e-4 !< gradient error tolerance
contains
    procedure :: new      => ftexp_shsrch_new
    procedure :: set_ptrs => ftexp_shsrch_set_ptrs
    procedure :: minimize => ftexp_shsrch_minimize
    procedure :: kill     => ftexp_shsrch_kill
end type ftexp_shsrch

contains

    !> Initialise  ftexp_shsrch
    subroutine ftexp_shsrch_new( self, ref, ptcl, trs, motion_correct_ftol, motion_correct_gtol )
        use simple_opt_factory, only: opt_factory
        class(ftexp_shsrch),        intent(inout) :: self
        class(ft_expanded), target, intent(in)    :: ref, ptcl
        real,                       intent(in)    :: trs
        real,             optional, intent(in)    :: motion_correct_ftol, motion_correct_gtol
        type(opt_factory) :: ofac
        real              :: lims(2,2)
        self%reference  => ref
        self%particle   => ptcl
        self%maxHWshift =  trs
        if( present(motion_correct_ftol) )then
            self%motion_correctftol = motion_correct_ftol
        else
            self%motion_correctftol = TOL
        end if
        if( present(motion_correct_gtol) )then
            self%motion_correctgtol = motion_correct_gtol
        else
            self%motion_correctgtol = TOL
        end if
        lims(1,1) = - self%maxHWshift
        lims(1,2) =   self%maxHWshift
        lims(2,1) = - self%maxHWshift
        lims(2,2) =   self%maxHWshift
        call self%ospec%specify('lbfgsb', 2, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=lims)
        call self%ospec%set_costfun_8(ftexp_shsrch_cost_8)
        call self%ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
        call self%ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%nlopt) )then
            call self%nlopt%kill
            deallocate(self%nlopt)
        end if
        call ofac%new(self%ospec, self%nlopt)
    end subroutine ftexp_shsrch_new

    subroutine ftexp_shsrch_set_ptrs( self, ref, ptcl )
        class(ftexp_shsrch),        intent(inout) :: self
        class(ft_expanded), target, intent(in)    :: ref, ptcl
        self%reference => ref
        self%particle  => ptcl
    end subroutine ftexp_shsrch_set_ptrs

    !> Cost function, double precision
    function ftexp_shsrch_cost_8( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8) :: cost
        select type(self)
            class is (ftexp_shsrch)
                cost = -self%reference%corr_shifted_8(self%particle, -vec)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_cost_8')
        end select
    end function ftexp_shsrch_cost_8

    !> Gradient function, double precision
    subroutine ftexp_shsrch_gcost_8( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                call self%reference%corr_gshifted_8(self%particle, -vec, grad)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_gcost_8')
        end select
    end subroutine ftexp_shsrch_gcost_8

    !> Gradient & cost function, double precision
    subroutine ftexp_shsrch_fdfcost_8( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                call self%reference%corr_fdfshifted_8(self%particle, -vec, f, grad)
                f = f * (-1.0_dp)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_fdfcost_8')
        end select
    end subroutine ftexp_shsrch_fdfcost_8

    !> Main search routine
    function ftexp_shsrch_minimize( self, prev_corr, prev_shift ) result( cxy )
        class(ftexp_shsrch), intent(inout) :: self
        real, optional,      intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3)
        self%ospec%limits(1,1) = - self%maxHWshift
        self%ospec%limits(1,2) =   self%maxHWshift
        self%ospec%limits(2,1) = - self%maxHWshift
        self%ospec%limits(2,2) =   self%maxHWshift
        if( present(prev_shift) )then
            self%ospec%limits(1,:) = self%ospec%limits(1,:) + prev_shift(1)
            self%ospec%limits(2,:) = self%ospec%limits(2,:) + prev_shift(2)
        endif
        if( present(prev_shift) ) self%ospec%x = prev_shift
        ! set initial solution to previous shift
        call self%nlopt%minimize(self%ospec, self, cxy(1))
        call self%reference%corr_normalize(self%particle, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = self%ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= TOL )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
    end function ftexp_shsrch_minimize

    subroutine ftexp_shsrch_kill( self )
        class(ftexp_shsrch), intent(inout) :: self
        call self%ospec%kill
        if( associated(self%nlopt) )then
            call self%nlopt%kill
            nullify(self%nlopt)
        end if
        if( associated(self%reference) ) self%reference => null()
        if( associated(self%particle)  ) self%particle  => null()
    end subroutine ftexp_shsrch_kill

    subroutine test_ftexp_shsrch
        use simple_image, only: image
        type(image)       :: img_ref, img_ptcl
        type(ft_expanded) :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch)  :: ftexp_srch
        real, parameter   :: TRS=5.0, lp=6., hp=100.
        real              :: cxy(3), x, y, lims(2,2)
        integer           :: i
        lims(:,1) = -TRS
        lims(:,2) = TRS
        call img_ref%new([100,100,1], 2.)
        call img_ptcl%new([100,100,1], 2.)
        call img_ref%square(20)
        call img_ref%fft()
        call ftexp_ref%new(img_ref, hp, lp, .true.)
        call ftexp_ptcl%new(img_ptcl, hp, lp, .true.)
        call ftexp_srch%new(ftexp_ref, ftexp_ptcl, trs)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            img_ptcl = img_ref
            call img_ptcl%shift([x,y,0.])
            call ftexp_ptcl%new(img_ptcl, hp, lp, .true.)
            cxy = ftexp_srch%minimize()
            if( cxy(1) < 0.995 )then
                THROW_HARD('shift alignment failed; test_ftexp_shsrch')
            endif
        end do
        write(logfhandle,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftexp_shsrch

end module simple_ftexp_shsrch
