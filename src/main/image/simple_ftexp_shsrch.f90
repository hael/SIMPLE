!@descr: shift search with L-BFGS-B using expanded Fourier transforms (used in motion_correct)
module simple_ftexp_shsrch
use simple_core_module_api
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded
implicit none

public :: ftexp_shsrch
private
#include "simple_local_flags.inc"

real,     parameter :: TOL    = 1e-4    !< tolerance parameter
integer,  parameter :: MAXITS = 30      !< maximum number of iterations

type :: ftexp_shsrch
    private
    type(opt_spec), public      :: ospec                     !< optimizer specification object
    class(optimizer),   pointer :: opt_obj      => null()    !< pointer to nonlinear optimizer
    class(ft_expanded), pointer :: reference    => null()    !< reference ft_exp
    class(ft_expanded), pointer :: particle     => null()    !< particle ft_exp
    real(dp)                    :: denominator        = 0.d0
    real                        :: maxHWshift         = 0.   !< maximum half-width of shift
    real                        :: motion_correctftol = 1e-4 !< function error tolerance
    real                        :: motion_correctgtol = 1e-4 !< gradient error tolerance
    real                        :: shsrch_tol         = TOL
    logical                     :: existence = .false.
  contains
    procedure          :: new            => ftexp_shsrch_new
    procedure          :: minimize       => ftexp_shsrch_minimize
    procedure          :: corr_shifted_8 => ftexp_shsrch_corr_shifted_8
    procedure          :: kill           => ftexp_shsrch_kill
    procedure          :: set_shsrch_tol
    procedure          :: set_factr_pgtol
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
        real              :: opt_lims(2,2)
        call self%kill()
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
        opt_lims(1,1) = - self%maxHWshift
        opt_lims(1,2) =   self%maxHWshift
        opt_lims(2,1) = - self%maxHWshift
        opt_lims(2,2) =   self%maxHWshift
        call self%ospec%specify('lbfgsb', 2, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=opt_lims)
        call self%ospec%set_costfun_8(ftexp_shsrch_cost_8)
        call self%ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
        call self%ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%opt_obj) )then
            call self%opt_obj%kill
            deallocate(self%opt_obj) ! because this extended type is allocated by opt_factory
            nullify(self%opt_obj)
        end if
        call ofac%new(self%ospec, self%opt_obj)
        self%existence = .true.
    end subroutine ftexp_shsrch_new

    !> Main search routine
    function ftexp_shsrch_minimize( self, prev_corr, prev_shift ) result( cxy )
        class(ftexp_shsrch), intent(inout) :: self
        real, optional,      intent(in)    :: prev_corr, prev_shift(2)
        real :: cxy(3)
        self%ospec%limits(1,1) = - self%maxHWshift
        self%ospec%limits(1,2) =   self%maxHWshift
        self%ospec%limits(2,1) = - self%maxHWshift
        self%ospec%limits(2,2) =   self%maxHWshift
        if( present(prev_shift) )then
            self%ospec%limits(1,:) = self%ospec%limits(1,:) + prev_shift(1)
            self%ospec%limits(2,:) = self%ospec%limits(2,:) + prev_shift(2)
        endif
        if( present(prev_shift) ) then
            self%ospec%x = prev_shift
        else
            self%ospec%x   = 0.
        end if
        self%ospec%x_8 = real(self%ospec%x,dp)
        ! self%kind_shift = self%reference%get_kind_shift()
        ! call self%set_dims_and_alloc()
        ! call self%calc_tmp_cmat12()
        ! the temp matrix is built on the particle only!
        call self%particle%alloc_and_calc_tmp_cmat12(self%reference, self%denominator )
        ! set initial solution to previous shift
        call self%opt_obj%minimize(self%ospec, self, cxy(1))
        call self%reference%corr_normalize(self%particle, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = self%ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= self%shsrch_tol )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
    end function ftexp_shsrch_minimize

    subroutine ftexp_shsrch_kill( self )
        class(ftexp_shsrch), intent(inout) :: self
        if ( self%existence ) then
            call self%ospec%kill
            if( associated( self%opt_obj ) )then
                call self%opt_obj%kill
                deallocate(self%opt_obj)
                nullify(self%opt_obj)
            end if
            if( associated(self%reference) ) self%reference => null()
            if( associated(self%particle)  )then
                call self%particle%dealloc_tmp_cmat12
                self%particle  => null()
            endif
            self%existence = .false.
        end if
    end subroutine ftexp_shsrch_kill

    pure subroutine set_shsrch_tol( self, shsrch_tol )
        class(ftexp_shsrch), intent(inout) :: self
        real,                intent(in)    :: shsrch_tol
        self%shsrch_tol = shsrch_tol
    end subroutine set_shsrch_tol

    pure subroutine set_factr_pgtol( self, factr, pgtol )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: factr, pgtol
        self%ospec%factr = factr
        self%ospec%pgtol = pgtol
    end subroutine set_factr_pgtol

    ! Correlation
    real(dp) function ftexp_shsrch_corr_shifted_8( self, shvec )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        call self%particle%alloc_and_calc_tmp_cmat12(self%reference, self%denominator)
        ftexp_shsrch_corr_shifted_8 = self%particle%corr_shifted_cost_8(shvec, self%denominator)
        call self%particle%corr_normalize(self%reference, ftexp_shsrch_corr_shifted_8)
    end function ftexp_shsrch_corr_shifted_8

    !> Cost function, double precision
    function ftexp_shsrch_cost_8( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8) :: cost
        select type(self)
            class is (ftexp_shsrch)
                cost = -ftexp_shsrch_corr_shifted_8(self, -vec)
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
                call self%particle%corr_gshifted_cost_8( -vec, self%denominator, grad )
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
                call self%particle%corr_fdfshifted_cost_8( -vec, self%denominator, f, grad )
                f = f * (-1.0_dp)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_fdfcost_8')
        end select
    end subroutine ftexp_shsrch_fdfcost_8

end module simple_ftexp_shsrch
