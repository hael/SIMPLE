! Cartesian-based shift search, gradient-based lbfgsb
module simple_cftcc_shsrch_grad
include 'simple_lib.f08'
use simple_opt_spec,         only: opt_spec
use simple_cartft_corrcalc,  only: cftcc_glob
use simple_optimizer,        only: optimizer
implicit none

public :: cftcc_shsrch_grad
private
#include "simple_local_flags.inc"

real(dp), parameter :: init_range = 2.0_dp  ! range for random initialization (negative to positive)

type :: cftcc_shsrch_grad
    private
    type(opt_spec)            :: ospec                 !< optimizer specification object
    class(optimizer), pointer :: nlopt       =>null()  !< optimizer object
    integer                   :: particle    = 0       !< particle index
    integer                   :: maxits      = 100     !< max # iterations
    logical                   :: shbarr      = .true.  !< shift barrier constraint or not
    real                      :: max_shift   = 0.      !< maximal shift
contains
    procedure :: new      => c_grad_shsrch_new
    procedure :: set_pind => c_grad_shsrch_set_pind
    procedure :: minimize => c_grad_shsrch_minimize
    procedure :: kill     => c_grad_shsrch_kill
end type cftcc_shsrch_grad

contains

    subroutine c_grad_shsrch_new( self, lims, lims_init, shbarrier, maxits )
        use simple_opt_factory, only: opt_factory
        class(cftcc_shsrch_grad),   intent(inout) :: self           !< instance
        real,                       intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional, intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional, intent(in)    :: maxits         !< maximum iterations
        type(opt_factory) :: opt_fact
        call self%kill
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! set costfun pointers
        self%ospec%costfun_8    => c_grad_shsrch_costfun
        self%ospec%gcostfun_8   => c_grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => c_grad_shsrch_fdfcostfun
    end subroutine c_grad_shsrch_new

    !> set indicies for shift search
    subroutine c_grad_shsrch_set_pind( self, ptcl )
        class(cftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(in)    :: ptcl
        self%particle  = ptcl
    end subroutine c_grad_shsrch_set_pind

    function c_grad_shsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        real(sp)                :: cost_sp
        select type(self)
            class is (cftcc_shsrch_grad)
                call cftcc_glob%corr_shifted(self%particle, real(vec), cost_sp)
                cost = real(-cost_sp,kind=dp)
                self%ospec%nevals = self%ospec%nevals + 1
            class default
                THROW_HARD('error in grad_shsrch_costfun: unknown type; grad_shsrch_costfun')
        end select
    end function c_grad_shsrch_costfun

    subroutine c_grad_shsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(sp)                :: corrs_grad(2)
        real(sp)                :: cost
        grad = 0.
        select type(self)
            class is (cftcc_shsrch_grad)
                call cftcc_glob%corr_shifted(self%particle, real(vec), cost, corrs_grad)
                cost = - cost
                grad = - corrs_grad
                self%ospec%nevals  = self%ospec%nevals  + 1
                self%ospec%ngevals = self%ospec%ngevals + 1
            class default
                THROW_HARD('error in grad_shsrch_gcostfun: unknown type; grad_shsrch_gcostfun')
        end select
    end subroutine c_grad_shsrch_gcostfun

    subroutine c_grad_shsrch_fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(sp)                :: corrs_grad(2), f_sp
        f    = 0.
        grad = 0.
        select type(self)
            class is (cftcc_shsrch_grad)
                call cftcc_glob%corr_shifted(self%particle, real(vec), f_sp, corrs_grad)
                f    = real(-f_sp,kind=dp)
                grad = - corrs_grad
                self%ospec%nevals  = self%ospec%nevals  + 1
                self%ospec%ngevals = self%ospec%ngevals + 1
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine c_grad_shsrch_fdfcostfun

    !> minimisation routine
    function c_grad_shsrch_minimize( self, nevals, shvec ) result( cxy )
        class(cftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(inout) :: nevals(2)
        real, optional,           intent(in)    :: shvec(2)
        real     :: cxy(3), lowest_cost
        real(dp) :: init_xy(2)
        self%ospec%x       = [0.,0.]
        self%ospec%x_8     = [0.d0,0.d0]
        self%ospec%nevals  = 0
        self%ospec%ngevals = 0
        if( present(shvec) )then
            self%ospec%x_8 = shvec
            self%ospec%x   = real(shvec)
        else
            init_xy(1)     = 2.*(ran3()-0.5) * init_range
            init_xy(2)     = 2.*(ran3()-0.5) * init_range
            self%ospec%x_8 = init_xy
            self%ospec%x   = real(init_xy)
        endif
        ! shift search
        call self%nlopt%minimize(self%ospec, self, lowest_cost)
        cxy(1)    = - real(lowest_cost) ! correlation
        cxy(2:)   =   self%ospec%x      ! shift
        nevals(1) = self%ospec%nevals
        nevals(2) = self%ospec%ngevals
    end function c_grad_shsrch_minimize

    subroutine c_grad_shsrch_kill( self )
        class(cftcc_shsrch_grad), intent(inout) :: self
        if( associated(self%nlopt) )then
            call self%ospec%kill
            call self%nlopt%kill
            nullify(self%nlopt)
        end if
    end subroutine c_grad_shsrch_kill

end module simple_cftcc_shsrch_grad
