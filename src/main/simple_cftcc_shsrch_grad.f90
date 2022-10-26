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

logical,  parameter :: perform_rndstart = .true.
real(dp), parameter :: init_range       = 2.0_dp  ! range for random initialization (negative to positive)
integer,  parameter :: coarse_num_steps = 5       ! no. of coarse search steps in x AND y (hence real no. is its square)

type :: cftcc_shsrch_grad
    private
    type(opt_spec)            :: ospec                  !< optimizer specification object
    class(optimizer), pointer :: nlopt        =>null()  !< optimizer object
    integer                   :: particle     = 0       !< particle index
    integer                   :: nrots        = 0       !< # rotations
    integer                   :: maxits       = 100     !< max # iterations
    logical                   :: shbarr       = .true.  !< shift barrier constraint or not
    integer                   :: max_evals    = 5       !< max # inplrot/shsrch cycles
    real                      :: max_shift    = 0.      !< maximal shift
    logical                   :: coarse_init  = .false. !< whether to perform an intial coarse search over the range
contains
    procedure :: new         => grad_carshsrch_new
    procedure :: set_pind    => grad_carshsrch_set_pind
    procedure :: minimize    => grad_carshsrch_minimize
    procedure :: kill        => grad_carshsrch_kill
    procedure :: coarse_search
end type cftcc_shsrch_grad

contains

    subroutine grad_carshsrch_new( self, lims, lims_init, shbarrier, maxits, coarse_init )
        use simple_opt_factory, only: opt_factory
        class(cftcc_shsrch_grad),   intent(inout) :: self           !< instance
        real,                       intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional, intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional, intent(in)    :: maxits         !< maximum iterations
        logical,          optional, intent(in)    :: coarse_init    !< coarse inital search
        type(opt_factory) :: opt_fact
        call self%kill
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%coarse_init = .false.
        if( present(coarse_init) ) self%coarse_init = coarse_init
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, ftol=1e-1, gtol=1e-3, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! set costfun pointers
        self%ospec%costfun_8    => grad_carshsrch_costfun
        self%ospec%gcostfun_8   => grad_carshsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_carshsrch_fdfcostfun
    end subroutine grad_carshsrch_new

    !> set indicies for shift search
    subroutine grad_carshsrch_set_pind( self, ptcl )
        class(cftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(in)    :: ptcl
        self%particle  = ptcl
    end subroutine grad_carshsrch_set_pind

    function grad_carshsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        select type(self)
            class is (cftcc_shsrch_grad)
                cost = - cftcc_glob%corr_shifted(self%particle, real(vec))
            class default
                THROW_HARD('error in grad_shsrch_costfun: unknown type; grad_shsrch_costfun')
        end select
    end function grad_carshsrch_costfun

    subroutine grad_carshsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(sp)                :: corrs_grad(2)
        real(sp)                :: cost
        grad = 0.
        select type(self)
            class is (cftcc_shsrch_grad)
                cost = - cftcc_glob%corr_shifted(self%particle, real(vec), corrs_grad)
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_gcostfun: unknown type; grad_shsrch_gcostfun')
        end select
    end subroutine grad_carshsrch_gcostfun

    subroutine grad_carshsrch_fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(sp)                :: corrs_grad(2)
        f    = 0.
        grad = 0.
        select type(self)
            class is (cftcc_shsrch_grad)
                f    = - cftcc_glob%corr_shifted(self%particle, real(vec), corrs_grad)
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine grad_carshsrch_fdfcostfun

    !> minimisation routine
    function grad_carshsrch_minimize( self ) result( cxy )
        class(cftcc_shsrch_grad), intent(inout) :: self
        real     :: cxy(3), lowest_cost
        real(dp) :: init_xy(2), coarse_cost
        self%ospec%x      = [0.,0.]
        self%ospec%x_8    = [0.d0,0.d0]
        if( perform_rndstart )then
            init_xy(1)     = 2.*(ran3()-0.5) * init_range
            init_xy(2)     = 2.*(ran3()-0.5) * init_range
            self%ospec%x_8 = init_xy
            self%ospec%x   = real(init_xy)
        endif
        if( self%coarse_init )then
            call self%coarse_search(coarse_cost, init_xy)
            self%ospec%x_8      = init_xy
            self%ospec%x        = real(init_xy)
        end if
        ! shift search
        call self%nlopt%minimize(self%ospec, self, lowest_cost)
        cxy(1)  = - real(lowest_cost)  ! correlation
        cxy(2:) =   self%ospec%x       ! shift
    end function grad_carshsrch_minimize

    subroutine coarse_search(self, lowest_cost, init_xy)
        class(cftcc_shsrch_grad), intent(inout) :: self
        real(dp),                 intent(out)   :: lowest_cost, init_xy(2)
        real(dp) :: x, y, cost, stepx, stepy
        integer  :: ix, iy
        lowest_cost = huge(lowest_cost)
        init_xy     = 0.d0
        if (coarse_num_steps .le. 1) return
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            x = self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,coarse_num_steps
                y = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                cost = -cftcc_glob%corr_shifted(self%particle, real([x,y]))
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    init_xy     = [x,y]
                end if
            enddo
        enddo
    end subroutine coarse_search

    subroutine grad_carshsrch_kill( self )
        class(cftcc_shsrch_grad), intent(inout) :: self
        if( associated(self%nlopt) )then
            call self%ospec%kill
            call self%nlopt%kill
            nullify(self%nlopt)
        end if
    end subroutine grad_carshsrch_kill

end module simple_cftcc_shsrch_grad
