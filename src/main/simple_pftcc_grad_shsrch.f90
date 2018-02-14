! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftcc_grad_shsrch
#include "simple_lib.f08"
use simple_opt_spec,          only: opt_spec
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_opt_factory,       only: opt_factory
use simple_optimizer,         only: optimizer
implicit none

public :: pftcc_grad_shsrch
private

type :: pftcc_grad_shsrch
    private
    type(opt_spec)                   :: ospec                  !< optimizer specification object
    class(optimizer), pointer        :: nlopt        =>null()  !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr    =>null()  !< pointer to pftcc object
    integer                          :: reference    = 0       !< reference pft
    integer                          :: particle     = 0       !< particle pft
    integer                          :: nrots        = 0       !< # rotations
    integer                          :: maxits       = 100     !< max # iterations
    logical                          :: shbarr       = .true.  !< shift barrier constraint or not
    integer                          :: cur_inpl_idx = 0       !< index of inplane angle for shift search
    real                             :: maxshift     = 0.      !< maximal shift
    integer                          :: max_evals    = 5       !< max # inplrot/shsrch cycles
  contains
    procedure          :: new         => grad_shsrch_new
    procedure          :: set_indices => grad_shsrch_set_indices
    procedure          :: minimize    => grad_shsrch_minimize
    procedure          :: kill        => grad_shsrch_kill
    procedure, private :: grad_shsrch_set_costfun
end type pftcc_grad_shsrch

contains

    !> Shift search constructor
    subroutine grad_shsrch_new( self, pftcc, lims, lims_init, shbarrier, maxits )
        class(pftcc_grad_shsrch),           intent(inout) :: self           !< instance
        class(polarft_corrcalc),    target, intent(in)    :: pftcc          !< correlator
        real,                               intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional,         intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional,         intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional,         intent(in)    :: maxits         !< maximum iterations
        type(opt_factory) :: opt_fact
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, ftol=1e-1, gtol=1e-3, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get # rotations
        self%nrots = pftcc%get_nrots()
        call self%grad_shsrch_set_costfun
    end subroutine grad_shsrch_new

    subroutine grad_shsrch_set_costfun( self )
        class(pftcc_grad_shsrch), intent(inout) :: self
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        self%ospec%opt_callback => grad_shsrch_optimize_angle_wrapper
    end subroutine grad_shsrch_set_costfun

    function grad_shsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        select type(self)
        class is (pftcc_grad_shsrch)
            cost = - self%pftcc_ptr%gencorr_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx)
        class default
            write (*,*) 'error in grad_shsrch_costfun: unknown type'
            stop
        end select
    end function grad_shsrch_costfun

    subroutine grad_shsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp)                :: corrs_grad(2)
        select type(self)
        class is (pftcc_grad_shsrch)
            call self%pftcc_ptr%gencorr_grad_only_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs_grad)
            grad = - corrs_grad
        class default
            write (*,*) 'error in grad_shsrch_gcostfun: unknown type'
            grad = 0.
            stop
        end select
    end subroutine grad_shsrch_gcostfun

    subroutine grad_shsrch_fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(dp)                :: corrs
        real(dp)                :: corrs_grad(2)
        select type(self)
        class is (pftcc_grad_shsrch)
            call self%pftcc_ptr%gencorr_grad_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
            f    = - corrs
            grad = - corrs_grad
        class default
            write (*,*) 'error in grad_shsrch_fdfcostfun: unknown type'
            f = 0.
            grad = 0.
            stop
        end select
    end subroutine grad_shsrch_fdfcostfun

    subroutine grad_shsrch_optimize_angle( self )
        class(pftcc_grad_shsrch), intent(inout) :: self
        real                                    :: corrs(self%nrots)
        integer                                 :: loc(1)
        call self%pftcc_ptr%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
        loc = maxloc(corrs)
        self%cur_inpl_idx = loc(1)
    end subroutine grad_shsrch_optimize_angle

    subroutine grad_shsrch_optimize_angle_wrapper( self )
        class(*), intent(inout) :: self
        select type(self)
        class is (pftcc_grad_shsrch)
            call grad_shsrch_optimize_angle(self)
        class default
            write (*,*) 'error in grad_shsrch_optimize_angle_wrapper: unknown type'
            stop
        end select
    end subroutine grad_shsrch_optimize_angle_wrapper

    !> shsrch_set_indices Set indicies for shift search
    !! \param ref reference
    !! \param ptcl particle index
    !! \param rot rotational index
    !! \param state current state
    !!
    subroutine grad_shsrch_set_indices( self, ref, ptcl )
        class(pftcc_grad_shsrch), intent(inout) :: self
        integer,                  intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine grad_shsrch_set_indices

    !> minimisation routine
    function grad_shsrch_minimize( self, irot ) result( cxy )
        class(pftcc_grad_shsrch), intent(inout) :: self
        integer,                  intent(out)   :: irot
        real    :: cost, corrs(self%nrots), cxy(3)
        real    :: lowest_cost, grad(2), f, lowest_cost_overall, lowest_shift(2)
        integer :: loc(1), i, irestart, lowest_rot, inpl_idx_zero_sh!, nevals, ngevals
        logical :: found_better
        found_better      = .false.
        call self%pftcc_ptr%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
        loc               = maxloc(corrs)
        self%cur_inpl_idx = loc(1)
        lowest_cost_overall = -corrs(self%cur_inpl_idx)
        ! shift search / in-plane rot update
        do i = 1,self%max_evals
            call self%nlopt%minimize(self%ospec, self, lowest_cost)
            call self%pftcc_ptr%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
            loc = maxloc(corrs)
            if( loc(1) == self%cur_inpl_idx ) exit
            self%cur_inpl_idx = loc(1)
        end do
        ! update best
        lowest_cost = -corrs(self%cur_inpl_idx)
        if( lowest_cost < lowest_cost_overall )then
            found_better        = .true.
            lowest_cost_overall = lowest_cost
            lowest_rot          = self%cur_inpl_idx
            lowest_shift        = self%ospec%x
        endif

        if( found_better )then
            irot    =   lowest_rot           ! in-plane index
            cxy(1)  = - lowest_cost_overall  ! correlation
            cxy(2:) =   lowest_shift         ! shift
            ! rotate the shift vector to the frame of reference
            cxy(2:) = matmul(cxy(2:), rotmat2d(self%pftcc_ptr%get_rot(irot)))
        else
            irot = 0 ! to communicate that a better solution was not found
        endif
    end function grad_shsrch_minimize

    subroutine grad_shsrch_kill( self )
        class(pftcc_grad_shsrch), intent(inout) :: self
        if( associated(self%nlopt) )then
            call self%nlopt%kill
            nullify(self%nlopt)
        end if
    end subroutine grad_shsrch_kill

    function grad_shsrch_get_nevals( self ) result( nevals )
        class(pftcc_grad_shsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function grad_shsrch_get_nevals

    subroutine grad_shsrch_get_peaks( self, peaks )
        class(pftcc_grad_shsrch), intent(inout) :: self
        real, allocatable,        intent(out)   :: peaks(:,:) !< output peak matrix
        allocate(peaks(1,2))
        peaks(1,:) = self%ospec%x
    end subroutine grad_shsrch_get_peaks

end module simple_pftcc_grad_shsrch
