module simple_pftcc_orisrch
#include "simple_lib.f08"
!$ use omp_lib
!$ use omp_lib_kinds
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_build,             only: build
use simple_ori,               only: ori
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad ! gradient-based angle and shift search
use simple_opt_spec,          only: opt_spec
use simple_optimizer,         only: optimizer
implicit none

public :: pftcc_orisrch
private

integer, parameter :: MAXITS = 60

type :: pftcc_orisrch
    private
    class(build),            pointer     :: bp          => null()  !< pointer to build
    class(polarft_corrcalc), pointer     :: pftcc_ptr   => null()  !< pointer to pftcc object
    class(optimizer),        pointer     :: nlopt                  !< optimizer object
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_obj        !< origin shift search objects, L-BFGS with gradient
    type(ori)                            :: e_trial                !< trial orientation
    type(opt_spec)                       :: ospec                  !< optimizer specification object
    real                                 :: cost_best              !< current best cost
    real                                 :: sh_best(2)  = 0        !< current best shift
    integer                              :: inpl_best   = 0        !< current best in-plane rotation
    integer                              :: particle    = 0        !< particle pft index
    integer                              :: nrots       = 0        !< # rotations
    integer                              :: nrestarts   = 3        !< simplex restarts (randomized bounds)
    logical                              :: exists=.false.
  contains
    procedure          :: new
    procedure          :: set_particle
    procedure          :: minimize
    procedure, private :: costfun
    procedure          :: kill
end type pftcc_orisrch

contains

    !> constructor
    subroutine new( self, pftcc, b, p, nrestarts )
        use simple_projector,   only: projector
        use simple_params,      only: params
        use simple_opt_factory, only: opt_factory
        class(pftcc_orisrch),            intent(inout) :: self      !< instance
        class(polarft_corrcalc), target, intent(in)    :: pftcc     !< correlator
        class(build),            target, intent(in)    :: b         !< builder
        class(params),                   intent(in)    :: p         !< parameters
        integer,               optional, intent(in)    :: nrestarts !< simplex restarts (randomized bounds)
        type(opt_factory) :: opt_fact
        integer :: ithr, i
        real    :: lims(2,2), lims_sh(2,2), lims_sh_init(2,2)
        ! set nrestarts
        call self%kill
        allocate(self%grad_shsrch_obj)
        self%nrestarts = 3
        if( present(nrestarts) ) self%nrestarts = nrestarts
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! set pointer to build
        self%bp => b
        ! get # rotations
        self%nrots = pftcc%get_nrots()
        ! create trial orientation
        call self%e_trial%new()
        ! create in-plane search object
        lims_sh(:,1)      = -p%trs
        lims_sh(:,2)      =  p%trs
        lims_sh_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_sh_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%grad_shsrch_obj%new(self%pftcc_ptr, lims_sh, lims_init=lims_sh_init,&
            &shbarrier=p%shbarrier, maxits=MAXITS)
        ! make optimizer spec
        lims = 0.
        lims(1,2) = 359.99
        lims(2,2) = 180.00
        call self%ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4, limits=lims,&
            &limits_init=lims, nrestarts=self%nrestarts, maxits=MAXITS)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! associate costfun
        self%ospec%costfun => costfun_wrapper
    end subroutine new

    !> set particle index for search
    subroutine set_particle( self, ptcl )
        class(pftcc_orisrch), intent(inout) :: self
        integer,              intent(in)    :: ptcl
        self%particle  = ptcl
    end subroutine set_particle

    !> minimisation
    function minimize( self, o_inout, angerr_deg, irot ) result( cxy )
        class(pftcc_orisrch), intent(inout) :: self
        class(ori),           intent(inout) :: o_inout
        real,                 intent(inout) :: angerr_deg
        integer,              intent(out)   :: irot
        type pftcc_ref
            complex, allocatable :: pft_ref_even(:,:), pft_ref_odd(:,:)
        end type pftcc_ref
        type(pftcc_ref), allocatable :: pftcc_refs(:)
        real    :: cost, cost_init, lims(2,2), lims_init(2,2), cxy(3)
        integer :: ithr
        ! copy nthr_glob pftcc references so we can put them back after minimization is done
        allocate(pftcc_refs(nthr_glob))
        do ithr=1,nthr_glob
            pftcc_refs(ithr)%pft_ref_even = self%pftcc_ptr%get_ref_pft(ithr, iseven=.true.)
            pftcc_refs(ithr)%pft_ref_odd  = self%pftcc_ptr%get_ref_pft(ithr, iseven=.false.)
        end do

        ! PREPARATION
        ! initialize with input projection direction
        self%ospec%x(1) = o_inout%e1get()
        self%ospec%x(2) = o_inout%e2get()
        ! set limits for simplex initialisation
        self%ospec%limits_init(1,1) = max(self%ospec%x(1) - angerr_deg, 0.)
        self%ospec%limits_init(1,2) = min(self%ospec%x(1) + angerr_deg, lims(1,2))
        self%ospec%limits_init(2,1) = max(self%ospec%x(2) - angerr_deg, 0.)
        self%ospec%limits_init(2,2) = min(self%ospec%x(2) + angerr_deg, lims(2,2))

        ! MINIMISATION
        self%ospec%nevals   = 0
        self%cost_best = self%costfun(self%ospec%x, self%ospec%ndim)
        cost_init      = self%cost_best ! for later comparison
        call self%nlopt%minimize(self%ospec, self, cost)

        ! OUTPUT
        if( cost <= cost_init )then
            ! set output
            cxy(1)  = -cost ! correlation
            ! set Euler
            call o_inout%set_euler([self%ospec%x(1),self%ospec%x(2),360. - self%pftcc_ptr%get_rot(self%inpl_best)])
            ! set shift (shift vector already rotated to the frame of reference in pftcc_shsrch_grad)
            cxy(2:) = self%sh_best
            call o_inout%set_shift(self%sh_best)
        else
            irot    = 0
            cxy(1)  = -cost_init ! correlation
            cxy(2:) = 0.
        endif

        ! put back references & deallocate
        do ithr=1,nthr_glob
            call self%pftcc_ptr%set_ref_pft(ithr, pftcc_refs(ithr)%pft_ref_even, iseven=.true.)
            call self%pftcc_ptr%set_ref_pft(ithr, pftcc_refs(ithr)%pft_ref_odd,  iseven=.false.)
            deallocate(pftcc_refs(ithr)%pft_ref_even, pftcc_refs(ithr)%pft_ref_odd)
        end do
        deallocate(pftcc_refs)
    end function minimize

    !> wrapper for cost function (gcc7+)
    function costfun_wrapper(self, vec, D) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        select type (self)
            class is (pftcc_orisrch)
                cost = self%costfun(vec, D)
            class default
                write (*,*) 'error in simple_pftcc_orisrch, costfun_wrapper: unknown type'
                stop
        end select
    end function costfun_wrapper

    !> cost function
    function costfun( self, vec, D ) result( cost )
        class(pftcc_orisrch), intent(inout) :: self
        integer,              intent(in)    :: D
        real,                 intent(in)    :: vec(D)
        real     :: cost
        integer  :: ithr, irot
        real     :: cxy(3)
        ! set Euler angle
        call self%e_trial%set_euler([vec(1),vec(2),0.])
        ! thread-safe extraction of projection (because pftcc is an OpenMP shared data structure)
        ithr = omp_get_thread_num() + 1
        if( self%pftcc_ptr%ptcl_iseven(self%particle) )then
            call self%bp%vol%fproject_polar(ithr, self%e_trial, self%pftcc_ptr, iseven=.true.)
        else
            call self%bp%vol_odd%fproject_polar(ithr, self%e_trial, self%pftcc_ptr, iseven=.false.)
        endif
        ! in-plane search with L-BFGS-B and callback for exhaustive in-plane rotation search
        call self%grad_shsrch_obj%set_indices(ithr, self%particle)
        cxy = self%grad_shsrch_obj%minimize(irot=irot)
        if( irot == 0 )then ! no better solution could be identified
            cost = self%cost_best
            return
        endif
        ! cost function
        cost = -cxy(1)
        if( cost < self%cost_best )then
            self%cost_best = cost
            self%sh_best   = cxy(2:3)
            self%inpl_best = irot
        endif
    end function costfun

    !> \brief  is a destructor
    subroutine kill( self )
        class(pftcc_orisrch), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate( self%grad_shsrch_obj )
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pftcc_orisrch
