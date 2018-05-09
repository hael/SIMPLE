module simple_pftcc_orisrch
#include "simple_lib.f08"
!$ use omp_lib
!$ use omp_lib_kinds
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_build,             only: build
use simple_ori,               only: ori
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad ! gradient-based angle and shift search
implicit none

public :: pftcc_orisrch
private

integer, parameter :: MAXITS = 60

type :: pftcc_orisrch
    private
    class(build),            pointer     :: bp          => null()  !< pointer to build
    class(polarft_corrcalc), pointer     :: pftcc_ptr   => null()  !< pointer to pftcc object
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_objs(:)    !< origin shift search objects, L-BFGS with gradient
    type(ori),               allocatable :: e_trials(:)            !< trial orientations (one per thread)
    real                                 :: cost_best              !< current best cost
    real                                 :: sh_best(2)  = 0        !< current best shift
    integer                              :: inpl_best   = 0        !< current best in-plane rotation
    integer                              :: particle    = 0        !< particle pft index
    integer                              :: nrots       = 0        !< # rotations
    integer                              :: nrestarts   = 3        !< simplex restarts (randomized bounds)
  contains
    procedure          :: new
    procedure          :: set_particle
    procedure          :: minimize
    procedure, private :: costfun
end type pftcc_orisrch

contains

    !> constructor
    subroutine new( self, pftcc, b, p, nrestarts )
        use simple_projector, only: projector
        use simple_params,    only: params
        class(pftcc_orisrch),            intent(inout) :: self      !< instance
        class(polarft_corrcalc), target, intent(in)    :: pftcc     !< correlator
        class(build),            target, intent(in)    :: b         !< builder
        class(params),                   intent(in)    :: p         !< parameters
        integer,               optional, intent(in)    :: nrestarts !< simplex restarts (randomized bounds)
        integer :: ithr, i
        real    :: lims_sh(2,2), lims_sh_init(2,2)
        ! kill allocatables
        if( allocated(self%e_trials) )then
            do i=1,size(self%e_trials)
                call self%e_trials(i)%kill
            end do
            deallocate(self%e_trials)
        endif
        if( allocated(self%grad_shsrch_objs) )then
            do i=1,size(self%grad_shsrch_objs)
                call self%grad_shsrch_objs(i)%kill
            end do
            deallocate(self%grad_shsrch_objs)
        endif
        ! set nrestarts
        self%nrestarts = 3
        if( present(nrestarts) ) self%nrestarts = nrestarts
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! set pointer to build
        self%bp => b
        ! get # rotations
        self%nrots = pftcc%get_nrots()
        ! create trial orientations and in-plane search objects (one per thread)
        allocate(self%grad_shsrch_objs(nthr_glob), self%e_trials(nthr_glob))
        lims_sh(:,1)      = -p%trs
        lims_sh(:,2)      =  p%trs
        lims_sh_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_sh_init(:,2) =  SHC_INPL_TRSHWDTH
        do ithr=1,nthr_glob
            call self%e_trials(ithr)%new_ori_clean
            call self%grad_shsrch_objs(ithr)%new(self%pftcc_ptr, lims_sh, lims_init=lims_sh_init,&
                &shbarrier=p%shbarrier, maxits=MAXITS)
        end do
    end subroutine new

    !> set particle index for search
    subroutine set_particle( self, ptcl )
        class(pftcc_orisrch), intent(inout) :: self
        integer,              intent(in)    :: ptcl
        self%particle  = ptcl
    end subroutine set_particle

    !> minimisation
    function minimize( self, o_inout, angerr_deg, irot ) result( cxy )
        use simple_ori,         only: ori
        use simple_opt_factory, only: opt_factory
        use simple_optimizer,   only: optimizer
        use simple_opt_spec,    only: opt_spec
        class(pftcc_orisrch), intent(inout) :: self
        class(ori),           intent(inout) :: o_inout
        real,                 intent(inout) :: angerr_deg
        integer,              intent(out)   :: irot
        type pftcc_ref
            complex, allocatable :: pft_ref_even(:,:), pft_ref_odd(:,:)
        end type pftcc_ref
        type(pftcc_ref),  allocatable :: pftcc_refs(:)
        real,             allocatable :: cxy(:)
        class(optimizer), pointer     :: nlopt
        type(opt_factory)             :: opt_fact
        type(opt_spec)                :: ospec
        real    :: cost, cost_init, lims(2,2), lims_init(2,2)
        integer :: ithr
        ! copy nthr_glob pftcc references so we can put them back after minimization is done
        allocate(pftcc_refs(nthr_glob), cxy(3))
        do ithr=1,nthr_glob
            pftcc_refs(ithr)%pft_ref_even = self%pftcc_ptr%get_ref_pft(ithr, iseven=.true.)
            pftcc_refs(ithr)%pft_ref_odd  = self%pftcc_ptr%get_ref_pft(ithr, iseven=.false.)
        end do

        ! PREPARATION
        ! associate costfun
        ospec%costfun => costfun_wrapper
        ! initialize with input projection direction
        ospec%x(1) = o_inout%e1get()
        ospec%x(2) = o_inout%e2get()
        ! set limits
        lims = 0.
        lims(1,2) = 359.99
        lims(2,2) = 180.00
        lims_init(1,1) = max(ospec%x(1) - angerr_deg, 0.)
        lims_init(1,2) = min(ospec%x(1) + angerr_deg, lims(1,2))
        lims_init(2,1) = max(ospec%x(2) - angerr_deg, 0.)
        lims_init(2,2) = min(ospec%x(2) + angerr_deg, lims(2,2))
        ! make optimizer spec
        call ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4, limits=lims,&
            &limits_init=lims_init, nrestarts=self%nrestarts, maxits=MAXITS)
        ! generate the optimizer object
        call opt_fact%new(ospec, nlopt)

        ! MINIMISATION
        ospec%nevals   = 0
        self%cost_best = self%costfun(ospec%x, ospec%ndim)
        cost_init      = self%cost_best ! for later comparison
        call nlopt%minimize(ospec, self, cost)

        ! OUTPUT
        if( cost <= cost_init )then
            ! set output
            cxy(1)  = -cost ! correlation
            ! set Euler
            call o_inout%set_euler([ospec%x(1),ospec%x(2),360. - self%pftcc_ptr%get_rot(self%inpl_best)])
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
        real              :: cost
        integer           :: ithr, irot
        real, allocatable :: cxy(:)
        ! thread-safe extraction of projection
        ithr = omp_get_thread_num() + 1
        call self%e_trials(ithr)%set_euler([vec(1),vec(2),0.])
        if( self%pftcc_ptr%ptcl_iseven(self%particle) )then
            call self%bp%vol%fproject_polar(ithr, self%e_trials(ithr), self%pftcc_ptr, iseven=.true.)
        else
            call self%bp%vol_odd%fproject_polar(ithr, self%e_trials(ithr), self%pftcc_ptr, iseven=.false.)
        endif
        ! in-plane search with L-BFGS-B and callback for exhaustive in-plane rotation search
        call self%grad_shsrch_objs(ithr)%set_indices(ithr, self%particle)
        cxy  = self%grad_shsrch_objs(ithr)%minimize(irot=irot)
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

end module simple_pftcc_orisrch
