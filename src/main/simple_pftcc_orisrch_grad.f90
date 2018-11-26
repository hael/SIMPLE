module simple_pftcc_orisrch_grad
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_polarft_corrcalc,  only: pftcc_glob
use simple_builder,           only: build_glob
use simple_optimizer,         only: optimizer
use simple_opt_spec,          only: opt_spec
implicit none

public :: pftcc_orisrch_grad
private

type :: pftcc_orisrch_grad
    private
    class(optimizer), pointer :: nlopt             !< optimizer object
    type(opt_spec)            :: ospec             !< optimizer spec
    integer                   :: particle = 0      !< particle pft index
    integer                   :: nrots    = 0      !< # rotations
    logical                   :: exists   =.false.
  contains
    procedure :: new
    procedure :: set_particle
    procedure :: minimize
    procedure :: kill
end type pftcc_orisrch_grad

contains

    !> constructor
    subroutine new( self )
        use simple_projector,   only: projector
        use simple_opt_factory, only: opt_factory
        class(pftcc_orisrch_grad), intent(inout) :: self !< instance
        type(opt_factory) :: ofac
        real              :: lims(5,2)
        call self%kill
        ! get # rotations
        self%nrots = pftcc_glob%get_nrots()
        ! just dummy limits for construction, these will be updated before minimization
        lims = 0.
        lims(:,2) = 360
        ! specify & set cost / gradient function
        call self%ospec%specify('lbfgsb', 5, ftol=1e-4, gtol=1e-4, limits=lims)
        call self%ospec%set_fdfcostfun_8(fdfcostfun)
        call self%ospec%set_costfun_8(costfun)
        call self%ospec%set_gcostfun_8(gcostfun)
        ! generate optimizer object with the factory
        call ofac%new(self%ospec, self%nlopt)
        ! flag existence
        self%exists = .true.
    end subroutine new

    !> set particle index for search
    subroutine set_particle( self, ptcl )
        class(pftcc_orisrch_grad), intent(inout) :: self
        integer,                   intent(in)    :: ptcl
        self%particle  = ptcl
    end subroutine set_particle

    !> minimisation
    function minimize( self, o_inout, angerr_deg, maxHWshift, found_better ) result( cxy )
        use simple_ori, only: ori
        class(pftcc_orisrch_grad), intent(inout) :: self
        class(ori),                intent(inout) :: o_inout
        real,                      intent(in)    :: angerr_deg, maxHWshift
        logical,                   intent(out)   :: found_better
        complex, allocatable :: pft_ref_even(:,:) !< for thread safe use of pftcc
        complex, allocatable :: pft_ref_odd(:,:)  !< -"-
        real     :: cxy(3), rotmat(2,2), cost, cost_init
        real(dp) :: f, grad(5), vec(5)
        integer  :: ithr
        ! copy the pftcc references so we can put them back after minimization is done
        ithr         = omp_get_thread_num() + 1
        pft_ref_even = pftcc_glob%get_ref_pft(ithr, iseven=.true.)
        pft_ref_odd  = pftcc_glob%get_ref_pft(ithr, iseven=.false.)

        ! PREPARATION
        ! initialize with input projection direction
        self%ospec%x(1) = o_inout%e1get()
        self%ospec%x(2) = o_inout%e2get()
        self%ospec%x(3) = o_inout%e3get()
        self%ospec%x(4) = 0. ! because shifts are updated with vector addition between successive rounds
        self%ospec%x(5) = 0. ! -"-
        ! set limits for L-BFGS-B initialisation
        self%ospec%limits(1,1) = self%ospec%x(1) - angerr_deg
        self%ospec%limits(1,2) = self%ospec%x(1) + angerr_deg
        self%ospec%limits(2,1) = self%ospec%x(2) - angerr_deg
        self%ospec%limits(2,2) = self%ospec%x(2) + angerr_deg
        self%ospec%limits(3,1) = self%ospec%x(3) - angerr_deg
        self%ospec%limits(3,2) = self%ospec%x(3) + angerr_deg
        self%ospec%limits(4,1) = -maxHWshift
        self%ospec%limits(4,2) =  maxHWshift
        self%ospec%limits(5,:) = self%ospec%limits(4,:)

        ! MINIMISATION
        self%ospec%nevals = 0
        vec = self%ospec%x
        call fdfcostfun(self, vec, f, grad, 5)
        cost_init = real(f) ! for later comparison
        call self%nlopt%minimize(self%ospec, self, cost)

        ! OUTPUT
        if( cost <= cost_init )then
            ! set corr
            cxy(1)  = -cost ! correlation
            ! set Euler
            call o_inout%set_euler(self%ospec%x(1:3))
            ! rotate the shift vector to the frame of reference
            cxy(2:) = self%ospec%x(4:5)
            call rotmat2d(self%ospec%x(3), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
            ! update ori
            call o_inout%set_shift(cxy(2:))
            ! indicate that better was found
            found_better = .true.
        else
            cxy(1)       = -cost_init ! correlation
            cxy(2:)      = 0.
            found_better = .false.
        endif

        ! put back references
        call pftcc_glob%set_ref_pft(ithr, pft_ref_even, iseven=.true.)
        call pftcc_glob%set_ref_pft(ithr, pft_ref_odd,  iseven=.false.)
    end function minimize

    subroutine gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp)                :: corr!, cost
        real(dp)                :: corr_grad(5)
        integer :: ithr, irot
        ! thread-safe extraction of projection and derivatives (because pftcc is an OpenMP shared data structure)
        ithr = omp_get_thread_num() + 1
        ! because we start from a continous solution
        irot = 1
        select type(self)
            class is (pftcc_orisrch_grad)
                if( pftcc_glob%ptcl_iseven(self%particle) )then
                    call build_glob%vol%fdf_project_polar(ithr, vec(1:3), pftcc_glob, iseven=.true.)
                else
                    call build_glob%vol_odd%fdf_project_polar(ithr, vec(1:3), pftcc_glob, iseven=.false.)
                endif
                call pftcc_glob%gencorr_cont_shift_grad_cc_for_rot_8(ithr, self%particle, vec(4:5), irot, corr, corr_grad)
                grad = - corr_grad
            class default
                write (*,*) 'error in pftcc_orisrch_grad :: gcostfun: unknown type'
                grad = 0.
                stop
        end select
    end subroutine gcostfun

    function costfun( self, vec, D ) result( cost )
        use simple_ori, only: ori
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: corr, cost
        !real(dp)                :: corr_grad(5)
        integer   :: ithr, irot
        type(ori) :: e
        ! thread-safe extraction of projection and derivatives (because pftcc is an OpenMP shared data structure)
        ithr = omp_get_thread_num() + 1
        ! because we start from a continous solution
        irot = 1
        select type(self)
            class is (pftcc_orisrch_grad)
                call e%new()
                call e%set_euler(real(vec(1:3)))
                if( pftcc_glob%ptcl_iseven(self%particle) )then
                    call build_glob%vol%fproject_polar(ithr, e, pftcc_glob, iseven=.true.)
                else
                    call build_glob%vol_odd%fproject_polar(ithr, e, pftcc_glob, iseven=.false.)
                endif
                corr = pftcc_glob%gencorr_cc_for_rot_8(ithr, self%particle, vec(4:5), irot)
                cost = -corr
                call e%kill()
            class default
                write (*,*) 'error in pftcc_orisrch_grad :: costfun: unknown type'
                cost = 0.
                stop
        end select
    end function costfun

    subroutine fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(dp)                :: corr
        real(dp)                :: corr_grad(5)
        integer :: ithr, irot
        ! thread-safe extraction of projection and derivatives (because pftcc is an OpenMP shared data structure)
        ithr = omp_get_thread_num() + 1
        ! because we start from a continous solution
        irot = 1
        select type(self)
            class is (pftcc_orisrch_grad)
                if( pftcc_glob%ptcl_iseven(self%particle) )then
                    call build_glob%vol%fdf_project_polar(ithr, vec(1:3), pftcc_glob, iseven=.true.)
                else
                    call build_glob%vol_odd%fdf_project_polar(ithr, vec(1:3), pftcc_glob, iseven=.false.)
                endif
                call pftcc_glob%gencorr_cont_shift_grad_cc_for_rot_8(ithr, self%particle, vec(4:5), irot, corr, corr_grad)
                f    = - corr
                grad = - corr_grad
            class default
                write (*,*) 'error in pftcc_orisrch_grad :: fdfcostfun: unknown type'
                f    = 0.
                grad = 0.
                stop
        end select
    end subroutine fdfcostfun

    ! destructor
    subroutine kill( self )
        class(pftcc_orisrch_grad), intent(inout) :: self !< instance
        if( self%exists )then
            if( associated(self%nlopt) )then
                call self%nlopt%kill
                deallocate(self%nlopt)
            end if
            call self%ospec%kill
            pftcc_glob => null()
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pftcc_orisrch_grad
