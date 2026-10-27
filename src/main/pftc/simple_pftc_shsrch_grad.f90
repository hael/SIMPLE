!@descr: rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftc_shsrch_grad
use iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_opt_spec,  only: opt_spec
use simple_optimizer, only: optimizer
use simple_builder,   only: builder
implicit none

public :: pftc_shsrch_grad
public :: fit_quadratic_step, grid_stencil_exists
public :: shift_ptcl2ref, shift_ref2ptcl
public :: FIT_EDGE, FIT_CURVATURE, FIT_STEP, FIT_LIMITS, FIT_VALUE, FIT_NREASONS
private
#include "simple_local_flags.inc"

integer,  parameter :: coarse_num_steps = 5       ! no. of coarse search steps in x AND y (hence real no. is its square)
integer,  parameter :: SHSRCH_SHIFT = 1, SHSRCH_JOINT = 2
! reasons the quadratic grid step falls back to the L-BFGS-B alternations
integer,  parameter :: FIT_EDGE = 1, FIT_CURVATURE = 2, FIT_STEP = 3, FIT_LIMITS = 4, FIT_VALUE = 5, FIT_NREASONS = 5
real(dp), parameter :: FIT_MIN_CURV_RATIO = 0.01d0 ! smallest accepted ratio of the Hessian's eigenvalues

type :: pftc_shsrch_grad
    private
    type(opt_spec)            :: ospec                  !< optimizer specification object
    class(optimizer), pointer :: opt_obj => null()      !< optimizer object
    class(builder),   pointer :: b_ptr   => null()      !< pointer to polarft_calc object for cost function evaluations
    integer                   :: reference    = 0       !< reference pft
    integer                   :: particle     = 0       !< particle pft
    integer                   :: nrots        = 0       !< # rotations
    integer                   :: maxits       = 100     !< max # iterations
    integer                   :: cur_inpl_idx = 0       !< index of inplane angle for shift search
    real                      :: cur_inpl_rotind = 0.   !< continuous in-plane rotation index
    integer                   :: max_evals    = 5       !< max # inplrot/shsrch cycles
    logical                   :: opt_angle    = .true.  !< alternate discrete in-plane and shift optimization
    integer                   :: search_mode = SHSRCH_SHIFT !< configured numerical algorithm
    real(dp)                  :: joint_initial_cost = 0.d0 !< first optimizer evaluation for monotonic acceptance
    logical                   :: joint_initial_cost_valid = .false.
    logical                   :: joint_cc = .false.      !< joint objective is cc (loss = -cc)
    logical                   :: coarse_init  = .false. !< whether to perform an intial coarse search over the range
    integer(int64)            :: profile_objective_evals = 0_int64
    integer(int64)            :: profile_gradient_evals  = 0_int64
    integer(int64)            :: profile_alternations    = 0_int64 !< L-BFGS-B / discrete in-plane alternations
    integer(int64)            :: profile_grid_points     = 0_int64 !< coarse-grid points scored
    integer(int64)            :: profile_minimize_calls  = 0_int64 !< calls of minimize
    ! quadratic step from the coarse grid (the multi-state probabilistic table)
    logical                   :: grid_fit     = .false. !< quadratic step from the coarse grid before the alternations
    real,         allocatable :: grid_corrs(:,:)        !< objective of every grid point over all rotations (nrots, points)
    integer                   :: grid_ix_star = 0       !< grid point of the lowest cost, x index
    integer                   :: grid_iy_star = 0       !< grid point of the lowest cost, y index
    integer(int64)            :: profile_fit(0:FIT_NREASONS) = 0_int64 !< (0) fits kept, (1:5) fallbacks per FIT_* reason
    ! particle-frame coarse grid (the multi-state probabilistic table)
    logical                   :: grid_ptcl    = .false. !< the grid is the lattice of particle-frame shifts, scored from a bank
    complex(sp),  allocatable :: grid_bank(:,:,:)       !< shifted particle memos (pftsz+1, nk, lattice points + 1 scratch)
    integer                   :: grid_bank_ptcl = 0     !< particle the bank was built for (0: none)
    real(dp)                  :: grid_u_star(2) = 0.d0  !< lattice shift of the grid minimum
    integer(int64)            :: profile_bank_builds = 0_int64 !< shifted particle memos built
    ! the last all-rotation objective of minimize and its shift: L-BFGS-B's callback scores every accepted point and
    ! minimize scores the point it ends on, the same vector when the shifts are bitwise equal
    real,         allocatable :: vals_cache(:)
    real                      :: vals_cache_xy(2) = 0.
    logical                   :: vals_cache_on    = .false. !< only within minimize: the pftc may change between calls
    logical                   :: vals_cache_valid = .false.
    ! the shift prior (shift_prior=yes, the probability tables): a Gaussian on the shift's distance from the start
    ! of the minimisation, folded into the objective as the factor exp(-prior_w |xy - prior_xy0|^2) on the
    ! correlation exp(-L), so the search and the row hold the MAP loss L + |xy - xy0|^2 / (2 sigma^2 T);
    ! prior_w = 1 / (2 sigma^2 T) with sigma^2 the previous iteration's shift variance per axis and T the
    ! particle's temperature (set_shift_prior). 0: no prior. The joint route (the polish of the committed pose)
    ! adds the same term to its raw Euclidean loss; cc has no loss scale for it
    real(dp)                  :: prior_w          = 0.d0
    real(dp)                  :: prior_xy0(2)     = 0.d0
contains
    procedure          :: new_alternating      => grad_shsrch_new_alternating
    procedure          :: new_fixed       => grad_shsrch_new_fixed
    procedure          :: new_joint       => grad_shsrch_new_joint
    procedure          :: set_indices     => grad_shsrch_set_indices
    procedure          :: minimize        => grad_shsrch_minimize
    procedure          :: minimize_joint  => grad_shsrch_minimize_joint
    procedure          :: select_best_discrete_angle
    procedure          :: kill            => grad_shsrch_kill
    procedure          :: does_opt_angle
    procedure          :: uses_joint_inplane
    procedure          :: set_limits
    procedure          :: set_shift_prior
    procedure          :: clear_shift_prior
    procedure          :: has_shift_prior
    procedure, private :: prior_factor
    procedure, private :: coarse_search
    procedure, private :: coarse_search_opt_angle
    procedure          :: reset_profile     => grad_shsrch_reset_profile
    procedure          :: get_profile       => grad_shsrch_get_profile
    procedure          :: get_fit_profile   => grad_shsrch_get_fit_profile
    procedure          :: get_bank_builds   => grad_shsrch_get_bank_builds
    procedure, private :: try_grid_fit
    procedure, private :: build_grid_bank
    procedure, private :: objfun_vals
end type pftc_shsrch_grad

contains

    subroutine grad_shsrch_new_alternating( self, build, lims, lims_init, maxits, coarse_init, grid_fit, grid_ptcl )
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(:,:)
        real,              optional, intent(in)    :: lims_init(:,:)
        integer,           optional, intent(in)    :: maxits
        logical,           optional, intent(in)    :: coarse_init
        logical,           optional, intent(in)    :: grid_fit   !< quadratic step from the coarse grid, needs coarse_init
        logical,           optional, intent(in)    :: grid_ptcl  !< particle-frame coarse grid, needs coarse_init
        call grad_shsrch_new_mode(self, build, lims, lims_init, maxits, .true., coarse_init, grid_fit, grid_ptcl)
    end subroutine grad_shsrch_new_alternating

    subroutine grad_shsrch_new_fixed( self, build, lims, lims_init, maxits, coarse_init )
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(:,:)
        real,              optional, intent(in)    :: lims_init(:,:)
        integer,           optional, intent(in)    :: maxits
        logical,           optional, intent(in)    :: coarse_init
        call grad_shsrch_new_mode(self, build, lims, lims_init, maxits, .false., coarse_init)
    end subroutine grad_shsrch_new_fixed

    subroutine grad_shsrch_new_mode( self, build, lims, lims_init, maxits, opt_angle, coarse_init, grid_fit, grid_ptcl )
        use simple_opt_factory, only: opt_factory
        class(pftc_shsrch_grad),     intent(inout) :: self           !< instance
        class(builder),      target, intent(in)    :: build          !< builder object for pftc access
        real,                        intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,              optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        integer,           optional, intent(in)    :: maxits         !< maximum iterations
        logical,                     intent(in)    :: opt_angle      !< alternate discrete in-plane and shift optimization
        logical,           optional, intent(in)    :: coarse_init    !< coarse inital search
        logical,           optional, intent(in)    :: grid_fit       !< quadratic step from the coarse grid
        logical,           optional, intent(in)    :: grid_ptcl      !< particle-frame coarse grid
        type(opt_factory) :: opt_fact
        integer :: kfromto(2)
        call self%kill
        self%cur_inpl_rotind = 0.
        ! set pointer to pftc instance for cost function evaluations
        self%b_ptr => build
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%opt_angle = opt_angle
        self%coarse_init = .false.
        if( present(coarse_init) ) self%coarse_init = coarse_init
        self%search_mode = SHSRCH_SHIFT
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%opt_obj)
        ! get # rotations
        self%nrots = self%b_ptr%pftc%get_nrots()
        ! the quadratic step reads the all-rotation grid of the alternating search
        self%grid_fit = .false.
        if( present(grid_fit) ) self%grid_fit = grid_fit .and. self%opt_angle .and. self%coarse_init
        if( self%grid_fit ) allocate(self%grid_corrs(self%nrots, coarse_num_steps*coarse_num_steps), source=0.)
        ! the bank holds one shifted particle memo per lattice point and a scratch entry for the fitted shift
        self%grid_ptcl = .false.
        if( present(grid_ptcl) ) self%grid_ptcl = grid_ptcl .and. self%opt_angle .and. self%coarse_init
        if( self%grid_ptcl )then
            kfromto = self%b_ptr%pftc%get_kfromto()
            allocate(self%grid_bank(self%b_ptr%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1, &
                &coarse_num_steps*coarse_num_steps+1))
        endif
        self%grid_bank_ptcl = 0
        allocate(self%vals_cache(self%nrots))
        ! set costfun pointers
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        if( self%opt_angle ) self%ospec%opt_callback => grad_shsrch_update_discrete_angle_wrapper
    end subroutine grad_shsrch_new_mode

    subroutine grad_shsrch_new_joint( self, build, lims, maxits )
        use simple_opt_factory, only: opt_factory
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(3,2)
        integer,                     intent(in)    :: maxits
        type(opt_factory) :: opt_fact
        call self%kill
        self%b_ptr => build
        if( .not. self%b_ptr%pftc%is_joint_grad_objfun() )then
            THROW_HARD('joint angle path requires a supported Euclidean or cc objective')
        endif
        self%joint_cc       = self%b_ptr%pftc%is_cc_objfun()
        self%nrots          = self%b_ptr%pftc%get_nrots()
        self%maxits         = maxits
        self%opt_angle      = .false.
        self%search_mode     = SHSRCH_JOINT
        self%coarse_init     = .false.
        call self%ospec%specify('lbfgsb', 3, factr=1.0d+7, pgtol=1.0d-5, limits=lims, &
            &max_step=0.01, maxits=self%maxits)
        call opt_fact%new(self%ospec, self%opt_obj)
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
    end subroutine grad_shsrch_new_joint

    pure logical function does_opt_angle( self )
        class(pftc_shsrch_grad), intent(in) :: self
        does_opt_angle = self%opt_angle
    end function does_opt_angle

    pure logical function uses_joint_inplane( self )
        class(pftc_shsrch_grad), intent(in) :: self
        uses_joint_inplane = self%search_mode == SHSRCH_JOINT
    end function uses_joint_inplane

    subroutine set_limits( self, lims )
        class(pftc_shsrch_grad), intent(inout) :: self              !< instance
        real,                    intent(in)    :: lims(self%ospec%ndim,2) !< new limits
        call self%ospec%set_limits(lims)
    end subroutine set_limits

    !> The shift prior of the next minimisations: sigma2 the shift variance per axis (pixels^2, from the previous
    !! iteration's shift increments), temp the particle's temperature (posterior_temperature, untempered: the prior
    !! is part of the loss the tempering acts on). Centred at each minimisation's start (minimize sets it).
    subroutine set_shift_prior( self, sigma2, temp )
        class(pftc_shsrch_grad), intent(inout) :: self
        real,                    intent(in)    :: sigma2, temp
        if( self%search_mode == SHSRCH_JOINT .and. self%joint_cc ) THROW_HARD('the shift prior is a Euclidean loss term; no cc joint route')
        if( .not. (sigma2 > 0. .and. temp > 0.) ) THROW_HARD('set_shift_prior needs a positive variance and temperature')
        self%prior_w = 1.d0 / (2.d0 * real(sigma2,dp) * real(temp,dp))
    end subroutine set_shift_prior

    subroutine clear_shift_prior( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        self%prior_w   = 0.d0
        self%prior_xy0 = 0.d0
    end subroutine clear_shift_prior

    pure logical function has_shift_prior( self )
        class(pftc_shsrch_grad), intent(in) :: self
        has_shift_prior = self%prior_w > 0.d0
    end function has_shift_prior

    !> the prior's factor on the correlation at shift xy: exp(-prior_w |xy - xy0|^2), 1 without a prior
    pure real(dp) function prior_factor( self, xy )
        class(pftc_shsrch_grad), intent(in) :: self
        real(dp),                intent(in) :: xy(2)
        prior_factor = 1.d0
        if( self%prior_w > 0.d0 ) prior_factor = exp(-self%prior_w * sum((xy - self%prior_xy0)**2))
    end function prior_factor

    ! Objective dispatch for the joint continuous route. Evaluators return a
    ! quantity to minimize: raw Euclidean loss or -cc.
    subroutine joint_grad_at_angle( self, vec, cost, grad )
        class(pftc_shsrch_grad), intent(inout) :: self
        real(dp),                intent(in)    :: vec(3)
        real(dp),                intent(out)   :: cost, grad(3)
        if( self%joint_cc )then
            call self%b_ptr%pftc%gen_corr_grad_at_angle(self%reference, self%particle, &
                &vec(1:2), vec(3), cost, grad)
        else
            call self%b_ptr%pftc%gen_raw_euclid_grad_at_angle(self%reference, self%particle, &
                &vec(1:2), vec(3), cost, grad)
            if( self%prior_w > 0.d0 )then
                ! the MAP loss: L + w |xy - xy0|^2, its shift gradient + 2w (xy - xy0)
                cost      = cost      + self%prior_w * sum((vec(1:2) - self%prior_xy0)**2)
                grad(1:2) = grad(1:2) + 2.d0 * self%prior_w * (vec(1:2) - self%prior_xy0)
            endif
        endif
    end subroutine joint_grad_at_angle

    ! Map a joint cost to the matcher score: exp(-loss) in [0,1] for
    ! raw Euclidean, or -cost clamped to [-1,1] for cc (benign series
    ! overshoot must not mint an invalid score).
    real function joint_cost_to_score( self, cost )
        class(pftc_shsrch_grad), intent(in) :: self
        real(dp),                intent(in) :: cost
        if( self%joint_cc )then
            joint_cost_to_score = real(min(1.d0, max(-1.d0, -cost)))
        else
            joint_cost_to_score = real(exp(-max(0.d0, cost)))
        endif
    end function joint_cost_to_score

    function grad_shsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                if( self%search_mode == SHSRCH_JOINT )then
                    block
                        real(dp) :: joint_grad(3)
                        call joint_grad_at_angle(self, vec(1:3), cost, joint_grad)
                    end block
                    if( .not. self%joint_initial_cost_valid )then
                        self%joint_initial_cost = cost
                        self%joint_initial_cost_valid = .true.
                    endif
                else
                    cost = - self%b_ptr%pftc%gen_corr_for_rot_8_aln(self%reference, self%particle, vec, self%cur_inpl_idx)
                    if( self%prior_w > 0.d0 ) cost = cost * self%prior_factor(vec(1:2))
                endif
            class default
                THROW_HARD('error in grad_shsrch_costfun: unknown type; grad_shsrch_costfun')
        end select
    end function grad_shsrch_costfun

    subroutine grad_shsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp)                :: corrs, corrs_grad(2), joint_grad(3), joint_loss, g
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_gradient_evals = self%profile_gradient_evals + 1_int64
                if( self%search_mode == SHSRCH_JOINT )then
                    call joint_grad_at_angle(self, vec(1:3), joint_loss, joint_grad)
                    grad = joint_grad
                else
                    call self%b_ptr%pftc%gen_corr_grad_for_rot_8_aln(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
                    grad = - corrs_grad
                    if( self%prior_w > 0.d0 )then
                        ! d/dxy of -(c g): -(dc) g + c g 2w (xy - xy0)
                        g    = self%prior_factor(vec(1:2))
                        grad = grad * g + corrs * g * 2.d0 * self%prior_w * (vec(1:2) - self%prior_xy0)
                    endif
                endif
            class default
                THROW_HARD('error in grad_shsrch_gcostfun: unknown type; grad_shsrch_gcostfun')
        end select
    end subroutine grad_shsrch_gcostfun

    subroutine grad_shsrch_fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(dp)                :: corrs, g
        real(dp)                :: corrs_grad(2), joint_grad(3)
        f    = 0.
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                self%profile_gradient_evals  = self%profile_gradient_evals  + 1_int64
                if( self%search_mode == SHSRCH_JOINT )then
                    call joint_grad_at_angle(self, vec(1:3), f, joint_grad)
                    grad = joint_grad
                    if( .not. self%joint_initial_cost_valid )then
                        self%joint_initial_cost = f
                        self%joint_initial_cost_valid = .true.
                    endif
                else
                    call self%b_ptr%pftc%gen_corr_grad_for_rot_8_aln(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
                    f    = - corrs
                    grad = - corrs_grad
                    if( self%prior_w > 0.d0 )then
                        g    = self%prior_factor(vec(1:2))
                        f    = f * g
                        grad = grad * g + corrs * g * 2.d0 * self%prior_w * (vec(1:2) - self%prior_xy0)
                    endif
                endif
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine grad_shsrch_fdfcostfun

    subroutine grad_shsrch_update_discrete_angle( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        real :: corr
        call self%select_best_discrete_angle(self%ospec%x(1:2), self%cur_inpl_idx, corr)
        self%cur_inpl_rotind = real(self%cur_inpl_idx)
    end subroutine grad_shsrch_update_discrete_angle

    !> Apply discrete all-angle selection once at a caller-supplied
    !! native-frame shift without searching alternative shifts.
    subroutine select_best_discrete_angle( self, xy, irot, corr )
        class(pftc_shsrch_grad), intent(inout) :: self
        real,                    intent(in)    :: xy(2)
        integer,                 intent(out)   :: irot
        real,                    intent(out)   :: corr
        real :: objective(self%nrots)

        call self%objfun_vals(xy, objective)
        irot = maxloc(objective, dim=1)
        corr = objective(irot)
    end subroutine select_best_discrete_angle

    !> the objective over all rotations at shift xy; within minimize the last vector is reused at a bitwise equal shift
    subroutine objfun_vals( self, xy, vals )
        class(pftc_shsrch_grad), intent(inout) :: self
        real,                    intent(in)    :: xy(2)
        real,                    intent(out)   :: vals(self%nrots)
        if( self%vals_cache_valid )then
            if( all(transfer(xy, 0, 2) == transfer(self%vals_cache_xy, 0, 2)) )then
                vals = self%vals_cache
                return
            endif
        endif
        self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
        call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, xy, vals)
        if( self%prior_w > 0.d0 ) vals = vals * real(self%prior_factor(real(xy,dp)))
        if( self%vals_cache_on )then
            self%vals_cache       = vals
            self%vals_cache_xy    = xy
            self%vals_cache_valid = .true.
        endif
    end subroutine objfun_vals

    real function grad_shsrch_get_angle( self ) result(angle)
        class(pftc_shsrch_grad), intent(in) :: self
        angle = (self%cur_inpl_rotind - 1.) * self%b_ptr%pftc%get_dang()
    end function grad_shsrch_get_angle

    subroutine grad_shsrch_update_discrete_angle_wrapper( self )
        class(*), intent(inout) :: self
        select type(self)
            class is (pftc_shsrch_grad)
                call grad_shsrch_update_discrete_angle(self)
            class DEFAULT
                THROW_HARD('unknown type in discrete angle update; simple_pftc_shsrch_grad')
        end select
    end subroutine grad_shsrch_update_discrete_angle_wrapper

    !> set indicies for shift search
    subroutine grad_shsrch_set_indices( self, ref, ptcl )
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine grad_shsrch_set_indices

    !> minimisation routine
    function grad_shsrch_minimize( self, irot, sh_rot, xy_in ) result( cxy )
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(inout) :: irot
        logical, optional,       intent(in)    :: sh_rot
        real,    optional,       intent(in)    :: xy_in(2)
        real     :: corrs(self%nrots), rotmat(2,2), cxy(3), lowest_shift(2), lowest_cost
        real(dp) :: init_xy(2), lowest_cost_overall, coarse_cost
        integer  :: loc, i, lowest_rot, init_rot
        logical  :: found_better, l_sh_rot, coarse_init_orig, l_fitted
        if( self%search_mode /= SHSRCH_SHIFT )then
            THROW_HARD('shift minimization requested from a joint search object')
        endif
        l_sh_rot = .true.
        if( present(sh_rot)  ) l_sh_rot = sh_rot
        if( present(xy_in)   )then
            coarse_init_orig = self%coarse_init
            self%coarse_init = .false.
            self%ospec%x     = xy_in
        else
            self%ospec%x     = [0.,0.]
        endif
        self%ospec%x_8 = dble(self%ospec%x)
        self%prior_xy0 = self%ospec%x_8 ! the prior is centred at the start: the seed, or zero
        found_better   = .false.
        self%profile_minimize_calls = self%profile_minimize_calls + 1_int64
        ! the objective's double-precision reference, copied once per minimisation
        call self%b_ptr%pftc%load_ref_aln_8(self%reference, self%particle)
        if( self%opt_angle )then
            self%vals_cache_on    = .true.
            self%vals_cache_valid = .false.
            call self%objfun_vals(self%ospec%x, corrs)
            self%cur_inpl_idx = maxloc(corrs, dim=1)
            self%cur_inpl_rotind = real(self%cur_inpl_idx)
            lowest_cost_overall = -corrs(self%cur_inpl_idx)
            l_fitted = .false.
            if( self%coarse_init )then
                call self%coarse_search_opt_angle(init_xy, init_rot)
                if( init_rot /= 0 )then
                    self%ospec%x_8      = init_xy
                    self%ospec%x        = real(init_xy)
                    self%cur_inpl_idx   = init_rot
                    self%cur_inpl_rotind = real(init_rot)
                    ! no fit: the alternations below run from the grid minimum as before
                    if( self%grid_fit ) call self%try_grid_fit(corrs, l_fitted)
                endif
            end if
            ! shift search / in-plane rot update
            if( .not. l_fitted )then
                do i = 1,self%max_evals
                    self%profile_alternations = self%profile_alternations + 1_int64
                    call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                    loc = self%cur_inpl_idx
                    call self%objfun_vals(self%ospec%x, corrs)
                    self%cur_inpl_idx = maxloc(corrs, dim=1)
                    self%cur_inpl_rotind = real(self%cur_inpl_idx)
                    if( self%cur_inpl_idx == loc ) exit
                end do
            endif
            ! update best
            lowest_cost = -corrs(self%cur_inpl_idx)
            if( lowest_cost < lowest_cost_overall )then
                found_better        = .true.
                lowest_cost_overall = lowest_cost
                lowest_rot          = self%cur_inpl_idx
                lowest_shift        = self%ospec%x
            endif
            if( found_better )then
                irot    =   lowest_rot                 ! in-plane index
                cxy(1) = - real(lowest_cost_overall)  ! correlation
                cxy(2:) =   real(lowest_shift)         ! shift
                if( l_sh_rot )then
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                    cxy(2:) = matmul(cxy(2:), rotmat)
                endif
            else
                irot = 0 ! to communicate that a better solution was not found
            endif
        else
            self%cur_inpl_idx   = irot
            self%cur_inpl_rotind = real(irot)
            self%profile_objective_evals = self%profile_objective_evals + 1_int64
            lowest_cost_overall = -self%b_ptr%pftc%gen_corr_for_rot_8_aln(self%reference, self%particle, self%ospec%x_8, self%cur_inpl_idx)
            if( self%coarse_init )then
                call self%coarse_search(coarse_cost, init_xy)
                if( coarse_cost < lowest_cost_overall )then
                    lowest_cost_overall = coarse_cost
                    self%ospec%x_8      = init_xy
                    self%ospec%x        = real(init_xy)
                endif
            end if
            ! shift search
            call self%opt_obj%minimize(self%ospec, self, lowest_cost)
            if( lowest_cost < lowest_cost_overall )then
                found_better        = .true.
                lowest_cost_overall = lowest_cost
                lowest_shift        = self%ospec%x
            endif
            if( found_better )then
                cxy(1)  = - real(lowest_cost_overall)  ! correlation
                cxy(2:) =   lowest_shift               ! shift
                if( l_sh_rot )then
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                    cxy(2:) = matmul(cxy(2:), rotmat)
                endif
            else
                irot = 0 ! to communicate that a better solution was not found
            endif
        end if
        self%vals_cache_on    = .false.
        self%vals_cache_valid = .false.
        if( present(xy_in) ) self%coarse_init = coarse_init_orig
    end function grad_shsrch_minimize

    !> Polish the caller's pose over (sx,sy,rotind_frac) within +/-2 cells of irot_in (euclid, cc);
    !! irot_in is authoritative, there is no all-angle rescan. An invalid or non-material solve returns the
    !! seed cell re-scored at xy_in (irot never 0, always committable); evaluation_valid/improved say which.
    function grad_shsrch_minimize_joint( self, irot, xy_in, sh_rot, rotind_frac, &
            &evaluation_valid, improved, initial_cost_out, irot_in ) result(cxy)
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(out)   :: irot
        real,                    intent(in)    :: xy_in(2)
        logical,                 intent(in)    :: sh_rot
        real(dp),                intent(out)   :: rotind_frac
        logical,       optional, intent(out)   :: evaluation_valid, improved
        real(dp),      optional, intent(out)   :: initial_cost_out
        integer,                 intent(in)    :: irot_in
        ! a final cost this far outside its physical range is a series artifact: invalid solve
        real(dp), parameter :: JOINT_NEG_COST_TOL = 1.d-2
        ! a gain must be material, not solver noise
        real(dp), parameter :: JOINT_IMPROVE_REL_TOL = 1.d-4
        ! bound-pinned solutions never displace the seed
        real(dp), parameter :: JOINT_BOUND_TOL = 1.d-3
        real :: cxy(3), rotmat(2,2), lowest_cost, seed_corr, joint_lims(3,2)
        real(dp) :: initial_cost, final_cost, improve_tol, coordinate_tol(3), brange
        integer  :: idim
        logical :: valid_result, improved_result, valid_coordinates

        if( self%search_mode /= SHSRCH_JOINT )then
            THROW_HARD('joint minimization requested from a non-joint search object')
        endif
        if( .not. self%b_ptr%pftc%is_joint_grad_objfun() )then
            THROW_HARD('joint minimization requires a supported Euclidean or cc objective')
        endif
        if( irot_in < 1 .or. irot_in > self%nrots )then
            THROW_HARD('minimize_joint requires the selected in-plane cell as irot_in seed')
        endif
        ! seed at the caller's authoritative assignment, never rescan; its
        ! score at xy_in is recovered from the solver's initial cost below
        irot      = irot_in
        seed_corr = 0.0
        self%cur_inpl_idx    = irot
        self%cur_inpl_rotind = real(irot)
        joint_lims = self%ospec%limits
        joint_lims(3,:) = [real(irot)-2., real(irot)+2.]
        call self%set_limits(joint_lims)
        self%ospec%x = [xy_in, real(irot)]
        self%ospec%x_8 = dble(self%ospec%x)
        self%prior_xy0 = self%ospec%x_8(1:2) ! the prior is centred at the start: the committed shift
        self%joint_initial_cost = 0.d0
        self%joint_initial_cost_valid = .false.
        call self%opt_obj%minimize(self%ospec, self, lowest_cost)
        initial_cost = self%joint_initial_cost
        final_cost = real(lowest_cost,dp)
        if( present(initial_cost_out) ) initial_cost_out = initial_cost
        if( self%joint_initial_cost_valid .and. ieee_is_finite(initial_cost) )then
            ! the incoming cell re-scored at xy_in, same mapping as the improved path
            seed_corr = joint_cost_to_score(self, initial_cost)
        endif
        if( self%joint_cc )then
            ! Correlation costs can approach or cross zero, where a
            ! purely relative guard is toothless; floor the material scale at
            ! a 1% score change.
            improve_tol = max(64.d0 * epsilon(1.d0), &
                &JOINT_IMPROVE_REL_TOL * max(abs(initial_cost), 1.d-2))
        else
            improve_tol = max(64.d0 * epsilon(1.d0) * max(1.d0, abs(initial_cost), abs(final_cost)), &
                &JOINT_IMPROVE_REL_TOL * abs(initial_cost))
        endif
        coordinate_tol = 64.d0 * real(epsilon(1.0),dp) * &
            &max(1.d0, max(abs(real(self%ospec%limits(:,1),dp)), &
            &abs(real(self%ospec%limits(:,2),dp))))
        valid_coordinates = all(ieee_is_finite(self%ospec%x))
        if( valid_coordinates )then
            valid_coordinates = all(real(self%ospec%x,dp) >= &
                &real(self%ospec%limits(:,1),dp)-coordinate_tol) .and. &
                &all(real(self%ospec%x,dp) <= &
                &real(self%ospec%limits(:,2),dp)+coordinate_tol)
        endif
        valid_result = self%joint_initial_cost_valid .and. ieee_is_finite(initial_cost) .and. &
            &ieee_is_finite(final_cost) .and. valid_coordinates
        if( valid_result )then
            if( self%joint_cc )then
                ! a correlation materially beyond +/-1 is an unphysical
                ! interpolation artifact; the pose found by descending into it
                ! cannot be trusted
                if( abs(final_cost) > 1.d0 + JOINT_NEG_COST_TOL ) valid_result = .false.
            else
                ! a materially negative loss invalidates the solve:
                ! exponentiating it would mint a score > 1 that corrupts
                ! probability-table distances and downstream weighting, and
                ! the pose itself is a descent into series artifact, not signal
                if( final_cost < -JOINT_NEG_COST_TOL ) valid_result = .false.
            endif
        endif
        improved_result = valid_result .and. final_cost < initial_cost - improve_tol
        ! demote bound-pinned solutions to the discrete floor: collapsed
        ! dimensions (doshift=no) have zero range and are skipped
        if( improved_result )then
            do idim = 1, 3
                brange = real(self%ospec%limits(idim,2), dp) - real(self%ospec%limits(idim,1), dp)
                if( brange <= 0.d0 ) cycle
                if( real(self%ospec%x(idim),dp) - real(self%ospec%limits(idim,1),dp) < JOINT_BOUND_TOL * brange .or. &
                    &real(self%ospec%limits(idim,2),dp) - real(self%ospec%x(idim),dp) < JOINT_BOUND_TOL * brange )then
                    improved_result = .false.
                    exit
                endif
            end do
        endif
        if( .not. valid_result )then
            ! return the incoming seed cell (irot = irot_in) and its score
            if( present(evaluation_valid) ) evaluation_valid = .false.
            if( present(improved) ) improved = .false.
            rotind_frac = real(irot,dp)
            cxy(1)  = seed_corr
            cxy(2:) = xy_in
            if( sh_rot )then
                call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                cxy(2:) = matmul(cxy(2:), rotmat)
            endif
            return
        endif
        if( present(evaluation_valid) ) evaluation_valid = valid_result
        if( present(improved) ) improved = improved_result
        if( .not. improved_result )then
            rotind_frac = real(irot,dp)
            cxy(1)  = seed_corr
            cxy(2:) = xy_in
            if( sh_rot )then
                call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                cxy(2:) = matmul(cxy(2:), rotmat)
            endif
            return
        endif
        self%cur_inpl_rotind = self%ospec%x(3)
        self%cur_inpl_idx = modulo(nint(self%cur_inpl_rotind)-1,self%nrots)+1
        irot = self%cur_inpl_idx
        ! clamp benign sub-tolerance series artifacts so the score keeps the
        ! matcher normalization contract ([0,1] Euclidean, [-1,1] cc)
        cxy(1) = joint_cost_to_score(self, final_cost)
        cxy(2:) = self%ospec%x(1:2)
        rotind_frac = self%cur_inpl_rotind
        if( sh_rot )then
            call rotmat2d(grad_shsrch_get_angle(self), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
        endif
    end function grad_shsrch_minimize_joint

    subroutine coarse_search(self, lowest_cost, init_xy)
        class(pftc_shsrch_grad), intent(inout) :: self
        real(dp),                intent(out)   :: lowest_cost, init_xy(2)
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
                y    = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                self%profile_grid_points     = self%profile_grid_points + 1_int64
                cost = -self%b_ptr%pftc%gen_corr_for_rot_8_aln(self%reference, self%particle, [x,y], self%cur_inpl_idx)
                if( self%prior_w > 0.d0 ) cost = cost * self%prior_factor([x,y])
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    init_xy     = [x,y]
                end if
            enddo
        enddo
    end subroutine coarse_search

    subroutine coarse_search_opt_angle(self, init_xy, irot)
        class(pftc_shsrch_grad), intent(inout) :: self
        real(dp),                intent(out)   :: init_xy(2)
        integer,                 intent(out)   :: irot
        real(dp) :: x, y, stepx,stepy
        real     :: corrs(self%nrots), lowest_cost, cost
        integer  :: loc, ix,iy
        init_xy = 0.d0
        irot    = 0
        self%grid_ix_star = 0
        self%grid_iy_star = 0
        if (coarse_num_steps .le. 1) return
        if( self%grid_ptcl .and. self%grid_bank_ptcl /= self%particle ) call self%build_grid_bank
        lowest_cost = huge(lowest_cost)
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            x = self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,coarse_num_steps
                y = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
                self%profile_grid_points     = self%profile_grid_points + 1_int64
                if( self%grid_ptcl )then
                    call self%b_ptr%pftc%gen_objfun_vals_bank(self%reference, self%particle, &
                        &self%grid_bank(:,:,(ix-1)*coarse_num_steps+iy), corrs)
                else
                    call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, real([x,y]), corrs)
                endif
                if( self%prior_w > 0.d0 ) corrs = corrs * real(self%prior_factor([x,y]))
                if( self%grid_fit ) self%grid_corrs(:,(ix-1)*coarse_num_steps+iy) = corrs
                loc  = maxloc(corrs,dim=1)
                cost = - corrs(loc)
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    irot        = loc
                    init_xy(1)  = x
                    init_xy(2)  = y
                    self%grid_ix_star = ix
                    self%grid_iy_star = iy
                end if
            end do
        end do
        if( self%grid_ptcl .and. irot /= 0 )then
            ! the lattice point is a particle-frame shift: the search starts from it in the reference frame at irot
            self%grid_u_star = init_xy
            init_xy = real(shift_ptcl2ref(real(init_xy), self%b_ptr%pftc%get_rot(irot)), dp)
            init_xy(1) = max(real(self%ospec%limits(1,1),dp), min(real(self%ospec%limits(1,2),dp), init_xy(1)))
            init_xy(2) = max(real(self%ospec%limits(2,1),dp), min(real(self%ospec%limits(2,2),dp), init_xy(2)))
        endif
    end subroutine coarse_search_opt_angle

    !> one shifted particle memo per lattice point of the coarse grid (the shift is the particle-frame shift)
    subroutine build_grid_bank( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        real(dp) :: stepx, stepy
        real     :: u(2)
        integer  :: ix, iy
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            u(1) = real(self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx)
            do iy = 1,coarse_num_steps
                u(2) = real(self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy)
                call self%b_ptr%pftc%memoize_shifted_ptcl(self%particle, u, self%grid_bank(:,:,(ix-1)*coarse_num_steps+iy))
            end do
        end do
        self%profile_bank_builds = self%profile_bank_builds + int(coarse_num_steps*coarse_num_steps, int64)
        self%grid_bank_ptcl      = self%particle
    end subroutine build_grid_bank

    !> a shift of the particle frame as the shift of the reference at in-plane angle rot (degrees);
    !! inverse of shift_ref2ptcl, the conversion minimize applies with sh_rot
    pure function shift_ptcl2ref( u, rot ) result( s )
        real, intent(in) :: u(2), rot
        real :: s(2), rotmat(2,2)
        call rotmat2d(rot, rotmat)
        s(1) = u(1) * rotmat(1,1) + u(2) * rotmat(1,2)
        s(2) = u(1) * rotmat(2,1) + u(2) * rotmat(2,2)
    end function shift_ptcl2ref

    !> a shift of the reference at in-plane angle rot (degrees) as the shift of the particle frame
    pure function shift_ref2ptcl( s, rot ) result( u )
        real, intent(in) :: s(2), rot
        real :: u(2), rotmat(2,2)
        call rotmat2d(rot, rotmat)
        u(1) = s(1) * rotmat(1,1) + s(2) * rotmat(2,1)
        u(2) = s(1) * rotmat(1,2) + s(2) * rotmat(2,2)
    end function shift_ref2ptcl

    !> a grid point has the full 3x3 stencil of the quadratic step unless it lies on the edge of the lattice
    pure logical function grid_stencil_exists( ix, iy )
        integer, intent(in) :: ix, iy
        grid_stencil_exists = ix >= 2 .and. ix <= coarse_num_steps-1 .and. iy >= 2 .and. iy <= coarse_num_steps-1
    end function grid_stencil_exists

    !> Stationary point (da,db), in grid spacings, of the least-squares quadratic over a 3x3 stencil v(a,b).
    !! reason: 0 = step, FIT_CURVATURE = H not positive definite or ill conditioned, FIT_STEP = leaves the centre cell.
    pure subroutine fit_quadratic_step( v, da, db, reason )
        real(dp), intent(in)  :: v(-1:1,-1:1)
        real(dp), intent(out) :: da, db
        integer,  intent(out) :: reason
        real(dp) :: S, Sa, Sb, Sab, Sa2, Sb2, g1, g2, h11, h22, h12, det, tr, disc, lmax, lmin
        integer  :: ia, ib
        da = 0.d0
        db = 0.d0
        S = 0.d0; Sa = 0.d0; Sb = 0.d0; Sab = 0.d0; Sa2 = 0.d0; Sb2 = 0.d0
        do ia = -1,1
            do ib = -1,1
                S   = S   + v(ia,ib)
                Sa  = Sa  + real(ia,dp) * v(ia,ib)
                Sb  = Sb  + real(ib,dp) * v(ia,ib)
                Sab = Sab + real(ia*ib,dp) * v(ia,ib)
                Sa2 = Sa2 + real(ia*ia,dp) * v(ia,ib)
                Sb2 = Sb2 + real(ib*ib,dp) * v(ia,ib)
            enddo
        enddo
        g1  = Sa / 6.d0
        g2  = Sb / 6.d0
        h11 = Sa2 - 2.d0 * S / 3.d0
        h22 = Sb2 - 2.d0 * S / 3.d0
        h12 = Sab / 4.d0
        det  = h11 * h22 - h12 * h12
        tr   = h11 + h22
        disc = sqrt(max(0.d0, tr * tr - 4.d0 * det))
        lmax = 0.5d0 * (tr + disc)
        lmin = 0.5d0 * (tr - disc)
        if( .not. (h11 > 0.d0 .and. det > 0.d0 .and. lmax > 0.d0 .and. lmin / lmax >= FIT_MIN_CURV_RATIO) )then
            reason = FIT_CURVATURE
            return
        endif
        da = -( h22 * g1 - h12 * g2) / det
        db = -(-h12 * g1 + h11 * g2) / det
        if( .not. (abs(da) <= 0.5d0 .and. abs(db) <= 0.5d0) )then
            reason = FIT_STEP
            return
        endif
        reason = 0
    end subroutine fit_quadratic_step

    !> Quadratic step from the coarse grid minimum: the 3x3 stencil of the grid loss at its rotation, -log(objective)
    !! for Euclid (the objective is exp(-loss)), -objective otherwise. A failed guard counts its reason and leaves the
    !! state untouched; see doc/policies/3D/refine3D_policy.md.
    subroutine try_grid_fit( self, corrs, fitted )
        class(pftc_shsrch_grad), intent(inout) :: self
        real,                    intent(out)   :: corrs(self%nrots)  !< objective over all rotations at the fitted shift
        logical,                 intent(out)   :: fitted
        real(dp) :: v(-1:1,-1:1), da, db, stepx, stepy, fx, fy
        real     :: fit_xy(2), best_grid_corr
        integer  :: ix, iy, ia, ib, irot_star, irot_fit, reason
        logical  :: l_euclid
        fitted = .false.
        ix = self%grid_ix_star
        iy = self%grid_iy_star
        irot_star = self%cur_inpl_idx
        if( .not. grid_stencil_exists(ix, iy) )then
            self%profile_fit(FIT_EDGE) = self%profile_fit(FIT_EDGE) + 1_int64
            return
        endif
        l_euclid = self%b_ptr%pftc%is_euclid_objfun()
        do ia = -1,1
            do ib = -1,1
                v(ia,ib) = real(self%grid_corrs(irot_star,(ix+ia-1)*coarse_num_steps+iy+ib), dp)
                if( l_euclid )then
                    v(ia,ib) = -log(max(v(ia,ib), 1.d-30))
                else
                    v(ia,ib) = -v(ia,ib)
                endif
            enddo
        enddo
        call fit_quadratic_step(v, da, db, reason)
        if( reason /= 0 )then
            self%profile_fit(reason) = self%profile_fit(reason) + 1_int64
            return
        endif
        stepx = real(self%ospec%limits(1,2) - self%ospec%limits(1,1), dp) / real(coarse_num_steps, dp)
        stepy = real(self%ospec%limits(2,2) - self%ospec%limits(2,1), dp) / real(coarse_num_steps, dp)
        best_grid_corr = self%grid_corrs(irot_star,(ix-1)*coarse_num_steps+iy)
        if( self%grid_ptcl )then
            ! the stencil lives in the particle frame: the fitted lattice shift is converted to the reference
            ! frame at the rotation it is kept at, which must keep it inside the box
            fx = self%grid_u_star(1) + da * stepx
            fy = self%grid_u_star(2) + db * stepy
            fit_xy = real([fx,fy])
            if( .not. in_limits(self, real(shift_ptcl2ref(fit_xy, self%b_ptr%pftc%get_rot(irot_star)),dp)) )then
                self%profile_fit(FIT_LIMITS) = self%profile_fit(FIT_LIMITS) + 1_int64
                return
            endif
            ! kept only if the exact value at the fitted shift beats the grid minimum
            call self%b_ptr%pftc%memoize_shifted_ptcl(self%particle, fit_xy, self%grid_bank(:,:,size(self%grid_bank,3)))
            self%profile_bank_builds     = self%profile_bank_builds + 1_int64
            self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
            call self%b_ptr%pftc%gen_objfun_vals_bank(self%reference, self%particle, &
                &self%grid_bank(:,:,size(self%grid_bank,3)), corrs)
            ! the prior's factor at the particle-frame shift, the frame the lattice's factors were taken in
            if( self%prior_w > 0.d0 ) corrs = corrs * real(self%prior_factor(real(fit_xy,dp)))
            irot_fit = maxloc(corrs, dim=1)
            if( .not. (corrs(irot_fit) > best_grid_corr) )then
                self%profile_fit(FIT_VALUE) = self%profile_fit(FIT_VALUE) + 1_int64
                return
            endif
            fit_xy = shift_ptcl2ref(fit_xy, self%b_ptr%pftc%get_rot(irot_fit))
            if( .not. in_limits(self, real(fit_xy,dp)) )then
                self%profile_fit(FIT_LIMITS) = self%profile_fit(FIT_LIMITS) + 1_int64
                return
            endif
        else
            fx = real(self%ospec%x_8(1), dp) + da * stepx
            fy = real(self%ospec%x_8(2), dp) + db * stepy
            if( .not. in_limits(self, [fx,fy]) )then
                self%profile_fit(FIT_LIMITS) = self%profile_fit(FIT_LIMITS) + 1_int64
                return
            endif
            ! kept only if the exact value at the fitted shift beats the grid minimum
            fit_xy = real([fx,fy])
            self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
            call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, fit_xy, corrs)
            if( self%prior_w > 0.d0 ) corrs = corrs * real(self%prior_factor(real(fit_xy,dp)))
            irot_fit = maxloc(corrs, dim=1)
            if( .not. (corrs(irot_fit) > best_grid_corr) )then
                self%profile_fit(FIT_VALUE) = self%profile_fit(FIT_VALUE) + 1_int64
                return
            endif
        endif
        self%ospec%x         = fit_xy
        self%ospec%x_8       = dble(fit_xy)
        self%cur_inpl_idx    = irot_fit
        self%cur_inpl_rotind = real(irot_fit)
        self%profile_fit(0)  = self%profile_fit(0) + 1_int64
        fitted = .true.
    end subroutine try_grid_fit

    !> shifted particle memos built by the particle-frame grid
    pure integer(int64) function grad_shsrch_get_bank_builds( self )
        class(pftc_shsrch_grad), intent(in) :: self
        grad_shsrch_get_bank_builds = self%profile_bank_builds
    end function grad_shsrch_get_bank_builds

    !> a reference-frame shift lies inside the search box
    logical function in_limits( self, xy )
        class(pftc_shsrch_grad), intent(in) :: self
        real(dp),                intent(in) :: xy(2)
        in_limits = xy(1) >= real(self%ospec%limits(1,1),dp) .and. xy(1) <= real(self%ospec%limits(1,2),dp) .and. &
                   &xy(2) >= real(self%ospec%limits(2,1),dp) .and. xy(2) <= real(self%ospec%limits(2,2),dp)
    end function in_limits

    !> fits kept (0) and fallbacks per FIT_* reason (1:5)
    subroutine grad_shsrch_get_fit_profile( self, counts )
        class(pftc_shsrch_grad), intent(in)  :: self
        integer(int64),          intent(out) :: counts(0:FIT_NREASONS)
        counts = self%profile_fit
    end subroutine grad_shsrch_get_fit_profile

    subroutine grad_shsrch_reset_profile( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        self%profile_objective_evals = 0_int64
        self%profile_gradient_evals  = 0_int64
        self%profile_alternations    = 0_int64
        self%profile_grid_points     = 0_int64
        self%profile_minimize_calls  = 0_int64
        self%profile_fit             = 0_int64
        self%profile_bank_builds     = 0_int64
    end subroutine grad_shsrch_reset_profile

    subroutine grad_shsrch_get_profile( self, objective_evals, gradient_evals, alternations, grid_points, minimize_calls )
        class(pftc_shsrch_grad), intent(in)  :: self
        integer(int64),          intent(out) :: objective_evals, gradient_evals, alternations, grid_points, minimize_calls
        objective_evals = self%profile_objective_evals
        gradient_evals  = self%profile_gradient_evals
        alternations    = self%profile_alternations
        grid_points     = self%profile_grid_points
        minimize_calls  = self%profile_minimize_calls
    end subroutine grad_shsrch_get_profile

    subroutine grad_shsrch_kill( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        if( associated(self%opt_obj) )then
            call self%opt_obj%kill
            nullify(self%opt_obj)
        end if
        call self%ospec%kill
        nullify(self%b_ptr)
        if( allocated(self%grid_corrs) ) deallocate(self%grid_corrs)
        if( allocated(self%grid_bank)  ) deallocate(self%grid_bank)
        if( allocated(self%vals_cache) ) deallocate(self%vals_cache)
        self%vals_cache_on    = .false.
        self%vals_cache_valid = .false.
        self%grid_fit       = .false.
        self%grid_ptcl      = .false.
        self%grid_bank_ptcl = 0
        self%search_mode = SHSRCH_SHIFT
        self%joint_initial_cost = 0.d0
        self%joint_initial_cost_valid = .false.
        self%joint_cc = .false.
        self%prior_w   = 0.d0
        self%prior_xy0 = 0.d0
        call self%reset_profile
    end subroutine grad_shsrch_kill

end module simple_pftc_shsrch_grad
