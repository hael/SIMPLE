!@descr: rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftc_shsrch_grad
use iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_opt_spec,  only: opt_spec
use simple_optimizer, only: optimizer
use simple_builder,   only: builder
implicit none

public :: pftc_shsrch_grad, bounded_shift_trial
private
#include "simple_local_flags.inc"

integer,  parameter :: coarse_num_steps = 5       ! no. of coarse search steps in x AND y (hence real no. is its square)
integer,  parameter :: direct_max_backtracks = 8
integer,  parameter :: SHSRCH_LEGACY = 1, SHSRCH_DIRECT = 2, SHSRCH_JOINT = 3

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
    integer                   :: search_mode = SHSRCH_LEGACY !< configured numerical algorithm
    real(dp)                  :: joint_initial_cost = 0.d0 !< first optimizer evaluation for monotonic acceptance
    logical                   :: joint_initial_cost_valid = .false.
    logical                   :: coarse_init  = .false. !< whether to perform an intial coarse search over the range
    logical                   :: raw_roundtrip_check = .false. !< validate vector/scalar raw loss when diagnostics are enabled
    logical                   :: raw_roundtrip_failed = .false. !< diagnostic mismatch observed during this search
    logical                   :: raw_roundtrip_offset_set = .false. !< first vector-scalar offset calibrated
    real(dp)                  :: raw_roundtrip_offset = 0.d0 !< diagnostic-only legacy score baseline offset
    integer(int64)            :: profile_objective_evals = 0_int64
    integer(int64)            :: profile_gradient_evals  = 0_int64
contains
    procedure          :: new_legacy      => grad_shsrch_new_legacy
    procedure          :: new_fixed       => grad_shsrch_new_fixed
    procedure          :: new_direct      => grad_shsrch_new_direct
    procedure          :: new_joint       => grad_shsrch_new_joint
    procedure          :: set_indices     => grad_shsrch_set_indices
    procedure          :: minimize        => grad_shsrch_minimize
    procedure          :: minimize_direct => grad_shsrch_minimize_direct
    procedure          :: minimize_joint  => grad_shsrch_minimize_joint
    procedure          :: minimize_joint_rounded => grad_shsrch_minimize_joint_rounded
    procedure          :: select_best_discrete_angle
    procedure          :: kill            => grad_shsrch_kill
    procedure          :: does_opt_angle
    procedure          :: uses_joint_inplane
    procedure          :: set_limits
    procedure          :: coarse_search
    procedure          :: coarse_search_opt_angle
    procedure          :: is_direct_shift_only
    procedure          :: set_diagnostic_mode
    procedure          :: diagnostic_failed => grad_shsrch_diagnostic_failed
    procedure          :: reset_profile     => grad_shsrch_reset_profile
    procedure          :: get_profile       => grad_shsrch_get_profile
end type pftc_shsrch_grad

contains

    subroutine set_diagnostic_mode( self, enabled )
        class(pftc_shsrch_grad), intent(inout) :: self
        logical,                 intent(in)    :: enabled
        self%raw_roundtrip_check = enabled
        self%raw_roundtrip_failed = .false.
        self%raw_roundtrip_offset_set = .false.
        self%raw_roundtrip_offset = 0.d0
    end subroutine set_diagnostic_mode

    logical function grad_shsrch_diagnostic_failed( self )
        class(pftc_shsrch_grad), intent(in) :: self
        grad_shsrch_diagnostic_failed = self%raw_roundtrip_failed
    end function grad_shsrch_diagnostic_failed

    subroutine grad_shsrch_new_legacy( self, build, lims, lims_init, maxits, coarse_init )
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(:,:)
        real,              optional, intent(in)    :: lims_init(:,:)
        integer,           optional, intent(in)    :: maxits
        logical,           optional, intent(in)    :: coarse_init
        call grad_shsrch_new_mode(self, build, lims, lims_init, maxits, .true., coarse_init)
    end subroutine grad_shsrch_new_legacy

    subroutine grad_shsrch_new_fixed( self, build, lims, lims_init, maxits, coarse_init )
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(:,:)
        real,              optional, intent(in)    :: lims_init(:,:)
        integer,           optional, intent(in)    :: maxits
        logical,           optional, intent(in)    :: coarse_init
        call grad_shsrch_new_mode(self, build, lims, lims_init, maxits, .false., coarse_init)
    end subroutine grad_shsrch_new_fixed

    subroutine grad_shsrch_new_mode( self, build, lims, lims_init, maxits, opt_angle, coarse_init )
        use simple_opt_factory, only: opt_factory
        class(pftc_shsrch_grad),     intent(inout) :: self           !< instance
        class(builder),      target, intent(in)    :: build          !< builder object for pftc access
        real,                        intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,              optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        integer,           optional, intent(in)    :: maxits         !< maximum iterations
        logical,                     intent(in)    :: opt_angle      !< alternate discrete in-plane and shift optimization
        logical,           optional, intent(in)    :: coarse_init    !< coarse inital search
        type(opt_factory) :: opt_fact
        call self%kill
        self%cur_inpl_rotind = 0.
        ! set pointer to pftc instance for cost function evaluations
        self%b_ptr => build
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%opt_angle = opt_angle
        self%coarse_init = .false.
        if( present(coarse_init) ) self%coarse_init = coarse_init
        self%search_mode = SHSRCH_LEGACY
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%opt_obj)
        ! get # rotations
        self%nrots = self%b_ptr%pftc%get_nrots()
        ! set costfun pointers
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        if( self%opt_angle ) self%ospec%opt_callback => grad_shsrch_update_discrete_angle_wrapper
    end subroutine grad_shsrch_new_mode

    subroutine grad_shsrch_new_direct( self, build, lims )
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(2,2)
        call self%kill
        self%b_ptr         => build
        self%nrots          = self%b_ptr%pftc%get_nrots()
        self%opt_angle      = .false.
        self%search_mode     = SHSRCH_DIRECT
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims, &
            &max_step=0.01, maxits=1)
    end subroutine grad_shsrch_new_direct

    subroutine grad_shsrch_new_joint( self, build, lims, maxits )
        use simple_opt_factory, only: opt_factory
        class(pftc_shsrch_grad),     intent(inout) :: self
        class(builder),      target, intent(in)    :: build
        real,                        intent(in)    :: lims(3,2)
        integer,                     intent(in)    :: maxits
        type(opt_factory) :: opt_fact
        call self%kill
        self%b_ptr => build
        if( .not. self%b_ptr%pftc%is_raw_euclid_objfun() )then
            THROW_HARD('joint angle path requires raw Euclidean objective; hybrid derivative is unavailable')
        endif
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

    pure logical function is_direct_shift_only( self )
        class(pftc_shsrch_grad), intent(in) :: self
        ! The streaming optimizer operates on a fixed discrete rotation and
        ! updates only s=(sx,sy); it must not allocate or fall through to the
        ! legacy L-BFGS-B implementation.
        is_direct_shift_only = self%search_mode == SHSRCH_DIRECT
    end function is_direct_shift_only

    subroutine set_limits( self, lims )
        class(pftc_shsrch_grad), intent(inout) :: self              !< instance
        real,                    intent(in)    :: lims(self%ospec%ndim,2) !< new limits
        call self%ospec%set_limits(lims)
    end subroutine set_limits

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
                        call self%b_ptr%pftc%gen_raw_euclid_grad_at_angle(self%reference, self%particle, &
                            &vec(1:2), vec(3), cost, joint_grad)
                    end block
                    if( .not. self%joint_initial_cost_valid )then
                        self%joint_initial_cost = cost
                        self%joint_initial_cost_valid = .true.
                    endif
                else
                    cost = - self%b_ptr%pftc%gen_corr_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx)
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
        real(dp)                :: corrs_grad(2), joint_grad(3), joint_loss
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_gradient_evals = self%profile_gradient_evals + 1_int64
                if( self%search_mode == SHSRCH_JOINT )then
                    call self%b_ptr%pftc%gen_raw_euclid_grad_at_angle(self%reference, self%particle, &
                        &vec(1:2), vec(3), joint_loss, joint_grad)
                    grad = joint_grad
                else
                    call self%b_ptr%pftc%gen_corr_grad_only_for_rot_8(self%reference, &
                        &self%particle, vec, self%cur_inpl_idx, corrs_grad)
                    grad = - corrs_grad
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
        real(dp)                :: corrs
        real(dp)                :: corrs_grad(2), joint_grad(3)
        f    = 0.
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                self%profile_gradient_evals  = self%profile_gradient_evals  + 1_int64
                if( self%search_mode == SHSRCH_JOINT )then
                    call self%b_ptr%pftc%gen_raw_euclid_grad_at_angle(self%reference, self%particle, &
                        &vec(1:2), vec(3), f, joint_grad)
                    grad = joint_grad
                    if( .not. self%joint_initial_cost_valid )then
                        self%joint_initial_cost = f
                        self%joint_initial_cost_valid = .true.
                    endif
                else
                    call self%b_ptr%pftc%gen_corr_grad_for_rot_8(self%reference, self%particle, &
                        &vec, self%cur_inpl_idx, corrs, corrs_grad)
                    f    = - corrs
                    grad = - corrs_grad
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

        self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
        call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, xy, objective)
        irot = maxloc(objective, dim=1)
        corr = objective(irot)
    end subroutine select_best_discrete_angle

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
        logical  :: found_better, l_sh_rot, coarse_init_orig
        if( self%search_mode /= SHSRCH_LEGACY )then
            THROW_HARD('legacy minimization requested from a specialized search object')
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
        found_better   = .false.
        if( self%opt_angle )then
            call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, self%ospec%x, corrs)
            self%cur_inpl_idx = maxloc(corrs, dim=1)
            self%cur_inpl_rotind = real(self%cur_inpl_idx)
            lowest_cost_overall = -corrs(self%cur_inpl_idx)
            if( self%coarse_init )then
                call self%coarse_search_opt_angle(init_xy, init_rot)
                if( init_rot /= 0 )then
                    self%ospec%x_8      = init_xy
                    self%ospec%x        = real(init_xy)
                    self%cur_inpl_idx   = init_rot
                    self%cur_inpl_rotind = real(init_rot)
                endif
            end if
            ! shift search / in-plane rot update
            do i = 1,self%max_evals
                call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                loc = self%cur_inpl_idx
                call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, self%ospec%x, corrs)
                self%cur_inpl_idx = maxloc(corrs, dim=1)
                self%cur_inpl_rotind = real(self%cur_inpl_idx)
                if( self%cur_inpl_idx == loc ) exit
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
            lowest_cost_overall = -self%b_ptr%pftc%gen_corr_for_rot_8(self%reference, &
                &self%particle, self%ospec%x_8, self%cur_inpl_idx)
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
        if( present(xy_in) ) self%coarse_init = coarse_init_orig
    end function grad_shsrch_minimize

    !> Classical Euclidean joint refinement over (sx,sy,rotind_frac).
    !! The selected reference/class is fixed. Without irot_in, a discrete
    !! all-angle scan chooses the best grid cell at xy_in. With irot_in, the
    !! caller's selected grid cell is authoritative and no scan is performed.
    !! The continuous solve starts from the selected cell with bounds of plus
    !! or minus two rotation indices. The PFTC objective is periodic, so those
    !! bounds may safely straddle the first or last grid index.
    !!
    !! A non-improving or numerically invalid solve returns its selected seed
    !! cell and score. evaluation_valid=.false. is a diagnostic signal only;
    !! irot is never 0 and callers always receive a committable discrete pose.
    !!
    !! irot_in switches to LOCAL mode for durable passes: the caller's selected
    !! assignment (already the product of exhaustive candidate construction) is
    !! authoritative, so no global all-angle reselection is performed -- under
    !! exactly degenerate in-plane branches (e.g. dihedral C2 about the symmetry
    !! axis) a second global selection at a slightly different shift hops
    !! branches on floating-point noise. The solve refines within +/-2 cells of
    !! irot_in; a non-improving solve returns the incoming cell re-scored at
    !! xy_in, so committing it retains the incoming pose.
    function grad_shsrch_minimize_joint( self, irot, xy_in, sh_rot, rotind_frac, &
            &evaluation_valid, improved, initial_cost_out, irot_in ) result(cxy)
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(out)   :: irot
        real,                    intent(in)    :: xy_in(2)
        logical,                 intent(in)    :: sh_rot
        real(dp),                intent(out)   :: rotind_frac
        logical,       optional, intent(out)   :: evaluation_valid, improved
        real(dp),      optional, intent(out)   :: initial_cost_out
        integer,       optional, intent(in)    :: irot_in
        ! The raw Euclidean loss is nonnegative by construction; the truncated
        ! coefficient series can undershoot slightly at fractional angles, but
        ! a final cost below this tolerance is an unphysical evaluator
        ! artifact and the pose found by descending into it cannot be trusted
        real(dp), parameter :: JOINT_NEG_COST_TOL = 1.d-2
        ! An improvement must be material to displace the exhaustive discrete
        ! floor: the roundoff-scale guard alone lets any solver twitch count as
        ! "improved" (L-BFGS-B's own convergence tolerance is orders looser),
        ! which commits fractional poses on noise and, under symmetric point
        ! groups with near-degenerate in-plane branches, flips assignments en
        ! masse between iterations
        real(dp), parameter :: JOINT_IMPROVE_REL_TOL = 1.d-4
        ! A solution pinned to a search bound is an uncontrolled excursion --
        ! for the rotation window it contradicts the exhaustive seed scan, and
        ! for shifts it is the corner-landing signature of a descent into
        ! series artifacts -- so it never displaces the floor either
        real(dp), parameter :: JOINT_BOUND_TOL = 1.d-3
        real :: cxy(3), rotmat(2,2), lowest_cost, seed_corr, joint_lims(3,2)
        real(dp) :: initial_cost, final_cost, improve_tol, coordinate_tol(3), brange
        integer  :: idim
        logical :: valid_result, improved_result, valid_coordinates, l_local

        if( self%search_mode /= SHSRCH_JOINT )then
            THROW_HARD('joint minimization requested from a non-joint search object')
        endif
        if( .not. self%b_ptr%pftc%is_raw_euclid_objfun() )then
            THROW_HARD('joint minimization requires raw Euclidean objective; hybrid derivative is unavailable')
        endif
        l_local = present(irot_in)
        if( l_local )then
            ! local mode: seed at the caller's authoritative assignment, no rescan;
            ! its score at xy_in is recovered from the solver's initial cost below
            irot      = irot_in
            seed_corr = 0.0
        else
            ! Select the grid seed once at the supplied shift before refinement.
            call self%select_best_discrete_angle(xy_in, irot, seed_corr)
        endif
        self%cur_inpl_idx    = irot
        self%cur_inpl_rotind = real(irot)
        joint_lims = self%ospec%limits
        joint_lims(3,:) = [real(irot)-2., real(irot)+2.]
        call self%set_limits(joint_lims)
        self%ospec%x = [xy_in, real(irot)]
        self%ospec%x_8 = dble(self%ospec%x)
        self%joint_initial_cost = 0.d0
        self%joint_initial_cost_valid = .false.
        call self%opt_obj%minimize(self%ospec, self, lowest_cost)
        initial_cost = self%joint_initial_cost
        final_cost = real(lowest_cost,dp)
        if( present(initial_cost_out) ) initial_cost_out = initial_cost
        if( l_local .and. self%joint_initial_cost_valid .and. ieee_is_finite(initial_cost) )then
            ! the incoming cell re-scored at xy_in, same mapping as the improved path
            seed_corr = real(exp(-max(0.d0, initial_cost)))
        endif
        improve_tol = max(64.d0 * epsilon(1.d0) * max(1.d0, abs(initial_cost), abs(final_cost)), &
            &JOINT_IMPROVE_REL_TOL * abs(initial_cost))
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
        ! a materially negative loss invalidates the solve: exponentiating it
        ! would mint a score > 1 that corrupts probability-table distances and
        ! downstream weighting, and the pose itself is a descent into series
        ! artifact, not signal
        if( valid_result .and. final_cost < -JOINT_NEG_COST_TOL ) valid_result = .false.
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
            ! Return the selected seed cell and its score; irot still holds the
            ! global scan result or the caller's authoritative local seed.
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
        ! clamp benign sub-tolerance series undershoot so the score keeps the
        ! legacy [0,1] normalization contract
        cxy(1) = real(exp(-max(0.d0, final_cost)))
        cxy(2:) = self%ospec%x(1:2)
        rotind_frac = self%cur_inpl_rotind
        if( sh_rot )then
            call rotmat2d(grad_shsrch_get_angle(self), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
        endif
    end function grad_shsrch_minimize_joint

    !> Jointly profile shift and continuous in-plane rotation while returning
    !! only the canonical integer rotation cell. Probabilistic searches use
    !! this form so their existing assignment artifacts remain integer based;
    !! the selected assignment is refined once more when durable e3 metadata
    !! is written.
    function grad_shsrch_minimize_joint_rounded( self, irot, sh_rot, xy_in, evaluation_valid ) result(cxy)
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(inout) :: irot
        logical,                 intent(in)    :: sh_rot
        real,          optional, intent(in)    :: xy_in(2)
        logical,       optional, intent(out)   :: evaluation_valid
        real :: cxy(3), xy_seed(2), rotmat(2,2)
        real(dp) :: rotind_frac
        integer :: selected_irot
        logical :: valid_result

        if( self%search_mode /= SHSRCH_JOINT )then
            THROW_HARD('rounded joint minimization requested from a non-joint search object')
        endif
        xy_seed   = 0.
        if( present(xy_in) ) xy_seed = xy_in
        ! The probability artifact drops rotind_frac, so retain its shift in
        ! the frame of the rounded integer cell.  This makes the later inverse
        ! transform recover the exact native shift used by the joint score.
        cxy = self%minimize_joint(selected_irot, xy_seed, .false., rotind_frac, &
            &evaluation_valid=valid_result)
        irot = selected_irot
        if( present(evaluation_valid) ) evaluation_valid = valid_result
        if( irot > 0 .and. sh_rot )then
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
        endif
    end function grad_shsrch_minimize_joint_rounded

    !> Bounded direct-gradient minimization for the streaming SGD path.
    !! The PFTC routine supplies grad(C) for the correlation objective.  We
    !! minimize the equivalent loss L=-C, hence grad(L)=-grad(C), and update
    !! s_{t+1}=Pi_bounds[s_t-eta_s grad(L)].  A trial is committed only when
    !! its evaluated loss is strictly lower; otherwise the original state is
    !! retained.  This is the continuous shift part of Design A; class and
    !! angle remain a discrete argmax in the surrounding search strategy.
    !! The current candidate shift is the initial state.  Each normalized step
    !! is projected into the legal shift box and accepted only when it lowers
    !! the objective; otherwise a short backtracking line search is used.
    function grad_shsrch_minimize_direct( self, irot, xy_in, step_size, max_steps,&
            &sh_rot, accepted_steps, objective_initial, objective_final, raw_euclid ) result( cxy )
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(inout) :: irot
        real,                    intent(in)    :: xy_in(2), step_size
        integer,                 intent(in)    :: max_steps
        logical,                 intent(in)    :: sh_rot
        integer,       optional, intent(out)   :: accepted_steps
        ! Optional diagnostics expose the minimized cost C=-score without
        ! changing the default update path: accepted steps require
        ! C_{t+1}<C_t, while the original state is retained on rejection.
        real(dp),      optional, intent(out)   :: objective_initial, objective_final
        logical,                 intent(in)    :: raw_euclid
        real :: cxy(3), rotmat(2,2)
        real(dp) :: current_xy(2), trial_xy(2), grad(2), corr_grad(2)
        real(dp) :: current_corr, current_cost, initial_cost, trial_corr, trial_cost
        real(dp) :: alpha, improve_tol
        real(sp), allocatable :: vector_losses(:)
        real(dp) :: vector_loss, roundtrip_tol
        integer :: istep, iback, naccepted
        logical :: accepted

        if( self%search_mode /= SHSRCH_DIRECT )then
            THROW_HARD('direct minimization requested from a non-direct search object')
        endif
        if( irot < 1 .or. irot > self%nrots )then
            THROW_HARD('direct shift minimizer received invalid rotation')
        endif
        if( step_size <= 0. )then
            THROW_HARD('direct shift minimizer step_size must be > 0')
        endif
        if( max_steps < 1 )then
            THROW_HARD('direct shift minimizer max_steps must be >= 1')
        endif
        naccepted = 0
        if( present(accepted_steps) ) accepted_steps = 0
        if( present(objective_initial) ) objective_initial = 0.d0
        if( present(objective_final) ) objective_final = 0.d0
        cxy = 0.
        self%cur_inpl_idx = irot
        current_xy = real(xy_in,dp)
        self%profile_objective_evals = self%profile_objective_evals + 1_int64
        if( raw_euclid )then
            call self%b_ptr%pftc%gen_raw_euclid_grad_for_rot_8(&
                &self%reference, self%particle, current_xy, self%cur_inpl_idx, current_corr, grad)
            if( self%raw_roundtrip_check )then
                ! Diagnostic-only invariant: the vector candidate loss and the
                ! scalar gradient loss must agree at the identical state.
                allocate(vector_losses(self%nrots))
                call self%b_ptr%pftc%gen_raw_euclid_vals(self%reference, self%particle, &
                    &real(current_xy,sp), vector_losses)
                vector_loss = real(vector_losses(self%cur_inpl_idx),dp)
                roundtrip_tol = 1.d-5 * (1.d0 + abs(current_corr))
                if( ieee_is_finite(vector_loss) )then
                    ! The legacy FFT score and the direct analytical loss have
                    ! the same shift-dependent landscape but differ by a
                    ! constant baseline in this SIMPLE representation.  The
                    ! diagnostic therefore calibrates that baseline once and
                    ! checks that it remains constant, rather than requiring
                    ! two independently normalized implementations to match
                    ! in absolute value.
                    if( .not. self%raw_roundtrip_offset_set )then
                        self%raw_roundtrip_offset = vector_loss - current_corr
                        self%raw_roundtrip_offset_set = .true.
                    endif
                    if( abs((vector_loss-current_corr)-self%raw_roundtrip_offset) > roundtrip_tol )then
                        write(*,'(A,2I8,5ES20.10)') &
                            'RAW ROUNDTRIP offset mismatch (ref,rot,xy0,xy1,delta,expected): ', &
                            self%reference, self%cur_inpl_idx, current_xy(1), current_xy(2), &
                            vector_loss-current_corr, self%raw_roundtrip_offset
                        self%raw_roundtrip_failed = .true.
                    endif
                else
                    self%raw_roundtrip_failed = .true.
                endif
                deallocate(vector_losses)
            endif
        else
            current_corr = self%b_ptr%pftc%gen_corr_for_rot_8(&
                &self%reference, self%particle, current_xy, self%cur_inpl_idx)
        endif
        if( .not. ieee_is_finite(current_corr) )then
            irot = 0
            return
        endif
        if( raw_euclid )then
            current_cost = current_corr
        else
            current_cost = -current_corr
        endif
        initial_cost = current_cost
        if( present(objective_initial) ) objective_initial = initial_cost

        do istep = 1, max_steps
            self%profile_objective_evals = self%profile_objective_evals + 1_int64
            self%profile_gradient_evals  = self%profile_gradient_evals  + 1_int64
            if( raw_euclid )then
                call self%b_ptr%pftc%gen_raw_euclid_grad_for_rot_8(&
                    &self%reference, self%particle, current_xy, self%cur_inpl_idx, current_corr, grad)
            else
                call self%b_ptr%pftc%gen_corr_grad_for_rot_8(self%reference, self%particle,&
                    &current_xy, self%cur_inpl_idx, current_corr, corr_grad)
                grad = -corr_grad
            endif
            if( .not. ieee_is_finite(current_corr) .or. any(.not. ieee_is_finite(grad)) )then
                exit
            endif
            if( sqrt(sum(grad * grad)) <= real(DTINY,dp) )then
                exit
            endif
            if( raw_euclid )then
                current_cost = current_corr
            else
                current_cost = -current_corr
            endif
            alpha = real(step_size,dp)
            accepted = .false.
            do iback = 0, direct_max_backtracks - 1
                trial_xy = bounded_shift_trial(current_xy, grad, alpha, real(self%ospec%limits,dp))
                if( maxval(abs(trial_xy - current_xy)) <= real(DTINY,dp) )then
                    alpha = 0.5_dp * alpha
                    cycle
                endif
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                if( raw_euclid )then
                    call self%b_ptr%pftc%gen_raw_euclid_grad_for_rot_8(&
                        &self%reference, self%particle, trial_xy, self%cur_inpl_idx, trial_corr, corr_grad)
                else
                    trial_corr = self%b_ptr%pftc%gen_corr_for_rot_8(&
                        &self%reference, self%particle, trial_xy, self%cur_inpl_idx)
                endif
                if( ieee_is_finite(trial_corr) )then
                    if( raw_euclid )then
                        trial_cost = trial_corr
                    else
                        trial_cost = -trial_corr
                    endif
                    improve_tol = 64.0_dp * epsilon(1.0_dp) *&
                        &max(1.0_dp, abs(current_cost), abs(trial_cost))
                    if( trial_cost < current_cost - improve_tol )then
                        current_xy   = trial_xy
                        current_cost = trial_cost
                        naccepted    = naccepted + 1
                        accepted     = .true.
                        exit
                    endif
                endif
                alpha = 0.5_dp * alpha
            end do
            if( .not. accepted )then
                exit
            endif
        end do

        if( present(accepted_steps) ) accepted_steps = naccepted
        if( present(objective_final) ) objective_final = current_cost
        improve_tol = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(initial_cost), abs(current_cost))
        if( naccepted < 1 .or. current_cost >= initial_cost - improve_tol )then
            irot = 0
            return
        endif
        cxy(1)  = -real(current_cost)
        cxy(2:) = real(current_xy)
        if( sh_rot )then
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            cxy(2:) = matmul(cxy(2:), rotmat)
        endif
    end function grad_shsrch_minimize_direct

    pure function bounded_shift_trial( xy, gradient, step_size, limits ) result( trial )
        real(dp), intent(in) :: xy(2), gradient(2), step_size, limits(2,2)
        real(dp) :: trial(2), grad_norm
        grad_norm = sqrt(sum(gradient * gradient))
        trial = xy
        if( step_size <= 0.0_dp .or. grad_norm <= real(DTINY,dp) ) return
        trial = xy - step_size * gradient / grad_norm
        trial = max(limits(:,1), min(limits(:,2), trial))
    end function bounded_shift_trial

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
                cost = -self%b_ptr%pftc%gen_corr_for_rot_8(self%reference, self%particle, [x,y], self%cur_inpl_idx)
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
        if (coarse_num_steps .le. 1) return
        lowest_cost = huge(lowest_cost)
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            x = self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,coarse_num_steps
                y = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                self%profile_objective_evals = self%profile_objective_evals + int(self%nrots,int64)
                call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, real([x,y]), corrs)
                loc  = maxloc(corrs,dim=1)
                cost = - corrs(loc)
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    irot        = loc
                    init_xy(1)  = x
                    init_xy(2)  = y
                end if
            end do
        end do
    end subroutine coarse_search_opt_angle

    subroutine grad_shsrch_reset_profile( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        self%profile_objective_evals = 0_int64
        self%profile_gradient_evals  = 0_int64
    end subroutine grad_shsrch_reset_profile

    subroutine grad_shsrch_get_profile( self, objective_evals, gradient_evals )
        class(pftc_shsrch_grad), intent(in) :: self
        integer(int64), intent(out) :: objective_evals, gradient_evals
        objective_evals = self%profile_objective_evals
        gradient_evals  = self%profile_gradient_evals
    end subroutine grad_shsrch_get_profile

    subroutine grad_shsrch_kill( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        if( associated(self%opt_obj) )then
            call self%opt_obj%kill
            nullify(self%opt_obj)
        end if
        call self%ospec%kill
        nullify(self%b_ptr)
        self%search_mode = SHSRCH_LEGACY
        self%joint_initial_cost = 0.d0
        self%joint_initial_cost_valid = .false.
        call self%reset_profile
    end subroutine grad_shsrch_kill

end module simple_pftc_shsrch_grad
