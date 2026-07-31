!@descr: rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftc_shsrch_grad
use iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_opt_spec,  only: opt_spec
use simple_optimizer, only: optimizer
use simple_builder,   only: builder
implicit none

public :: pftc_shsrch_grad, bounded_shift_trial, parabolic_peak_offset
private
#include "simple_local_flags.inc"

integer,  parameter :: coarse_num_steps = 5       ! no. of coarse search steps in x AND y (hence real no. is its square)
integer,  parameter :: direct_max_backtracks = 8

type :: pftc_shsrch_grad
    private
    type(opt_spec)            :: ospec                  !< optimizer specification object
    class(optimizer), pointer :: opt_obj => null()      !< optimizer object
    class(builder),   pointer :: b_ptr   => null()      !< pointer to polarft_calc object for cost function evaluations
    integer                   :: reference    = 0       !< reference pft
    integer                   :: particle     = 0       !< particle pft
    integer                   :: nrots        = 0       !< # rotations
    integer                   :: maxits       = 100     !< max # iterations
    logical                   :: shbarr       = .true.  !< shift barrier constraint or not
    integer                   :: cur_inpl_idx = 0       !< index of inplane angle for shift search
    real                      :: cur_inpl_ang = 0.      !< continuous angle in grid-index units
    integer                   :: max_evals    = 5       !< max # inplrot/shsrch cycles
    logical                   :: opt_angle    = .true.  !< optimise in-plane angle with callback flag
    logical                   :: coarse_init  = .false. !< whether to perform an intial coarse search over the range
    logical                   :: direct_only  = .false. !< skip allocating an L-BFGS-B object for direct-gradient use
    integer(int64)            :: profile_objective_evals = 0_int64
    integer(int64)            :: profile_gradient_evals  = 0_int64
contains
    procedure          :: new         => grad_shsrch_new
    procedure          :: set_indices => grad_shsrch_set_indices
    procedure          :: minimize    => grad_shsrch_minimize
    procedure          :: minimize_direct => grad_shsrch_minimize_direct
    procedure          :: kill        => grad_shsrch_kill
    procedure          :: does_opt_angle
    procedure          :: set_limits
    procedure          :: coarse_search
    procedure          :: coarse_search_opt_angle
    procedure          :: is_direct_shift_only
    procedure          :: reset_profile => grad_shsrch_reset_profile
    procedure          :: get_profile   => grad_shsrch_get_profile
end type pftc_shsrch_grad

contains

    pure real function parabolic_peak_offset( vals, j ) result(offset)
        real,    intent(in) :: vals(:)
        integer, intent(in) :: j
        integer :: jm, jp, n
        real    :: denom, scale
        n = size(vals)
        offset = 0.
        if( n < 3 .or. j < 1 .or. j > n ) return
        jm = modulo(j - 2, n) + 1
        jp = modulo(j,     n) + 1
        denom = vals(jm) - 2.*vals(j) + vals(jp)
        scale  = max(1., abs(vals(jm)), abs(vals(j)), abs(vals(jp)))
        if( denom >= -10.*epsilon(1.)*scale ) return
        offset = 0.5 * (vals(jm) - vals(jp)) / denom
        if( .not. ieee_is_finite(offset) .or. abs(offset) > 0.5 ) offset = 0.
    end function parabolic_peak_offset


    subroutine grad_shsrch_new( self, build, lims, lims_init, shbarrier, maxits, opt_angle, coarse_init, direct_only )
        use simple_opt_factory, only: opt_factory
        class(pftc_shsrch_grad),     intent(inout) :: self           !< instance
        class(builder),      target, intent(in)    :: build          !< builder object for pftc access 
        real,                        intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,              optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*),  optional, intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,           optional, intent(in)    :: maxits         !< maximum iterations
        logical,           optional, intent(in)    :: opt_angle      !< optimise in-plane angle with callback flag
        logical,           optional, intent(in)    :: coarse_init    !< coarse inital search
        logical,           optional, intent(in)    :: direct_only    !< configure for the bounded direct-gradient path only
        type(opt_factory) :: opt_fact
        call self%kill
        ! set pointer to pftc instance for cost function evaluations
        self%b_ptr => build
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%opt_angle = .true.
        if( present(opt_angle) ) self%opt_angle = opt_angle
        self%coarse_init = .false.
        if( present(coarse_init) ) self%coarse_init = coarse_init
        self%direct_only = .false.
        if( present(direct_only) ) self%direct_only = direct_only
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        if( .not. self%direct_only ) call opt_fact%new(self%ospec, self%opt_obj)
        ! get # rotations
        self%nrots = self%b_ptr%pftc%get_nrots()
        ! set costfun pointers
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        if( self%opt_angle ) self%ospec%opt_callback => grad_shsrch_optimize_angle_wrapper
    end subroutine grad_shsrch_new

    pure logical function does_opt_angle( self )
        class(pftc_shsrch_grad), intent(in) :: self
        does_opt_angle = self%opt_angle
    end function does_opt_angle

    pure logical function is_direct_shift_only( self )
        class(pftc_shsrch_grad), intent(in) :: self
        ! The streaming optimizer operates on a fixed discrete rotation and
        ! updates only s=(sx,sy); it must neither run the angle callback nor
        ! allocate/fall through to the legacy L-BFGS-B implementation.
        is_direct_shift_only = self%direct_only .and. (.not. self%opt_angle)
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
                cost = - self%b_ptr%pftc%gen_corr_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx)
            class default
                THROW_HARD('error in grad_shsrch_costfun: unknown type; grad_shsrch_costfun')
        end select
    end function grad_shsrch_costfun

    subroutine grad_shsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp)                :: corrs_grad(2)
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_gradient_evals = self%profile_gradient_evals + 1_int64
                call self%b_ptr%pftc%gen_corr_grad_only_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs_grad)
                grad = - corrs_grad
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
        real(dp)                :: corrs_grad(2)
        f    = 0.
        grad = 0.
        select type(self)
            class is (pftc_shsrch_grad)
                self%profile_objective_evals = self%profile_objective_evals + 1_int64
                self%profile_gradient_evals  = self%profile_gradient_evals  + 1_int64
                call self%b_ptr%pftc%gen_corr_grad_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
                f    = - corrs
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine grad_shsrch_fdfcostfun

    subroutine grad_shsrch_optimize_angle( self )
        class(pftc_shsrch_grad), intent(inout) :: self
        real :: objective(self%nrots), raw_losses(self%nrots)
        integer :: jmax
        if( self%b_ptr%pftc%is_euclid_objfun() )then
            call self%b_ptr%pftc%gen_raw_euclid_vals(self%reference, self%particle, &
                &real(self%ospec%x,sp), raw_losses)
            objective = -raw_losses
            jmax = maxloc(objective, dim=1)
            self%cur_inpl_ang = real(jmax) + parabolic_peak_offset(objective, jmax)
        else
            call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, self%ospec%x, objective)
            jmax = maxloc(objective, dim=1)
            self%cur_inpl_ang = real(jmax)
        endif
        self%cur_inpl_idx = modulo(nint(self%cur_inpl_ang) - 1, self%nrots) + 1
    end subroutine grad_shsrch_optimize_angle

    real function grad_shsrch_get_angle( self ) result(angle)
        class(pftc_shsrch_grad), intent(in) :: self
        integer :: iang
        iang = modulo(nint(self%cur_inpl_ang) - 1, self%nrots) + 1
        angle = self%b_ptr%pftc%get_rot(iang) + &
            (self%cur_inpl_ang - real(iang)) * self%b_ptr%pftc%get_dang()
    end function grad_shsrch_get_angle


    subroutine grad_shsrch_optimize_angle_wrapper( self )
        class(*), intent(inout) :: self
        select type(self)
            class is (pftc_shsrch_grad)
                call grad_shsrch_optimize_angle(self)
            class DEFAULT
                THROW_HARD('error in grad_shsrch_optimize_angle_wrapper: unknown type; simple_pftc_shsrch_grad')
        end select
    end subroutine grad_shsrch_optimize_angle_wrapper

    !> set indicies for shift search
    subroutine grad_shsrch_set_indices( self, ref, ptcl )
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine grad_shsrch_set_indices

    !> minimisation routine
    function grad_shsrch_minimize( self, irot, sh_rot, xy_in, theta ) result( cxy )
        class(pftc_shsrch_grad), intent(inout) :: self
        integer,                 intent(inout) :: irot
        logical, optional,       intent(in)    :: sh_rot
        real,    optional,       intent(in)    :: xy_in(2)
        real,    optional,       intent(out)   :: theta
        real     :: corrs(self%nrots), rotmat(2,2), cxy(3), lowest_shift(2), lowest_cost
        real(dp) :: init_xy(2), lowest_cost_overall, coarse_cost, initial_cost
        integer  :: loc, i, lowest_rot, init_rot
        logical  :: found_better, l_sh_rot, coarse_init_orig
        if( self%direct_only )then
            THROW_HARD('L-BFGS-B minimize requested from a direct-only shift search object')
        endif
        if( present(theta) ) theta = 0.
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
            call grad_shsrch_optimize_angle(self)
            lowest_cost_overall = -corrs(self%cur_inpl_idx)
            initial_cost        = lowest_cost_overall
            if( self%coarse_init )then
                call self%coarse_search_opt_angle(init_xy, init_rot)
                if( init_rot /= 0 )then
                    self%ospec%x_8      = init_xy
                    self%ospec%x        = real(init_xy)
                    self%cur_inpl_idx   = init_rot
                    self%cur_inpl_ang   = real(init_rot)
                endif
            end if
            ! shift search / in-plane rot update
            do i = 1,self%max_evals
                call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                loc = self%cur_inpl_idx
                call self%b_ptr%pftc%gen_objfun_vals(self%reference, self%particle, self%ospec%x, corrs)
                call grad_shsrch_optimize_angle(self)
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
                cxy(1)  = - real(lowest_cost_overall)  ! correlation
                cxy(2:) =   real(lowest_shift)         ! shift
                if( present(theta) ) theta = self%cur_inpl_ang
                if( l_sh_rot )then
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(grad_shsrch_get_angle(self), rotmat)
                    cxy(2:) = matmul(cxy(2:), rotmat)
                endif
            else
                irot = 0 ! to communicate that a better solution was not found
            endif
        else
            self%cur_inpl_idx   = irot
            self%cur_inpl_ang   = real(irot)
            self%profile_objective_evals = self%profile_objective_evals + 1_int64
            lowest_cost_overall = -self%b_ptr%pftc%gen_corr_for_rot_8(self%reference, self%particle, self%ospec%x_8, self%cur_inpl_idx)
            initial_cost        = lowest_cost_overall
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
                if( present(theta) ) theta = real(irot)
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
        logical,       optional, intent(in)    :: sh_rot
        integer,       optional, intent(out)   :: accepted_steps
        ! Optional diagnostics expose the minimized cost C=-score without
        ! changing the default update path: accepted steps require
        ! C_{t+1}<C_t, while the original state is retained on rejection.
        real(dp),      optional, intent(out)   :: objective_initial, objective_final
        logical,       optional, intent(in)    :: raw_euclid
        real :: cxy(3), rotmat(2,2)
        real(dp) :: current_xy(2), trial_xy(2), grad(2), corr_grad(2)
        real(dp) :: current_corr, current_cost, initial_cost, trial_corr, trial_cost
        real(dp) :: alpha, improve_tol
        integer :: istep, iback, naccepted
        logical :: accepted, l_sh_rot, l_raw_euclid

        if( self%opt_angle )then
            THROW_HARD('direct joint shift minimizer requires a fixed in-plane angle')
        endif
        if( irot < 1 .or. irot > self%nrots )then
            THROW_HARD('direct joint shift minimizer received invalid rotation')
        endif
        if( step_size <= 0. )then
            THROW_HARD('direct joint shift minimizer step_size must be > 0')
        endif
        if( max_steps < 1 )then
            THROW_HARD('direct joint shift minimizer max_steps must be >= 1')
        endif
        l_sh_rot = .true.
        if( present(sh_rot) ) l_sh_rot = sh_rot
        l_raw_euclid = .false.
        if( present(raw_euclid) ) l_raw_euclid = raw_euclid
        naccepted = 0
        if( present(accepted_steps) ) accepted_steps = 0
        if( present(objective_initial) ) objective_initial = 0.d0
        if( present(objective_final) ) objective_final = 0.d0
        cxy = 0.
        self%cur_inpl_idx = irot
        current_xy = real(xy_in,dp)
        self%profile_objective_evals = self%profile_objective_evals + 1_int64
        if( l_raw_euclid )then
            call self%b_ptr%pftc%gen_raw_euclid_grad_for_rot_8(&
                &self%reference, self%particle, current_xy, self%cur_inpl_idx, current_corr, grad)
        else
            current_corr = self%b_ptr%pftc%gen_corr_for_rot_8(&
                &self%reference, self%particle, current_xy, self%cur_inpl_idx)
        endif
        if( .not. ieee_is_finite(current_corr) )then
            irot = 0
            return
        endif
        if( l_raw_euclid )then
            current_cost = current_corr
        else
            current_cost = -current_corr
        endif
        initial_cost = current_cost
        if( present(objective_initial) ) objective_initial = initial_cost

        do istep = 1, max_steps
            self%profile_objective_evals = self%profile_objective_evals + 1_int64
            self%profile_gradient_evals  = self%profile_gradient_evals  + 1_int64
            if( l_raw_euclid )then
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
            if( l_raw_euclid )then
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
                if( l_raw_euclid )then
                    call self%b_ptr%pftc%gen_raw_euclid_grad_for_rot_8(&
                        &self%reference, self%particle, trial_xy, self%cur_inpl_idx, trial_corr, corr_grad)
                else
                    trial_corr = self%b_ptr%pftc%gen_corr_for_rot_8(&
                        &self%reference, self%particle, trial_xy, self%cur_inpl_idx)
                endif
                if( ieee_is_finite(trial_corr) )then
                    if( l_raw_euclid )then
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
        if( l_sh_rot )then
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
        self%direct_only = .false.
        call self%reset_profile
    end subroutine grad_shsrch_kill

end module simple_pftc_shsrch_grad
