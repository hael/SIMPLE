!@descr: prior terms, the observability test and the LM problems of Cartesian pose refinement
submodule (simple_cartft_pose_opt) simple_cartft_pose_opt_system
implicit none

! lower edges of RATIO_LABELS' bins 2..NRATIO_BINS
real(dp), parameter :: RATIO_EDGES(NRATIO_BINS-1) = [0.75_dp, 1.25_dp, 1.75_dp, 2.5_dp]

contains

    !> Combine enabled priors once, in physical pose coordinates and normalized loss units.
    pure module subroutine pose_prior_terms( shift_weight, rotation_weight, rotation, shift, seed_rotation, seed_shift, &
        &penalty, gradient, hessian )
        real(dp), intent(in)  :: shift_weight, rotation_weight, rotation(3,3), shift(2), seed_rotation(3,3), seed_shift(2)
        real(dp), intent(out) :: penalty, gradient(5), hessian(5,5)
        real(dp) :: shift_penalty, rotation_penalty, sg(2), sh(2,2), rg(3), rh(3,3)
        gradient = 0._dp
        hessian = 0._dp
        shift_penalty = 0._dp
        rotation_penalty = 0._dp
        if( shift_weight > 0._dp )then
            call shift_prior_terms(shift_weight, shift, seed_shift, shift_penalty, sg, sh)
            gradient(4:5) = sg
            hessian(4:5,4:5) = sh
        endif
        if( rotation_weight > 0._dp )then
            call rotation_prior_terms(rotation_weight, rotation, seed_rotation, rotation_penalty, rg, rh)
            gradient(1:3) = rg
            hessian(1:3,1:3) = rh
        endif
        penalty = shift_penalty + rotation_penalty
    end subroutine pose_prior_terms

    !> Matrix-Fisher residual with right-increment gradient and Gauss-Newton curvature.
    pure module subroutine rotation_prior_terms( weight, rotation, seed, penalty, gradient, hessian )
        real(dp), intent(in)  :: weight, rotation(3,3), seed(3,3)
        real(dp), intent(out) :: penalty, gradient(3), hessian(3,3)
        real(dp) :: residual(3,3), jacobian(3,3,3), generators(3,3,3)
        integer :: i, j
        generators = 0._dp
        generators(3,2,1) = 1._dp
        generators(2,3,1) = -1._dp
        generators(1,3,2) = 1._dp
        generators(3,1,2) = -1._dp
        generators(2,1,3) = 1._dp
        generators(1,2,3) = -1._dp
        residual = rotation - seed
        penalty = 0.5_dp*weight*sum(residual**2)
        do i = 1, 3
            jacobian(:,:,i) = matmul(rotation, generators(:,:,i))
            gradient(i) = weight*sum(residual*jacobian(:,:,i))
        end do
        do j = 1, 3
            do i = 1, 3
                hessian(i,j) = weight*sum(jacobian(:,:,i)*jacobian(:,:,j))
            end do
        end do
    end subroutine rotation_prior_terms

    !> Gaussian shift penalty in normalized data-loss units, anchored to the transaction seed.
    pure module subroutine shift_prior_terms( weight, shift, seed, penalty, gradient, hessian )
        real(dp), intent(in)  :: weight, shift(2), seed(2)
        real(dp), intent(out) :: penalty, gradient(2), hessian(2,2)
        real(dp) :: delta(2)
        delta = shift - seed
        penalty = weight*dot_product(delta, delta)
        gradient = 2._dp*weight*delta
        hessian = 0._dp
        hessian(1,1) = 2._dp*weight
        hessian(2,2) = 2._dp*weight
    end subroutine shift_prior_terms

    !> The data-only observability test of the joint solve: the scaled, masked system factors with the
    !! relative pivot floor. It sees no prior, so a prior can never rescue an unobservable pose.
    pure module function pose_data_observable( scaled_gradient, scaled_hessian ) result( observable )
        use simple_linalg, only: spd_solve_rel
        real(dp), intent(in) :: scaled_gradient(:), scaled_hessian(:,:)
        logical :: observable
        real(dp) :: ignored_step(size(scaled_gradient))
        call spd_solve_rel(scaled_hessian, -scaled_gradient, ignored_step, observable)
    end function pose_data_observable

    ! THE JOINT SOLVE (refine_joint)

    !> objective, gradient and Gauss-Newton matrix at the trial pose, and its shell residual
    module subroutine pose_problem_eval( self, f, g, curv )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(out)   :: f, g(:), curv(:,:)
        call self%calc%pose_normal_terms(self%state, self%iseven, self%iptcl, self%trial_rotmat, self%trial_shift, &
            &f, g, curv, self%trial_residual)
        self%npasses = self%npasses + 1
    end subroutine pose_problem_eval

    !> the rotation and shift norm caps of one step
    module subroutine pose_problem_limit_step( self, step, bounded )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(inout) :: step(:)
        logical,                    intent(out)   :: bounded
        real(dp) :: rotation_norm, shift_norm
        bounded       = .false.
        rotation_norm = sqrt(dot_product(step(1:3), step(1:3)))
        if( rotation_norm > self%rotation_scale )then
            step(1:3) = step(1:3)*(self%rotation_scale/rotation_norm)
            bounded = .true.
        endif
        shift_norm = sqrt(dot_product(step(4:5), step(4:5)))
        if( shift_norm > self%shift_step_bound )then
            step(4:5) = step(4:5)*(self%shift_step_bound/shift_norm)
            bounded = .true.
        endif
        if( all(ieee_is_finite(step)) )then
            self%max_rotation_step = max(self%max_rotation_step, sqrt(dot_product(step(1:3), step(1:3))))
            self%max_shift_step    = max(self%max_shift_step,    sqrt(dot_product(step(4:5), step(4:5))))
        endif
    end subroutine pose_problem_limit_step

    !> the right SO(3) increment and the shift increment; infeasible outside the cumulative guard
    module function pose_problem_propose( self, step ) result( feasible )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: step(:)
        logical :: feasible
        real(dp) :: sine_half, cumulative_rotation, cumulative_shift
        self%trial_rotmat = right_increment_rotation(self%rotmat, step(1:3))
        self%trial_shift  = self%shift + step(4:5)
        feasible = .true.
        if( .not. self%cumulative_guard ) return
        ! the capture basin of the anchor
        sine_half           = sqrt(sum((self%trial_rotmat - self%anchor_rotmat)**2))/(2._dp*sqrt(2._dp))
        cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
        cumulative_shift    = sqrt(sum((self%trial_shift - self%anchor_shift)**2))
        feasible = .not. (cumulative_rotation > self%max_total_rotation + 10._dp*epsilon(1._dp) .or. &
            &cumulative_shift > self%max_total_shift + 10._dp*epsilon(1._dp))
    end function pose_problem_propose

    !> adopt the trial pose and its shell residual; bin the ratio for the bench
    module function pose_problem_accept( self, ratio ) result( changed )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: ratio
        logical :: changed
        self%nratio(ratio_bin(ratio)) = self%nratio(ratio_bin(ratio)) + 1
        self%rotmat = self%trial_rotmat
        self%shift  = self%trial_shift
        if( allocated(self%residual) ) self%residual = self%trial_residual
        changed = .false.
    end function pose_problem_accept

    !> the priors anchored at the transaction seed, at the current or the trial pose
    module function pose_problem_penalty( self, at_trial, penalty, g, curv ) result( has_penalty )
        class(cartft_pose_problem), intent(inout) :: self
        logical,                    intent(in)    :: at_trial
        real(dp),                   intent(out)   :: penalty, g(:), curv(:,:)
        logical :: has_penalty
        has_penalty = self%shift_prior_weight > 0._dp .or. self%rotation_prior_weight > 0._dp
        if( .not. has_penalty )then
            penalty = 0._dp
            g       = 0._dp
            curv    = 0._dp
        else if( at_trial )then
            call pose_prior_terms(self%shift_prior_weight, self%rotation_prior_weight, self%trial_rotmat, &
                &self%trial_shift, self%anchor_rotmat, self%anchor_shift, penalty, g, curv)
        else
            call pose_prior_terms(self%shift_prior_weight, self%rotation_prior_weight, self%rotmat, self%shift, &
                &self%anchor_rotmat, self%anchor_shift, penalty, g, curv)
        endif
    end function pose_problem_penalty

    module function pose_problem_observable( self, scaled_g, scaled_curv ) result( observable )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: scaled_g(:), scaled_curv(:,:)
        logical :: observable
        observable = pose_data_observable(scaled_g, scaled_curv)
    end function pose_problem_observable

    !> neither the rotation (radians) nor the shift (cropped pixels) of an accepted step resolves
    module function pose_problem_negligible_step( self, step ) result( negligible )
        class(cartft_pose_problem), intent(inout) :: self
        real(dp),                   intent(in)    :: step(:)
        logical :: negligible
        negligible = max(sqrt(dot_product(step(1:3), step(1:3))), sqrt(dot_product(step(4:5), step(4:5)))) &
            &< NEGLIGIBLE_STEP
    end function pose_problem_negligible_step

    ! THE SHIFT-ONLY SOLVE (refine_shift)

    module subroutine shift_problem_eval( self, f, g, curv )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(out)   :: f, g(:), curv(:,:)
        call self%calc%shift_normal_terms(self%state, self%iseven, self%iptcl, self%rotmat, self%trial_shift, f, g, curv)
        self%npasses = self%npasses + 1
    end subroutine shift_problem_eval

    module subroutine shift_problem_limit_step( self, step, bounded )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(inout) :: step(:)
        logical,                     intent(out)   :: bounded
        real(dp) :: step_norm
        step_norm = sqrt(dot_product(step, step))
        bounded   = step_norm > self%shift_step_bound
        if( bounded ) step = step*(self%shift_step_bound/step_norm)
        if( all(ieee_is_finite(step)) ) self%max_shift_step = max(self%max_shift_step, sqrt(dot_product(step, step)))
    end subroutine shift_problem_limit_step

    !> infeasible beyond trs of the input shift
    module function shift_problem_propose( self, step ) result( feasible )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: step(:)
        logical :: feasible
        self%trial_shift = self%shift + step
        feasible = sqrt(sum((self%trial_shift - self%anchor_shift)**2)) <= self%max_total_shift + 10._dp*epsilon(1._dp)
    end function shift_problem_propose

    module function shift_problem_accept( self, ratio ) result( changed )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: ratio
        logical :: changed
        self%nratio(ratio_bin(ratio)) = self%nratio(ratio_bin(ratio)) + 1
        self%shift = self%trial_shift
        changed    = .false.
    end function shift_problem_accept

    module function shift_problem_negligible_step( self, step ) result( negligible )
        class(cartft_shift_problem), intent(inout) :: self
        real(dp),                    intent(in)    :: step(:)
        logical :: negligible
        negligible = sqrt(dot_product(step, step)) < NEGLIGIBLE_STEP
    end function shift_problem_negligible_step

    !> The bench histogram bin of an accepted step's actual over predicted reduction.
    pure integer function ratio_bin( ratio )
        real(dp), intent(in) :: ratio
        ratio_bin = 1 + count(ratio >= RATIO_EDGES)
    end function ratio_bin

end submodule simple_cartft_pose_opt_system
