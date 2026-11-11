!@descr: private prior terms and linear systems for Cartesian pose refinement
submodule (simple_cartft_pose_opt) simple_cartft_pose_opt_system
implicit none

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

    !> Scale and test the data-only pose system before constructing a damped proposal.
    pure module subroutine prepare_pose_data_system( gradient, hessian, rotation_scale, active, &
        &scaled_gradient, scaled_hessian, identifiable )
        real(dp), intent(in)  :: gradient(5), hessian(5,5), rotation_scale
        logical,  intent(in)  :: active(5)
        real(dp), intent(out) :: scaled_gradient(5), scaled_hessian(5,5)
        logical,  intent(out) :: identifiable
        real(dp) :: ignored_step(5)
        call scale_pose_system(gradient, hessian, rotation_scale, active, scaled_gradient, scaled_hessian)
        call solve_pose_cholesky(scaled_hessian, -scaled_gradient, ignored_step, identifiable)
    end subroutine prepare_pose_data_system

    pure module subroutine scale_pose_system( gradient, hessian, rotation_scale, active, scaled_gradient, scaled_hessian )
        real(dp), intent(in)  :: gradient(5), hessian(5,5), rotation_scale
        logical, intent(in)   :: active(5)
        real(dp), intent(out) :: scaled_gradient(5), scaled_hessian(5,5)
        real(dp) :: coordinate_scale(5)
        integer :: axis, jaxis
        coordinate_scale = [rotation_scale, rotation_scale, rotation_scale, 1._dp, 1._dp]
        do axis = 1, 5
            scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
            do jaxis = 1, 5
                scaled_hessian(axis,jaxis) = coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
            end do
        end do
        call apply_pose_parameter_mask(scaled_gradient, scaled_hessian, active)
    end subroutine scale_pose_system

    !> One scaled, damped and independently bounded pose proposal.
    pure module subroutine build_pose_lm_system( scaled_gradient, scaled_hessian, rotation_scale, mu, &
        &damping_diagonal, solve_matrix, scaled_step, physical_step, stationary, reliable, bounded, shift_step_bound )
        real(dp), intent(in)  :: scaled_gradient(5), scaled_hessian(5,5), rotation_scale, mu
        real(dp), intent(out) :: damping_diagonal(5), solve_matrix(5,5)
        real(dp), intent(out) :: scaled_step(5), physical_step(5)
        logical,  intent(out) :: stationary, reliable, bounded
        real(dp), intent(in)  :: shift_step_bound
        real(dp) :: coordinate_scale(5), hessian_scale
        real(dp) :: rotation_norm, shift_norm
        integer  :: axis
        coordinate_scale = [rotation_scale, rotation_scale, rotation_scale, 1._dp, 1._dp]
        stationary       = sqrt(dot_product(scaled_gradient, scaled_gradient)) < 1.e-8_dp
        damping_diagonal = 0._dp
        solve_matrix     = scaled_hessian
        scaled_step      = 0._dp
        physical_step    = 0._dp
        reliable         = .false.
        bounded          = .false.
        if( stationary ) return
        hessian_scale = max(maxval(abs(scaled_hessian)), NUMERIC_FLOOR)
        do axis = 1, 5
            damping_diagonal(axis) = max(scaled_hessian(axis,axis), sqrt(epsilon(1._dp))*hessian_scale, NUMERIC_FLOOR)
            solve_matrix(axis,axis) = solve_matrix(axis,axis) + mu*damping_diagonal(axis)
        end do
        call solve_pose_cholesky(solve_matrix, -scaled_gradient, scaled_step, reliable)
        if( .not. reliable ) return
        physical_step = coordinate_scale*scaled_step
        rotation_norm = sqrt(dot_product(physical_step(1:3), physical_step(1:3)))
        if( rotation_norm > rotation_scale )then
            physical_step(1:3) = physical_step(1:3)*(rotation_scale/rotation_norm)
            bounded = .true.
        endif
        shift_norm = sqrt(dot_product(physical_step(4:5), physical_step(4:5)))
        if( shift_norm > shift_step_bound )then
            physical_step(4:5) = physical_step(4:5)*(shift_step_bound/shift_norm)
            bounded = .true.
        endif
    end subroutine build_pose_lm_system

    !> Freeze inactive pose coordinates while retaining one five-vector LM path.
    pure module subroutine apply_pose_parameter_mask( gradient, hessian, active )
        real(dp), intent(inout) :: gradient(5), hessian(5,5)
        logical,  intent(in)    :: active(5)
        integer :: axis
        do axis = 1, 5
            if( active(axis) ) cycle
            gradient(axis)     = 0._dp
            hessian(axis,:)    = 0._dp
            hessian(:,axis)    = 0._dp
            hessian(axis,axis) = 1._dp
        end do
    end subroutine apply_pose_parameter_mask

    !> Cholesky solve with a relative pivot test for a 5x5 symmetric positive-definite block.
    pure module subroutine solve_pose_cholesky( matrix, rhs, solution, reliable )
        real(dp), intent(in)  :: matrix(5,5), rhs(5)
        real(dp), intent(out) :: solution(5)
        logical,  intent(out) :: reliable
        real(dp) :: lower(5,5), intermediate(5), pivot, pivot_floor, matrix_scale
        integer  :: i, j
        solution     = 0._dp
        intermediate = 0._dp
        lower        = 0._dp
        reliable     = .false.
        if( any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs)) ) return
        matrix_scale = maxval(abs(matrix))
        if( matrix_scale <= NUMERIC_FLOOR ) return
        pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
        do i = 1, 5
            do j = 1, i - 1
                lower(i,j) = (matrix(i,j) - dot_product(lower(i,1:j-1), lower(j,1:j-1)))/lower(j,j)
            end do
            pivot = matrix(i,i) - dot_product(lower(i,1:i-1), lower(i,1:i-1))
            if( .not. ieee_is_finite(pivot) .or. pivot <= pivot_floor ) return
            lower(i,i) = sqrt(pivot)
        end do
        do i = 1, 5
            intermediate(i) = (rhs(i) - dot_product(lower(i,1:i-1), intermediate(1:i-1)))/lower(i,i)
        end do
        do i = 5, 1, -1
            solution(i) = (intermediate(i) - dot_product(lower(i+1:5,i), solution(i+1:5)))/lower(i,i)
        end do
        reliable = all(ieee_is_finite(solution))
    end subroutine solve_pose_cholesky


end submodule simple_cartft_pose_opt_system
