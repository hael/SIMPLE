!@descr: standalone Cartesian local-pose strategy for already aligned particles
module simple_strategy3D_pose_cont
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder, only: builder
use simple_cartesian_pose_refiner, only: cartesian_pose_data, LM_ACCEPTED_IMPROVEMENT
use simple_core_module_api, only: ctfparams, dp, simple_exception
use simple_image, only: image
use simple_linalg, only: arg
use simple_ori, only: ori
use simple_ori_utils, only: dm2euler
use simple_oris, only: oris
use simple_parameters, only: parameters
use simple_pose_cont_refine3D_adapter, only: pose_cont_config, pose_cont_limits, &
    &pose_cont_pose, pose_cont_reference_workspace, pose_cont_transaction_result, &
    &prepare_pose_cont_observation, shift_crop_to_native, shift_native_to_crop, &
    &POSE_CONT_ROUTE_JOINT, POSE_CONT_ROUTE_SHIFT_THEN_JOINT
use simple_strategy3D, only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
use simple_type_defs, only: OBJFUN_EUCLID
implicit none
private

#include "simple_local_flags.inc"

public :: strategy3D_pose_cont, pose_cont_seed_is_valid

!> Cartesian-only local search starting from an initialized ptcl3D pose.
!!
!! The inherited strategy3D_srch object is intentionally unused: initializing
!! it would construct PFTC search state and violate this strategy's ownership
!! contract. The matcher supplies its shared Cartesian reference workspace and
!! a thread-local particle image buffer through bind_context before srch.
type, extends(strategy3D) :: strategy3D_pose_cont
    private
    class(parameters), pointer :: p_ptr => null()
    class(builder), pointer :: b_ptr => null()
    class(pose_cont_reference_workspace), pointer :: refs_ptr => null()
    class(image), pointer :: work_img_ptr => null()
    integer :: iptcl_batch = 0
    type(pose_cont_config) :: config
    type(pose_cont_limits) :: limits
    type(ori) :: output_ori
    type(pose_cont_transaction_result) :: result
    logical :: exists = .false.
    logical :: context_bound = .false.
    logical :: searched = .false.
contains
    procedure :: new => new_pose_cont
    procedure :: bind_context => bind_pose_cont_context
    procedure :: srch => srch_pose_cont
    procedure :: oris_assign => oris_assign_pose_cont
    procedure :: kill => kill_pose_cont
end type strategy3D_pose_cont

contains

    !> Identify an initialized local-search seed without treating the valid
    !! identity rotation as an uninitialized sentinel.
    pure logical function pose_cont_seed_is_valid(seed) result(valid)
        class(ori), intent(in) :: seed
        real :: euler(3), shift(2)
        integer :: eo

        euler = seed%get_euler()
        shift = seed%get_2Dshift()
        eo = seed%get_eo()
        valid = seed%get_state() > 0 .and. seed%get_proj() > 0 .and. &
            &(eo == 0 .or. eo == 1) .and. &
            &all(ieee_is_finite(euler)) .and. all(ieee_is_finite(shift))
    end function pose_cont_seed_is_valid

    !> Initialize particle identity and local-LM policy without touching PFTC.
    subroutine new_pose_cont(self, params, spec, build)
        class(strategy3D_pose_cont), intent(inout) :: self
        class(parameters), intent(in) :: params
        class(strategy3D_spec), intent(inout) :: spec
        class(builder), intent(in) :: build
        real(dp) :: crop_scale

        call self%kill
        if (trim(params%oritype) /= 'ptcl3D') &
            &THROW_HARD('strategy3D_pose_cont requires oritype=ptcl3D')
        if (params%cc_objfun /= OBJFUN_EUCLID) &
            &THROW_HARD('strategy3D_pose_cont requires objfun=euclid')
        if (trim(params%inpl_cont) /= 'no') &
            &THROW_HARD('strategy3D_pose_cont cannot execute with inpl_cont=yes')
        if (.not. associated(build%spproj_field)) &
            &THROW_HARD('strategy3D_pose_cont requires an active ptcl3D project field')
        if (spec%iptcl < 1 .or. spec%iptcl > build%spproj_field%get_noris()) &
            &THROW_HARD('strategy3D_pose_cont particle index is outside ptcl3D')

        self%spec = spec
        select case (trim(params%pose_cont_route))
        case ('shift_then_joint')
            self%config%route = POSE_CONT_ROUTE_SHIFT_THEN_JOINT
        case ('joint')
            self%config%route = POSE_CONT_ROUTE_JOINT
        case default
            THROW_HARD('unsupported pose_cont_route in strategy3D_pose_cont')
        end select
        crop_scale = real(params%box_crop, dp)/real(params%box, dp)
        self%limits = pose_cont_limits(shift_step_bound=crop_scale, &
            &max_total_shift=5._dp*crop_scale)
        self%exists = .true.
    end subroutine new_pose_cont

    !> Supply matcher-owned immutable references and this particle's batch data.
    !! The image buffer must be unique to the calling OpenMP thread.
    subroutine bind_pose_cont_context(self, params, build, refs, work_img, iptcl_batch)
        class(strategy3D_pose_cont), intent(inout) :: self
        class(parameters), target, intent(in) :: params
        class(builder), target, intent(in) :: build
        class(pose_cont_reference_workspace), target, intent(in) :: refs
        class(image), target, intent(inout) :: work_img
        integer, intent(in) :: iptcl_batch
        integer :: ldim(3)

        if (.not. self%exists) THROW_HARD('strategy3D_pose_cont must be initialized before binding')
        if (iptcl_batch < 1 .or. .not. allocated(build%imgbatch)) &
            &THROW_HARD('strategy3D_pose_cont requires a loaded particle batch')
        if (iptcl_batch > size(build%imgbatch)) &
            &THROW_HARD('strategy3D_pose_cont batch index is outside the loaded batch')
        ldim = work_img%get_ldim()
        if (any(ldim /= [params%box_crop, params%box_crop, 1])) &
            &THROW_HARD('strategy3D_pose_cont received an incompatible image workspace')
        self%p_ptr => params
        self%b_ptr => build
        self%refs_ptr => refs
        self%work_img_ptr => work_img
        self%iptcl_batch = iptcl_batch
        self%context_bound = .true.
    end subroutine bind_pose_cont_context

    !> Refine one stored pose against the matching Cartesian state/half volume.
    subroutine srch_pose_cont(self, os, ithr)
        class(strategy3D_pose_cont), intent(inout) :: self
        class(oris), intent(inout) :: os
        integer, intent(in) :: ithr
        type(ctfparams) :: ctfparms, cropped_ctfparms
        type(cartesian_pose_data) :: data
        type(pose_cont_pose) :: seed
        type(pose_cont_pose) :: terminal_pose
        type(ori) :: input_ori
        complex, allocatable :: observed(:, :)
        real, allocatable :: sigma2(:), sigma_contrib(:)
        real :: euler(3), shift_native(2)
        integer :: eo, state
        logical :: even

        if (.not. self%exists .or. .not. self%context_bound) &
            &THROW_HARD('strategy3D_pose_cont requires initialized and bound context')
        if (ithr < 1) THROW_HARD('strategy3D_pose_cont received an invalid thread index')
        if (os%get_state(self%spec%iptcl) <= 0) then
            call os%reject(self%spec%iptcl)
            return
        end if
        ! Routing validates ptcl3D initialization once for the complete project.
        ! Zero Euler angles are a valid identity pose and cannot identify an
        ! uninitialized individual particle.

        ! Stage 1: preserve the stored pose as the transaction seed and rollback.
        call os%get_ori(self%spec%iptcl, input_ori)
        self%output_ori = input_ori
        state = input_ori%get_state()
        eo = input_ori%get_eo()
        select case (eo)
        case (0)
            even = .true.
        case (1)
            even = .false.
        case default
            THROW_HARD('strategy3D_pose_cont requires an even/odd assignment')
        end select
        if (.not. self%refs_ptr%is_ready(state, even)) &
            &THROW_HARD('strategy3D_pose_cont reference state/half is unavailable')

        ! Stage 2: prepare the already-loaded raw particle on the Cartesian grid.
        ctfparms = self%b_ptr%spproj%get_ctfparams('ptcl3D', self%spec%iptcl)
        call prepare_pose_cont_observation(self%b_ptr%imgbatch(self%iptcl_batch), &
            &self%b_ptr%lmsk, self%work_img_ptr, self%p_ptr%msk_crop, &
            &self%p_ptr%smpd_crop, ctfparms, observed, cropped_ctfparms)
        call build_particle_sigma(self, sigma2)
        call self%refs_ptr%prepare_particle(state, even, observed, cropped_ctfparms, &
            &sigma2, self%p_ptr%kfromto, data)

        ! Stage 3: run only the selected Cartesian LM route from the stored pose.
        seed%rotmat = real(input_ori%get_mat(), dp)
        seed%shift = real(shift_native_to_crop(input_ori%get_2Dshift(), &
            &self%p_ptr%box, self%p_ptr%box_crop), dp)
        call self%refs_ptr%refine_particle(state, even, seed, data, self%config, &
            &self%limits, self%result)
        terminal_pose = seed

        ! Stage 4: stage only an accepted improvement for project assignment.
        if (self%result%status == LM_ACCEPTED_IMPROVEMENT) then
            terminal_pose = self%result%pose
            euler = real(dm2euler(self%result%pose%rotmat))
            shift_native = shift_crop_to_native(real(self%result%pose%shift), &
                &self%p_ptr%box, self%p_ptr%box_crop)
            call self%output_ori%set_euler(euler)
            call self%output_ori%set_shift(shift_native)
        end if
        call self%refs_ptr%sigma_contribution(state, even, terminal_pose, data, sigma_contrib)
        call self%b_ptr%esig%set_particle_contribution(self%spec%iptcl, sigma_contrib)
        self%searched = .true.
        call self%oris_assign
        call input_ori%kill
    end subroutine srch_pose_cont

    !> Commit an accepted Cartesian pose and record seed-to-terminal motion.
    !! A rejected transaction keeps the stored seed and records zero motion.
    !! The Cartesian objective intentionally never replaces the PFTC `corr`.
    subroutine oris_assign_pose_cont(self)
        class(strategy3D_pose_cont), intent(inout) :: self
        type(ori) :: input_ori, symmetry_equivalent_ori
        real :: shift_increment(2), euler_distance, inplane_distance, mi_proj

        if (.not. self%searched) THROW_HARD('strategy3D_pose_cont has no search result to assign')
        call self%b_ptr%spproj_field%get_ori(self%spec%iptcl, input_ori)

        ! Stage 1: commit only a valid improvement; output_ori remains the
        ! original seed for every rejected or invalid transaction.
        if (self%result%status == LM_ACCEPTED_IMPROVEMENT) then
            call self%b_ptr%spproj_field%set_ori(self%spec%iptcl, self%output_ori)
        end if

        ! Stage 2: replace inherited search statistics with measurements of
        ! this local Cartesian transaction. State and projection identity are
        ! fixed, while the Euler/shift endpoint may change after acceptance.
        call self%b_ptr%pgrpsyms%sym_dists(input_ori, self%output_ori, &
            &symmetry_equivalent_ori, euler_distance, inplane_distance)
        shift_increment = self%output_ori%get_2Dshift() - input_ori%get_2Dshift()
        mi_proj = merge(1., 0., euler_distance <= self%p_ptr%angthres_mi_proj)
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'shincarg', arg(shift_increment))
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'dist', euler_distance)
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'dist_inpl', inplane_distance)
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'mi_proj', mi_proj)
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'mi_state', 1.)

        ! The strategy exhausts its one seed-centered local search domain; this
        ! is local coverage, not a claim that the global orientation space was
        ! scanned.
        call self%b_ptr%spproj_field%set(self%spec%iptcl, 'frac', 100.)
        call symmetry_equivalent_ori%kill
        call input_ori%kill
    end subroutine oris_assign_pose_cont

    !> Copy this particle's active sigma shells into a compact zero-based vector.
    subroutine build_particle_sigma(self, sigma2)
        class(strategy3D_pose_cont), intent(in) :: self
        real, allocatable, intent(out) :: sigma2(:)
        integer :: shell_from, shell_to

        if (.not. allocated(self%b_ptr%esig%sigma2_noise)) &
            &THROW_HARD('strategy3D_pose_cont requires allocated sigma2 noise')
        shell_from = self%p_ptr%kfromto(1)
        shell_to = self%p_ptr%kfromto(2)
        if (shell_from < lbound(self%b_ptr%esig%sigma2_noise, 1) .or. &
            &shell_to > ubound(self%b_ptr%esig%sigma2_noise, 1)) &
            &THROW_HARD('strategy3D_pose_cont shell range exceeds sigma2 noise')
        if (self%spec%iptcl < lbound(self%b_ptr%esig%sigma2_noise, 2) .or. &
            &self%spec%iptcl > ubound(self%b_ptr%esig%sigma2_noise, 2)) &
            &THROW_HARD('strategy3D_pose_cont particle index exceeds sigma2 noise')
        allocate (sigma2(0:shell_to), source=1.)
        sigma2(shell_from:shell_to) = &
            &self%b_ptr%esig%sigma2_noise(shell_from:shell_to, self%spec%iptcl)
    end subroutine build_particle_sigma

    !> Release local state while leaving matcher-owned shared data untouched.
    subroutine kill_pose_cont(self)
        class(strategy3D_pose_cont), intent(inout) :: self

        call self%output_ori%kill
        nullify (self%p_ptr, self%b_ptr, self%refs_ptr, self%work_img_ptr)
        nullify (self%spec%eulprob_obj_part)
        self%iptcl_batch = 0
        self%config = pose_cont_config()
        self%limits = pose_cont_limits()
        self%result = pose_cont_transaction_result()
        self%exists = .false.
        self%context_bound = .false.
        self%searched = .false.
    end subroutine kill_pose_cont

end module simple_strategy3D_pose_cont
