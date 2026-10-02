!@descr: Cartesian continuous-pose refinement of one particle (refine=cont), the 3D strategy of a Cartesian pass
! The Cartesian child of strategy3D (plan section 6.3): it extends the representation-neutral
! base directly and carries no polar search state. One particle per object: seed from the stored
! ptcl3D pose, one transaction of the pose optimizer (cartft_pose_opt) against the particle's
! state and half reference in build%cftc, on the particle slot the batch preparation filled
! (prep_cart_batch; the observation is centred on the stored shift, so the solve is for the shift
! increment). oris_assign commits pose and score together: the pose only when the solve was
! accepted, corr_cart (the Cartesian score, cc or exp(-L)) always, at the seed when the solve was
! rejected (C14); corr is never written. The improved flag records acceptance; a Cartesian pass
! attempts every particle it samples, so attempted needs no field (C15, O8). Under objfun=euclid
! the sigma owner records the residual at the committed pose (C5). The convergence fields (dist,
! dist_inpl, shincarg, mi_proj, mi_state, frac) are written from the seed-to-result motion.
module simple_strategy3D_cont
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder,          only: builder
use simple_core_module_api,  only: dp, logfhandle, simple_exception
use simple_linalg,           only: arg
use simple_ori,              only: ori
use simple_oris,             only: oris
use simple_parameters,       only: parameters
use simple_strategy3D,       only: strategy3D
use simple_strategy3D_srch,  only: strategy3D_spec
use simple_type_defs,        only: OBJFUN_CC, OBJFUN_EUCLID
use simple_cartft_pose_opt,  only: cartft_pose_opt, CARTFT_ACCEPTED, CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION
implicit none

public :: strategy3D_cont, cont_seed_is_valid, check_cont_seeds
private

#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_cont
    private
    class(parameters), pointer :: p_ptr => null()
    class(builder),    pointer :: b_ptr => null()
    type(cartft_pose_opt) :: opt                    !< the optimizer of this particle's search
    type(ori)             :: seed_ori               !< the stored pose, the transaction seed
    type(ori)             :: result_ori             !< the seed, or the accepted pose
    integer               :: status = CARTFT_NOT_ATTEMPTED
    real                  :: score  = 0.            !< corr_cart: at the result, at the seed when rejected
    logical               :: searched = .false.
    logical               :: exists   = .false.
contains
    procedure :: new         => new_cont
    procedure :: srch        => srch_cont
    procedure :: oris_assign => oris_assign_cont
    procedure :: kill        => kill_cont
    procedure, private :: bind
end type strategy3D_cont

contains

    !> An initialized seed of a Cartesian pass: positive state, a half-set, a projection
    !! direction from the discrete workflow the poses came from, finite Euler angles and shift.
    !! The identity rotation is a valid pose. The caller turns false into a stop.
    pure logical function cont_seed_is_valid( seed ) result( valid )
        class(ori), intent(in) :: seed
        real    :: euler(3), shift(2)
        integer :: eo
        euler = seed%get_euler()
        shift = seed%get_2Dshift()
        eo    = seed%get_eo()
        valid = seed%get_state() > 0 .and. seed%get_proj() > 0 .and. (eo == 0 .or. eo == 1) .and. &
            &all(ieee_is_finite(euler)) .and. all(ieee_is_finite(shift))
    end function cont_seed_is_valid

    !> The entry check of a Cartesian pass (C11): every active particle (state > 0) of the field
    !! holds a valid seed (cont_seed_is_valid), so the pass continues from poses a discrete
    !! workflow left in the project. THROW_HARD naming the requirement otherwise, or when no
    !! particle is active.
    subroutine check_cont_seeds( os )
        class(oris), intent(in) :: os
        type(ori) :: seed
        integer   :: iptcl, nactive, ninvalid
        nactive  = 0
        ninvalid = 0
        do iptcl = 1, os%get_noris()
            if( os%get_state(iptcl) <= 0 ) cycle
            nactive = nactive + 1
            call os%get_ori(iptcl, seed)
            if( .not. cont_seed_is_valid(seed) ) ninvalid = ninvalid + 1
        end do
        call seed%kill
        if( nactive == 0 ) THROW_HARD('refine=cont found no active particle in ptcl3D')
        if( ninvalid > 0 )then
            write(logfhandle,'(a,i0,a,i0)') 'particles without a valid 3D pose: ', ninvalid, ' of ', nactive
            THROW_HARD('refine=cont continues from the 3D poses of a discrete workflow (state, half-set, projection direction, finite angles and shifts); run one first')
        endif
    end subroutine check_cont_seeds

    !> Particle identity and the optimizer of the pass's policy: the total bounds trs and
    !! athres_cont (O6), the route of cont_route (C16; joint in production).
    subroutine new_cont( self, params, spec, build )
        class(strategy3D_cont), intent(inout) :: self
        class(parameters),      intent(in)    :: params
        class(strategy3D_spec), intent(inout) :: spec
        class(builder),         intent(in)    :: build
        call self%kill
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('strategy3D_cont requires oritype=ptcl3D')
        if( params%cc_objfun /= OBJFUN_EUCLID .and. params%cc_objfun /= OBJFUN_CC ) &
            &THROW_HARD('strategy3D_cont supports only objfun=euclid or objfun=cc')
        if( trim(params%inpl_cont) /= 'no' ) THROW_HARD('strategy3D_cont cannot execute with inpl_cont=yes')
        if( trim(params%projrec) == 'yes' ) THROW_HARD('strategy3D_cont does not support projrec=yes')
        if( .not. associated(build%spproj_field) ) THROW_HARD('strategy3D_cont requires an active ptcl3D project field')
        if( spec%iptcl < 1 .or. spec%iptcl > build%spproj_field%get_noris() ) &
            &THROW_HARD('strategy3D_cont particle index is outside ptcl3D')
        self%spec = spec
        call self%bind(params, build)
        call self%opt%new(params%box, params%box_crop, params%trs, params%athres_cont, shift_first=params%l_cont_shift_first)
        self%exists = .true.
    end subroutine new_cont

    subroutine bind( self, params, build )
        class(strategy3D_cont),    intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        self%p_ptr => params
        self%b_ptr => build
    end subroutine bind

    !> Refine the stored pose of the particle on its slot of the prepared batch. A state-zero
    !! particle is rejected; an invalid seed or a particle the batch does not hold is a caller
    !! contract violation.
    subroutine srch_cont( self, os, ithr )
        class(strategy3D_cont), intent(inout) :: self
        class(oris),            intent(inout) :: os
        integer,                intent(in)    :: ithr
        real(dp) :: objective_before, objective_after, rotmat(3,3), shift(2)
        integer  :: islot, state
        logical  :: iseven
        if( .not. self%exists ) THROW_HARD('strategy3D_cont used before new')
        if( ithr < 1 ) THROW_HARD('strategy3D_cont received an invalid thread index')
        if( os%get_state(self%spec%iptcl) <= 0 )then
            call os%reject(self%spec%iptcl)
            return
        endif
        call os%get_ori(self%spec%iptcl, self%seed_ori)
        if( .not. cont_seed_is_valid(self%seed_ori) )then
            write(logfhandle,'(a,i0)') 'invalid refine=cont seed particle: ', self%spec%iptcl
            THROW_HARD('a Cartesian pass requires finite ptcl3D poses with state, half and projection')
        endif
        state  = self%seed_ori%get_state()
        iseven = self%seed_ori%get_eo() == 0
        if( .not. self%b_ptr%cftc%ref_exists(state, iseven) ) THROW_HARD('strategy3D_cont reference state/half is unavailable')
        islot = self%b_ptr%cftc%get_ptcl_slot(self%spec%iptcl)
        if( islot < 1 ) THROW_HARD('strategy3D_cont particle is not in the prepared batch')
        ! one transaction; result_ori changes only when the solve is accepted
        self%result_ori = self%seed_ori
        call self%opt%refine(self%b_ptr%cftc, islot, self%result_ori)
        self%status = self%opt%get_status()
        call self%opt%get_objectives(objective_before, objective_after)
        ! the score at the committed pose: the result, or the seed after a rejection
        self%score = 0.
        if( self%status /= CARTFT_INVALID_PREPARATION .and. objective_before >= 0._dp ) &
            &self%score = real(self%b_ptr%cftc%score(islot, objective_after))
        ! sigma2 belongs to the Euclidean objective only (C5); the committed pose relative to the
        ! slot's observation, centred on the stored shift
        if( self%p_ptr%cc_objfun == OBJFUN_EUCLID .and. self%status /= CARTFT_INVALID_PREPARATION )then
            rotmat = real(self%result_ori%get_mat(), dp)
            shift  = real((self%result_ori%get_2Dshift() - self%seed_ori%get_2Dshift()) * &
                &real(self%p_ptr%box_crop)/real(self%p_ptr%box), dp)
            call self%b_ptr%esig%calc_sigma2(self%b_ptr%cftc, islot, self%spec%iptcl, state, iseven, rotmat, shift)
        endif
        self%searched = .true.
        call self%oris_assign
    end subroutine srch_cont

    !> Commit the pose (when accepted) and corr_cart together with the improved flag and, in a
    !! refine=cont pass, the convergence fields of the seed-to-result motion; a polish pass
    !! leaves those of the discrete search. corr stays as the last polar pass left it.
    subroutine oris_assign_cont( self )
        class(strategy3D_cont), intent(inout) :: self
        type(ori) :: symmetry_equivalent_ori
        real      :: euler_distance, inplane_distance
        if( .not. self%searched ) THROW_HARD('strategy3D_cont has no search result to assign')
        associate( iptcl => self%spec%iptcl, field => self%b_ptr%spproj_field )
            if( self%status == CARTFT_ACCEPTED ) call field%set_ori(iptcl, self%result_ori)
            call field%set(iptcl, 'corr_cart', self%score)
            call field%set(iptcl, 'pose_cont_improved', merge(1., 0., self%status == CARTFT_ACCEPTED))
            ! the polish pass leaves the convergence fields of the discrete search it follows, so
            ! the main-loop rule keeps measuring that search (C8, C19; the committed record is the
            ! seed's, which carries them)
            if( .not. self%p_ptr%l_cont_polish )then
                call self%b_ptr%pgrpsyms%sym_dists(self%seed_ori, self%result_ori, symmetry_equivalent_ori, &
                    &euler_distance, inplane_distance)
                call field%set(iptcl, 'dist',      euler_distance)
                call field%set(iptcl, 'dist_inpl', inplane_distance)
                call field%set(iptcl, 'shincarg',  arg(self%result_ori%get_2Dshift() - self%seed_ori%get_2Dshift()))
                call field%set(iptcl, 'mi_proj',   merge(1., 0., euler_distance <= self%p_ptr%angthres_mi_proj))
                call field%set(iptcl, 'mi_state',  1.)
                ! the transaction exhausts its one seed-centred local domain
                call field%set(iptcl, 'frac',      100.)
            endif
        end associate
        call symmetry_equivalent_ori%kill
    end subroutine oris_assign_cont

    subroutine kill_cont( self )
        class(strategy3D_cont), intent(inout) :: self
        call self%opt%kill
        call self%seed_ori%kill
        call self%result_ori%kill
        nullify(self%p_ptr, self%b_ptr, self%spec%eulprob_obj_part)
        self%status   = CARTFT_NOT_ATTEMPTED
        self%score    = 0.
        self%searched = .false.
        self%exists   = .false.
    end subroutine kill_cont

end module simple_strategy3D_cont
