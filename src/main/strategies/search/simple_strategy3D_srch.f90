!@descr: common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_strategy3D_alloc, only: prep_strategy3D_thread, s3D, &
    &seed_continuous_inplane_candidate
use simple_builder,          only: builder
use simple_eul_prob_tab,     only: eul_prob_tab
use simple_parameters,       only: parameters
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
implicit none

integer, parameter :: CONT_ROUTE_NOT_ATTEMPTED = 0
integer, parameter :: CONT_ROUTE_IMPROVED       = 1
integer, parameter :: CONT_ROUTE_NO_IMPROVEMENT = 2
integer, parameter :: CONT_ROUTE_FALLBACK       = 3

public :: strategy3D_srch, strategy3D_spec
private
#include "simple_local_flags.inc"

type strategy3D_spec
    type(eul_prob_tab), pointer :: eulprob_obj_part
    integer :: iptcl=0, iptcl_map=0
end type strategy3D_spec

type strategy3D_srch
    character(len=:), allocatable :: refine                   !< 3D refinement flag
    type(pftc_shsrch_grad)     :: grad_shsrch_obj             !< origin shift search object, L-BFGS with gradient
    type(pftc_shsrch_grad)     :: grad_shsrch_obj2            !<
    type(pftc_shsrch_grad)     :: grad_shsrch_first_obj       !< origin shift search object, L-BFGS with gradient, used for initial shift search on previous ref
    type(pftc_shsrch_grad)     :: joint_inpl_optimizer        !< selected-reference joint (sx,sy,rotind_frac) refinement
    type(ori)                  :: o_prev                      !< previous orientation
    type(oris)                 :: opeaks                      !< peak orientations to consider for refinement
    class(parameters), pointer :: p_ptr        => null()      !< global parameters
    class(builder),    pointer :: b_ptr        => null()      !< build handle for access to spproj/eulspace/symmetry
    integer                    :: iptcl           = 0         !< global particle index
    integer                    :: iptcl_map       = 0         !< map particle index
    integer                    :: ithr            = 0         !< thread index
    integer                    :: nrefs           = 0         !< total # references (nstates*nprojs)
    integer                    :: nrefs_sub       = 0         !< total # references (nstates*nprojs), subspace
    integer                    :: npeaks          = 0         !< # peak subspace orientations to consider
    integer                    :: npeaks_inpl     = 0         !< # # multi-neighborhood peaks to refine with L-BFGS
    integer                    :: nstates         = 0         !< # states
    integer                    :: nprojs          = 0         !< # projections
    integer                    :: nprojs_sub      = 0         !< # projections, subspace
    integer                    :: nrots           = 0         !< # in-plane rotations
    integer                    :: nsym            = 0         !< symmetry order
    integer                    :: nnn             = 0         !< # nearest neighbors
    integer                    :: nbetter         = 0         !< # better orientations identified
    integer                    :: nrefs_eval      = 0         !< # references evaluated
    integer                    :: nsolns          = 0         !< # distinct refs with a stored solution this search
    integer                    :: ntrs_eval       = 0         !< # shifts evaluated
    integer                    :: prev_roind      = 0         !< previous in-plane rotation index
    integer                    :: prev_state      = 0         !< previous state index
    integer                    :: prev_ref        = 0         !< previous reference index
    integer                    :: prev_proj       = 0         !< previous projection direction index
    real                       :: athres          = 10.       !< angular treshold (refine=neighc) for neighborhood continuous Cartesian search
    real                       :: prev_corr       = 1.        !< previous best correlation
    real                       :: prev_shvec(2)   = 0.        !< previous origin shift vector
    real                       :: xy_first(2)     = 0.        !< initial shifts identified by searching the previous best reference
    real                       :: xy_first_rot(2) = 0.        !< initial shifts identified by searching the previous best reference, rotated
    real                       :: prev_rotind_frac = 0.       !< durable continuous restart coordinate for the previous pose
    logical                    :: xy_first_valid  = .false.   !< pre-selection shift seed was evaluated for this particle
    logical                    :: prev_rotind_frac_valid = .false. !< previous coordinate is finite and coupled to prev_roind
    logical                    :: l_neigh         = .false.   !< neighbourhood refinement flag
    logical                    :: l_greedy        = .false.   !< greedy        refinement flag
    logical                    :: l_sh_first      = .false.   !< use previous-ref shift seed before orientation/state scoring
    logical                    :: continuous_active = .false. !< selected candidates receive continuous in-plane refinement
    integer                    :: continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED !< final selected-reference route outcome
    logical                    :: doshift         = .true.    !< 2 indicate whether 2 serch shifts
    logical                    :: exists          = .false.   !< 2 indicate existence
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: prep4prob
    procedure :: inpl_srch_first
    procedure :: inpl_srch
    procedure :: inpl_srch_peaks
    procedure :: refine_selected_continuously
    procedure :: store_solution
    procedure :: store_continuous_solution
    procedure :: retain_selected_seed_shift
    procedure :: uses_continuous_refinement
    procedure :: bypasses_legacy_post_refinement
    procedure :: requires_legacy_fallback
    procedure :: get_continuous_route_status
    procedure :: resolve_restart_inpl_coordinate
    procedure :: prepare_restart_inplane_seed
    procedure :: kill
end type strategy3D_srch

contains

    subroutine new( self, params, spec, build )
        class(strategy3D_srch),    intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(strategy3D_spec),    intent(in)    :: spec
        class(builder),    target, intent(in)    :: build
        real :: lims(2,2), lims_init(2,2), joint_lims(3,2)
        ! set pointers
        self%p_ptr         => params
        self%b_ptr         => build
        ! set constants
        self%iptcl         = spec%iptcl
        self%iptcl_map     = spec%iptcl_map
        self%nstates       = self%p_ptr%nstates
        self%nprojs        = self%p_ptr%nspace
        self%nprojs_sub    = self%p_ptr%nspace_sub
        self%nrefs         = self%nprojs     * self%nstates
        self%nrefs_sub     = self%nprojs_sub * self%nstates
        self%npeaks        = self%p_ptr%npeaks
        self%npeaks_inpl   = self%p_ptr%npeaks_inpl
        self%athres        = self%p_ptr%athres
        self%nbetter       = 0
        self%nrefs_eval    = 0
        self%ntrs_eval     = 0
        self%nsym          = self%b_ptr%pgrpsyms%get_nsym()
        self%doshift       = self%p_ptr%l_doshift
        self%l_sh_first    = self%doshift .and. self%nstates <= 1
        self%l_neigh       = self%p_ptr%l_neigh
        self%refine        = trim(self%p_ptr%refine)
        self%l_greedy      = str_has_substr(self%p_ptr%refine, 'greedy')
        self%continuous_active = trim(self%p_ptr%inpl_cont) == 'yes' .and. &
            &self%b_ptr%pftc%is_raw_euclid_objfun() .and. &
            &(.not. self%p_ptr%l_objfun_den) .and. &
            &(.not. str_has_substr(self%refine, 'prob'))
        self%continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED
        lims(:,1)          = -self%p_ptr%trs
        lims(:,2)          =  self%p_ptr%trs
        lims_init(:,1)     = -SHC_INPL_TRSHWDTH
        lims_init(:,2)     =  SHC_INPL_TRSHWDTH
        call self%o_prev%new(is_ptcl=.true.)
        call self%opeaks%new(self%npeaks, is_ptcl=.true.)
        ! create in-plane search objects
        self%nrots = self%b_ptr%pftc%get_nrots()
        if( .not. str_has_substr(self%refine, 'prob') )then
            call self%grad_shsrch_obj%new(self%b_ptr, lims, lims_init=lims_init, &
            &maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
            call self%grad_shsrch_obj2%new(self%b_ptr, lims, lims_init=lims_init, maxits=self%p_ptr%maxits_sh,&
            &opt_angle=.false.)
            call self%grad_shsrch_first_obj%new(self%b_ptr, lims, lims_init=lims_init, &
            &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            if( self%continuous_active )then
                joint_lims(1:2,:) = lims
                joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
                call self%joint_inpl_optimizer%new_joint(self%b_ptr, joint_lims, self%p_ptr%maxits_sh)
            endif
        endif
        self%exists = .true.
    end subroutine new

    subroutine prep4prob( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%b_ptr%spproj_field%get_ori(self%iptcl, self%o_prev)         ! previous ori
        self%prev_state = self%o_prev%get_state()                             ! state index
        self%prev_roind = self%b_ptr%pftc%get_roind(360.-self%o_prev%e3get()) ! in-plane angle index
        call self%prepare_restart_inplane_seed()
        self%prev_shvec = self%o_prev%get_2Dshift()                           ! shift vector
        self%prev_proj  = self%b_ptr%eulspace%find_closest_proj(self%o_prev)  ! previous projection direction
        self%prev_ref   = (self%prev_state-1)*self%nprojs + self%prev_proj
        if( trim(self%p_ptr%multivol_mode) /= 'input_oris_fixed' )then
            call self%b_ptr%spproj_field%set(self%iptcl, 'proj', real(self%prev_proj))
        endif
        call prep_strategy3D_thread(self%ithr)
        self%nsolns     = 0
        self%xy_first_valid = .false.
        self%nrefs_eval = 0
        self%ntrs_eval  = 0
        self%prev_corr  = 0.
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4prob')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4prob')
        endif
    end subroutine prep4prob

    subroutine prep4srch( self )
        class(strategy3D_srch), intent(inout) :: self
        real      :: corrs(self%nrots), corr
        integer   :: ipeak, iref_sub, iref_sub_rand, prev_proj_sub
        ! previous parameters
        call self%b_ptr%spproj_field%get_ori(self%iptcl, self%o_prev)         ! previous ori
        self%prev_state = self%o_prev%get_state()                             ! state index
        self%prev_roind = self%b_ptr%pftc%get_roind(360.-self%o_prev%e3get()) ! in-plane angle index
        call self%prepare_restart_inplane_seed()
        self%prev_shvec = self%o_prev%get_2Dshift()                           ! shift vector
        self%prev_proj  = self%b_ptr%eulspace%find_closest_proj(self%o_prev)  ! previous projection direction
        ! map previous projection direction (always use full-space ref index)
        self%prev_ref = (self%prev_state-1)*self%nprojs + self%prev_proj
        call self%b_ptr%spproj_field%set(self%iptcl, 'proj', self%prev_proj)
        ! copy self%o_prev to opeaks to transfer paticle-dependent parameters
        do ipeak = 1, self%npeaks
            call self%opeaks%set_ori(ipeak, self%o_prev)
        end do
        ! init threaded search arrays
        call prep_strategy3D_thread(self%ithr)
        self%nsolns = 0
        self%xy_first_valid = .false.
        ! search order
        ! -- > full space
        call s3D%rts(     self%ithr)%ne_ran_iarr(s3D%srch_order(:,self%ithr))
        call s3D%rts_inpl(self%ithr)%ne_ran_iarr(s3D%inpl_order(:,self%ithr))
        call put_last(self%prev_ref, s3D%srch_order(:,self%ithr))
        ! --> subspace
        if( self%l_neigh )then
            call s3D%rts_sub(self%ithr)%ne_ran_iarr(s3D%srch_order_sub(:,self%ithr))
            ! --> do the mapping
            do iref_sub = 1, self%nrefs_sub
                iref_sub_rand = s3D%srch_order_sub(iref_sub,self%ithr)
                ! index for b_ptr%subspace_inds needs to be a projection direction
                prev_proj_sub = self%b_ptr%subspace_inds(&
                &iref_sub_rand - (self%prev_state - 1) * self%p_ptr%nspace_sub)
                ! but then we need to turn it back into a reference index in the full search space
                s3D%srch_order_sub(iref_sub,self%ithr) = (self%prev_state - 1) * self%nprojs + prev_proj_sub
            enddo
            ! --> check if we have prev_ref in subspace, put last
            if( any(s3D%srch_order_sub(:,self%ithr) == self%prev_ref ) )then
                call put_last(self%prev_ref, s3D%srch_order_sub(:,self%ithr))
            endif
        endif
        ! sanity check
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4srch')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4srch')
        endif
        ! prep corr
        call self%b_ptr%pftc%gen_objfun_vals(self%prev_ref, self%iptcl, [0.,0.], corrs)
        corr = max(0.,maxval(corrs))
        self%prev_corr = corr
    end subroutine prep4srch

    subroutine inpl_srch_first( self )
        class(strategy3D_srch), intent(inout) :: self
        real    :: cxy(3), rotmat(2,2)
        integer :: irot
        if( .not. self%l_sh_first ) return
        ! BFGS over shifts with in-plane rot exhaustive callback
        irot = self%prev_roind
        call self%grad_shsrch_first_obj%set_indices(self%prev_ref, self%iptcl)
        cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
        endif
        self%xy_first_valid = .true.
    end subroutine inpl_srch_first

    subroutine inpl_srch( self, ref, irot_in, force_legacy )
        class(strategy3D_srch), intent(inout) :: self
        integer, optional,      intent(in)    :: ref
        integer, optional,      intent(inout) :: irot_in
        logical, optional,      intent(in)    :: force_legacy
        real    :: cxy(3)
        integer :: iref, irot, loc(1)
        logical :: run_legacy
        if( .not. self%doshift ) return
        run_legacy = .false.
        if( present(force_legacy) ) run_legacy = force_legacy
        if( self%bypasses_legacy_post_refinement() .and. .not. run_legacy ) return
        if( present(ref) )then
            iref = ref
        else
            loc  = maxloc(s3D%proj_space_corrs(:,self%ithr))
            iref = loc(1)
        endif
        irot = s3D%proj_space_inplinds(iref,self%ithr)
        if( present(irot_in) ) irot = irot_in
        ! BFGS over shifts with in-plane rot exhaustive callback
        call self%grad_shsrch_obj%set_indices(iref, self%iptcl)
        if( self%l_sh_first )then
            cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
        else
            cxy = self%grad_shsrch_obj%minimize(irot=irot)
        endif
        if( irot > 0 ) call self%store_solution(iref, irot, cxy(1), sh=cxy(2:3))
        if( present(irot_in) ) irot_in = irot
    end subroutine inpl_srch

    subroutine inpl_srch_peaks( self, npeaks_inpl )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: npeaks_inpl
        real    :: cxy(3)
        integer :: refs(npeaks_inpl), irot, ipeak
        if( .not. self%doshift ) return
        if( self%bypasses_legacy_post_refinement() ) return
        ! BFGS over shifts with in-plane rot exhaustive callback
        refs = maxnloc(s3D%proj_space_corrs(:,self%ithr), npeaks_inpl)
        do ipeak = 1, npeaks_inpl
            call self%grad_shsrch_obj%set_indices(refs(ipeak), self%iptcl)
            if( self%l_sh_first )then
                cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
            else
                cxy = self%grad_shsrch_obj%minimize(irot=irot)
            endif
            if( irot > 0 )then
                ! irot > 0 guarantees improvement found, update solution
                call self%store_solution(refs(ipeak), irot, cxy(1), sh=cxy(2:3))
            endif
        enddo
    end subroutine inpl_srch_peaks

    subroutine store_solution( self, ref, inpl_ind, corr, sh )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr
        real,       optional,   intent(in)    :: sh(2)
        if( s3D%proj_space_corrs(ref, self%ithr) <= -huge(1.0)/2.0 )then
            self%nsolns = self%nsolns + 1
        elseif( corr <= s3D%proj_space_corrs(ref, self%ithr) )then
            return
        endif
        if( present(sh) ) s3D%proj_space_shift(:,ref,self%ithr) = sh
        s3D%proj_space_inplinds(ref,self%ithr) = inpl_ind
        call seed_continuous_inplane_candidate(inpl_ind, &
            &s3D%proj_space_inplcoords(ref,self%ithr), &
            &s3D%proj_space_inplvalid(ref,self%ithr))
        s3D%proj_space_corrs(   ref,self%ithr) = corr
    end subroutine store_solution

    !> Jointly refine shifts and in-plane rotation for the selected 3D
    !! reference. The state and projection direction encoded by ref remain
    !! fixed. A valid non-improving result retains the discrete score and angle
    !! while coupling them to the shift used during candidate scoring. A
    !! numerically invalid result invokes the original legacy post-selection
    !! optimizer for this same reference.
    subroutine refine_selected_continuously( self, ref )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref
        real     :: cxy(3), joint_lims(3,2), rotmat(2,2), xy_native(2), xy_seed_rot(2)
        real(dp) :: rotind_frac, rotind_seed, initial_cost
        integer  :: irot, selected_irot
        logical  :: evaluation_valid, improved, preserve_restart_seed

        if( .not. self%continuous_active ) return
        if( ref < 1 .or. ref > self%nrefs ) return
        irot = s3D%proj_space_inplinds(ref,self%ithr)
        if( irot < 1 .or. irot > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_FALLBACK
            call self%inpl_srch(ref=ref, force_legacy=.true.)
            return
        endif

        if( self%xy_first_valid )then
            ! The global discrete search was evaluated with this native-frame
            ! pre-selection shift seed. Start joint refinement at the same
            ! point rather than at the zero-valued storage placeholder.
            xy_native = self%xy_first
        else
            ! Stored 3D-search shifts are in the selected reference frame.
            ! The PFTC joint objective consumes the native particle frame.
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            xy_native = matmul(s3D%proj_space_shift(:,ref,self%ithr), transpose(rotmat))
        endif
        if( .not. all(ieee_is_finite(xy_native)) .or. &
            &.not. ieee_is_finite(s3D%proj_space_inplcoords(ref,self%ithr)) )then
            self%continuous_route_outcome = CONT_ROUTE_FALLBACK
            call self%inpl_srch(ref=ref, force_legacy=.true.)
            return
        endif
        joint_lims(1,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        joint_lims(2,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        if( .not. self%doshift )then
            joint_lims(1,:) = xy_native(1)
            joint_lims(2,:) = xy_native(2)
        endif
        joint_lims(3,:) = [real(irot)-2., real(irot)+2.]
        rotind_seed = real(s3D%proj_space_inplcoords(ref,self%ithr),dp)
        preserve_restart_seed = self%prev_rotind_frac_valid .and. &
            &ref == self%prev_ref .and. irot == self%prev_roind .and. &
            &abs(real(self%prev_rotind_frac,dp) - real(irot,dp)) > 1.d-6
        if( preserve_restart_seed ) rotind_seed = real(self%prev_rotind_frac,dp)
        selected_irot = irot

        call self%joint_inpl_optimizer%set_indices(ref, self%iptcl)
        call self%joint_inpl_optimizer%set_limits(joint_lims)
        cxy = self%joint_inpl_optimizer%minimize_joint(irot, &
            &real(rotind_seed), xy_native, &
            &sh_rot=.true., rotind_frac=rotind_frac, &
            &evaluation_valid=evaluation_valid, improved=improved, &
            &initial_cost_out=initial_cost)
        if( self%requires_legacy_fallback(evaluation_valid) )then
            self%continuous_route_outcome = CONT_ROUTE_FALLBACK
            call self%inpl_srch(ref=ref, force_legacy=.true.)
            return
        endif
        if( .not. improved )then
            self%continuous_route_outcome = CONT_ROUTE_NO_IMPROVEMENT
            if( preserve_restart_seed )then
                call rotmat2d(real((rotind_seed - 1.d0) * &
                    &real(self%b_ptr%pftc%get_dang(),dp)), rotmat)
                xy_seed_rot = matmul(xy_native, rotmat)
                call self%store_continuous_solution(ref, self%prev_roind, &
                    &rotind_seed, real(exp(-initial_cost)), xy_seed_rot)
            elseif( self%xy_first_valid )then
                ! The discrete winner was scored at xy_first in the native
                ! particle frame. Retain that same seed in the selected-angle
                ! reference frame so the stored score and shift describe one
                ! pose even when the joint optimizer finds no improvement.
                ! minimize_joint returns irot=0 for this outcome, so use the
                ! selected grid index saved before invoking the optimizer.
                call rotmat2d(self%b_ptr%pftc%get_rot(selected_irot), rotmat)
                xy_seed_rot = matmul(xy_native, rotmat)
                call self%retain_selected_seed_shift(ref, xy_seed_rot)
            endif
            return
        endif
        if( irot <= 0 )then
            self%continuous_route_outcome = CONT_ROUTE_FALLBACK
            call self%inpl_srch(ref=ref, force_legacy=.true.)
            return
        endif
        call self%store_continuous_solution(ref, irot, rotind_frac, cxy(1), cxy(2:3))
        self%continuous_route_outcome = CONT_ROUTE_IMPROVED
    end subroutine refine_selected_continuously

    !> Commit an already accepted joint result without changing its selected
    !! state or projection reference.
    subroutine store_continuous_solution( self, ref, inpl_ind, inpl_coord, corr, sh )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real(dp),               intent(in)    :: inpl_coord
        real,                   intent(in)    :: corr, sh(2)

        s3D%proj_space_shift(:,ref,self%ithr)    = sh
        s3D%proj_space_inplinds(ref,self%ithr)   = inpl_ind
        s3D%proj_space_inplcoords(ref,self%ithr) = real(inpl_coord)
        s3D%proj_space_inplvalid(ref,self%ithr)  = .true.
        s3D%proj_space_corrs(ref,self%ithr)      = corr
    end subroutine store_continuous_solution

    !> Retain the shift used to score an already-selected grid candidate.
    !! The score, integer angle, fractional grid seed, and validity flag must
    !! remain unchanged because no continuous improvement was accepted.
    subroutine retain_selected_seed_shift( self, ref, sh )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref
        real,                   intent(in)    :: sh(2)

        if( ref < 1 .or. ref > size(s3D%proj_space_corrs,1) )then
            THROW_HARD('selected reference outside stored search space; retain_selected_seed_shift')
        endif
        if( self%ithr < 1 .or. self%ithr > size(s3D%proj_space_corrs,2) )then
            THROW_HARD('thread index outside stored search space; retain_selected_seed_shift')
        endif
        if( .not. all(ieee_is_finite(sh)) )then
            THROW_HARD('non-finite selected seed shift; retain_selected_seed_shift')
        endif
        s3D%proj_space_shift(:,ref,self%ithr) = sh
    end subroutine retain_selected_seed_shift

    pure logical function uses_continuous_refinement( self )
        class(strategy3D_srch), intent(in) :: self
        uses_continuous_refinement = self%continuous_active
    end function uses_continuous_refinement

    pure logical function bypasses_legacy_post_refinement( self )
        class(strategy3D_srch), intent(in) :: self
        bypasses_legacy_post_refinement = self%continuous_active
    end function bypasses_legacy_post_refinement

    pure logical function requires_legacy_fallback( self, evaluation_valid )
        class(strategy3D_srch), intent(in) :: self
        logical,                intent(in) :: evaluation_valid
        requires_legacy_fallback = self%continuous_active .and. .not. evaluation_valid
    end function requires_legacy_fallback

    pure subroutine get_continuous_route_status( self, attempted, improved, no_improvement, fallback )
        class(strategy3D_srch), intent(in)  :: self
        logical,                intent(out) :: attempted, improved, no_improvement, fallback

        attempted      = self%continuous_route_outcome /= CONT_ROUTE_NOT_ATTEMPTED
        improved       = self%continuous_route_outcome == CONT_ROUTE_IMPROVED
        no_improvement = self%continuous_route_outcome == CONT_ROUTE_NO_IMPROVEMENT
        fallback       = self%continuous_route_outcome == CONT_ROUTE_FALLBACK
    end subroutine get_continuous_route_status

    !> Recover the continuous rotation-index coordinate from the durable Euler
    !! angle on the current PFTC grid and unwrap it to the periodic copy nearest
    !! the current e3-derived integer index. The persisted integer `inpl` is
    !! not used because refine3D may rebuild a different grid on restart.
    pure subroutine resolve_restart_inpl_coordinate( self, e3, current_inpl, dang, coordinate, valid )
        class(strategy3D_srch), intent(in)  :: self
        real,                   intent(in)  :: e3
        integer,                intent(in)  :: current_inpl
        real(dp),               intent(in)  :: dang
        real(dp),               intent(out) :: coordinate
        logical,                intent(out) :: valid
        real(dp) :: raw_coordinate, period
        integer :: nearest_inpl

        coordinate = 0.d0
        valid = .false.
        if( self%nrots < 1 .or. current_inpl < 1 .or. current_inpl > self%nrots ) return
        if( dang <= 0.d0 .or. .not. ieee_is_finite(e3) ) return
        period = real(self%nrots,dp)
        raw_coordinate = modulo(360.d0 - real(e3,dp), 360.d0) / dang + 1.d0
        nearest_inpl = modulo(nint(raw_coordinate)-1,self%nrots)+1
        if( nearest_inpl /= current_inpl ) return
        coordinate = real(current_inpl,dp) + modulo(raw_coordinate - real(current_inpl,dp) + &
            &0.5d0*period, period) - 0.5d0*period
        valid = ieee_is_finite(coordinate)
    end subroutine resolve_restart_inpl_coordinate

    subroutine prepare_restart_inplane_seed( self )
        class(strategy3D_srch), intent(inout) :: self
        real(dp) :: coordinate, dang
        self%prev_rotind_frac = 0.
        self%prev_rotind_frac_valid = .false.
        if( .not. self%continuous_active ) return
        dang = real(self%b_ptr%pftc%get_dang(),dp)
        call self%resolve_restart_inpl_coordinate(self%o_prev%e3get(), self%prev_roind, dang, &
            &coordinate, self%prev_rotind_frac_valid)
        if( .not. self%prev_rotind_frac_valid ) return
        self%prev_rotind_frac = real(coordinate)
    end subroutine prepare_restart_inplane_seed

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%grad_shsrch_obj2%kill
        call self%grad_shsrch_first_obj%kill
        call self%joint_inpl_optimizer%kill
        call self%opeaks%kill
        call self%o_prev%kill
        if( allocated(self%refine) ) deallocate(self%refine)
        self%continuous_active = .false.
        self%continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED
        self%xy_first_valid = .false.
        self%prev_rotind_frac = 0.
        self%prev_rotind_frac_valid = .false.
        nullify(self%b_ptr)
    end subroutine kill

end module simple_strategy3D_srch
