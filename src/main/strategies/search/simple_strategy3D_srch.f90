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
integer, parameter :: CONT_ROUTE_INVALID        = 3

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
    procedure :: refine_assignment_continuously
    procedure :: store_solution
    procedure :: store_continuous_solution
    procedure :: store_discrete_seed_solution
    procedure :: uses_continuous_refinement
    procedure :: joint_evaluation_invalid
    procedure :: get_continuous_route_status
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
        if( trim(self%p_ptr%inpl_cont) == 'yes' .and. &
            &(.not. self%b_ptr%pftc%is_joint_grad_objfun()) )then
            THROW_HARD('inpl_cont=yes requires a supported Euclidean, hybrid, or cc joint objective')
        endif
        self%continuous_active = trim(self%p_ptr%inpl_cont) == 'yes' .and. &
            &self%p_ptr%l_doshift .and. &
            &self%b_ptr%pftc%is_joint_grad_objfun()
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
            ! selection parity: candidate scoring, shift-seed estimation and
            ! peak refinement always use the legacy objects; the joint
            ! optimizer only polishes the committed assignment
            call self%grad_shsrch_obj%new_legacy(self%b_ptr, lims, lims_init=lims_init, &
                &maxits=self%p_ptr%maxits_sh)
            call self%grad_shsrch_first_obj%new_legacy(self%b_ptr, lims, lims_init=lims_init, &
                &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
            if( self%continuous_active )then
                joint_lims(1:2,:) = lims
                joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
                call self%joint_inpl_optimizer%new_joint(self%b_ptr, joint_lims, self%p_ptr%maxits_sh)
            endif
            call self%grad_shsrch_obj2%new_fixed(self%b_ptr, lims, lims_init=lims_init, &
                &maxits=self%p_ptr%maxits_sh)
        else if( self%continuous_active )then
            joint_lims(1:2,:) = lims
            joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
            call self%joint_inpl_optimizer%new_joint(self%b_ptr, joint_lims, self%p_ptr%maxits_sh)
        endif
        self%exists = .true.
    end subroutine new

    subroutine prep4prob( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%b_ptr%spproj_field%get_ori(self%iptcl, self%o_prev)         ! previous ori
        self%prev_state = self%o_prev%get_state()                             ! state index
        self%prev_roind = self%b_ptr%pftc%get_roind(360.-self%o_prev%e3get()) ! in-plane angle index
        self%prev_shvec = self%o_prev%get_2Dshift()                           ! shift vector
        self%prev_proj  = self%b_ptr%eulspace%find_closest_proj(self%o_prev)  ! previous projection direction
        self%prev_ref   = (self%prev_state-1)*self%nprojs + self%prev_proj
        if( trim(self%p_ptr%multivol_mode) /= 'input_oris_fixed' )then
            call self%b_ptr%spproj_field%set(self%iptcl, 'proj', real(self%prev_proj))
        endif
        call prep_strategy3D_thread(self%ithr)
        self%nsolns     = 0
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
        irot = self%prev_roind
        call self%grad_shsrch_first_obj%set_indices(self%prev_ref, self%iptcl)
        ! selection parity: the pre-selection shift seed is estimated with the
        ! legacy optimizer regardless of inpl_cont
        cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
        endif
    end subroutine inpl_srch_first

    subroutine inpl_srch( self, ref, irot_in )
        class(strategy3D_srch), intent(inout) :: self
        integer, optional,      intent(in)    :: ref
        integer, optional,      intent(inout) :: irot_in
        real    :: cxy(3)
        integer :: iref, irot, loc(1)
        if( .not. self%doshift ) return
        ! selection parity: legacy post-selection refinement runs unchanged
        ! under both inpl_cont values; the joint solve polishes afterwards
        ! (extract_peak_ori), never instead
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
        ! selection parity: multi-peak shift refinement shapes the final
        ! selection, so it runs as legacy under the continuous route too
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
        if( allocated(s3D%proj_space_inplcoords) .neqv. &
            &allocated(s3D%proj_space_inplvalid) )then
            THROW_HARD('incomplete continuous in-plane candidate storage')
        endif
        if( allocated(s3D%proj_space_inplcoords) )then
            call seed_continuous_inplane_candidate(inpl_ind, &
                &s3D%proj_space_inplcoords(ref,self%ithr), &
                &s3D%proj_space_inplvalid(ref,self%ithr))
        endif
        s3D%proj_space_corrs(   ref,self%ithr) = corr
    end subroutine store_solution

    !> Jointly refine shifts and in-plane rotation for the selected 3D
    !! reference. The state, projection direction, and incoming discrete
    !! in-plane cell remain authoritative. The solve locally refines that cell
    !! and its native-frame shift. A valid non-improving result retains the
    !! discrete pose with a consistent re-scored objective; a numerically
    !! invalid result leaves the selected candidate untouched.
    subroutine refine_selected_continuously( self, ref )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref
        real     :: cxy(3), joint_lims(3,2), rotmat(2,2), xy_native(2)
        real(dp) :: rotind_frac
        integer  :: irot, incoming_irot
        logical  :: evaluation_valid, improved

        if( .not. self%continuous_active ) return
        if( ref < 1 .or. ref > self%nrefs ) return
        incoming_irot = s3D%proj_space_inplinds(ref,self%ithr)
        if( incoming_irot < 1 .or. incoming_irot > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif

        ! The committed candidate storage is authoritative: legacy
        ! post-selection refinement has already replaced any pre-selection
        ! seed (xy_first), so the polish always starts from the stored
        ! selected-candidate shift.  Stored 3D-search shifts are in the
        ! selected reference frame; the PFTC joint objective consumes the
        ! native particle frame.
        call rotmat2d(self%b_ptr%pftc%get_rot(incoming_irot), rotmat)
        xy_native = matmul(s3D%proj_space_shift(:,ref,self%ithr), transpose(rotmat))
        if( .not. all(ieee_is_finite(xy_native)) )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        joint_lims(1,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        joint_lims(2,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        if( .not. self%doshift )then
            joint_lims(1,:) = xy_native(1)
            joint_lims(2,:) = xy_native(2)
        endif
        joint_lims(3,:) = [1.-2., real(self%nrots)+2.]

        call self%joint_inpl_optimizer%set_indices(ref, self%iptcl)
        call self%joint_inpl_optimizer%set_limits(joint_lims)
        ! LOCAL refinement around the discrete search's selection (irot_in): the
        ! durable pass never performs another global all-angle reselection
        cxy = self%joint_inpl_optimizer%minimize_joint(irot, xy_native, &
            &sh_rot=.true., rotind_frac=rotind_frac, &
            &evaluation_valid=evaluation_valid, improved=improved, irot_in=incoming_irot)
        if( self%joint_evaluation_invalid(evaluation_valid) )then
            ! retain the incoming assignment untouched
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( .not. improved )then
            ! the seed is the incoming cell re-scored, so this commit retains the
            ! incoming pose with a consistent objective value
            self%continuous_route_outcome = CONT_ROUTE_NO_IMPROVEMENT
            call self%store_discrete_seed_solution(ref, irot, cxy(1), cxy(2:3))
            return
        endif
        if( irot <= 0 )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        call self%store_continuous_solution(ref, irot, rotind_frac, cxy(1), cxy(2:3))
        self%continuous_route_outcome = CONT_ROUTE_IMPROVED
    end subroutine refine_selected_continuously

    !> Run the durable joint refinement after a probabilistic assignment has
    !! selected its state, projection, integer in-plane cell, and shift.  The
    !! probability table deliberately carries no fractional angle. The incoming
    !! assignment is authoritative: the solve refines LOCALLY around its in-plane
    !! cell and shift (no global all-angle reselection -- that would hop between
    !! degenerate in-plane branches under dihedral symmetry). An invalid solve
    !! retains the incoming assignment untouched; a valid non-improving solve
    !! retains the pose with its re-scored objective value.
    subroutine refine_assignment_continuously( self, ref, inpl, corr, sh, inpl_coord, inpl_valid )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref
        integer,                intent(inout) :: inpl
        real,                   intent(inout) :: corr, sh(2)
        real,                   intent(out)   :: inpl_coord
        logical,                intent(out)   :: inpl_valid
        real     :: cxy(3), joint_lims(3,2), rotmat(2,2), xy_native(2)
        real(dp) :: rotind_frac
        integer  :: refined_inpl
        logical  :: evaluation_valid, improved

        inpl_coord = real(inpl)
        inpl_valid = .false.
        if( .not. self%continuous_active ) return
        if( ref < 1 .or. ref > self%nrefs )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( inpl < 1 .or. inpl > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( .not. all(ieee_is_finite(sh)) .or. .not. ieee_is_finite(corr) )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif

        ! Probabilistic assignment shifts are stored in the selected reference
        ! frame; the joint PFTC objective consumes the native particle frame.
        call rotmat2d(self%b_ptr%pftc%get_rot(inpl), rotmat)
        xy_native = matmul(sh, transpose(rotmat))
        joint_lims(1,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        joint_lims(2,:) = [-self%p_ptr%trs, self%p_ptr%trs]
        if( .not. self%doshift )then
            joint_lims(1,:) = xy_native(1)
            joint_lims(2,:) = xy_native(2)
        endif
        joint_lims(3,:) = [1.-2., real(self%nrots)+2.]

        call self%joint_inpl_optimizer%set_indices(ref, self%iptcl)
        call self%joint_inpl_optimizer%set_limits(joint_lims)
        ! LOCAL refinement: the probabilistic assignment is authoritative -- it came
        ! from exhaustive candidate construction -- so no global all-angle
        ! reselection happens here (irot_in), which would hop between degenerate
        ! in-plane branches on floating-point noise under dihedral symmetry
        cxy = self%joint_inpl_optimizer%minimize_joint(refined_inpl, xy_native, &
            &sh_rot=.true., rotind_frac=rotind_frac, evaluation_valid=evaluation_valid, &
            &improved=improved, irot_in=inpl)
        if( .not. evaluation_valid )then
            ! retain the incoming assignment untouched
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( .not. improved )then
            ! the seed is the incoming cell re-scored at its own shift, so this
            ! commit retains the incoming pose with a consistent objective value
            inpl       = refined_inpl
            corr       = cxy(1)
            sh         = cxy(2:3)
            inpl_coord = real(refined_inpl)
            self%continuous_route_outcome = CONT_ROUTE_NO_IMPROVEMENT
            return
        endif
        if( refined_inpl < 1 .or. refined_inpl > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif

        inpl       = refined_inpl
        corr       = cxy(1)
        sh         = cxy(2:3)
        inpl_coord = real(rotind_frac)
        inpl_valid = .true.
        self%continuous_route_outcome = CONT_ROUTE_IMPROVED
    end subroutine refine_assignment_continuously

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

    !> Commit the discrete seed returned when a valid joint solve does not
    !! materially improve its starting pose.
    subroutine store_discrete_seed_solution( self, ref, inpl_ind, corr, sh )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr, sh(2)

        s3D%proj_space_shift(:,ref,self%ithr)    = sh
        s3D%proj_space_inplinds(ref,self%ithr)   = inpl_ind
        s3D%proj_space_inplcoords(ref,self%ithr) = real(inpl_ind)
        s3D%proj_space_inplvalid(ref,self%ithr)  = .false.
        s3D%proj_space_corrs(ref,self%ithr)      = corr
    end subroutine store_discrete_seed_solution

    pure logical function uses_continuous_refinement( self )
        class(strategy3D_srch), intent(in) :: self
        uses_continuous_refinement = self%continuous_active
    end function uses_continuous_refinement

    pure logical function joint_evaluation_invalid( self, evaluation_valid )
        class(strategy3D_srch), intent(in) :: self
        logical,                intent(in) :: evaluation_valid
        joint_evaluation_invalid = self%continuous_active .and. .not. evaluation_valid
    end function joint_evaluation_invalid

    pure subroutine get_continuous_route_status( self, attempted, improved, no_improvement, invalid )
        class(strategy3D_srch), intent(in)  :: self
        logical,                intent(out) :: attempted, improved, no_improvement, invalid

        attempted      = self%continuous_route_outcome /= CONT_ROUTE_NOT_ATTEMPTED
        improved       = self%continuous_route_outcome == CONT_ROUTE_IMPROVED
        no_improvement = self%continuous_route_outcome == CONT_ROUTE_NO_IMPROVEMENT
        invalid        = self%continuous_route_outcome == CONT_ROUTE_INVALID
    end subroutine get_continuous_route_status

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
        nullify(self%b_ptr)
    end subroutine kill

end module simple_strategy3D_srch
