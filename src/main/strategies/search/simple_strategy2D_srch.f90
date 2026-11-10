!@descr: common strategy2D methods and type specification for polymorphic strategy2D object creation are delegated to this class
module simple_strategy2D_srch
use simple_core_module_api, only: dp
use simple_pftc_srch_api
use simple_strategy2D_alloc, only: prep_strategy2D_thread, s2D, is_fresh_2D_start
use simple_eul_prob_tab2D,   only: eul_prob_tab2D
use simple_eul_prob_tab_utils, only: posterior_temperature, record_shift_state, shift_box_radius, prob_rnd_scope
use simple_rnd,              only: rnd_keyed_int, rnd_draw_key, RND_FAM_ASSIGN
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_shiftmap_srch,    only: shiftmap_srch
use simple_builder,          only: builder
implicit none

public :: strategy2D_srch, strategy2D_spec
private
#include "simple_local_flags.inc"

! continuous joint-route outcomes, mirroring strategy3D_srch
integer, parameter :: CONT_ROUTE_NOT_ATTEMPTED  = 0
integer, parameter :: CONT_ROUTE_IMPROVED       = 1
integer, parameter :: CONT_ROUTE_NO_IMPROVEMENT = 2
integer, parameter :: CONT_ROUTE_INVALID        = 3

type strategy2D_spec
    type(eul_prob_tab2D), pointer :: eulprob_obj_part2D => null()
    real    :: stoch_bound = 0.
    integer :: iptcl       = 0  ! global particle index
    integer :: iptcl_batch = 0  ! maps to index in batch
    integer :: iptcl_map   = 0  ! index in all contiguous batches
end type strategy2D_spec

type strategy2D_srch
    class(parameters), pointer :: p_ptr => null()     ! pointer to parameters
    class(builder),    pointer :: b_ptr => null()     !< build handle for access to pftc
    type(pftc_shsrch_grad) :: grad_shsrch_obj         !< origin shift search object, L-BFGS with gradient
    type(pftc_shsrch_grad) :: grad_shsrch_obj2        !< fixed-angle origin shift search object
    type(pftc_shsrch_grad) :: grad_shsrch_first_obj   !< initial shift search on the previous reference
    type(pftc_shsrch_grad) :: joint_inpl_optimizer    !< continuous joint (sx,sy,rotind_frac) refinement
    type(shiftmap_srch), allocatable :: smap         !< the polar shift map under the key (l_map_run)
    integer                 :: nrefs           =  0   !< number of references
    integer                 :: nrots           =  0   !< number of in-plane rotations in polar representation
    integer                 :: nrefs_eval      =  0   !< nr of references evaluated
    integer                 :: nsolns          =  0   !< # distinct refs with a stored solution this search
    integer                 :: prev_class      =  0   !< previous class index
    integer                 :: prev_class_mi   =  0   !< previous class used for overlap reporting
    integer                 :: best_class      =  0   !< best class index found by search
    integer                 :: best_rot        =  0   !< best in-plane rotation found by search
    integer                 :: prev_rot        =  0   !< previous in-plane rotation found by search
    integer                 :: iptcl           =  0   !< global particle index
    integer                 :: iptcl_batch     =  0   !< index in pre-allocated batch array
    integer                 :: iptcl_map       =  0   !< index in all contiguous batches combined
    integer                 :: ithr            =  0   !< current thread
    real                    :: prev_shvec(2)   =  0.  !< previous origin shift vector
    real                    :: best_shvec(2)   =  0.  !< best shift vector found by search
    real                    :: xy_first(2)     =  0.  !< initial shifts identified by searching the previous best reference
    real                    :: xy_first_rot(2) =  0.  !< initial shifts identified by searching the previous best reference, rotated
    real                    :: prev_corr       = -1.  !< previous best correlation
    real                    :: best_corr       = -1.  !< best corr found by search
    real                    :: trs             =  0.  !< shift boundary
    logical                 :: l_sh_first      = .false. !< Whether to search the shifts on previous best reference
    logical                 :: l_fresh_start   = .false. !< Whether previous alignment parameters are intentionally ignored
    logical                 :: l_map_run       = .false. !< the classic search centres with the polar shift map (shift_sigma)
    logical                 :: l_map           = .false. !< ... and searches this particle's shifts this iteration
    logical                 :: l_peaks_mapped  = .false. !< the peaks' screen polished the classes this search
    logical                 :: continuous_active = .false. !< selected candidates receive continuous refinement
    logical                 :: has_continuous_e3 = .false. !< whether continuous_e3 is valid
    real                    :: continuous_e3 = 0. !< continuous e3 angle for the selected candidate
    integer                 :: continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED !< joint-route outcome for the selected candidate
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch_first
    procedure :: inpl_srch
    procedure :: refine_selected_continuously
    procedure :: inpl_srch_peaks
    procedure, private :: map_srch
    procedure :: store_solution
    procedure :: assign_ori
    procedure :: uses_continuous_refinement
    procedure :: continuous_pose_committed
    procedure :: get_continuous_route_status
    procedure :: kill
end type strategy2D_srch

contains

    subroutine new( self, params, spec, build )
        class(strategy2D_srch),    intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(strategy2D_spec),    intent(in)    :: spec
        class(builder),    target, intent(in)    :: build
        real    :: lims(2,2), lims_init(2,2), joint_lims(3,2)
        logical :: continuous_eligible
        call self%kill
        ! set pointer to parameters
        self%p_ptr => params
        ! set pointer to builder
        self%b_ptr => build
        ! set constants
        self%iptcl       = spec%iptcl
        self%iptcl_batch = spec%iptcl_batch
        self%iptcl_map   = spec%iptcl_map
        self%nrefs       = self%p_ptr%ncls
        self%nrots       = self%b_ptr%pftc%get_nrots()
        if( self%nrots < 1 ) THROW_HARD('strategy2D_srch constructed before PFTC rotations were initialized')
        self%nrefs_eval  = 0
        continuous_eligible = self%b_ptr%pftc%is_joint_grad_objfun() .and. trim(self%p_ptr%tseries) /= 'yes'
        if( trim(self%p_ptr%inpl_cont) == 'yes' .and. trim(self%p_ptr%tseries) /= 'yes' .and. &
            &(.not. continuous_eligible) )then
            THROW_HARD('inpl_cont=yes requires a supported Euclidean or cc joint objective')
        endif
        self%continuous_active = continuous_eligible .and. self%p_ptr%l_doshift .and. &
            &trim(self%p_ptr%inpl_cont) == 'yes'
        self%has_continuous_e3 = .false.
        self%continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED
        ! the polar shift map under the key: the classic modes' own path (the probability tables' consumer reads
        ! its poses from the table)
        self%l_map_run = self%p_ptr%l_centring .and. (.not. self%p_ptr%l_prob_align_mode) .and. &
            &trim(self%p_ptr%tseries) /= 'yes'
        if( self%l_map_run )then
            if( .not. allocated(self%smap) ) allocate(self%smap)
            call self%smap%new(params, build)
        endif
        ! construct composites
        self%trs        =  self%p_ptr%trs
        lims(:,1)       = -self%p_ptr%trs
        lims(:,2)       =  self%p_ptr%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        if( trim(self%p_ptr%tseries).eq.'yes' )then
            ! shift only search
            call self%grad_shsrch_obj%new_fixed(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh)
            call self%grad_shsrch_first_obj%new_fixed(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
        else
            ! selection parity: candidate scoring and shift-seed estimation
            ! always use the legacy objects; the joint optimizer only
            ! polishes the committed assignment
            call self%grad_shsrch_obj%new_alternating(self%b_ptr, lims, lims_init=lims_init,&
                &maxits=self%p_ptr%maxits_sh)
            if( .not. self%l_map_run ) call self%grad_shsrch_first_obj%new_alternating(self%b_ptr, lims,&
                &lims_init=lims_init, maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
        endif
        call self%grad_shsrch_obj2%new_fixed(self%b_ptr, lims, lims_init=lims_init,&
        &maxits=self%p_ptr%maxits_sh)
        if( self%continuous_active )then
            joint_lims(1:2,:) = lims
            joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
            call self%joint_inpl_optimizer%new_joint(self%b_ptr, joint_lims, self%p_ptr%maxits_sh)
        endif
    end subroutine new

    subroutine prep4srch( self, os )
        class(strategy2D_srch), intent(inout) :: self
        class(oris),            intent(inout) :: os
        real    :: corrs(self%nrots)
        logical :: has_been_searched
        self%nrefs_eval = 0
        self%nsolns     = 0
        self%has_continuous_e3 = .false.
        self%continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED
        self%ithr       = omp_get_thread_num() + 1
        ! the shift prior on the polish of the committed class (the prob modes), at the particle's temperature
        if( self%p_ptr%l_shift_prior .and. self%continuous_active .and. s2D%shift_prior_sigma2 > 0. ) &
            &call self%joint_inpl_optimizer%set_shift_prior(s2D%shift_prior_sigma2, &
            &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(self%iptcl), self%b_ptr%pftc%get_pftsz()))
        ! find previous discrete alignment parameters
        self%l_fresh_start = is_fresh_2D_start(self%p_ptr, self%p_ptr%which_iter, os%isthere(self%iptcl,'class'))
        if( self%l_fresh_start )then
            self%prev_class    = 0
            self%prev_class_mi = 0
            self%prev_rot      = 1
            self%prev_shvec    = 0.
        else
            self%prev_class_mi = nint(os%get(self%iptcl,'class'))            ! class index before any fallback
            self%prev_class    = self%prev_class_mi
            self%prev_rot      = self%b_ptr%pftc%get_roind(360.-os%e3get(self%iptcl)) ! in-plane angle index
            if( self%prev_rot < 1 .or. self%prev_rot > self%nrots ) &
                THROW_HARD('Invalid previous in-plane rotation index')
            self%prev_shvec = os%get_2Dshift(self%iptcl)                  ! shift vector
        endif
        self%best_shvec = 0.
        if( self%prev_class > 0 )then
            if( s2D%cls_pops(self%prev_class) > 0 )then
                ! all done
            else
                ! for limiting cases
                call pick_populated_class
            endif
        else
            ! initialization
            call pick_populated_class
        endif
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better)
        call self%b_ptr%pftc%gen_objfun_vals(self%prev_class, self%iptcl, [0.,0.], corrs)
        if( self%p_ptr%cc_objfun == OBJFUN_CC )then
            self%prev_corr  = max(0., corrs(self%prev_rot))
        else
            self%prev_corr  = corrs(self%prev_rot)
        endif
        self%best_corr = self%prev_corr
        ! whether to search shifts first
        has_been_searched = (.not. self%l_fresh_start) .and. os%has_been_searched(self%iptcl)
        self%l_sh_first   = s2D%do_inplsrch(self%iptcl_batch) .and. self%p_ptr%l_doshift &
            &.and. has_been_searched
        ! under the key no seed at the previous class: the candidates are scored at the stored shift and the map
        ! locates the shift (ruling 20)
        self%l_map          = self%l_map_run .and. s2D%do_inplsrch(self%iptcl_batch)
        self%l_peaks_mapped = .false.
        if( self%l_map )then
            self%l_sh_first = .false.
            call self%smap%prep(self%iptcl)
        endif
        self%xy_first     =  0.
        self%xy_first_rot =  0.
        ! init per-thread class-space arrays
        call prep_strategy2D_thread(self%ithr)
        ! under the key the previous class competes in the screen at its stored-shift score
        if( self%l_map ) call self%store_solution(self%prev_class, self%prev_rot, self%prev_corr)
    contains

        !> a populated class drawn uniformly until one is found: on the probabilistic routes keyed on the particle
        !! (family ASSIGN with the pick's ordinal from 1; the assignment draw is ordinal 0), so the pick depends on
        !! the particle and the iteration alone; on the classic routes the intrinsic stream, as before
        subroutine pick_populated_class
            integer :: ipick
            ipick = 0
            do
                ipick = ipick + 1
                if( self%p_ptr%l_prob_align_mode )then
                    self%prev_class = rnd_keyed_int(rnd_draw_key(prob_rnd_scope(self%p_ptr%which_iter, .true.), &
                        &self%iptcl, RND_FAM_ASSIGN, ordinal=ipick), self%nrefs)
                else
                    self%prev_class = irnd_uni(self%nrefs)
                endif
                if( s2D%cls_pops(self%prev_class) > 0 ) exit
            enddo
        end subroutine pick_populated_class

    end subroutine prep4srch

    subroutine inpl_srch_first( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3), rotmat(2,2)
        integer :: irot
        self%best_shvec = [0.,0.]
        if( .not. self%l_sh_first ) return
        irot = 0
        call self%grad_shsrch_first_obj%set_indices(self%prev_class, self%iptcl)
        if( .not.self%grad_shsrch_first_obj%does_opt_angle() )then
            ! shift-only optimization
            irot = self%prev_rot
        endif
        cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
            ! update best
            self%best_corr  = cxy(1)
            self%best_rot   = irot
            self%best_shvec = self%xy_first_rot
        endif
    end subroutine inpl_srch_first

    subroutine inpl_srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3)
        integer :: irot
        irot = 0
        if( self%l_map )then
            ! the screen of the scored classes and the committed pose (the peaks' screen already ran when it did)
            if( .not. self%l_peaks_mapped ) call self%map_srch(1, .false.)
            call commit_mapped_class
        else if( s2D%do_inplsrch(self%iptcl_batch) )then
            ! selection parity: legacy post-selection refinement runs
            ! unchanged under both inpl_cont values, so the committed
            ! discrete pose is identical to the inpl_cont=no route
            self%best_shvec = [0.,0.]
            call self%grad_shsrch_obj%set_indices(self%best_class, self%iptcl)
            if( .not.self%grad_shsrch_obj%does_opt_angle() )then
                ! shift-only optimization
                irot = self%best_rot
            endif
            if( self%l_sh_first )then
                cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
            else
                cxy = self%grad_shsrch_obj%minimize(irot=irot)
            endif
            if( irot > 0 )then
                self%best_corr  = cxy(1)
                self%best_rot   = irot
                self%best_shvec = cxy(2:3)
            endif
        endif
        ! polish-only: exactly one local joint solve of the committed pose,
        ! after legacy refinement, independent of the shift-search schedule
        if( self%continuous_active ) call self%refine_selected_continuously

    contains

        !> the committed pose under the key: the class the peaks' sampling chose, else the best class entry after the
        !! screen (the mode's choice or a screened class that beat it); its mapped pose when the polish (or its cell)
        !! beat its stored-shift score, else the stored shift
        subroutine commit_mapped_class
            real    :: corr, sh(2), sdev
            integer :: irot_m, iref
            logical :: found, valid
            if( .not. self%l_peaks_mapped )then
                do iref = 1, self%nrefs
                    if( s2D%cls_pops(iref) <= 0 .or. s2D%class_space_inplinds(iref,self%ithr) <= 0 ) cycle
                    if( s2D%class_space_corrs(iref,self%ithr) > self%best_corr )then
                        self%best_class = iref
                        self%best_rot   = s2D%class_space_inplinds(iref,self%ithr)
                        self%best_corr  = s2D%class_space_corrs(iref,self%ithr)
                    endif
                enddo
            endif
            self%best_shvec = 0.
            call self%smap%get_result(self%best_class, found, valid, irot_m, corr, sh, sdev)
            if( found .and. valid )then
                self%best_rot   = irot_m
                self%best_corr  = corr
                self%best_shvec = sh
            endif
        end subroutine commit_mapped_class

    end subroutine inpl_srch

    !> The screen and the polish of the polar shift map (ruling 20) over the classes scored at the stored shift, the
    !! previous one included, as the tables screen their rows: one polished for the committed pose, npeaks for the
    !! peaks, whose entries then are the only valid ones, as in the legacy peaks. A polished pose that beat its
    !! class's stored-shift score replaces the class's entry
    subroutine map_srch( self, npeaks, l_peaks )
        class(strategy2D_srch), intent(inout) :: self
        integer,                intent(in)    :: npeaks
        logical,                intent(in)    :: l_peaks
        integer, allocatable :: refs(:), irots0(:)
        real,    allocatable :: corrs0(:)
        real    :: corr, sh(2), sdev
        integer :: iref, irot_m, i, n
        logical :: found, valid
        n = count(s2D%cls_pops > 0 .and. s2D%class_space_inplinds(:,self%ithr) > 0)
        allocate(refs(n), irots0(n), corrs0(n))
        n = 0
        do iref = 1, self%nrefs
            if( s2D%cls_pops(iref) <= 0 .or. s2D%class_space_inplinds(iref,self%ithr) <= 0 ) cycle
            n = n + 1
            refs(n)   = iref
            irots0(n) = s2D%class_space_inplinds(iref,self%ithr)
            corrs0(n) = s2D%class_space_corrs(iref,self%ithr)
        enddo
        if( n == 0 )then
            ! the committed class alone, also when the strategy stored no class entry
            if( allocated(refs) ) deallocate(refs, irots0, corrs0)
            refs   = [self%best_class]
            irots0 = [self%best_rot]
            corrs0 = [self%best_corr]
        endif
        call self%smap%srch(refs, corrs0, npeaks)
        if( l_peaks )then
            s2D%class_space_corrs(   :,self%ithr) = -1.
            s2D%class_space_inplinds(:,self%ithr) = 0
        endif
        do i = 1, size(refs)
            call self%smap%get_result(refs(i), found, valid, irot_m, corr, sh, sdev)
            if( .not. found ) cycle
            if( valid )then
                call self%store_solution(refs(i), irot_m, corr)
            else if( l_peaks )then
                call self%store_solution(refs(i), irots0(i), corrs0(i))
            endif
        enddo
        deallocate(refs, irots0, corrs0)
    end subroutine map_srch

    !> Jointly refine the shift and in-plane angle of one selected candidate.
    !! Probabilistic selection is complete before this routine is entered and is
    !! authoritative: the solve refines LOCALLY around the incoming integer angle
    !! and shift, with no global all-angle reselection -- rotationally
    !! self-similar class averages have near-degenerate in-plane branches that a
    !! global reselection at a slightly different shift hops on floating-point
    !! noise (see the identical 3D treatment in refine_assignment_continuously).
    !! An invalid solve retains the incoming assignment untouched; a valid
    !! non-improving solve retains the pose with its re-scored objective value.
    subroutine refine_selected_continuously( self, l_polish_shifts )
        class(strategy2D_srch), intent(inout) :: self
        logical, optional,      intent(in)    :: l_polish_shifts
        real    :: rotmat(2,2), xy_native(2)
        logical :: l_sh
        if( .not. self%continuous_active ) return
        ! shifts are polished only when the legacy schedule searched them this
        ! iteration; routes whose legacy behavior never refines shifts (e.g.
        ! refine=inpl_smpl) pass l_polish_shifts=.false. for strict parity
        l_sh = self%p_ptr%l_doshift .and. s2D%do_inplsrch(self%iptcl_batch)
        if( present(l_polish_shifts) ) l_sh = l_sh .and. l_polish_shifts
        if( self%best_class < 1 .or. self%best_class > self%nrefs )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( self%best_rot < 1 .or. self%best_rot > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        call rotmat2d(self%b_ptr%pftc%get_rot(self%best_rot), rotmat)
        xy_native = matmul(self%best_shvec, transpose(rotmat))
        block
        real     :: cxy(3), joint_lims(3,2), rbox
        real(dp) :: rotind_frac
        integer  :: irot
        logical  :: evaluation_valid, improved
        ! the joint domain contains the particle's box (its disk under the key) and the incoming shift, so the solve
        ! starts at, and is judged against, the pose it polishes
        rbox = self%trs
        if( self%p_ptr%l_centring ) rbox = shift_box_radius(self%b_ptr%spproj_field%get(self%iptcl, 'shsig'),&
            &self%p_ptr%shift_sigma, self%p_ptr%smpd, self%p_ptr%box, self%trs)
        rbox = max(rbox, maxval(abs(xy_native)))
        joint_lims(1,:) = [-rbox, rbox]
        joint_lims(2,:) = [-rbox, rbox]
        if( .not. l_sh )then
            ! pin the shifts and refine the angle only
            joint_lims(1,:) = xy_native(1)
            joint_lims(2,:) = xy_native(2)
        endif
        joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
        call self%joint_inpl_optimizer%set_indices(self%best_class, self%iptcl)
        call self%joint_inpl_optimizer%set_limits(joint_lims)
        cxy = self%joint_inpl_optimizer%minimize_joint(irot, xy_native, &
            &sh_rot=.true., rotind_frac=rotind_frac, &
            &evaluation_valid=evaluation_valid, improved=improved, irot_in=self%best_rot)
        if( .not. evaluation_valid )then
            ! retain the incoming assignment untouched
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        if( irot < 1 .or. irot > self%nrots )then
            self%continuous_route_outcome = CONT_ROUTE_INVALID
            return
        endif
        self%best_rot   = irot
        self%best_corr  = cxy(1)
        self%best_shvec = cxy(2:3)
        if( improved )then
            call store_continuous_e3(self, rotind_frac)
            self%continuous_route_outcome = CONT_ROUTE_IMPROVED
        else
            ! the seed is the incoming cell re-scored at its own shift, so this
            ! commit retains the incoming pose with a consistent objective value;
            ! has_continuous_e3 stays false so assign_ori writes the grid angle
            self%continuous_route_outcome = CONT_ROUTE_NO_IMPROVEMENT
        endif
        end block
    end subroutine refine_selected_continuously

    subroutine store_continuous_e3( self, rotind_frac )
        class(strategy2D_srch), intent(inout) :: self
        real(dp),                intent(in)    :: rotind_frac
        real(dp) :: e3
        e3 = (rotind_frac - 1.d0) * real(self%b_ptr%pftc%get_dang(),dp)
        self%continuous_e3 = real(modulo(360.d0 - e3,360.d0))
        self%has_continuous_e3 = .true.
    end subroutine store_continuous_e3

    pure logical function uses_continuous_refinement( self )
        class(strategy2D_srch), intent(in) :: self
        uses_continuous_refinement = self%continuous_active
    end function uses_continuous_refinement

    !> whether the joint route ran to a committed pose (improved or grid-seed
    !! commit) for the currently selected candidate
    pure logical function continuous_pose_committed( self )
        class(strategy2D_srch), intent(in) :: self
        continuous_pose_committed = &
            &self%continuous_route_outcome == CONT_ROUTE_IMPROVED .or. &
            &self%continuous_route_outcome == CONT_ROUTE_NO_IMPROVEMENT
    end function continuous_pose_committed

    pure subroutine get_continuous_route_status( self, attempted, improved, no_improvement, invalid )
        class(strategy2D_srch), intent(in)  :: self
        logical,                intent(out) :: attempted, improved, no_improvement, invalid
        attempted      = self%continuous_route_outcome /= CONT_ROUTE_NOT_ATTEMPTED
        improved       = self%continuous_route_outcome == CONT_ROUTE_IMPROVED
        no_improvement = self%continuous_route_outcome == CONT_ROUTE_NO_IMPROVEMENT
        invalid        = self%continuous_route_outcome == CONT_ROUTE_INVALID
    end subroutine get_continuous_route_status

    subroutine inpl_srch_peaks( self, npeaks_inpl )
        class(strategy2D_srch), intent(inout) :: self
        integer,                intent(in)    :: npeaks_inpl
        real    :: sorted_cls_corrs(self%nrefs), cxy(3)
        integer :: sorted_cls_inds(self%nrefs), saved_inpl_inds(self%nrefs), iref, isample, inpl_ind
        if( .not. s2D%do_inplsrch(self%iptcl_batch) ) return
        if( self%l_map )then
            call self%map_srch(npeaks_inpl, .true.)
            self%l_peaks_mapped = .true.
            return
        endif
        sorted_cls_corrs         = s2D%class_space_corrs(   :, self%ithr)
        saved_inpl_inds          = s2D%class_space_inplinds(:, self%ithr)
        sorted_cls_inds          = (/(iref,iref=1,self%nrefs)/)
        call hpsort(sorted_cls_corrs, sorted_cls_inds)
        ! reset class-space arrays so only shift-refined entries will be valid
        s2D%class_space_corrs(   :,self%ithr) = -1.
        s2D%class_space_inplinds(:,self%ithr) = 0
        do isample = self%nrefs-npeaks_inpl+1,self%nrefs
            iref     = sorted_cls_inds(isample)
            inpl_ind = saved_inpl_inds(iref)
            if( s2D%cls_pops(iref) == 0 ) cycle
            if( inpl_ind == 0 ) cycle
            call self%grad_shsrch_obj2%set_indices(iref, self%iptcl)
            if( self%l_sh_first )then
                cxy = self%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%xy_first)
                if( inpl_ind == 0 )then
                    inpl_ind = saved_inpl_inds(iref)
                    cxy(1)   = real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, self%iptcl, real(self%xy_first,dp), inpl_ind))
                    cxy(2:3) = self%xy_first_rot
                endif
            else
                cxy = self%grad_shsrch_obj2%minimize(irot=inpl_ind)
                if( inpl_ind == 0 )then
                    inpl_ind = saved_inpl_inds(iref)
                    cxy      = [real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, self%iptcl, inpl_ind)), 0.,0.]
                endif
            endif
            call self%store_solution(iref, inpl_ind, cxy(1))
        enddo
    end subroutine inpl_srch_peaks

    subroutine store_solution( self, ref, inpl_ind, corr )
        class(strategy2D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr
        if( inpl_ind < 1 .or. inpl_ind > self%nrots )then
            return
        endif
        if( s2D%class_space_corrs(ref, self%ithr) <= -huge(1.0)/2.0 )then
            self%nsolns = self%nsolns + 1
        elseif( self%continuous_pose_committed() .and. ref == self%best_class .and. &
            &inpl_ind == self%best_rot )then
            ! the joint continuous solve owns the final pose for this class
            ! (improved or grid-seed commit): its discrete index must be
            ! stored even when its score does not exceed the stored
            ! peak-search score (obtained at a shift that is discarded),
            ! otherwise assign_ori couples the committed joint shift with a
            ! stale integer angle and drops the fractional e3
        elseif( corr <= s2D%class_space_corrs(ref, self%ithr) )then
            return
        endif
        s2D%class_space_inplinds(ref,self%ithr) = inpl_ind
        s2D%class_space_corrs(   ref,self%ithr) = corr
    end subroutine store_solution

    subroutine assign_ori( self, os )
        class(strategy2D_srch), intent(in)    :: self
        class(oris),            intent(inout) :: os
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, frac
        real :: best_corr_local, shvec_incr(2), sh_m(2), corr_m, sdev
        integer :: iref, best_class_local, best_rot_local, irot_m
        logical :: found_valid, found, valid
        integer :: ithr
        ithr = self%ithr
        found_valid = .false.
        best_corr_local = -huge(1.0)
        if( self%prev_rot   <= 0 ) THROW_HARD('Previous in-plane rotation index is invalid, cannot assign orientation.')
        if( self%prev_class <= 0 ) THROW_HARD('Previous in-plane class index is invalid, cannot assign orientation.')
        best_class_local = self%prev_class
        best_rot_local   = self%prev_rot
        do iref = 1, self%nrefs
            if( s2D%cls_pops(iref) <= 0 ) cycle
            if( s2D%class_space_inplinds(iref,ithr) <= 0 ) cycle
            if( s2D%class_space_corrs(iref,ithr) > best_corr_local )then
                best_corr_local  = s2D%class_space_corrs(iref,ithr)
                best_class_local = iref
                best_rot_local   = s2D%class_space_inplinds(iref,ithr)
                found_valid      = .true.
            endif
        enddo
        if( .not. found_valid )then
            best_class_local = self%prev_class
            best_rot_local   = self%prev_rot
            best_corr_local  = self%prev_corr
        endif
        ! Keep the legacy integer rotation and continuous metadata coupled.
        ! Only the explicit Euclidean opt-in gateway may replace the grid angle
        ! with a continuous value; all other paths reconstruct e3 from inpl.
        if( self%has_continuous_e3 .and. best_class_local == self%best_class .and. &
            &best_rot_local == self%best_rot )then
            e3 = self%continuous_e3
        else
            e3 = 360. - self%b_ptr%pftc%get_rot(best_rot_local)
        endif
        ! calculate in-plane rot dist (radians) from the DISCRETE cells:
        ! statistics parity -- search-control statistics keep legacy grid
        ! semantics regardless of fractional pose precision, so convergence
        ! thresholds calibrated against the legacy route stay meaningful
        if( self%l_fresh_start )then
            dist = 0.
        else
            u(1) = 0.
            u(2) = 1.
            call rotmat2d(360. - self%b_ptr%pftc%get_rot(best_rot_local), mat)
            x1   = matmul(u,mat)
            call rotmat2d(360. - self%b_ptr%pftc%get_rot(self%prev_rot), mat)
            x2   = matmul(u,mat)
            dist = myacos(dot_product(x1,x2))
        endif
        ! calculate overlap between distributions
        mi_class = 0.
        if( (.not. self%l_fresh_start) .and. self%prev_class_mi == best_class_local ) mi_class = 1.
        ! search space explored
        frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        ! under the key the written class carries its own mapped shift (a peaks' class need not be the sampled one)
        shvec_incr = self%best_shvec
        sdev       = -1.
        if( self%l_map )then
            call self%smap%get_result(best_class_local, found, valid, irot_m, corr_m, sh_m, sdev)
            if( best_class_local /= self%best_class )then
                shvec_incr = 0.
                if( found .and. valid ) shvec_incr = sh_m
            endif
        endif
        ! update parameters
        call os%e3set(self%iptcl, e3)
        call os%set_shift(self%iptcl, self%prev_shvec + shvec_incr)
        call os%set(self%iptcl, 'shincarg',   arg(shvec_incr))
        call os%set(self%iptcl, 'inpl',       real(best_rot_local))
        call os%set(self%iptcl, 'class',      real(best_class_local))
        call os%set(self%iptcl, 'corr',       best_corr_local)
        call os%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call os%set(self%iptcl, 'mi_class',   mi_class)
        call os%set(self%iptcl, 'frac',       frac)
        ! the dynamic centring's state: this iteration's move and the prior width of the next search
        if( self%l_map ) call record_shift_state(os, self%iptcl, norm2(shvec_incr), sdev, self%p_ptr%shift_sigma,&
            &self%p_ptr%smpd, self%p_ptr%box, self%p_ptr%trs, self%b_ptr%pftc%get_kfromto())
    end subroutine assign_ori

    subroutine kill( self )
        class(strategy2D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%grad_shsrch_obj2%kill
        call self%grad_shsrch_first_obj%kill
        call self%joint_inpl_optimizer%kill
        if( allocated(self%smap) )then
            call self%smap%kill
            deallocate(self%smap)
        endif
        self%l_map_run         = .false.
        self%l_map             = .false.
        self%continuous_active = .false.
        self%has_continuous_e3 = .false.
        self%continuous_route_outcome = CONT_ROUTE_NOT_ATTEMPTED
    end subroutine kill

end module simple_strategy2D_srch
