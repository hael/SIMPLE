!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use omp_lib,                   only: omp_get_max_threads
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab_utils, only: update_screen_best, &
    &eulprob_dist_switch, record_seed_shift, &
    &prob_candidate, prob_dense_entry, prob_refined_entry, prob_fill_profile, &
    &prob_posterior_stat, posterior_temperature, posterior_row_stats, draw_posterior_row, &
    &prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN, top_posterior_extras, temperature_or_argmin, &
    &prob_insert_table, prob_insert_cone, screen_pool_size, shift_needs_disk, disk_radius_px, insert_pose_loss, &
    &prob_rnd_scope, calc_num2sample, shift_prior_sigma2, prob_readout, prob_readout_fname, write_prob_readout, &
    &report_prob_readout, PROB_READOUT_3D_DENSE
use simple_rnd,                only: rnd_key, rnd_keyed, rnd_draw_key, RND_FAM_ASSIGN
use simple_type_defs,          only: OBJFUN_EUCLID
use simple_defs,               only: PROB_INSERT_WMIN
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_polarft_calc,       only: SHERRSQ, SHIFTMAP_SCREEN, SHIFTMAP_WINDOW, NSHIFTMAPS
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"


type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
    real,           allocatable :: seed_shifts(:,:) !< per-particle seeded shift       (2,nptcls)
    real,           allocatable :: temps(:)         !< per-particle posterior temperature pi wsqsum / (2 pftsz) (nptcls; the fill writes it, the selection draws at it)
    logical,        allocatable :: seed_has_sh(:)   !< per-particle seeded shift flag  (nptcls)
    integer                     :: seed_nrots = 0   !< rotation grid of the table (checked when the parts are read)
    integer,        allocatable :: pinds(:)        !< particle indices for processing
    integer,        allocatable :: ssinds(:)       !< non-empty state indices
    integer,        allocatable :: state_to_active_rank(:) !< requested state -> compact active-state rank
    logical,        allocatable :: state_exists(:)
    integer                     :: post_unit      = -1      !< the pose-posterior diagnostic file (prob_post_diag=yes)
    logical                     :: post_is_open   = .false.
    logical                     :: post_slots_identical = .false. !< the two reference slots hold the same sections (header)
    logical                     :: post_has_halves = .false. !< the real halves are loaded beside the references (header)
    integer                     :: post_ref_offset = 0       !< the diagnostic model's first reference minus one: the full production model's size nspace*nstates (not nrefs, which counts the active states only)
    integer                     :: nptcls          !< size of pinds array
    integer                     :: nstates         !< states number
    integer                     :: nrefs           !< reference number
    type(prob_fill_profile)     :: fill_profile    !< accumulated cost record of fill_tab_range (bench file)
    type(prob_insert_table)     :: insert_tab      !< the weighted insertion's extra poses of the assignment (topk > 1)
    type(prob_insert_cone)      :: insert_cone     !< the weighted insertion's cone about the committed pose (set_insert_cone; unset: no extras)
    type(rnd_key)               :: rnd_scope       !< the keyed draws' scope of the iteration (prob_rnd_scope), taken at every fill and assignment
    integer                     :: projs_ns_model = 0      !< the iteration's refined-reference budget (calc_num2sample), set by the fill's caller
    real                        :: sigma2_prior_model = 0. !< the iteration's shift prior variance per axis (shift_prior_sigma2), set by the fill's caller
    logical                     :: l_model_set = .false.   !< set_iteration_model called
    ! the batch's rows, kept by the fill for the batch's own selection
    type(prob_dense_entry),   allocatable :: batch_rows(:,:)    !< (nrefs, i_first:i_last) the batch's rows
    type(prob_refined_entry), allocatable :: batch_refined(:,:) !< (max(1,projs_ns), i_first:i_last) its refined entries
    integer,                  allocatable :: batch_nrefined(:)  !< (i_first:i_last)
    ! the assignment across the batches (open_assignment .. close_assignment): the readout and the counts
    type(prob_neff_readout)     :: neff
    logical                     :: l_assign_open = .false.
    logical                     :: l_extras      = .false.  !< the weighted insertion's extras are selected (topk > 1, a cone)
    integer                     :: kmax = 1                 !< poses per particle of the weighted insertion
    integer                     :: nmode = 0, ntemp_missing = 0, nbad = 0, nflat = 0
    contains
    ! CONSTRUCTOR
    procedure :: new_pipeline
    procedure, private :: new_common
    procedure, private :: initialize_storage
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: set_iteration_model
    procedure :: set_model_from_field
    procedure :: fill_tab_range
    procedure, private :: keep_batch_rows
    procedure :: open_assignment
    procedure :: assign_batch
    procedure :: close_assignment
    procedure, private :: assign_rows
    procedure, private :: get_row
    procedure, private :: max_refined
    procedure, private :: write_posterior_diag
    procedure :: get_fill_profile
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: set_insert_cone
    ! REFERENCE INDEX MAPPING
    procedure :: ref_state
    procedure :: ref_proj
    procedure :: ref_full
    procedure :: full_to_compact_ref
    procedure :: assign_candidate
    ! DESTRUCTOR
    procedure :: kill
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new_common( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer, parameter :: MIN_POP = 5   ! ignoring cavgs with less than 5 particles
        integer :: istate, si
        call self%kill
        self%p_ptr => params
        self%b_ptr  => build
        ! the posterior is defined by the Euclidean objective (the particle temperature is its noise scale)
        if( params%cc_objfun /= OBJFUN_EUCLID )&
            &THROW_HARD('the probability tables draw from the tempered posterior: refine=prob|prob_neigh needs objfun=euclid')
        self%nptcls       = size(pinds)
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        ! weighted: a state is populated by the effective sample size of its frozen weights
        if( self%p_ptr%l_m_estimator_flex ) call weighted_states_exist
        self%nstates      = count(self%state_exists .eqv. .true.)
        self%nrefs = self%p_ptr%nspace * self%nstates
        allocate(self%ssinds(self%nstates), self%state_to_active_rank(self%p_ptr%nstates), source=0)
        si = 0
        do istate = 1, self%p_ptr%nstates
            if( .not. self%state_exists(istate) )cycle
            si                                = si + 1
            self%ssinds(si)                   = istate
            self%state_to_active_rank(istate) = si
        enddo
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
        allocate(self%temps(self%nptcls), source=0.)

    contains

        subroutine weighted_states_exist
            use simple_state_weight_set, only: state_weight_set
            type(state_weight_set) :: wset
            call wset%new(self%b_ptr%spproj, self%b_ptr%spproj_field)
            if( wset%get_nstates() /= self%p_ptr%nstates ) THROW_HARD('state weight set and nstates disagree; eul_prob_tab')
            do istate = 1, self%p_ptr%nstates
                self%state_exists(istate) = wset%get_ess(istate) > real(MIN_POP)
            enddo
            call wset%kill
        end subroutine weighted_states_exist

    end subroutine new_common

    subroutine initialize_storage( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i
        real    :: x
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            self%assgn_map(i)%pind   = self%pinds(i)
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
        end do
        !$omp end parallel do
    end subroutine initialize_storage

    !> the table of a part (the matcher's, in shared memory and in a distributed part alike): the assignment map and,
    !! batch by batch, the rows the fill keeps for the batch's own selection (fill_tab_range, assign_batch)
    subroutine new_pipeline( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls))
        call self%initialize_storage
        self%seed_nrots = build%pftc%get_nrots()
    end subroutine new_pipeline

    !> the iteration's model from the complete incoming field (section 3.2 of the pipeline note): the refined-reference
    !! budget, the largest calc_num2sample over the active states, and the shift prior's width (shift_prior=yes)
    subroutine set_model_from_field( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: si, n, projs_ns
        real    :: sigma2_prior
        projs_ns = 0
        do si = 1, self%nstates
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, &
                &state=self%ssinds(si))
            projs_ns = max(projs_ns, n)
        enddo
        sigma2_prior = 0.
        if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
        call self%set_iteration_model(projs_ns, sigma2_prior)
    end subroutine set_model_from_field

    !> The iteration's model the fill reads instead of the particle field: the refined-reference budget projs_ns
    !! (the largest calc_num2sample over the active states) and the shift prior's variance per axis
    !! (shift_prior_sigma2), both computed once per iteration by the fill's caller from the complete incoming
    !! field, so that a batch's rows never depend on what earlier batches committed
    subroutine set_iteration_model( self, projs_ns, sigma2_prior )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: projs_ns
        real,                intent(in)    :: sigma2_prior
        self%projs_ns_model     = projs_ns
        self%sigma2_prior_model = sigma2_prior
        self%l_model_set        = .true.
    end subroutine set_iteration_model

    subroutine fill_tab_range( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        integer, allocatable   :: locn(:,:)
        type(pftc_shsrch_grad), target :: grad_shsrch_obj(nthr_glob)  !< the shift search object, L-BFGS with gradient
        type(pftc_shsrch_grad), target :: grad_shsrch_disk(nthr_glob) !< the disk's (shift_sigma > 0): its limits bound the polish box
        type(prob_dense_entry),   allocatable :: rows_batch(:,:)     !< the batch's rows (nrefs,i_first:i_last)
        type(prob_refined_entry), allocatable :: refined_batch(:,:)  !< refined entries per particle (projs_ns,i_first:i_last)
        integer,                  allocatable :: nrefined_batch(:)
        type(prob_posterior_stat), allocatable :: post_batch(:)      !< pose-posterior statistics per particle (prob_post_diag=yes)
        real,                      allocatable :: post_w(:,:)        !< per thread: the posterior over the row, sorted
        integer,                   allocatable :: post_inds(:,:)     !< per thread: its reference order
        real,                      allocatable :: post_d(:,:)        !< per thread: the deterministic row, own slot
        real,                      allocatable :: post_wd(:,:)       !< per thread: its posterior, sorted
        integer,                   allocatable :: post_indsd(:,:)    !< per thread: its reference order
        real,                      allocatable :: post_h(:,:)        !< per thread: the deterministic row, own real half
        real,                      allocatable :: post_wh(:,:)       !< per thread: its posterior, sorted
        integer,                   allocatable :: post_indsh(:,:)    !< per thread: its reference order
        real,                      allocatable :: post_x(:,:)        !< per thread: the deterministic row, other real half
        real,                      allocatable :: post_wx(:,:)       !< per thread: its posterior, sorted
        integer,                   allocatable :: post_indsx(:,:)    !< per thread: its reference order
        logical :: l_post
        integer :: i, iptcl, projs_ns, ithr
        integer :: inpl_refs(self%nrefs,nthr_glob)
        real    :: lims(2,2), lims_init(2,2)
        real    :: dists_refs(self%nrefs,nthr_glob)
        logical :: l_sh_first, l_grid
        real    :: sigma2_prior  !< the shift prior's variance per axis (shift_prior=yes)
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo (single state)
        integer(int64) :: seed_builds, bank_ithr
        integer :: kfromto_loc(2)
        integer(int64) :: fit_ithr(0:5)
        integer(int64) :: rank_ithr(0:11,nthr_glob)  !< per thread: the refined argmin's unrefined rank (bench, prob_fill_profile%rank_counts)
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr
        integer(int64) :: ndisk_thr(nthr_glob), nscreens_thr(nthr_glob)
        real(dp)       :: prior_w, sigma_px, radius_disk
        real           :: lims_disk(2,2)
        real           :: step_map(NSHIFTMAPS)
        real           :: scr_xy(2,self%nrefs,nthr_glob), scr_loss(self%nrefs,nthr_glob) !< per thread: the screen's coarse optima and losses
        real           :: scr_sdev(self%nrefs,nthr_glob)  !< per thread: the screen's posterior widths
        integer        :: scr_irot(self%nrefs,nthr_glob)
        integer        :: kscreen
        logical        :: l_map
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first )then
            THROW_HARD('invalid particle range in eul_prob_tab%fill_tab_range')
        endif
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        rank_ithr = 0_int64
        sigma2_prior = 0.
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        l_sh_first      = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1
        ! the particle-frame grid and its quadratic step in the multi-state table, the particle-frame seed in the
        ! seeded single-state table
        l_grid          = self%p_ptr%nstates > 1
        seed_builds     = 0_int64
        ndisk_thr       = 0_int64
        nscreens_thr    = 0_int64
        l_map           = self%p_ptr%l_centring
        prior_w         = 0.d0
        kscreen         = 0
        step_map        = 0.
        if( .not. self%l_model_set ) THROW_HARD('the iteration model is set by the caller before the fill (set_iteration_model); eul_prob_tab%fill_tab_range')
        ! the keyed draws of this fill: the in-plane draw of every scored reference on the calculator
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .false.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
        projs_ns = self%projs_ns_model
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), source=0)
        allocate(rows_batch(self%nrefs,i_first:i_last), refined_batch(max(1,projs_ns),i_first:i_last))
        allocate(nrefined_batch(i_first:i_last), source=0)
        ! the pose posterior is defined by the Euclidean objective only
        l_post = self%p_ptr%l_prob_post_diag .and. self%p_ptr%cc_objfun == OBJFUN_EUCLID
        if( l_post )then
            allocate(post_batch(i_first:i_last), post_w(self%nrefs,nthr_glob), post_inds(self%nrefs,nthr_glob), &
                &post_d(self%nrefs,nthr_glob), post_wd(self%nrefs,nthr_glob), post_indsd(self%nrefs,nthr_glob), &
                &post_h(self%nrefs,nthr_glob), post_wh(self%nrefs,nthr_glob), post_indsh(self%nrefs,nthr_glob), &
                &post_x(self%nrefs,nthr_glob), post_wx(self%nrefs,nthr_glob), post_indsx(self%nrefs,nthr_glob))
            ! solve3D's LP-set stages put one merged map in both of the matcher's slots (recorded in the header);
            ! the real halves are the diagnostic model loaded after the matcher's references, when there is one
            self%post_slots_identical = self%b_ptr%pftc%ref_slots_identical(self%ref_full(1))
            ! the calculator holds the full production model (every requested state) and, in diagnostic mode,
            ! the halves' model after it, so the second block starts at nspace*nstates whatever the active
            ! states, and it is there exactly when the calculator holds twice that
            self%post_ref_offset = self%p_ptr%nspace * self%p_ptr%nstates
            self%post_has_halves = self%b_ptr%pftc%get_nrefs() == 2 * self%post_ref_offset
        endif
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new_alternating(self%b_ptr, lims, lims_init=lims_init, &
                    &maxits=self%p_ptr%maxits_sh, coarse_init=.true., grid_fit=l_grid, grid_ptcl=l_grid)
            end do
            ! the shift prior's width of the iteration (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = self%sigma2_prior_model
            if( l_sh_first )then
                kfromto_loc = self%b_ptr%pftc%get_kfromto()
                allocate(seed_bank(self%b_ptr%pftc%get_pftsz()+1, kfromto_loc(2)-kfromto_loc(1)+1, nthr_glob))
            endif
            ! the polar shift map (shift_sigma > 0): the coarse map over a particle's box at its likely references, the
            ! polish at full resolution from its optimum. The box is the disk of three sigma under the prior for a
            ! particle not yet centred in this workflow or whose prior width exceeds a third of the window (the
            ! dynamic centring: shift_needs_disk), the window of trs for the others
            if( l_map )then
                if( .not. (self%p_ptr%prob_beta > 0.) ) THROW_HARD('the shift map draws at prob_beta x T: prob_beta must be positive')
                if( l_post ) THROW_HARD('prob_post_diag=yes is not implemented for the polar shift map (shift_sigma)')
                sigma_px    = real(self%p_ptr%shift_sigma, dp) / real(self%p_ptr%smpd, dp)
                radius_disk = real(disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd, self%p_ptr%box), dp)
                call self%b_ptr%pftc%setup_shiftmaps(sigma_px, radius_disk, real(self%p_ptr%trs, dp), prior_w, kscreen)
                ! each map's radial step: the box the polish starts in and stays in, about the coarse optimum
                step_map(SHIFTMAP_SCREEN) = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 2, 1) -&
                    &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 1, 1)))
                step_map(SHIFTMAP_WINDOW) = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 2, 1) -&
                    &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 1, 1)))
                ! the disk's own search object: the polish box clamps to the object's limits
                lims_disk(:,1) = -real(radius_disk)
                lims_disk(:,2) =  real(radius_disk)
                do ithr = 1,nthr_glob
                    call grad_shsrch_disk(ithr)%new_alternating(self%b_ptr, lims_disk, lims_init=lims_init, &
                        &maxits=self%p_ptr%maxits_sh, coarse_init=.true., grid_fit=l_grid, grid_ptcl=l_grid)
                end do
            endif
            ! fill the table
            t_region = omp_get_wtime()
            ! A particle's cost depends on its shift searches, so the schedule is dynamic. Each thread owns
            ! its search object and a particle is processed whole, writing only its own row; its draws are
            ! keyed on the particle, so the thread it meets changes nothing.
            !$omp parallel do default(shared) private(i,iptcl,ithr) proc_bind(close) schedule(dynamic)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_with_shift(i, iptcl, ithr)
            enddo
            !$omp end parallel do
            self%fill_profile%wall_region = self%fill_profile%wall_region + (omp_get_wtime() - t_region)
        else
            ! fill the table
            t_region = omp_get_wtime()
            !$omp parallel do default(shared) private(i,iptcl,ithr) proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_without_shift(i, iptcl, ithr)
            enddo
            !$omp end parallel do
            self%fill_profile%wall_region = self%fill_profile%wall_region + (omp_get_wtime() - t_region)
        endif
        ! cost record; the search objects' counters are read before kill resets them
        self%fill_profile%ts_seed   = self%fill_profile%ts_seed   + sum(ts_seed)
        self%fill_profile%ts_score  = self%fill_profile%ts_score  + sum(ts_score)
        self%fill_profile%ts_refine = self%fill_profile%ts_refine + sum(ts_refine)
        self%fill_profile%nthr      = nthr_glob
        self%fill_profile%nrefs     = self%nrefs
        self%fill_profile%nptcls    = self%fill_profile%nptcls + int(i_last - i_first + 1, int64)
        if( self%p_ptr%l_doshift )then
            self%fill_profile%projs_ns = max(self%fill_profile%projs_ns, projs_ns)
            self%fill_profile%nrefined = self%fill_profile%nrefined + int(i_last - i_first + 1, int64) * int(projs_ns, int64)
            self%fill_profile%rank_counts = self%fill_profile%rank_counts + sum(rank_ithr, dim=2)
        endif
        do ithr = 1,nthr_glob
            call add_counters(grad_shsrch_obj(ithr))
            if( l_map ) call add_counters(grad_shsrch_disk(ithr))
        end do
        self%fill_profile%bank_builds = self%fill_profile%bank_builds + seed_builds
        self%fill_profile%screens     = self%fill_profile%screens    + sum(nscreens_thr)
        self%fill_profile%disk_ptcls  = self%fill_profile%disk_ptcls + sum(ndisk_thr)
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
            if( l_map ) call grad_shsrch_disk(ithr)%kill
        end do
        if( allocated(seed_bank) ) deallocate(seed_bank)
        call self%keep_batch_rows(rows_batch, refined_batch, nrefined_batch)
        if( l_post )then
            call self%write_posterior_diag(i_first, i_last, post_batch)
            deallocate(post_batch, post_w, post_inds, post_d, post_wd, post_indsd, post_h, post_wh, post_indsh, &
                &post_x, post_wx, post_indsx)
        endif

    contains

        !> one search object's counters into the cost record (read before kill resets them)
        subroutine add_counters( gso )
            type(pftc_shsrch_grad), intent(inout) :: gso
            call gso%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
            self%fill_profile%objective_evals = self%fill_profile%objective_evals + obj_ithr
            self%fill_profile%gradient_evals  = self%fill_profile%gradient_evals  + grad_ithr
            self%fill_profile%alternations    = self%fill_profile%alternations    + alt_ithr
            self%fill_profile%grid_points     = self%fill_profile%grid_points     + grid_ithr
            self%fill_profile%minimize_calls  = self%fill_profile%minimize_calls  + min_ithr
            call gso%get_fit_profile(fit_ithr)
            bank_ithr = gso%get_bank_builds()
            self%fill_profile%bank_builds = self%fill_profile%bank_builds + bank_ithr
            self%fill_profile%fit_counts  = self%fill_profile%fit_counts  + fit_ithr
        end subroutine add_counters

        subroutine process_particle_with_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            type(ori) :: o_prev_loc
            integer   :: prev_state, prev_proj, ri_loc, irot_loc, irot_seed
            real      :: shift_seed(3), dist_loc, shsig_loc
            real(dp)  :: t0
            logical   :: l_bank_seed, l_disk
            if( l_map )then
                ! the disk for a particle not yet centred in this workflow (under the map's prior) or whose prior width
                ! exceeds a third of the window (the data's alone); the window for the others
                shsig_loc = self%b_ptr%spproj_field%get(iptcl_loc, 'shsig')
                l_disk    = shift_needs_disk(shsig_loc, self%p_ptr%trs)
                call process_particle_with_map(i_loc, iptcl_loc, ithr_loc, l_disk, shsig_loc < 0.)
                return
            endif
            call get_particle_context(iptcl_loc, o_prev_loc, prev_state, prev_proj)
            t0 = omp_get_wtime()
            ! the shift prior at the particle's own temperature: the seed estimate and every refinement of the particle
            if( self%p_ptr%l_shift_prior ) call grad_shsrch_obj(ithr_loc)%set_shift_prior(sigma2_prior,&
                &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz()))
            call estimate_shift_seed(iptcl_loc, ithr_loc, prev_state, prev_proj, o_prev_loc, shift_seed, irot_seed)
            ! the seed is one particle-frame shift from here on: the estimate converted at its rotation (a zero seed
            ! without the single-state seeding)
            if( irot_seed > 0 ) shift_seed(2:3) = shift_ref2ptcl(shift_seed(2:3), self%b_ptr%pftc%get_rot(irot_seed))
            l_bank_seed = sum(shift_seed(2:3)**2) > SHERRSQ
            if( l_bank_seed )then
                call self%b_ptr%pftc%memoize_shifted_ptcl(iptcl_loc, shift_seed(2:3), seed_bank(:,:,ithr_loc))
                !$omp atomic
                seed_builds = seed_builds + 1_int64
            endif
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            self%seed_shifts(:,i_loc) = shift_seed(2:3)
            self%seed_has_sh(i_loc)   = l_sh_first
            t0 = omp_get_wtime()
            do ri_loc = 1,self%nrefs
                if( l_bank_seed )then
                    call score_ref_bank(ri_loc, iptcl_loc, ithr_loc, dist_loc, irot_loc)
                else
                    call score_ref(ri_loc, iptcl_loc, shift_seed(2:3), dist_loc, irot_loc)
                endif
                rows_batch(ri_loc,i_loc) = prob_dense_entry(dist_loc, irot_loc)
                dists_refs(ri_loc,ithr_loc) = dist_loc
                inpl_refs(ri_loc,ithr_loc)  = irot_loc
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            call refine_best_refs(i_loc, iptcl_loc, ithr_loc, shift_seed)
            ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            if( l_post ) call posterior_diag(i_loc, iptcl_loc, ithr_loc, shift_seed(2:3), l_bank_seed)
            call o_prev_loc%kill
        end subroutine process_particle_with_shift

        !> The map's path (shift_sigma > 0): every reference scored at the stored shift (the sections carry it); the
        !! coarse map over the box at the references holding SCREEN_QUANTILE of the particle's posterior at the origin
        !! (a pick far off leaves it flat and keeps every reference), at least twice the refined count, each keeping
        !! its coarse optimum; the projs_ns best of the screen polished at full resolution in a box of one coarse step
        !! about that optimum (minimize_in_cell). The box is the disk (l_disk), under the map's prior at the particle's
        !! first centring (l_prior) and without it after, or the window without it. Rows the polished losses, as the
        !! window path's; the refined entries carry the posterior loss for the insertion
        subroutine process_particle_with_map( i_loc, iptcl_loc, ithr_loc, l_disk, l_prior )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in) :: l_disk, l_prior
            type(pftc_shsrch_grad), pointer :: gso
            integer, allocatable :: pool(:)
            real(dp) :: t0, temp_i, temp_s, xcell(2), prior_w_map
            real     :: loss_cell, loss_pose, dist_loc, cxy(3), width_best(max(1,projs_ns))
            integer  :: imap, ri_loc, irot_loc, irho_loc, ialpha_loc, j_loc, k_loc, npool, nbest
            logical  :: l_polished
            imap        = SHIFTMAP_WINDOW
            prior_w_map = 0.d0
            gso         => grad_shsrch_obj(ithr_loc)
            if( l_disk )then
                imap        = SHIFTMAP_SCREEN
                if( l_prior ) prior_w_map = prior_w
                gso         => grad_shsrch_disk(ithr_loc)
                ndisk_thr(ithr_loc) = ndisk_thr(ithr_loc) + 1_int64
            endif
            t0 = omp_get_wtime()
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            temp_i = real(self%p_ptr%prob_beta, dp) * real(self%temps(i_loc), dp)
            temp_s = real(self%p_ptr%prob_beta, dp) * real(posterior_temperature(&
                &self%b_ptr%pftc%wsqsum_ptcl_band(iptcl_loc, kscreen), self%b_ptr%pftc%get_pftsz()), dp)
            self%seed_shifts(:,i_loc) = 0.
            self%seed_has_sh(i_loc)   = .true.
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            do ri_loc = 1,self%nrefs
                call score_ref(ri_loc, iptcl_loc, [0.,0.], dist_loc, irot_loc)
                rows_batch(ri_loc,i_loc)    = prob_dense_entry(dist_loc, irot_loc)
                dists_refs(ri_loc,ithr_loc) = dist_loc
                inpl_refs(ri_loc,ithr_loc)  = irot_loc
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            ! the screen: the coarse map at the origin's best references, each keeping its coarse optimum
            pool  = minnloc(dists_refs(:,ithr_loc), self%nrefs)
            npool = min(self%nrefs, max(2 * projs_ns, screen_pool_size(dists_refs(pool,ithr_loc), self%nrefs, temp_i)))
            scr_loss(:,ithr_loc) = huge(1.0)
            width_best = huge(1.)
            do k_loc = 1, npool
                ri_loc = pool(k_loc)
                if( ri_loc < 1 ) cycle
                call self%b_ptr%pftc%gen_shiftmap(imap, self%ref_full(ri_loc), iptcl_loc)
                call self%b_ptr%pftc%shiftmap_mode(imap, temp_s, prior_w_map, &
                    &irot_loc, irho_loc, ialpha_loc, loss_cell, &
                    &temp_phys=temp_s / real(self%p_ptr%prob_beta, dp), sdev=scr_sdev(ri_loc,ithr_loc), &
                    &width_cutoff=maxval(width_best))
                call update_screen_best(loss_cell, width_best)
                scr_irot(ri_loc,ithr_loc) = irot_loc
                scr_xy(:,ri_loc,ithr_loc) = real(self%b_ptr%pftc%shiftmap_shift(imap, irho_loc, ialpha_loc))
                scr_loss(ri_loc,ithr_loc) = loss_cell
                nscreens_thr(ithr_loc) = nscreens_thr(ithr_loc) + 1_int64
            enddo
            deallocate(pool)
            ! the polish of the screen's best: the box's optimum if it beat the origin, minimize's baseline; else the
            ! coarse optimum, scored at full resolution, if it beats the origin's row; else that row stands
            nbest = min(projs_ns, npool)
            locn(:,ithr_loc) = 0
            locn(1:nbest,ithr_loc) = minnloc(scr_loss(:,ithr_loc), nbest)
            do j_loc = 1, nbest
                ri_loc = locn(j_loc,ithr_loc)
                if( ri_loc < 1 ) cycle
                if( scr_loss(ri_loc,ithr_loc) >= huge(1.0)/2. ) cycle
                call gso%set_indices(self%ref_full(ri_loc), iptcl_loc)
                irot_loc = scr_irot(ri_loc,ithr_loc)
                call gso%minimize_in_cell(irot_loc, scr_xy(:,ri_loc,ithr_loc), step_map(imap), cxy, l_polished)
                dist_loc = eulprob_dist_switch(cxy(1), self%p_ptr%cc_objfun)
                if( .not. l_polished .and. .not. (dist_loc < rows_batch(ri_loc,i_loc)%dist) ) cycle
                rows_batch(ri_loc,i_loc) = prob_dense_entry(dist_loc, irot_loc)
                xcell = real(cxy(2:3), dp)
                nrefined_batch(i_loc) = nrefined_batch(i_loc) + 1
                ! the pose's loss for the insertion: the data term with the map's prior on the increment at the physical
                ! temperature
                loss_pose = insert_pose_loss(dist_loc, xcell, prior_w, self%temps(i_loc))
                refined_batch(nrefined_batch(i_loc),i_loc) = prob_refined_entry(i_loc, ri_loc, real(xcell(1)), real(xcell(2)), loss_pose,&
                    &scr_sdev(ri_loc,ithr_loc))
            enddo
            ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
        end subroutine process_particle_with_map

        subroutine process_particle_without_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            integer  :: ri_loc, irot_loc
            real     :: dist_loc
            real(dp) :: t0
            t0 = omp_get_wtime()
            do ri_loc = 1,self%nrefs
                call score_ref(ri_loc, iptcl_loc, [0.,0.], dist_loc, irot_loc)
                rows_batch(ri_loc,i_loc) = prob_dense_entry(dist_loc, irot_loc)
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            if( l_post ) call posterior_diag(i_loc, iptcl_loc, ithr_loc, [0.,0.], .false.)
        end subroutine process_particle_without_shift

        !> The pose posterior over the particle's finished production row (posterior-drawn in-plane,
        !! the best projs_ns shift-refined) at the temperature of its weighted power, with the angular
        !! spread of the references holding 99% of the mass about the mode (the mode's state, across the
        !! symmetry operations). Then the controls: the identical-pose null test, the kernel re-scored at the production mode's
        !! exact (reference, shift, grid rotation) against the stored loss, with the kernel the stored value
        !! came from (the shift search's shifted-reference kernel at the frame-converted stored shift, or the
        !! seed bank); the deterministic row on the matcher's references (best grid rotation at the
        !! particle-frame seed, through the bank when the particle has one, no sampling, no refinement), the
        !! search test; and, when the real halves are loaded beside the references, the same deterministic
        !! row on the particle's own real half and on the other one, the reference test, with the
        !! self-reinforcement test at the pose the particle went into its half with.
        subroutine posterior_diag( i_loc, iptcl_loc, ithr_loc, shift_xy, l_bank )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            real,    intent(in) :: shift_xy(2)
            logical, intent(in) :: l_bank
            type(ori) :: o_prev_loc
            logical   :: l_refined
            real     :: temp, pmax, neff, euls_min(3), euls_d(3), euls_h(3), euls_x(3), euls_sym(3), euldist, inpldist
            real     :: sh_mode(2), sh_kernel(2), dist_loc, dist_at, dist_own, dist_other
            integer  :: imin, k99, nvalid, ri_loc, istate_min, iproj_min, irot_loc, irot_min, j_loc
            integer  :: imin_d, k99_d, nvalid_d, imin_h, k99_h, nvalid_h, imin_x, k99_x, nvalid_x
            integer  :: prev_state, prev_proj, irot_prev, irot_ctrl, iref_prev, nrots_loc
            temp = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            call posterior_row_stats(self%nrefs, rows_batch(:,i_loc)%dist, temp, post_w(:,ithr_loc), &
                &post_inds(:,ithr_loc), imin, pmax, neff, k99, nvalid)
            post_batch(i_loc)%iptcl  = iptcl_loc
            post_batch(i_loc)%nvalid = nvalid
            post_batch(i_loc)%temp   = temp
            post_batch(i_loc)%pmax   = pmax
            post_batch(i_loc)%neff   = neff
            post_batch(i_loc)%k99    = k99
            if( nvalid == 0 ) return
            istate_min = self%ref_state(imin)
            iproj_min  = self%ref_proj(imin)
            irot_min   = rows_batch(imin,i_loc)%inpl
            post_batch(i_loc)%istate = istate_min
            post_batch(i_loc)%iproj  = iproj_min
            post_batch(i_loc)%lmin   = rows_batch(imin,i_loc)%dist
            euls_min = self%b_ptr%eulspace%get_euler(iproj_min)
            post_batch(i_loc)%spread = spread_about(istate_min, euls_min, k99, post_w(:,ithr_loc), post_inds(:,ithr_loc))
            ! the null test, with the kernel the stored value came from. A refined mode: the shift search's
            ! shifted-reference kernel (gen_euclids) at the stored shift, which the search returns in the
            ! particle frame (shift_ref2ptcl of the kernel-frame minimizer), converted back with shift_ptcl2ref;
            ! null_b is the stored shift fed to the kernel as it is, the wrong frame. An unrefined mode with a
            ! seed bank: the bank kernel at the stored rotation; null_b the shifted-reference kernel at the
            ! frame-converted seed, which differs by the angular Nyquist bin the latter drops. Neither: the
            ! unshifted kernel for both.
            l_refined = .false.
            sh_mode   = shift_xy
            do j_loc = 1,nrefined_batch(i_loc)
                if( refined_batch(j_loc,i_loc)%ri == imin )then
                    sh_mode   = [refined_batch(j_loc,i_loc)%x, refined_batch(j_loc,i_loc)%y]
                    l_refined = .true.
                    exit
                endif
            end do
            if( irot_min >= 1 )then
                if( l_refined )then
                    post_batch(i_loc)%mode_src = 1
                    sh_kernel = shift_ptcl2ref(sh_mode, self%b_ptr%pftc%get_rot(irot_min))
                    call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, sh_kernel, irot_min, dist_at)
                    post_batch(i_loc)%null_a = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                    call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, sh_mode, irot_min, dist_at)
                    post_batch(i_loc)%null_b = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                else if( l_bank )then
                    post_batch(i_loc)%mode_src = 2
                    call self%b_ptr%pftc%gen_euclid_val_at_bank(self%ref_full(imin), iptcl_loc, seed_bank(:,:,ithr_loc), irot_min, dist_at)
                    post_batch(i_loc)%null_a = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                    sh_kernel = shift_ptcl2ref(shift_xy, self%b_ptr%pftc%get_rot(irot_min))
                    call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, sh_kernel, irot_min, dist_at)
                    post_batch(i_loc)%null_b = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                else
                    post_batch(i_loc)%mode_src = 3
                    call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, [0.,0.], irot_min, dist_at)
                    post_batch(i_loc)%null_a = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                    post_batch(i_loc)%null_b = post_batch(i_loc)%null_a
                endif
            endif
            ! the deterministic row on the matcher's references: the search test. The production search scored
            ! every reference at the particle-frame seed through the bank (a reference shift is one frame per
            ! rotation, so the seed cannot be fed to the shifted-reference kernel); without a bank the seed is
            ! below the shift tolerance and the unshifted kernel is the production's
            do ri_loc = 1,self%nrefs
                call deterministic_best(self%ref_full(ri_loc), iptcl_loc, ithr_loc, l_bank, .false., dist_loc, irot_loc)
                post_d(ri_loc,ithr_loc) = dist_loc
            end do
            call posterior_row_stats(self%nrefs, post_d(:,ithr_loc), temp, post_wd(:,ithr_loc), &
                &post_indsd(:,ithr_loc), imin_d, pmax, neff, k99_d, nvalid_d)
            if( nvalid_d == 0 ) return
            post_batch(i_loc)%iproj_d = self%ref_proj(imin_d)
            post_batch(i_loc)%lmin_d  = post_d(imin_d,ithr_loc)
            post_batch(i_loc)%pmax_d  = pmax
            post_batch(i_loc)%neff_d  = neff
            post_batch(i_loc)%k99_d   = k99_d
            euls_d = self%b_ptr%eulspace%get_euler(self%ref_proj(imin_d))
            post_batch(i_loc)%spread_d = spread_about(self%ref_state(imin_d), euls_d, k99_d, post_wd(:,ithr_loc), post_indsd(:,ithr_loc))
            ! sym_dists returns degrees
            call self%b_ptr%pgrpsyms%sym_dists(euls_min, euls_d, euls_sym, euldist, inpldist)
            post_batch(i_loc)%dist_d = euldist
            post_batch(i_loc)%rank_d = rank_in_row(post_d(:,ithr_loc), imin)
            if( .not. self%post_has_halves ) return
            ! the deterministic row on the particle's own real half: the references of the diagnostic model
            ! follow the matcher's, direction for direction
            do ri_loc = 1,self%nrefs
                call deterministic_best(self%ref_full(ri_loc) + self%post_ref_offset, iptcl_loc, ithr_loc, l_bank, .false., dist_loc, irot_loc)
                post_h(ri_loc,ithr_loc) = dist_loc
            end do
            call posterior_row_stats(self%nrefs, post_h(:,ithr_loc), temp, post_wh(:,ithr_loc), &
                &post_indsh(:,ithr_loc), imin_h, pmax, neff, k99_h, nvalid_h)
            if( nvalid_h == 0 ) return
            post_batch(i_loc)%iproj_h = self%ref_proj(imin_h)
            post_batch(i_loc)%lmin_h  = post_h(imin_h,ithr_loc)
            post_batch(i_loc)%pmax_h  = pmax
            post_batch(i_loc)%neff_h  = neff
            post_batch(i_loc)%k99_h   = k99_h
            euls_h = self%b_ptr%eulspace%get_euler(self%ref_proj(imin_h))
            post_batch(i_loc)%spread_h = spread_about(self%ref_state(imin_h), euls_h, k99_h, post_wh(:,ithr_loc), post_indsh(:,ithr_loc))
            call self%b_ptr%pgrpsyms%sym_dists(euls_d, euls_h, euls_sym, euldist, inpldist)
            post_batch(i_loc)%dist_h = euldist
            post_batch(i_loc)%rank_h = rank_in_row(post_h(:,ithr_loc), imin_d)
            ! the same row on the other real half: the reference test
            do ri_loc = 1,self%nrefs
                call deterministic_best(self%ref_full(ri_loc) + self%post_ref_offset, iptcl_loc, ithr_loc, l_bank, .true., dist_loc, irot_loc)
                post_x(ri_loc,ithr_loc) = dist_loc
            end do
            call posterior_row_stats(self%nrefs, post_x(:,ithr_loc), temp, post_wx(:,ithr_loc), &
                &post_indsx(:,ithr_loc), imin_x, pmax, neff, k99_x, nvalid_x)
            if( nvalid_x == 0 ) return
            post_batch(i_loc)%iproj_x = self%ref_proj(imin_x)
            post_batch(i_loc)%lmin_x  = post_x(imin_x,ithr_loc)
            post_batch(i_loc)%pmax_x  = pmax
            post_batch(i_loc)%neff_x  = neff
            post_batch(i_loc)%k99_x   = k99_x
            euls_x = self%b_ptr%eulspace%get_euler(self%ref_proj(imin_x))
            post_batch(i_loc)%spread_x = spread_about(self%ref_state(imin_x), euls_x, k99_x, post_wx(:,ithr_loc), post_indsx(:,ithr_loc))
            call self%b_ptr%pgrpsyms%sym_dists(euls_h, euls_x, euls_sym, euldist, inpldist)
            post_batch(i_loc)%dist_x = euldist
            post_batch(i_loc)%rank_x = rank_in_row(post_x(:,ithr_loc), imin_h)
            ! the self-reinforcement test at the pose the particle went into its half with: its previous
            ! direction on the current grid, its previous in-plane rotation, no residual shift (the image is
            ! shifted by its stored shift when it is prepared); the control a half turn away in-plane
            call get_particle_context(iptcl_loc, o_prev_loc, prev_state, prev_proj)
            if( prev_state >= 1 .and. prev_state <= self%p_ptr%nstates .and. &
                &prev_proj >= 1 .and. prev_proj <= self%p_ptr%nspace )then
                if( self%state_exists(prev_state) )then
                    nrots_loc = self%b_ptr%pftc%get_nrots()
                    iref_prev = (prev_state-1)*self%p_ptr%nspace + prev_proj + self%post_ref_offset
                    irot_prev = self%b_ptr%pftc%get_roind(360.-o_prev_loc%e3get())
                    if( irot_prev >= 1 .and. irot_prev <= nrots_loc )then
                        irot_ctrl = irot_prev + nrots_loc/2
                        if( irot_ctrl > nrots_loc ) irot_ctrl = irot_ctrl - nrots_loc
                        call self%b_ptr%pftc%gen_euclid_val_at(iref_prev, iptcl_loc, [0.,0.], irot_prev, dist_own)
                        call self%b_ptr%pftc%gen_euclid_val_at(iref_prev, iptcl_loc, [0.,0.], irot_prev, dist_other, cross=.true.)
                        post_batch(i_loc)%iproj_p   = prev_proj
                        post_batch(i_loc)%self_nats = temp * (dist_other - dist_own)
                        call self%b_ptr%pftc%gen_euclid_val_at(iref_prev, iptcl_loc, [0.,0.], irot_ctrl, dist_own)
                        call self%b_ptr%pftc%gen_euclid_val_at(iref_prev, iptcl_loc, [0.,0.], irot_ctrl, dist_other, cross=.true.)
                        post_batch(i_loc)%ctrl_nats = temp * (dist_other - dist_own)
                    endif
                endif
            endif
            call o_prev_loc%kill
        end subroutine posterior_diag

        !> the Euclidean loss of reference iref (the calculator's index) at its best grid rotation, at the
        !! particle-frame seed through the bank when the particle has one, unshifted otherwise; own slot or,
        !! with cross, the other half's
        subroutine deterministic_best( iref, iptcl_loc, ithr_loc, l_bank_loc, cross, dist_loc, irot_loc )
            integer, intent(in)  :: iref, iptcl_loc, ithr_loc
            logical, intent(in)  :: l_bank_loc, cross
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            if( l_bank_loc )then
                call self%b_ptr%pftc%gen_best_euclid_val_bank(iref, iptcl_loc, seed_bank(:,:,ithr_loc), dist_loc, irot_loc, cross=cross)
            else if( cross )then
                call self%b_ptr%pftc%gen_best_euclid_val_cross(iref, iptcl_loc, [0.,0.], dist_loc, irot_loc)
            else
                call self%b_ptr%pftc%gen_best_objfun_val(iref, iptcl_loc, [0.,0.], dist_loc, irot_loc)
            endif
        end subroutine deterministic_best

        !> the rank of reference iref in a row of losses (1: the lowest)
        integer function rank_in_row( row, iref ) result( rank )
            real,    intent(in) :: row(:)
            integer, intent(in) :: iref
            integer :: ri_loc
            rank = 1
            do ri_loc = 1,size(row)
                if( row(ri_loc) < row(iref) ) rank = rank + 1
            end do
        end function rank_in_row

        !> the posterior-weighted mean angular distance (degrees) of the k99 references holding 99% of
        !! the mass from the mode's direction, in the mode's state, across the symmetry operations;
        !! w ascending with its reference order inds, as posterior_row_stats leaves them
        real function spread_about( istate_mode, euls_mode, k99_loc, w, inds ) result( spread )
            integer, intent(in) :: istate_mode, k99_loc
            real,    intent(in) :: euls_mode(3), w(:)
            integer, intent(in) :: inds(:)
            real(dp) :: wsum, dsum
            real     :: euls(3), euls_sym(3), euldist, inpldist
            integer  :: j_loc, ri_loc
            spread = 0.
            wsum   = 0.d0
            dsum   = 0.d0
            do j_loc = self%nrefs, self%nrefs - k99_loc + 1, -1
                ri_loc = inds(j_loc)
                if( self%ref_state(ri_loc) /= istate_mode ) cycle
                euls = self%b_ptr%eulspace%get_euler(self%ref_proj(ri_loc))
                ! sym_dists returns degrees
                call self%b_ptr%pgrpsyms%sym_dists(euls_mode, euls, euls_sym, euldist, inpldist)
                wsum = wsum + real(w(j_loc),dp)
                dsum = dsum + real(w(j_loc),dp) * real(euldist,dp)
            end do
            if( wsum > 0.d0 ) spread = real(dsum / wsum)
        end function spread_about

        subroutine get_particle_context( iptcl_loc, o_prev_loc, prev_state, prev_proj )
            integer,   intent(in)    :: iptcl_loc
            type(ori), intent(inout) :: o_prev_loc
            integer,   intent(out)   :: prev_state, prev_proj
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            prev_state = o_prev_loc%get_state()
            prev_proj  = self%b_ptr%eulspace%find_closest_proj(o_prev_loc)
        end subroutine get_particle_context

        subroutine estimate_shift_seed( iptcl_loc, ithr_loc, prev_state, prev_proj, o_prev_loc, shift_seed, irot_seed )
            integer,   intent(in)    :: iptcl_loc, ithr_loc, prev_state, prev_proj
            type(ori), intent(inout) :: o_prev_loc
            real,      intent(out)   :: shift_seed(3)
            integer,   intent(out)   :: irot_seed !< in-plane rotation the estimate ends at (0: zero shift)
            integer :: irot, iref
            shift_seed = 0.
            irot_seed  = 0
            if( .not. l_sh_first ) return
            if( prev_state < 1 .or. prev_state > self%p_ptr%nstates ) return
            if( prev_proj  < 1 .or. prev_proj  > self%p_ptr%nspace  ) return
            if( .not. self%state_exists(prev_state) ) return
            iref = (prev_state-1)*self%p_ptr%nspace + prev_proj
            irot = self%b_ptr%pftc%get_roind(360.-o_prev_loc%e3get())
            call grad_shsrch_obj(ithr_loc)%set_indices(iref, iptcl_loc)
            shift_seed = grad_shsrch_obj(ithr_loc)%minimize(irot=irot, sh_rot=.false.)
            if( irot == 0 ) shift_seed(2:3) = 0.
            irot_seed = irot
        end subroutine estimate_shift_seed

        !> one reference's row entry: the in-plane rotation drawn from the row posterior at the particle's
        !! temperature over every grid rotation (gen_posterior_val; the uniform draw within an angular window it
        !! replaced was the search's exploration, now the tempered posterior's)
        subroutine score_ref( ri_loc, iptcl_loc, shift_xy, dist_loc, irot_loc )
            integer, intent(in)  :: ri_loc, iptcl_loc
            real,    intent(in)  :: shift_xy(2)
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            call self%b_ptr%pftc%gen_posterior_val(self%ref_full(ri_loc), iptcl_loc, shift_xy,&
                &self%b_ptr%pftc%ptcl_posterior_temp(iptcl_loc), dist_loc, irot_loc)
        end subroutine score_ref

        subroutine score_ref_bank( ri_loc, iptcl_loc, ithr_loc, dist_loc, irot_loc )
            integer, intent(in)  :: ri_loc, iptcl_loc, ithr_loc
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            call self%b_ptr%pftc%gen_posterior_val_bank(self%ref_full(ri_loc), iptcl_loc, seed_bank(:,:,ithr_loc),&
                &self%b_ptr%pftc%ptcl_posterior_temp(iptcl_loc), dist_loc, irot_loc)
        end subroutine score_ref_bank

        subroutine refine_best_refs( i_loc, iptcl_loc, ithr_loc, shift_seed )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            real,    intent(in) :: shift_seed(3)
            integer :: j_loc, ri_loc, istate_loc, iproj_loc, irot_loc, jbest
            real    :: refined_shift(3), xy_start(2), dbest
            locn(:,ithr_loc) = minnloc(dists_refs(:,ithr_loc), projs_ns) ! ascending: j is the unrefined rank
            do j_loc = 1,projs_ns
                ri_loc     = locn(j_loc,ithr_loc)
                istate_loc = self%ref_state(ri_loc)
                iproj_loc  = self%ref_proj(ri_loc)
                call grad_shsrch_obj(ithr_loc)%set_indices(self%ref_full(ri_loc), iptcl_loc)
                irot_loc = inpl_refs(ri_loc,ithr_loc)
                if( l_sh_first )then
                    ! the particle-frame seed as the shift of this reference at the rotation it was scored at, in the box
                    xy_start = shift_ptcl2ref(shift_seed(2:3), self%b_ptr%pftc%get_rot(irot_loc))
                    xy_start(1) = max(lims(1,1), min(lims(1,2), xy_start(1)))
                    xy_start(2) = max(lims(2,1), min(lims(2,2), xy_start(2)))
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, &
                        &sh_rot=.true., xy_in=xy_start)
                else
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true.)
                endif
                if( irot_loc > 0 )then
                    rows_batch(ri_loc,i_loc) = prob_dense_entry(eulprob_dist_switch(refined_shift(1),&
                        &self%p_ptr%cc_objfun), irot_loc)
                    nrefined_batch(i_loc) = nrefined_batch(i_loc) + 1
                    refined_batch(nrefined_batch(i_loc),i_loc) = prob_refined_entry(i_loc, ri_loc,&
                        &refined_shift(2), refined_shift(3), rows_batch(ri_loc,i_loc)%dist, -1.)
                endif
            enddo
            ! the bench record of whether the budget binds: the unrefined rank of the entry lowest after refinement
            jbest = 0
            dbest = huge(dbest)
            do j_loc = 1,projs_ns
                if( rows_batch(locn(j_loc,ithr_loc),i_loc)%dist < dbest )then
                    dbest = rows_batch(locn(j_loc,ithr_loc),i_loc)%dist
                    jbest = j_loc
                endif
            enddo
            if( jbest < 1 .or. minval(rows_batch(:,i_loc)%dist) < dbest )then
                rank_ithr(11,ithr_loc) = rank_ithr(11,ithr_loc) + 1_int64
            else
                if( jbest == 1 ) rank_ithr(0,ithr_loc) = rank_ithr(0,ithr_loc) + 1_int64
                j_loc = min(10, (jbest - 1) * 10 / projs_ns + 1)
                rank_ithr(j_loc,ithr_loc) = rank_ithr(j_loc,ithr_loc) + 1_int64
            endif
        end subroutine refine_best_refs

    end subroutine fill_tab_range

    !> One line per particle of the batch to the part's pose-posterior file, prob_posterior_iterNNN_partNN.txt
    !! (no iteration tag when which_iter is unset), opened at the first batch and closed by kill. The
    !! columns are the fields of prob_posterior_stat, after header lines saying whether the matcher's two
    !! reference slots hold the same sections and whether the real halves are loaded beside them (the
    !! _h, _x and self-test columns are zero otherwise); nvalid 0 marks a particle without a valid distance.
    subroutine write_posterior_diag( self, i_first, i_last, post_batch )
        class(eul_prob_tab),       intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        type(prob_posterior_stat), intent(in)    :: post_batch(i_first:i_last)
        type(string) :: fname
        integer :: i, io_stat
        if( .not. self%post_is_open )then
            fname = string('prob_posterior_')
            if( self%p_ptr%which_iter > 0 ) fname = fname//'iter'//int2str_pad(self%p_ptr%which_iter,3)//'_'
            fname = fname//'part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//'.txt'
            call fopen(self%post_unit, fname, status='replace', action='write', iostat=io_stat)
            call fileiochk('eul_prob_tab%write_posterior_diag; open '//fname%to_char(), io_stat)
            write(self%post_unit,'(A)') '# iptcl istate iproj nvalid lmin temp pmax neff k99 spread_deg '//&
                &'null_a_nats null_b_nats '//&
                &'iproj_d lmin_d pmax_d neff_d k99_d spread_d_deg dist_d_deg rank_d '//&
                &'iproj_h lmin_h pmax_h neff_h k99_h spread_h_deg dist_h_deg rank_h '//&
                &'iproj_x lmin_x pmax_x neff_x k99_x spread_x_deg dist_x_deg rank_x '//&
                &'iproj_p self_nats ctrl_nats mode_src'
            if( self%post_slots_identical )then
                write(self%post_unit,'(A)') '# slots_identical yes'
            else
                write(self%post_unit,'(A)') '# slots_identical no'
            endif
            if( self%post_has_halves )then
                write(self%post_unit,'(A)') '# halves yes'
            else
                write(self%post_unit,'(A)') '# halves no'
            endif
            self%post_is_open = .true.
            call fname%kill
        endif
        do i = i_first,i_last
            write(self%post_unit,'(I0,1X,I0,1X,I0,1X,I0,1X,ES13.6,1X,ES13.6,1X,F8.5,1X,F10.3,1X,I0,1X,F8.3,'//&
                &'1X,ES13.6,1X,ES13.6,'//&
                &'1X,I0,1X,ES13.6,1X,F8.5,1X,F10.3,1X,I0,1X,F8.3,1X,F8.3,1X,I0,'//&
                &'1X,I0,1X,ES13.6,1X,F8.5,1X,F10.3,1X,I0,1X,F8.3,1X,F8.3,1X,I0,'//&
                &'1X,I0,1X,ES13.6,1X,F8.5,1X,F10.3,1X,I0,1X,F8.3,1X,F8.3,1X,I0,'//&
                &'1X,I0,1X,ES13.6,1X,ES13.6,1X,I0)') &
                &post_batch(i)%iptcl, post_batch(i)%istate, post_batch(i)%iproj, post_batch(i)%nvalid, &
                &post_batch(i)%lmin, post_batch(i)%temp, post_batch(i)%pmax, post_batch(i)%neff, &
                &post_batch(i)%k99, post_batch(i)%spread, &
                &post_batch(i)%null_a, post_batch(i)%null_b, &
                &post_batch(i)%iproj_d, post_batch(i)%lmin_d, post_batch(i)%pmax_d, post_batch(i)%neff_d, &
                &post_batch(i)%k99_d, post_batch(i)%spread_d, post_batch(i)%dist_d, post_batch(i)%rank_d, &
                &post_batch(i)%iproj_h, post_batch(i)%lmin_h, post_batch(i)%pmax_h, post_batch(i)%neff_h, &
                &post_batch(i)%k99_h, post_batch(i)%spread_h, post_batch(i)%dist_h, post_batch(i)%rank_h, &
                &post_batch(i)%iproj_x, post_batch(i)%lmin_x, post_batch(i)%pmax_x, post_batch(i)%neff_x, &
                &post_batch(i)%k99_x, post_batch(i)%spread_x, post_batch(i)%dist_x, post_batch(i)%rank_x, &
                &post_batch(i)%iproj_p, post_batch(i)%self_nats, post_batch(i)%ctrl_nats, post_batch(i)%mode_src
        end do
    end subroutine write_posterior_diag

    !> the cost record of fill_tab_range since new, read before kill
    function get_fill_profile( self ) result( profile )
        class(eul_prob_tab), intent(in) :: self
        type(prob_fill_profile) :: profile
        profile = self%fill_profile
        if( associated(self%b_ptr) )then
            call self%b_ptr%pftc%get_shiftmap_profile(profile%map_times, profile%map_counts, profile%map_geometry)
        endif
    end function get_fill_profile

    !> the batch's rows and refined entries, kept for assign_batch; the fill's arrays move in
    subroutine keep_batch_rows( self, rows, refined, nrefined )
        class(eul_prob_tab),                   intent(inout) :: self
        type(prob_dense_entry),   allocatable, intent(inout) :: rows(:,:)
        type(prob_refined_entry), allocatable, intent(inout) :: refined(:,:)
        integer,                  allocatable, intent(inout) :: nrefined(:)
        if( allocated(self%batch_rows) ) THROW_HARD('the previous batch is not assigned yet; eul_prob_tab%keep_batch_rows')
        call move_alloc(rows,     self%batch_rows)
        call move_alloc(refined,  self%batch_refined)
        call move_alloc(nrefined, self%batch_nrefined)
    end subroutine keep_batch_rows

    !> the row of particle i over the compact references and its refined entries (nrefined of them), from the batch
    !! the fill kept
    subroutine get_row( self, i, row, refined, nrefined )
        class(eul_prob_tab),      intent(in)    :: self
        integer,                  intent(in)    :: i
        type(prob_dense_entry),   intent(inout) :: row(:)
        type(prob_refined_entry), intent(inout) :: refined(:)
        integer,                  intent(out)   :: nrefined
        integer :: j
        row(1:self%nrefs) = self%batch_rows(:,i)
        nrefined = self%batch_nrefined(i)
        do j = 1,nrefined
            refined(j) = self%batch_refined(j,i)
        enddo
    end subroutine get_row

    !> the most refined entries any row of the batch holds: the fill's budget
    integer function max_refined( self )
        class(eul_prob_tab), intent(in) :: self
        max_refined = 0
        if( allocated(self%batch_refined) ) max_refined = size(self%batch_refined,1)
    end function max_refined

    !> the entry of compact reference ri of a row, with the shift of its search if the row refined it
    function entry_candidate( self, ri, row, refined, nrefined ) result( candidate )
        class(eul_prob_tab),      intent(in) :: self
        integer,                  intent(in) :: ri, nrefined
        type(prob_dense_entry),   intent(in) :: row(:)
        type(prob_refined_entry), intent(in) :: refined(:)
        type(prob_candidate) :: candidate
        integer :: j
        candidate%iref = self%ref_full(ri)
        candidate%inpl = row(ri)%inpl
        candidate%dist = row(ri)%dist
        do j = 1,nrefined
            if( refined(j)%ri == ri )then
                candidate%x      = refined(j)%x
                candidate%y      = refined(j)%y
                candidate%sdev   = refined(j)%sdev
                candidate%has_sh = .true.
                return
            endif
        enddo
    end function entry_candidate

    !> the assignment opened (the batches' assign_batch run between it and close_assignment): the insertion table, the cone's geometry, the n_eff readout over the rows of the global
    !! sample (nglob its size, rank0 the rank before this table's first row: a part's offset), the draws' scope, the
    !! counts; by default the table's own rows are the sample
    subroutine open_assignment( self, nglob, rank0 )
        class(eul_prob_tab), intent(inout) :: self
        integer, optional,   intent(in)    :: nglob, rank0
        integer :: nglob_here, rank0_here
        if( self%l_assign_open ) THROW_HARD('the assignment is open; eul_prob_tab%open_assignment')
        nglob_here = self%nptcls
        rank0_here = 0
        if( present(nglob) ) nglob_here = nglob
        if( present(rank0) ) rank0_here = rank0
        self%kmax = max(1, self%p_ptr%topk)
        call self%insert_tab%new(self%nptcls, self%pinds, self%kmax)
        self%l_extras = self%kmax > 1 .and. self%insert_cone%is_open()
        if( self%l_extras )then
            if( self%b_ptr%eulspace%get_noris() /= self%p_ptr%nspace ) THROW_HARD('the direction space is not built; the weighted insertion needs it; open_assignment')
            if( self%seed_nrots < 1 ) THROW_HARD('the table carries no rotation grid; the weighted insertion needs it; open_assignment')
            call self%insert_cone%set_space(self%b_ptr%eulspace, self%b_ptr%pgrpsyms, self%seed_nrots)
        endif
        call self%neff%new(nglob_here, PROB_NEFF_ROWS, PROB_NEFF_NMIN, rank0=rank0_here)
        ! the draw within the state's block, keyed on the particle
        self%rnd_scope     = prob_rnd_scope(self%p_ptr%which_iter, .false.)
        self%nmode         = 0
        self%ntemp_missing = 0
        self%nbad          = 0
        self%nflat         = 0
        self%l_assign_open = .true.
    end subroutine open_assignment

    !> the batch's selection (the pipeline): the rows the fill kept are assigned and released
    subroutine assign_batch( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        if( .not. allocated(self%batch_rows) ) THROW_HARD('no rows kept for this batch; eul_prob_tab%assign_batch')
        if( lbound(self%batch_rows,2) /= i_first .or. ubound(self%batch_rows,2) /= i_last )&
            &THROW_HARD('the rows kept are not this batch''s; eul_prob_tab%assign_batch')
        call self%assign_rows(i_first, i_last)
        deallocate(self%batch_rows, self%batch_refined, self%batch_nrefined)
    end subroutine assign_batch

    !> the assignment closed: its readouts as a record, published by a distributed part (the master merges the parts'
    !! records in part order and reports) or reported here by the project owner; the invalid and flat rows are
    !! counted per row and reported here
    subroutine close_assignment( self )
        class(eul_prob_tab), intent(inout) :: self
        type(prob_readout) :: rec
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab%close_assignment')
        rec%route         = PROB_READOUT_3D_DENSE
        rec%nptcls        = self%nptcls
        rec%nrefs         = self%nrefs
        rec%nmode         = self%nmode
        rec%ntemp_missing = self%ntemp_missing
        rec%kmax          = self%kmax
        rec%l_extras      = self%l_extras
        rec%nbad  = self%nbad
        rec%nflat = self%nflat
        if( self%l_extras ) call self%insert_tab%counts(rec%ins_nptcls, rec%ins_npose, rec%ins_wc_sum, rec%ins_nelig, &
            &rec%ins_ncapped, rec%kmax)
        rec%neff_nrows    = self%neff%nkept()
        rec%neff_sum_draw = self%neff%sum_neff(self%p_ptr%prob_beta)
        rec%neff_sum_phys = self%neff%sum_neff(1.)
        call self%neff%kill
        if( self%p_ptr%l_distr_worker )then
            call write_prob_readout(rec, prob_readout_fname(self%p_ptr%part, self%p_ptr%numlen))
        else
            call report_prob_readout(rec, self%p_ptr%prob_beta)
        endif
        self%l_assign_open = .false.
    end subroutine close_assignment

    !> the selection over the rows i_first..i_last of the batch: per particle
    !! the state of its closest reference, the direction drawn in that state's block, the drawn entry committed with
    !! its refined shift (the seed shift otherwise), then the extras of the weighted insertion over the same rows
    subroutine assign_rows( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        type(prob_dense_entry),   allocatable :: row(:)
        type(prob_refined_entry), allocatable :: refined(:)
        type(prob_candidate) :: candidate
        logical, allocatable :: l_refined(:)
        real,    allocatable :: w(:)
        integer :: i, iref, iref_best, si, ifirst, ilast, nblock, which, assigned_iref, nrefined, maxref
        real    :: score, score_best, dist, temp, invalid_dist
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab%assign_rows')
        invalid_dist = 0.1 * huge(invalid_dist)
        maxref = self%max_refined()
        allocate(row(self%nrefs), refined(max(1,maxref)), w(self%nrefs), l_refined(i_first:i_last))
        l_refined = .false.
        do i = i_first,i_last
            call self%get_row(i, row, refined, nrefined)
            ! the state label: the particle's closest reference over all states
            iref_best  = 1
            score_best = valid_score(row(1)%dist)
            do iref = 2,self%nrefs
                score = valid_score(row(iref)%dist)
                if( score < score_best )then
                    score_best = score
                    iref_best  = iref
                endif
            enddo
            call count_row_validity
            si     = self%state_to_active_rank(self%ref_state(iref_best))
            ifirst = (si - 1) * self%p_ptr%nspace + 1
            ilast  = si * self%p_ptr%nspace
            nblock = ilast - ifirst + 1
            ! the n_eff readout takes the rows it wants in particle order
            if( self%neff%wanted(i) ) call self%neff%offer(nblock, row(ifirst:ilast)%dist, self%temps(i))
            temp = self%p_ptr%prob_beta * self%temps(i)
            if( .not. (temp > 0.) )then
                ! no temperature from the fill: the argmin, as the particle has no scale to draw at
                self%ntemp_missing = self%ntemp_missing + 1
                temp = huge(temp)
            endif
            call draw_posterior_row(nblock, row(ifirst:ilast)%dist, temp, w, &
                &rnd_keyed(rnd_draw_key(self%rnd_scope, self%pinds(i), RND_FAM_ASSIGN)), which, dist)
            if( which < 1 )then
                assigned_iref = iref_best ! no valid distance in the block: the balanced rule's reference
            else
                assigned_iref = ifirst - 1 + which
            endif
            if( assigned_iref == iref_best ) self%nmode = self%nmode + 1
            candidate = entry_candidate(self, assigned_iref, row, refined, nrefined)
            ! the extras of the weighted insertion are weighed against the committed entry's distance, so the
            ! committed entry must be a refined one too (an unrefined entry's distance is at the seed shift)
            if( self%l_extras ) l_refined(i) = which > 0 .and. candidate%has_sh
            call self%assign_candidate(i, candidate)
            call record_seed_shift(self%assgn_map(i), self%seed_shifts(:,i), self%seed_has_sh(i), self%p_ptr%l_doshift)
            self%assgn_map(i)%frac = 100.
        enddo
        if( self%l_extras ) call select_extras
        deallocate(row, refined, w, l_refined)

    contains

        !> the distance, or the invalid marker when it is not finite or at the marker
        real function valid_score( d ) result( score )
            real, intent(in) :: d
            score = d
            if( .not.(ieee_is_finite(score) .and. score < invalid_dist) ) score = invalid_dist
        end function valid_score

        !> the rows the readout warns of: a row without a valid distance, a flat one
        subroutine count_row_validity
            real    :: d, min_dist, max_dist
            integer :: nvalid, iref_loc
            nvalid   = 0
            min_dist = huge(min_dist)
            max_dist = -huge(max_dist)
            do iref_loc = 1,self%nrefs
                d = row(iref_loc)%dist
                if( ieee_is_finite(d) .and. d < invalid_dist )then
                    min_dist = min(min_dist, d)
                    max_dist = max(max_dist, d)
                    nvalid   = nvalid + 1
                endif
            enddo
            if( nvalid == 0 )then
                self%nbad = self%nbad + 1
            else if( max_dist - min_dist <= 0. )then
                self%nflat = self%nflat + 1
            endif
        end subroutine count_row_validity

        !> the extras of the weighted insertion: for every particle whose committed entry is refined, the refined
        !! entries of its state's block whose pose lies within the insertion cone of the committed one (the rotation
        !! modulo the point group and the shift together moving the rim by less than the map resolves), the kmax - 1 of
        !! highest posterior at the particle's own temperature over the refined poses' own losses (the row's distance
        !! under the window search, the insertion loss with the map's prior under the map: insert_pose_loss), with the
        !! weights that posterior renormalised over the kept set (top_posterior_extras): the draw decides the
        !! hypothesis, the extras only spread the particle over the poses the map resolves alike. Deterministic given
        !! the committed poses, so the particles run in parallel after the draws, each thread on its own row and
        !! eligibility mask; the table takes the result in particle order
        subroutine select_extras
            type(prob_candidate) :: cands(self%kmax)
            type(prob_dense_entry),   allocatable :: rowt(:,:)
            type(prob_refined_entry), allocatable :: refinedt(:,:)
            integer, allocatable :: kinds(:,:), nkept(:), nelig(:), nreft(:)
            real,    allocatable :: kw(:,:), pose_loss(:,:)
            logical, allocatable :: eligible(:,:)
            integer :: i_loc, si_loc, ifirst_loc, ilast_loc, istate_loc, j, k, ithr, nthr_team
            nthr_team = omp_get_max_threads() ! the per-thread rows and masks are sized by the team the region runs
            allocate(kinds(self%kmax,i_first:i_last), nkept(i_first:i_last), nelig(i_first:i_last), source=0)
            allocate(kw(self%kmax,i_first:i_last), source=0.)
            allocate(eligible(self%p_ptr%nspace,nthr_team), source=.false.)
            allocate(pose_loss(self%p_ptr%nspace,nthr_team), source=0.)
            allocate(rowt(self%nrefs,nthr_team), refinedt(max(1,maxref),nthr_team), nreft(nthr_team))
            !$omp parallel do default(shared) proc_bind(close) schedule(static) &
            !$omp private(i_loc,si_loc,ifirst_loc,ilast_loc,istate_loc,j,k,ithr)
            do i_loc = i_first,i_last
                nkept(i_loc) = 1
                kw(1,i_loc)  = 1.
                if( .not. l_refined(i_loc) ) cycle
                ithr       = omp_get_thread_num() + 1
                call self%get_row(i_loc, rowt(:,ithr), refinedt(:,ithr), nreft(ithr))
                istate_loc = self%assgn_map(i_loc)%istate
                si_loc     = self%state_to_active_rank(istate_loc)
                ifirst_loc = (si_loc - 1) * self%p_ptr%nspace + 1
                ilast_loc  = si_loc * self%p_ptr%nspace
                eligible(:,ithr)  = .false.
                pose_loss(:,ithr) = rowt(ifirst_loc:ilast_loc,ithr)%dist
                do j = 1,nreft(ithr)
                    if( refinedt(j,ithr)%ri < ifirst_loc .or. refinedt(j,ithr)%ri > ilast_loc ) cycle
                    k = refinedt(j,ithr)%ri - ifirst_loc + 1
                    pose_loss(k,ithr) = refinedt(j,ithr)%loss
                    eligible(k,ithr)  = self%insert_cone%within(istate_loc, self%assgn_map(i_loc)%iproj, &
                        &self%assgn_map(i_loc)%inpl, [self%assgn_map(i_loc)%x, self%assgn_map(i_loc)%y], &
                        &k, rowt(refinedt(j,ithr)%ri,ithr)%inpl, [refinedt(j,ithr)%x, refinedt(j,ithr)%y])
                enddo
                call top_posterior_extras(self%p_ptr%nspace, pose_loss(:,ithr), &
                    &temperature_or_argmin(self%temps(i_loc)), eligible(:,ithr), self%assgn_map(i_loc)%iproj, self%kmax, &
                    &PROB_INSERT_WMIN, nkept(i_loc), kinds(:,i_loc), kw(:,i_loc), nelig(i_loc))
            enddo
            !$omp end parallel do
            do i_loc = i_first,i_last
                si_loc     = self%state_to_active_rank(self%assgn_map(i_loc)%istate)
                ifirst_loc = (si_loc - 1) * self%p_ptr%nspace + 1
                if( nkept(i_loc) > 1 )then
                    call self%get_row(i_loc, rowt(:,1), refinedt(:,1), nreft(1))
                    do k = 2,nkept(i_loc)
                        cands(k-1) = entry_candidate(self, ifirst_loc - 1 + kinds(k,i_loc), rowt(:,1), refinedt(:,1), nreft(1))
                    enddo
                endif
                call self%insert_tab%add(i_loc, kw(1,i_loc), nkept(i_loc) - 1, cands(1:nkept(i_loc)-1), &
                    &kw(2:nkept(i_loc),i_loc), nelig(i_loc))
            enddo
            deallocate(kinds, nkept, nelig, kw, eligible, pose_loss, rowt, refinedt, nreft)
        end subroutine select_extras

    end subroutine assign_rows

    ! ptcl -> state (using assigned iproj or previous iproj) assignment
    pure integer function ref_state( self, iref ) result( state )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        integer :: active_rank
        active_rank = (iref - 1) / self%p_ptr%nspace + 1
        state = self%ssinds(active_rank)
    end function ref_state

    pure integer function ref_proj( self, iref ) result( proj )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        proj = modulo(iref - 1, self%p_ptr%nspace) + 1
    end function ref_proj

    pure integer function ref_full( self, iref ) result( full_ref )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        full_ref = (self%ref_state(iref) - 1) * self%p_ptr%nspace + self%ref_proj(iref)
    end function ref_full

    pure integer function full_to_compact_ref( self, full_ref ) result( iref )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: full_ref
        integer :: state, proj, active_rank
        iref = 0
        if( full_ref < 1 .or. full_ref > self%p_ptr%nspace * self%p_ptr%nstates ) return
        state = (full_ref - 1) / self%p_ptr%nspace + 1
        proj  = modulo(full_ref - 1, self%p_ptr%nspace) + 1
        active_rank = self%state_to_active_rank(state)
        if( active_rank < 1 ) return
        iref = (active_rank - 1) * self%p_ptr%nspace + proj
    end function full_to_compact_ref

    subroutine assign_candidate( self, iptcl_map, candidate )
        class(eul_prob_tab),   intent(inout) :: self
        integer,               intent(in)    :: iptcl_map
        type(prob_candidate),  intent(in)    :: candidate
        integer :: state, proj
        if( candidate%iref < 1 .or. candidate%iref > self%p_ptr%nspace * self%p_ptr%nstates )then
            THROW_HARD('candidate reference is out of range in assign_candidate')
        endif
        state = (candidate%iref - 1) / self%p_ptr%nspace + 1
        proj  = modulo(candidate%iref - 1, self%p_ptr%nspace) + 1
        self%assgn_map(iptcl_map)%pind   = self%pinds(iptcl_map)
        self%assgn_map(iptcl_map)%istate = state
        self%assgn_map(iptcl_map)%iproj  = proj
        self%assgn_map(iptcl_map)%inpl   = candidate%inpl
        self%assgn_map(iptcl_map)%dist   = candidate%dist
        self%assgn_map(iptcl_map)%x      = candidate%x
        self%assgn_map(iptcl_map)%y      = candidate%y
        self%assgn_map(iptcl_map)%sdev   = candidate%sdev
        self%assgn_map(iptcl_map)%has_sh = candidate%has_sh
    end subroutine assign_candidate

    ! FILE IO

    ! write a global assignment map to binary file

    !> the weighted insertion's cone about the committed pose (per state the resolution that bounds the rim
    !! displacement); unset, no extra pose is selected. The geometry (direction space, point group, rotation grid)
    !! is this table's and is set at the assignment
    subroutine set_insert_cone( self, cone )
        class(eul_prob_tab),     intent(inout) :: self
        class(prob_insert_cone), intent(in)    :: cone
        self%insert_cone = cone
    end subroutine set_insert_cone

    ! read from the global assignment map to local partition for shift search and further refinement

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        self%l_model_set        = .false.
        self%projs_ns_model     = 0
        self%sigma2_prior_model = 0.
        if( self%post_is_open ) call fclose(self%post_unit)
        self%post_unit      = -1
        self%post_is_open   = .false.
        self%post_slots_identical = .false.
        self%post_has_halves = .false.
        self%post_ref_offset = 0
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%batch_rows)     ) deallocate(self%batch_rows)
        if( allocated(self%batch_refined)  ) deallocate(self%batch_refined)
        if( allocated(self%batch_nrefined) ) deallocate(self%batch_nrefined)
        call self%neff%kill
        self%l_assign_open = .false.
        self%l_extras      = .false.
        self%kmax          = 1
        self%nmode         = 0
        self%ntemp_missing = 0
        self%nbad          = 0
        self%nflat         = 0
        call self%insert_tab%kill
        call self%insert_cone%kill
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
        if( allocated(self%temps)        ) deallocate(self%temps)
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%ssinds)       ) deallocate(self%ssinds)
        if( allocated(self%state_to_active_rank) ) deallocate(self%state_to_active_rank)
        if( allocated(self%state_exists) ) deallocate(self%state_exists)
        self%seed_nrots = 0
        self%fill_profile = prob_fill_profile()
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

end module simple_eul_prob_tab
