!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab_utils, only: build_pind_lookup, calc_athres, calc_num2sample, &
    &eulprob_dist_switch, read_seed_shift_table, record_seed_shift, sample_likelihood_dist, &
    &write_seed_shift_table, prob_candidate, prob_dense_entry, prob_refined_entry, prob_fill_profile, &
    &prob_posterior_stat, posterior_temperature, posterior_row_stats
use simple_type_defs,          only: OBJFUN_EUCLID
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_polarft_calc,       only: SHERRSQ
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"

! the format word of a dense part file, in the header slot where a sparse part file holds its candidate count
integer(int64), parameter :: PROB_TAB_DENSE_FORMAT = -3_int64

! outcomes of reading one dense part file, raised after the parallel read
integer, parameter :: PART_READ_OK = 0, PART_READ_OPEN = 1, PART_READ_IO = 2, PART_READ_NREFS = 3
integer, parameter :: PART_READ_FORMAT = 4, PART_READ_COUNTS = 5, PART_READ_PTCL = 6, PART_READ_REFINED = 7
integer, parameter :: PART_READ_SIZE = 8, PART_READ_STATES = 9
! particles per block of the master's loops over references (a block's columns stay in cache)
integer, parameter :: PTCL_BLOCK = 256
! bytes of dense rows a part read takes in one statement
integer(int64), parameter :: ROWS_READ_BYTES = 8_int64 * 1024_int64 * 1024_int64

! the refined entries of one part file, master indices
type :: refined_list
    type(prob_refined_entry), allocatable :: e(:)
end type refined_list

type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    real,                     allocatable :: dists(:,:)         !< dense table distances (nptcls,active refs), master
    integer,                  allocatable :: inpls(:,:)         !< dense table in-plane indices (nptcls,active refs), master
    type(prob_refined_entry), allocatable :: refined(:)         !< refined entries in particle order, master
    integer,                  allocatable :: refined_offsets(:) !< first refined entry of each particle (nptcls+1)
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
    real,           allocatable :: seed_shifts(:,:) !< per-particle seeded shift       (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)   !< per-particle seeded shift flag  (nptcls)
    integer                     :: seed_nrots = 0   !< rotation grid of the table (checked when the parts are read)
    integer,        allocatable :: pinds(:)        !< particle indices for processing
    integer,        allocatable :: ssinds(:)       !< non-empty state indices
    integer,        allocatable :: state_to_active_rank(:) !< requested state -> compact active-state rank
    logical,        allocatable :: state_exists(:)
    integer                     :: table_unit    = -1
    logical                     :: table_is_open = .false.
    integer(int64)              :: table_addr    = 1_int64  !< next write position of the part file
    integer(int64)              :: rows_addr     = 1_int64  !< first row of a dense part file
    integer(int64)              :: table_nrefined = 0_int64 !< refined entries written to a dense part file
    integer                     :: table_nfilled  = 0       !< particles whose rows are written
    integer                     :: post_unit      = -1      !< the pose-posterior diagnostic file (prob_post_diag=yes)
    logical                     :: post_is_open   = .false.
    logical                     :: post_slots_identical = .false. !< the two reference slots hold the same sections (header)
    logical                     :: post_has_halves = .false. !< the real halves are loaded beside the references (header)
    integer                     :: nptcls          !< size of pinds array
    integer                     :: nstates         !< states number
    integer                     :: nrefs           !< reference number
    type(prob_fill_profile)     :: fill_profile    !< accumulated cost record of fill_tab_range (bench file)
    contains
    ! CONSTRUCTOR
    procedure :: new => new_global
    procedure :: new_worker
    procedure :: new_compact_global
    procedure :: new_assignment
    procedure, private :: new_common
    procedure, private :: initialize_storage
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: fill_tab_range
    procedure :: write_tab
    procedure :: begin_write
    procedure :: write_rows
    procedure, private :: write_dense_head
    procedure, private :: write_posterior_diag
    procedure :: get_fill_profile
    procedure :: read_assignment
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: read_tabs_to_glob
    procedure, private :: read_dense_part
    procedure :: dense_candidate
    procedure :: ref_assign
    procedure :: write_assignment
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

    subroutine new_global( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        ! every row is written by the read; no refined entries until then
        allocate(self%assgn_map(self%nptcls), self%dists(self%nptcls,self%nrefs), self%inpls(self%nptcls,self%nrefs))
        allocate(self%refined(0))
        allocate(self%refined_offsets(self%nptcls+1), source=1)
        call self%initialize_storage
    end subroutine new_global

    subroutine new_worker( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
    end subroutine new_worker

    subroutine new_compact_global( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls))
        call self%initialize_storage
    end subroutine new_compact_global

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

    subroutine new_assignment( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: i, iptcl
        real    :: x
        call self%kill
        self%p_ptr  => params
        self%b_ptr  => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        do i = 1,self%nptcls
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
        enddo
    end subroutine new_assignment

    subroutine fill_tab_range( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        integer, allocatable   :: locn(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
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
        integer :: i, si, iptcl, n, projs_ns, ithr, inds_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), istate
        integer :: inpl_refs(self%nrefs,nthr_glob)
        real    :: lims(2,2), lims_init(2,2), inpl_athres(self%p_ptr%nstates)
        real    :: dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), dists_refs(self%nrefs,nthr_glob)
        logical :: l_sh_first, l_grid
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo (single state)
        integer(int64) :: seed_builds, bank_ithr
        integer :: kfromto_loc(2)
        integer(int64) :: fit_ithr(0:5)
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first )then
            THROW_HARD('invalid particle range in eul_prob_tab%fill_tab_range')
        endif
        if( .not. self%table_is_open ) THROW_HARD('eul_prob_tab%fill_tab_range; no open part file')
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        l_sh_first      = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1
        ! the particle-frame grid and its quadratic step in the multi-state table, the particle-frame seed in the
        ! seeded single-state table
        l_grid          = self%p_ptr%nstates > 1
        seed_builds     = 0_int64
        call seed_rnd
        projs_ns = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, state=istate)
            projs_ns            = max(projs_ns, n)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
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
            self%post_has_halves      = self%b_ptr%pftc%get_nrefs() >= 2 * self%nrefs
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
            if( l_sh_first )then
                kfromto_loc = self%b_ptr%pftc%get_kfromto()
                allocate(seed_bank(self%b_ptr%pftc%get_pftsz()+1, kfromto_loc(2)-kfromto_loc(1)+1, nthr_glob))
            endif
            ! fill the table
            t_region = omp_get_wtime()
            ! A particle's cost depends on its shift searches, so the schedule is dynamic. Each thread owns
            ! its search object and a particle is processed whole, writing only its own row: the order
            ! changes only which thread, and so which random draws, a particle meets.
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
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
            self%fill_profile%objective_evals = self%fill_profile%objective_evals + obj_ithr
            self%fill_profile%gradient_evals  = self%fill_profile%gradient_evals  + grad_ithr
            self%fill_profile%alternations    = self%fill_profile%alternations    + alt_ithr
            self%fill_profile%grid_points     = self%fill_profile%grid_points     + grid_ithr
            self%fill_profile%minimize_calls  = self%fill_profile%minimize_calls  + min_ithr
            call grad_shsrch_obj(ithr)%get_fit_profile(fit_ithr)
            bank_ithr = grad_shsrch_obj(ithr)%get_bank_builds()
            self%fill_profile%bank_builds = self%fill_profile%bank_builds + bank_ithr
            self%fill_profile%fit_counts      = self%fill_profile%fit_counts      + fit_ithr
        end do
        self%fill_profile%bank_builds = self%fill_profile%bank_builds + seed_builds
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        if( allocated(seed_bank) ) deallocate(seed_bank)
        call self%write_rows(i_first, rows_batch, pack_refined())
        if( l_post )then
            call self%write_posterior_diag(i_first, i_last, post_batch)
            deallocate(post_batch, post_w, post_inds, post_d, post_wd, post_indsd, post_h, post_wh, post_indsh, &
                &post_x, post_wx, post_indsx)
        endif

    contains

        !> the batch's refined entries in particle order
        function pack_refined() result( refined )
            type(prob_refined_entry), allocatable :: refined(:)
            integer :: i_loc, n
            allocate(refined(sum(nrefined_batch)))
            n = 0
            do i_loc = i_first,i_last
                refined(n+1:n+nrefined_batch(i_loc)) = refined_batch(1:nrefined_batch(i_loc),i_loc)
                n = n + nrefined_batch(i_loc)
            enddo
        end function pack_refined

        subroutine process_particle_with_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            type(ori) :: o_prev_loc
            integer   :: prev_state, prev_proj, ri_loc, irot_loc, irot_seed
            real      :: shift_seed(3), dist_loc
            real(dp)  :: t0
            logical   :: l_bank_seed
            call get_particle_context(iptcl_loc, o_prev_loc, prev_state, prev_proj)
            t0 = omp_get_wtime()
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
                    call score_ref(ri_loc, iptcl_loc, shift_seed(2:3), ithr_loc, dist_loc, irot_loc)
                endif
                rows_batch(ri_loc,i_loc) = prob_dense_entry(dist_loc, irot_loc)
                dists_refs(ri_loc,ithr_loc) = dist_loc
                inpl_refs(ri_loc,ithr_loc)  = irot_loc
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            call refine_best_refs(i_loc, iptcl_loc, ithr_loc, shift_seed)
            ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            if( l_post ) call posterior_diag(i_loc, iptcl_loc, ithr_loc, shift_seed(2:3))
            call o_prev_loc%kill
        end subroutine process_particle_with_shift

        subroutine process_particle_without_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            integer  :: ri_loc, irot_loc
            real     :: dist_loc
            real(dp) :: t0
            t0 = omp_get_wtime()
            do ri_loc = 1,self%nrefs
                call score_ref(ri_loc, iptcl_loc, [0.,0.], ithr_loc, dist_loc, irot_loc)
                rows_batch(ri_loc,i_loc) = prob_dense_entry(dist_loc, irot_loc)
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            if( l_post ) call posterior_diag(i_loc, iptcl_loc, ithr_loc, [0.,0.])
        end subroutine process_particle_without_shift

        !> The pose posterior over the particle's finished production row (likelihood-sampled in-plane,
        !! the best projs_ns shift-refined) at the temperature of its weighted power, with the angular
        !! spread of the references holding 99% of the mass about the mode (the mode's state, across the
        !! symmetry operations). Then the controls (doc/implementation_notes/planned/pose_marginalization.md
        !! sections 5 and 11): the identical-pose null test, the kernel re-scored at the production mode's
        !! exact (reference, shift, grid rotation) against the stored loss; the deterministic row on the
        !! matcher's references (best grid rotation at the seed shift, no sampling, no refinement), the
        !! search test; and, when the real halves are loaded beside the references, the same deterministic
        !! row on the particle's own real half and on the other one, the reference test, with the
        !! self-reinforcement test at the pose the particle went into its half with.
        subroutine posterior_diag( i_loc, iptcl_loc, ithr_loc, shift_xy )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            real,    intent(in) :: shift_xy(2)
            type(ori) :: o_prev_loc
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
            ! the null test: the production mode's shift as stored (the refined shift when it was refined,
            ! rotated to the frame of reference by the search; the seed shift otherwise)
            sh_mode = shift_xy
            do j_loc = 1,nrefined_batch(i_loc)
                if( refined_batch(j_loc,i_loc)%ri == imin )then
                    sh_mode = [refined_batch(j_loc,i_loc)%x, refined_batch(j_loc,i_loc)%y]
                    exit
                endif
            end do
            if( irot_min >= 1 )then
                sh_kernel = shift_ref2ptcl(sh_mode, self%b_ptr%pftc%get_rot(irot_min))
                call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, sh_kernel, irot_min, dist_at)
                post_batch(i_loc)%null_a = temp * (dist_at - rows_batch(imin,i_loc)%dist)
                call self%b_ptr%pftc%gen_euclid_val_at(self%ref_full(imin), iptcl_loc, sh_mode, irot_min, dist_at)
                post_batch(i_loc)%null_b = temp * (dist_at - rows_batch(imin,i_loc)%dist)
            endif
            ! the deterministic row on the own slot: the search test
            do ri_loc = 1,self%nrefs
                call self%b_ptr%pftc%gen_best_objfun_val(self%ref_full(ri_loc), iptcl_loc, shift_xy, dist_loc, irot_loc)
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
                call self%b_ptr%pftc%gen_best_objfun_val(self%ref_full(ri_loc) + self%nrefs, iptcl_loc, shift_xy, dist_loc, irot_loc)
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
                call self%b_ptr%pftc%gen_best_euclid_val_cross(self%ref_full(ri_loc) + self%nrefs, iptcl_loc, shift_xy, dist_loc, irot_loc)
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
                    iref_prev = (prev_state-1)*self%p_ptr%nspace + prev_proj + self%nrefs
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

        subroutine score_ref( ri_loc, iptcl_loc, shift_xy, ithr_loc, dist_loc, irot_loc )
            integer, intent(in)  :: ri_loc, iptcl_loc, ithr_loc
            real,    intent(in)  :: shift_xy(2)
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            integer :: istate_loc, iproj_loc, full_ref
            real    :: corr_loc
            istate_loc = self%ref_state(ri_loc)
            iproj_loc  = self%ref_proj(ri_loc)
            full_ref   = self%ref_full(ri_loc)
            call self%b_ptr%pftc%gen_likelihood_val(full_ref, iptcl_loc, shift_xy,&
                &inpl_likelihood_nsample(istate_loc), dist_loc, corr_loc, irot_loc,&
                &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
        end subroutine score_ref

        subroutine score_ref_bank( ri_loc, iptcl_loc, ithr_loc, dist_loc, irot_loc )
            integer, intent(in)  :: ri_loc, iptcl_loc, ithr_loc
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            real :: corr_loc
            call self%b_ptr%pftc%gen_likelihood_val_bank(self%ref_full(ri_loc), iptcl_loc, seed_bank(:,:,ithr_loc),&
                &inpl_likelihood_nsample(self%ref_state(ri_loc)), dist_loc, corr_loc, irot_loc,&
                &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
        end subroutine score_ref_bank

        integer function inpl_likelihood_nsample( istate_loc ) result(nsample)
            integer, intent(in) :: istate_loc
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, inpl_athres(istate_loc))
            nsample = min(self%b_ptr%pftc%get_nrots(), max(1, int(athres_ub * real(self%b_ptr%pftc%get_nrots()) / 180.)))
        end function inpl_likelihood_nsample

        subroutine refine_best_refs( i_loc, iptcl_loc, ithr_loc, shift_seed )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            real,    intent(in) :: shift_seed(3)
            integer :: j_loc, ri_loc, istate_loc, iproj_loc, irot_loc
            real    :: refined_shift(3), xy_start(2)
            locn(:,ithr_loc) = minnloc(dists_refs(:,ithr_loc), projs_ns)
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
                        &refined_shift(2), refined_shift(3))
                endif
            enddo
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
                &'iproj_p self_nats ctrl_nats'
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
                &'1X,I0,1X,ES13.6,1X,ES13.6)') &
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
                &post_batch(i)%iproj_p, post_batch(i)%self_nats, post_batch(i)%ctrl_nats
        end do
    end subroutine write_posterior_diag

    !> the cost record of fill_tab_range since new, read before kill
    function get_fill_profile( self ) result( profile )
        class(eul_prob_tab), intent(in) :: self
        type(prob_fill_profile) :: profile
        profile = self%fill_profile
    end function get_fill_profile

    ! ptcl -> (proj, state) assignment using calibrated likelihood weights.
    ! Every reference sorts the particles by distance; its head is its closest available particle.
    ! One state: draw a reference among the heads, assign its head, advance the heads that pointed
    ! at that particle. Several states: a particle's state is the state of its closest reference over
    ! all states (ties: the lowest reference index), then the same draw runs within each state over
    ! that state's references. A particle-indexed list of the references whose head it is
    ! (head_first, ref_next) finds the heads to advance without rescanning every reference.
    subroutine ref_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer, allocatable :: stab_inds(:,:), inds_sorted(:), iref_dist_inds(:)
        integer, allocatable :: greedy_state(:), head_first(:), ref_next(:)
        real,    allocatable :: score_work(:,:)
        real,    allocatable :: iref_dist(:), dists_sorted(:), state_projs_athres(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: i, iref, assigned_iref, assigned_ptcl, istate, si, active_idx, ithr, nleft, ref_offset
        integer :: alloc_stat, nthr_work
        real    :: projs_athres, dist_tmp, corr_tmp
        character(len=256) :: alloc_msg
        nthr_work = 1
        !$ nthr_work = omp_get_max_threads()
        allocate(stab_inds(self%nptcls, self%nrefs), score_work(self%nptcls,nthr_work),&
            &inds_sorted(self%nrefs), iref_dist_inds(self%nrefs), greedy_state(self%nptcls),&
            &iref_dist(self%nrefs), dists_sorted(self%nrefs), state_projs_athres(self%p_ptr%nstates),&
            &ptcl_avail(self%nptcls), head_first(self%nptcls), ref_next(self%nrefs), stat=alloc_stat, errmsg=alloc_msg)
        if( alloc_stat /= 0 )then
            write(logfhandle,*) 'eul_prob_tab%ref_assign allocation failed for nptcls/nrefs: ', self%nptcls, self%nrefs
            write(logfhandle,*) trim(alloc_msg)
            THROW_HARD('failed allocating probability assignment work arrays')
        endif
        call count_invalid_and_flat
        ! Sort each reference using one scratch score column per worker thread.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i,ithr)
        do iref = 1, self%nrefs
            ithr = omp_get_thread_num() + 1
            do i = 1,self%nptcls
                stab_inds(i,iref) = i
                score_work(i,ithr) = ref_score(iref, i)
            enddo
            call hpsort(score_work(:,ithr), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        projs_athres = 0.
        state_projs_athres = 0.
        do si = 1, self%nstates
            istate = self%ssinds(si)
            state_projs_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist',&
                &self%p_ptr%prob_athres, state=istate)
            projs_athres = max(projs_athres, state_projs_athres(istate))
        enddo
        greedy_state = 0
        if( self%nstates > 1 )then
            call assign_state_labels()
            do si = 1,self%nstates
                call assign_refs_for_state(si)
            enddo
        else
            call reset_ref_frontier(1, self%nrefs)
            ptcl_avail = .true.
            nleft      = self%nptcls
            call build_head_lists(1, self%nrefs)
            do while( nleft > 0 )
                call sample_likelihood_dist(self%nrefs, frontier_ref_dist, likelihood_nsample(self%nrefs, projs_athres),&
                    &dist_tmp, corr_tmp, assigned_iref, dists_sorted, inds_sorted)
                call assign_current_ref()
                call advance_heads_at(assigned_ptcl)
            enddo
        endif
        deallocate(stab_inds, inds_sorted, iref_dist_inds, greedy_state, score_work,&
            &iref_dist, dists_sorted, state_projs_athres, ptcl_avail, head_first, ref_next)

    contains

        !> the particles without a valid distance and those whose valid distances are all equal: warnings, and a stop
        !! when every particle is one of them. References outside, so that each reference's column streams.
        subroutine count_invalid_and_flat()
            real,    allocatable :: min_dist(:), max_dist(:)
            integer, allocatable :: nvalid(:)
            real    :: dist_val, invalid_dist
            integer :: ib, i_loc, iref_loc, nflat, nbad
            invalid_dist = 0.1 * huge(invalid_dist)
            allocate(min_dist(self%nptcls), source=huge(invalid_dist))
            allocate(max_dist(self%nptcls), source=-huge(invalid_dist))
            allocate(nvalid(self%nptcls), source=0)
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ib,i_loc,iref_loc,dist_val)
            do ib = 1,self%nptcls,PTCL_BLOCK
                do iref_loc = 1,self%nrefs
                    do i_loc = ib,min(self%nptcls,ib+PTCL_BLOCK-1)
                        dist_val = self%dists(i_loc,iref_loc)
                        if( ieee_is_finite(dist_val) .and. dist_val < invalid_dist )then
                            min_dist(i_loc) = min(min_dist(i_loc), dist_val)
                            max_dist(i_loc) = max(max_dist(i_loc), dist_val)
                            nvalid(i_loc)   = nvalid(i_loc) + 1
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
            nbad  = 0
            nflat = 0
            do i_loc = 1,self%nptcls
                if( nvalid(i_loc) == 0 )then
                    nbad = nbad + 1
                else if( max_dist(i_loc) - min_dist(i_loc) <= 0. )then
                    nflat = nflat + 1
                endif
            enddo
            if( nbad == self%nptcls ) THROW_HARD('all probability-table reference distances are invalid')
            if( nflat == self%nptcls ) THROW_HARD('all probability-table reference distances are flat')
            if( nbad > 0 ) write(logfhandle,*) 'WARNING: particles with invalid probability-table distances: ', nbad
            if( nflat > 0 ) write(logfhandle,*) 'WARNING: particles with flat probability-table reference distances: ', nflat
            deallocate(min_dist, max_dist, nvalid)
        end subroutine count_invalid_and_flat

        real function ref_score( iref_loc, iptcl_loc ) result(score)
            integer, intent(in) :: iref_loc, iptcl_loc
            real :: dist_val, invalid_dist
            invalid_dist = 0.1 * huge(invalid_dist)
            dist_val = self%dists(iptcl_loc,iref_loc)
            score = 1.
            if( .not.(ieee_is_finite(dist_val) .and. dist_val < invalid_dist) )then
                score = invalid_dist
                return
            endif
            score = dist_val
        end function ref_score

        integer function likelihood_nsample( n, athres_ub_in ) result(nsample)
            integer, intent(in) :: n
            real,    intent(in) :: athres_ub_in
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, athres_ub_in)
            nsample = min(n, max(1, int(athres_ub * real(n) / 180.)))
        end function likelihood_nsample

        real function frontier_ref_dist( iref_loc ) result(dist)
            integer, intent(in) :: iref_loc
            dist = iref_dist(iref_loc)
        end function frontier_ref_dist

        !> the references ifirst..ilast with their heads at their closest particle
        subroutine reset_ref_frontier( ifirst, ilast )
            integer, intent(in) :: ifirst, ilast
            integer :: iref_loc
            do iref_loc = ifirst,ilast
                iref_dist(iref_loc) = ref_score(iref_loc, stab_inds(1,iref_loc))
            enddo
            iref_dist_inds(ifirst:ilast) = 1
        end subroutine reset_ref_frontier

        !> move the head of iref_current to its closest available particle (huge distance when none is left)
        subroutine advance_ref_head( iref_current )
            integer, intent(in) :: iref_current
            integer :: iptcl_current
            do while( iref_dist_inds(iref_current) <= self%nptcls )
                iptcl_current = stab_inds(iref_dist_inds(iref_current), iref_current)
                if( ptcl_avail(iptcl_current) )then
                    iref_dist(iref_current) = ref_score(iref_current, iptcl_current)
                    return
                endif
                iref_dist_inds(iref_current) = iref_dist_inds(iref_current) + 1
            enddo
            iref_dist(iref_current) = huge(iref_dist(iref_current))
        end subroutine advance_ref_head

        !> enter reference iref_loc in the list of its head particle
        subroutine push_head( iref_loc )
            integer, intent(in) :: iref_loc
            integer :: iptcl_head
            if( iref_dist_inds(iref_loc) > self%nptcls ) return
            iptcl_head             = stab_inds(iref_dist_inds(iref_loc), iref_loc)
            ref_next(iref_loc)     = head_first(iptcl_head)
            head_first(iptcl_head) = iref_loc
        end subroutine push_head

        subroutine build_head_lists( ifirst, ilast )
            integer, intent(in) :: ifirst, ilast
            integer :: iref_loc
            head_first = 0
            do iref_loc = ifirst,ilast
                call push_head(iref_loc)
            enddo
        end subroutine build_head_lists

        !> advance the heads that pointed at the particle just assigned: no other head changes
        subroutine advance_heads_at( iptcl_assigned )
            integer, intent(in) :: iptcl_assigned
            integer :: iref_loc, iref_next
            iref_loc = head_first(iptcl_assigned)
            head_first(iptcl_assigned) = 0
            do while( iref_loc > 0 )
                iref_next = ref_next(iref_loc)
                call advance_ref_head(iref_loc)
                call push_head(iref_loc)
                iref_loc = iref_next
            enddo
        end subroutine advance_heads_at

        subroutine assign_current_ref()
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)     = .false.
            nleft = nleft - 1
            call self%assign_candidate(assigned_ptcl, self%dense_candidate(assigned_iref,assigned_ptcl))
            call record_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift)
            self%assgn_map(assigned_ptcl)%frac = 100.
        end subroutine assign_current_ref

        !> Multi-state state labelling is deterministic (argmin distance) for both weighting schemes;
        !! probabilistic exploration is confined to the within-state projection assignment
        !! (assign_refs_for_state). A particle's state is the state of its closest reference over
        !! all states; an exact tie between references goes to the lowest reference index. This is
        !! the label the global greedy frontier (the closest available reference-particle pair first,
        !! the lowest reference index among equal distances) gives every particle.
        subroutine assign_state_labels()
            real,    allocatable :: score_best(:)
            integer, allocatable :: iref_best(:)
            integer :: ib, i_loc, iref_loc
            real    :: score
            allocate(score_best(self%nptcls), iref_best(self%nptcls))
            ! references outside, ascending, with the strict < per particle
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ib,i_loc,iref_loc,score)
            do ib = 1,self%nptcls,PTCL_BLOCK
                do i_loc = ib,min(self%nptcls,ib+PTCL_BLOCK-1)
                    iref_best(i_loc)  = 1
                    score_best(i_loc) = ref_score(1, i_loc)
                enddo
                do iref_loc = 2,self%nrefs
                    do i_loc = ib,min(self%nptcls,ib+PTCL_BLOCK-1)
                        score = ref_score(iref_loc, i_loc)
                        if( score < score_best(i_loc) )then
                            score_best(i_loc) = score
                            iref_best(i_loc)  = iref_loc
                        endif
                    enddo
                enddo
                do i_loc = ib,min(self%nptcls,ib+PTCL_BLOCK-1)
                    greedy_state(i_loc) = self%ref_state(iref_best(i_loc))
                enddo
            enddo
            !$omp end parallel do
            deallocate(score_best, iref_best)
        end subroutine assign_state_labels

        !> the projection assignment of the particles labelled with the state of rank si, over that
        !! state's references, the contiguous block (si-1)*nspace+1..si*nspace; each of them keeps a
        !! head while a particle of the state is left
        subroutine assign_refs_for_state( si_loc )
            integer, intent(in) :: si_loc
            integer :: state_filter, ifirst, ilast, nblock, iref_loc
            state_filter = self%ssinds(si_loc)
            ifirst       = (si_loc - 1) * self%p_ptr%nspace + 1
            ilast        = si_loc * self%p_ptr%nspace
            nblock       = ilast - ifirst + 1
            ref_offset   = ifirst - 1
            call reset_ref_frontier(ifirst, ilast)
            ptcl_avail = greedy_state == state_filter
            nleft      = count(ptcl_avail)
            do iref_loc = ifirst,ilast
                call advance_ref_head(iref_loc)
            enddo
            call build_head_lists(ifirst, ilast)
            do while( nleft > 0 )
                if( nblock == 1 )then
                    active_idx = 1
                else
                    call sample_likelihood_dist(nblock, block_ref_dist,&
                        &likelihood_nsample(nblock, state_projs_athres(state_filter)),&
                        &dist_tmp, corr_tmp, active_idx, dists_sorted, inds_sorted)
                endif
                assigned_iref = ref_offset + active_idx
                call assign_current_ref()
                call advance_heads_at(assigned_ptcl)
            enddo
        end subroutine assign_refs_for_state

        real function block_ref_dist( active_loc ) result(dist)
            integer, intent(in) :: active_loc
            dist = iref_dist(ref_offset + active_loc)
        end function block_ref_dist

    end subroutine ref_assign

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
        self%assgn_map(iptcl_map)%has_sh = candidate%has_sh
    end subroutine assign_candidate

    !> the dense table's entry of compact reference ri and particle i, with the shift of its search if refined
    function dense_candidate( self, ri, i ) result( candidate )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: ri, i
        type(prob_candidate) :: candidate
        integer :: j
        candidate%iref = self%ref_full(ri)
        candidate%inpl = self%inpls(i,ri)
        candidate%dist = self%dists(i,ri)
        do j = self%refined_offsets(i), self%refined_offsets(i+1) - 1
            if( self%refined(j)%ri == ri )then
                candidate%x      = self%refined(j)%x
                candidate%y      = self%refined(j)%y
                candidate%has_sh = .true.
                return
            endif
        enddo
    end function dense_candidate

    ! FILE IO

    !> Opens the dense part file: header (nrefs, nptcls, format word, refined entries), active states (count, list), particle indices,
    !! seed table, one row of (distance, in-plane index) per particle and reference, then the refined entries
    subroutine begin_write( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(prob_dense_entry) :: entry
        integer :: io_stat
        if( self%table_is_open ) THROW_HARD('probability table part file is already open')
        call fopen(self%table_unit,binfname,access='STREAM',action='WRITE',status='REPLACE',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab; begin_write; file: '//binfname%to_char(),io_stat)
        self%table_is_open  = .true.
        self%table_nrefined = 0
        self%table_nfilled  = 0
        call self%write_dense_head
        self%table_addr = self%rows_addr + int(self%nptcls,int64) * int(self%nrefs,int64) * int(storage_size(entry)/8,int64)
    end subroutine begin_write

    !> header, active states, particle indices and seed table of the dense part file; the rows follow at rows_addr
    subroutine write_dense_head( self )
        class(eul_prob_tab), intent(inout) :: self
        integer(int64) :: file_header(4), addr
        file_header = [int(self%nrefs,int64), int(self%nptcls,int64), PROB_TAB_DENSE_FORMAT, self%table_nrefined]
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%nstates, self%ssinds
        addr = addr + sizeof(self%nstates) + sizeof(self%ssinds)
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        self%rows_addr = addr
    end subroutine write_dense_head

    !> the rows of the particles i_first.. (consecutive) and their refined entries, to the open dense part file
    subroutine write_rows( self, i_first, rows, refined )
        class(eul_prob_tab),      intent(inout) :: self
        integer,                  intent(in)    :: i_first
        type(prob_dense_entry),   intent(in)    :: rows(:,:)
        type(prob_refined_entry), intent(in)    :: refined(:)
        integer(int64) :: addr
        if( .not. self%table_is_open ) THROW_HARD('eul_prob_tab%write_rows; no open part file')
        if( size(rows,1) /= self%nrefs .or. i_first < 1 .or. i_first + size(rows,2) - 1 > self%nptcls )then
            THROW_HARD('eul_prob_tab%write_rows; rows out of range')
        endif
        addr = self%rows_addr + int(i_first-1,int64) * int(self%nrefs,int64) * int(storage_size(rows)/8,int64)
        write(self%table_unit,pos=addr) rows
        if( size(refined) > 0 )then
            write(self%table_unit,pos=self%table_addr) refined
            self%table_addr = self%table_addr + int(size(refined),int64) * int(storage_size(refined)/8,int64)
        endif
        self%table_nrefined = self%table_nrefined + int(size(refined),int64)
        self%table_nfilled  = self%table_nfilled  + size(rows,2)
    end subroutine write_rows

    !> completes the dense part file: the header with the count of refined entries, and the seed table
    subroutine write_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        if( .not. self%table_is_open ) call self%begin_write(binfname)
        if( self%table_nfilled /= self%nptcls ) THROW_HARD('eul_prob_tab%write_tab; rows of the part are missing')
        call self%write_dense_head
        call fclose(self%table_unit)
        self%table_unit    = -1
        self%table_is_open = .false.
        self%table_addr    = 1
    end subroutine write_tab

    !> Reads the nparts dense part files into the master table, in parallel: a part's rows go to its particles'
    !! columns, its refined entries to the particles' lists, in part order. Every particle of the master must be in
    !! exactly one part, and the parts must agree on the rotation grid.
    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: fbody
        integer,             intent(in)    :: nparts, numlen
        type(refined_list), allocatable :: part_refined(:)
        type(string), allocatable :: part_fnames(:)
        integer, allocatable :: next(:), nparts_of(:), seed_nrots_part(:), part_status(:)
        character(len=:), allocatable :: err_msg
        integer(int64) :: nrefined_tot
        integer :: ipart, i, j
        if( .not. allocated(self%dists) ) THROW_HARD('dense master table is not allocated; read_tabs_to_glob')
        if( size(self%dists,1) /= self%nptcls .or. size(self%dists,2) /= self%nrefs )then
            THROW_HARD('dense master row dimensions disagree; read_tabs_to_glob')
        endif
        allocate(part_refined(nparts), seed_nrots_part(nparts))
        allocate(nparts_of(self%nptcls), next(self%nptcls), part_status(nparts), source=0)
        allocate(part_fnames(nparts))
        do ipart = 1,nparts
            part_fnames(ipart) = fbody//int2str_pad(ipart,numlen)//'.dat'
            if( part_fnames(ipart)%strlen() > XLONGSTRLEN ) THROW_HARD('probability table part file name too long')
            if( .not. file_exists(part_fnames(ipart)) )then
                THROW_HARD('missing probability table part file: '//part_fnames(ipart)%to_char())
            endif
        enddo
        !$omp parallel do default(shared) private(ipart) proc_bind(close) schedule(dynamic)
        do ipart = 1,nparts
            call self%read_dense_part(part_fnames(ipart), part_refined(ipart)%e,&
                &seed_nrots_part(ipart), nparts_of, next, part_status(ipart))
        enddo
        !$omp end parallel do
        do ipart = 1,nparts
            if( part_status(ipart) == PART_READ_OK ) cycle
            err_msg = 'dense probability table part file '//part_fnames(ipart)%to_char()//': '//&
                &trim(part_read_error(part_status(ipart)))
            THROW_HARD(err_msg)
        enddo
        if( any(nparts_of /= 1) ) THROW_HARD('a particle of the dense table is in no part file or in several')
        ! a part without particles never set its rotation grid (0)
        do ipart = 1,nparts
            if( seed_nrots_part(ipart) == 0 ) cycle
            if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_part(ipart)
            if( self%seed_nrots /= seed_nrots_part(ipart) ) THROW_HARD('seed_nrots mismatch between dense part files')
        enddo
        ! the refined entries indexed by particle (next holds the read's counts per particle): offsets by a prefix sum,
        ! then every part places its entries, in file order; a particle's entries all come from its one part
        nrefined_tot = 0_int64
        do ipart = 1,nparts
            nrefined_tot = nrefined_tot + int(size(part_refined(ipart)%e),int64)
        enddo
        if( nrefined_tot >= int(huge(1),int64) ) THROW_HARD('refined entries exceed default-integer bounds; read_tabs_to_glob')
        self%refined_offsets(1) = 1
        do i = 1,self%nptcls
            self%refined_offsets(i+1) = self%refined_offsets(i) + next(i)
        enddo
        next = self%refined_offsets(1:self%nptcls)
        deallocate(self%refined)
        allocate(self%refined(self%refined_offsets(self%nptcls+1)-1))
        !$omp parallel do default(shared) private(ipart,j,i) proc_bind(close) schedule(dynamic)
        do ipart = 1,nparts
            do j = 1,size(part_refined(ipart)%e)
                i = part_refined(ipart)%e(j)%i
                self%refined(next(i)) = part_refined(ipart)%e(j)
                next(i) = next(i) + 1
            enddo
        enddo
        !$omp end parallel do
        deallocate(part_refined, part_fnames, next, nparts_of, seed_nrots_part, part_status)

    contains

        function part_read_error( status ) result( reason )
            integer, intent(in) :: status
            character(len=64)   :: reason
            select case(status)
                case(PART_READ_OPEN);    reason = 'cannot be opened'
                case(PART_READ_IO);      reason = 'read error'
                case(PART_READ_SIZE);    reason = 'file size differs from the size its header implies'
                case(PART_READ_STATES);  reason = 'active states differ from the master'
                case(PART_READ_NREFS);   reason = 'reference count differs from the master'
                case(PART_READ_FORMAT);  reason = 'not a dense part file'
                case(PART_READ_COUNTS);  reason = 'invalid particle or refined-entry count'
                case(PART_READ_PTCL);    reason = 'a particle is invalid or in several parts'
                case(PART_READ_REFINED); reason = 'invalid refined entry'
                case DEFAULT;            reason = 'unknown error'
            end select
        end function part_read_error

    end subroutine read_tabs_to_glob

    !> one dense part file: rows and seed table to the master's particles, out the refined entries (master indices),
    !! the part's rotation grid (0 without particles) and a PART_READ_* status; nparts_of counts the parts holding each
    !! master particle, nrefined_of its refined entries. Runs inside the parallel read: no stop, no string temporaries.
    subroutine read_dense_part( self, binfname, refined, seed_nrots_loc, nparts_of, nrefined_of, status )
        class(eul_prob_tab),                   intent(inout) :: self
        class(string),                         intent(in)    :: binfname
        type(prob_refined_entry), allocatable, intent(out)   :: refined(:)
        integer,                               intent(out)   :: seed_nrots_loc
        integer,                               intent(inout) :: nparts_of(:), nrefined_of(:)
        integer,                               intent(out)   :: status
        type(prob_refined_entry), allocatable :: refined_loc(:)
        type(prob_dense_entry),   allocatable :: rows_loc(:,:)
        type(prob_refined_entry) :: refined_entry
        type(prob_dense_entry)   :: entry
        real,    allocatable :: seed_shifts_loc(:,:)
        logical, allocatable :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:), pinds_loc(:), glob_of(:), ssinds_loc(:)
        integer(int64) :: file_header(4), addr, rows_addr, row_bytes, file_bytes, ptcl_bytes
        logical :: flag
        integer :: funit, io_stat, nptcls_loc, nrefined_loc, nstates_loc, i_loc, i_glob, pind, max_pind, j, n
        integer :: previous_parts, nchunk, first
        seed_nrots_loc = 0
        allocate(refined(0))
        call fopen(funit,binfname,access='STREAM',form='UNFORMATTED',action='READ',status='OLD',iostat=io_stat)
        if( io_stat /= 0 )then
            status = PART_READ_OPEN
            return
        endif
        status = read_part()
        call fclose(funit)

    contains

        integer function read_part() result( status_loc )
            read(funit,pos=1,iostat=io_stat) file_header
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            if( int(file_header(1)) /= self%nrefs )then
                status_loc = PART_READ_NREFS
                return
            endif
            if( file_header(3) /= PROB_TAB_DENSE_FORMAT )then
                status_loc = PART_READ_FORMAT
                return
            endif
            if( file_header(2) < 0_int64 .or. file_header(2) > int(huge(nptcls_loc),int64) .or.&
                &file_header(4) < 0_int64 .or. file_header(4) > int(huge(nrefined_loc),int64) )then
                status_loc = PART_READ_COUNTS
                return
            endif
            nptcls_loc   = int(file_header(2))
            nrefined_loc = int(file_header(4))
            ! a row's reference index counts the active references: the part's states must be the master's
            addr = sizeof(file_header) + 1
            read(funit,pos=addr,iostat=io_stat) nstates_loc
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            if( nstates_loc /= self%nstates )then
                status_loc = PART_READ_STATES
                return
            endif
            ! every read below lies inside the size the header implies
            row_bytes  = int(self%nrefs,int64) * int(storage_size(entry)/8,int64)
            ptcl_bytes = int(storage_size(pind)/8 + 2*storage_size(1.)/8 + storage_size(flag)/8,int64) + row_bytes
            inquire(unit=funit, size=file_bytes)
            if( file_bytes /= sizeof(file_header) + int((self%nstates + 2)*storage_size(seed_nrots_loc)/8,int64) +&
                &int(nptcls_loc,int64) * ptcl_bytes + int(nrefined_loc,int64) * int(storage_size(refined_entry)/8,int64) )then
                status_loc = PART_READ_SIZE
                return
            endif
            allocate(ssinds_loc(nstates_loc))
            read(funit,pos=addr+sizeof(nstates_loc),iostat=io_stat) ssinds_loc
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            if( any(ssinds_loc /= self%ssinds) )then
                status_loc = PART_READ_STATES
                return
            endif
            addr = addr + sizeof(nstates_loc) + sizeof(ssinds_loc)
            allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
            read(funit,pos=addr,iostat=io_stat) pinds_loc
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            addr = addr + sizeof(pinds_loc)
            call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
            if( nptcls_loc < 1 ) seed_nrots_loc = 0
            rows_addr = addr
            call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
            allocate(glob_of(nptcls_loc), source=0)
            do i_loc = 1,nptcls_loc
                pind = pinds_loc(i_loc)
                if( pind < 1 .or. pind > max_pind ) cycle
                i_glob = pind2glob(pind)
                if( i_glob < 1 ) cycle
                if( i_glob > self%nptcls )then
                    status_loc = PART_READ_PTCL
                    return
                endif
                glob_of(i_loc) = i_glob
                !$omp atomic capture
                previous_parts = nparts_of(i_glob)
                nparts_of(i_glob) = nparts_of(i_glob) + 1
                !$omp end atomic
                if( previous_parts /= 0 )then
                    status_loc = PART_READ_PTCL
                    return
                endif
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            enddo
            ! consecutive particles' rows in one statement, scattered to the master's columns
            nchunk = int(max(1_int64, min(int(max(1,nptcls_loc),int64), ROWS_READ_BYTES / row_bytes)))
            allocate(rows_loc(self%nrefs,nchunk))
            do first = 1,nptcls_loc,nchunk
                n = min(nchunk, nptcls_loc - first + 1)
                read(funit,pos=rows_addr+int(first-1,int64)*row_bytes,iostat=io_stat) rows_loc(:,1:n)
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
                do j = 1,n
                    i_glob = glob_of(first+j-1)
                    if( i_glob < 1 ) cycle
                    self%dists(i_glob,:) = rows_loc(:,j)%dist
                    self%inpls(i_glob,:) = rows_loc(:,j)%inpl
                enddo
            enddo
            allocate(refined_loc(nrefined_loc))
            if( nrefined_loc > 0 )then
                read(funit,pos=rows_addr+int(nptcls_loc,int64)*row_bytes,iostat=io_stat) refined_loc
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
            endif
            n = 0
            do j = 1,nrefined_loc
                i_loc = refined_loc(j)%i
                if( i_loc < 1 .or. i_loc > nptcls_loc .or. refined_loc(j)%ri < 1 .or. refined_loc(j)%ri > self%nrefs )then
                    status_loc = PART_READ_REFINED
                    return
                endif
                if( glob_of(i_loc) < 1 ) cycle
                n = n + 1
                refined_loc(n)   = refined_loc(j)
                refined_loc(n)%i = glob_of(i_loc)
                !$omp atomic update
                nrefined_of(glob_of(i_loc)) = nrefined_of(glob_of(i_loc)) + 1
            enddo
            refined    = refined_loc(1:n)
            status_loc = PART_READ_OK
        end function read_part

    end subroutine read_dense_part

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: assgn_glob(:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit,pos=headsz + 1) assgn_glob
        call fclose(funit)
        pinds_glob = assgn_glob(:)%pind
        call build_pind_lookup(pinds_glob, self%pinds, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(assgn_glob, pind2glob, pinds_glob)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, self%nptcls
            pind = self%pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 ) self%assgn_map(i_loc) = assgn_glob(i_glob)
        end do
        !$omp end parallel do
        deallocate(assgn_glob, pind2glob, pinds_glob)
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        if( self%table_is_open ) call fclose(self%table_unit)
        self%table_unit     = -1
        self%table_is_open  = .false.
        if( self%post_is_open ) call fclose(self%post_unit)
        self%post_unit      = -1
        self%post_is_open   = .false.
        self%post_slots_identical = .false.
        self%post_has_halves = .false.
        self%table_addr     = 1
        self%rows_addr      = 1
        self%table_nrefined = 0
        self%table_nfilled  = 0
        if( allocated(self%dists)           ) deallocate(self%dists)
        if( allocated(self%inpls)           ) deallocate(self%inpls)
        if( allocated(self%refined)         ) deallocate(self%refined)
        if( allocated(self%refined_offsets) ) deallocate(self%refined_offsets)
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
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
