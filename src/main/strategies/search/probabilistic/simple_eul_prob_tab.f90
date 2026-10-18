!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab_utils, only: build_pind_lookup, calc_athres, calc_num2sample, &
    &eulprob_dist_switch, read_seed_shift_table, record_seed_shift, sample_likelihood_dist, &
    &write_seed_shift_table, prob_candidate, prob_dense_entry, prob_refined_entry, prob_fill_profile, &
    &PROB_TAB_DENSE_FORMAT
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_polarft_calc,       only: SHERRSQ
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"

! the refined entries of one part file, master indices
type :: refined_list
    type(prob_refined_entry), allocatable :: e(:)
end type refined_list

type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(prob_dense_entry),   allocatable :: rows(:,:)          !< dense table (active refs,nptcls), master
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
        allocate(self%assgn_map(self%nptcls), self%rows(self%nrefs,self%nptcls))
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
            ! A particle's cost depends on its shift searches (how many alternations they take, and
            ! with the grid keys how many of them the fit settles), so static blocks leave threads
            ! idle at the tail (load balance 0.90 against 0.98 dynamic, measured with the grid
            ! keys). Each thread owns its search object and candidate buffer and a particle is
            ! processed whole, so the dynamic order changes only which thread, and so which random
            ! draws, a particle meets.
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
        end subroutine process_particle_without_shift

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
        integer, allocatable :: greedy_state(:), score_mode(:), head_first(:), ref_next(:)
        real,    allocatable :: score_work(:,:), score_min(:), score_spread(:)
        real,    allocatable :: iref_dist(:), dists_sorted(:), state_projs_athres(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: i, iref, assigned_iref, assigned_ptcl, istate, si, active_idx, ithr, nleft, ref_offset
        integer :: alloc_stat
        real    :: projs_athres, dist_tmp, corr_tmp
        character(len=256) :: alloc_msg
        allocate(stab_inds(self%nptcls, self%nrefs), score_work(self%nptcls,nthr_glob),&
            &score_min(self%nptcls), score_spread(self%nptcls), score_mode(self%nptcls),&
            &inds_sorted(self%nrefs), iref_dist_inds(self%nrefs), greedy_state(self%nptcls),&
            &iref_dist(self%nrefs), dists_sorted(self%nrefs), state_projs_athres(self%p_ptr%nstates),&
            &ptcl_avail(self%nptcls), head_first(self%nptcls), ref_next(self%nrefs), stat=alloc_stat, errmsg=alloc_msg)
        if( alloc_stat /= 0 )then
            write(logfhandle,*) 'eul_prob_tab%ref_assign allocation failed for nptcls/nrefs: ', self%nptcls, self%nrefs
            write(logfhandle,*) trim(alloc_msg)
            THROW_HARD('failed allocating probability assignment work arrays')
        endif
        call prepare_ref_score_vectors
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
        deallocate(stab_inds, inds_sorted, iref_dist_inds, greedy_state, score_work, score_min, score_spread, score_mode,&
            &iref_dist, dists_sorted, state_projs_athres, ptcl_avail, head_first, ref_next)

    contains

        subroutine prepare_ref_score_vectors()
            real    :: min_dist, max_dist, dist_val, spread, invalid_dist
            integer :: i_loc, iref_loc, nvalid, nflat, nbad
            invalid_dist = 0.1 * huge(invalid_dist)
            nflat = 0
            nbad  = 0
            score_min    = 0.
            score_spread = 0.
            score_mode   = 0
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i_loc,iref_loc,min_dist,max_dist,dist_val,spread,nvalid) reduction(+:nflat,nbad)
            do i_loc = 1, self%nptcls
                min_dist = huge(min_dist)
                max_dist = -huge(max_dist)
                nvalid   = 0
                do iref_loc = 1,self%nrefs
                    dist_val = self%rows(iref_loc,i_loc)%dist
                    if( ieee_is_finite(dist_val) .and. dist_val < invalid_dist )then
                        min_dist = min(min_dist, dist_val)
                        max_dist = max(max_dist, dist_val)
                        nvalid   = nvalid + 1
                    endif
                enddo
                if( nvalid == 0 )then
                    nbad = nbad + 1
                    cycle
                endif
                spread = max_dist - min_dist
                score_min(i_loc)    = min_dist
                score_spread(i_loc) = spread
                if( spread <= 0. )then
                    score_mode(i_loc) = 1
                    nflat = nflat + 1
                else
                    score_mode(i_loc) = 2
                endif
            enddo
            !$omp end parallel do
            if( nbad == self%nptcls ) THROW_HARD('all probability-table reference distances are invalid')
            if( nflat == self%nptcls ) THROW_HARD('all probability-table reference distances are flat')
            if( nbad > 0 ) write(logfhandle,*) 'WARNING: particles with invalid probability-table distances: ', nbad
            if( nflat > 0 ) write(logfhandle,*) 'WARNING: particles with flat probability-table reference distances: ', nflat
        end subroutine prepare_ref_score_vectors

        real function ref_score( iref_loc, iptcl_loc ) result(score)
            integer, intent(in) :: iref_loc, iptcl_loc
            real :: dist_val, invalid_dist
            invalid_dist = 0.1 * huge(invalid_dist)
            dist_val = self%rows(iref_loc,iptcl_loc)%dist
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
            integer :: i_loc, iref_loc, iref_best
            real    :: score, score_best
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i_loc,iref_loc,iref_best,score,score_best)
            do i_loc = 1,self%nptcls
                iref_best  = 1
                score_best = ref_score(1, i_loc)
                do iref_loc = 2,self%nrefs
                    score = ref_score(iref_loc, i_loc)
                    if( score < score_best )then
                        score_best = score
                        iref_best  = iref_loc
                    endif
                enddo
                greedy_state(i_loc) = self%ref_state(iref_best)
            enddo
            !$omp end parallel do
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
        candidate%inpl = self%rows(ri,i)%inpl
        candidate%dist = self%rows(ri,i)%dist
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

    !> Opens the dense part file: header (nrefs, nptcls, format word, refined entries), particle indices, seed table,
    !! one row of (distance, in-plane index) per particle and reference, then the refined entries
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

    !> header, particle indices and seed table of the dense part file; the rows follow at rows_addr
    subroutine write_dense_head( self )
        class(eul_prob_tab), intent(inout) :: self
        integer(int64) :: file_header(4), addr
        file_header = [int(self%nrefs,int64), int(self%nptcls,int64), PROB_TAB_DENSE_FORMAT, self%table_nrefined]
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
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
        integer, allocatable :: next(:), nparts_of(:), seed_nrots_part(:)
        integer :: ipart, i, j
        allocate(part_refined(nparts), seed_nrots_part(nparts))
        allocate(nparts_of(self%nptcls), source=0)
        ! a part writes only its own particles' columns; each thread has its own unit and buffers
        !$omp parallel do default(shared) private(ipart) proc_bind(close) schedule(dynamic)
        do ipart = 1,nparts
            call self%read_dense_part(fbody//int2str_pad(ipart,numlen)//'.dat', part_refined(ipart)%e,&
                &seed_nrots_part(ipart), nparts_of)
        enddo
        !$omp end parallel do
        if( any(nparts_of /= 1) ) THROW_HARD('a particle of the dense table is in no part file or in several')
        ! a part without particles never set its rotation grid (0)
        do ipart = 1,nparts
            if( seed_nrots_part(ipart) == 0 ) cycle
            if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_part(ipart)
            if( self%seed_nrots /= seed_nrots_part(ipart) ) THROW_HARD('seed_nrots mismatch between dense part files')
        enddo
        ! the refined entries indexed by particle, in part order
        allocate(next(self%nptcls), source=0)
        do ipart = 1,nparts
            do j = 1,size(part_refined(ipart)%e)
                i = part_refined(ipart)%e(j)%i
                next(i) = next(i) + 1
            enddo
        enddo
        self%refined_offsets(1) = 1
        do i = 1,self%nptcls
            self%refined_offsets(i+1) = self%refined_offsets(i) + next(i)
        enddo
        next = self%refined_offsets(1:self%nptcls)
        deallocate(self%refined)
        allocate(self%refined(self%refined_offsets(self%nptcls+1)-1))
        do ipart = 1,nparts
            do j = 1,size(part_refined(ipart)%e)
                i = part_refined(ipart)%e(j)%i
                self%refined(next(i)) = part_refined(ipart)%e(j)
                next(i) = next(i) + 1
            enddo
        enddo
        deallocate(part_refined, next, nparts_of, seed_nrots_part)
    end subroutine read_tabs_to_glob

    !> one dense part file: rows and seed table to the master's particles, out the refined entries (master indices)
    !! and the part's rotation grid (0 without particles); nparts_of counts the parts holding each master particle
    subroutine read_dense_part( self, binfname, refined, seed_nrots_loc, nparts_of )
        class(eul_prob_tab),                   intent(inout) :: self
        class(string),                         intent(in)    :: binfname
        type(prob_refined_entry), allocatable, intent(out)   :: refined(:)
        integer,                               intent(out)   :: seed_nrots_loc
        integer,                               intent(inout) :: nparts_of(:)
        type(prob_refined_entry), allocatable :: refined_loc(:)
        type(prob_dense_entry) :: entry
        real,    allocatable :: seed_shifts_loc(:,:)
        logical, allocatable :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:), pinds_loc(:), glob_of(:)
        integer(int64) :: file_header(4), addr, rows_addr, row_bytes
        integer :: funit, io_stat, nptcls_loc, nrefined_loc, i_loc, i_glob, pind, max_pind, j, n
        if( .not. file_exists(binfname) ) THROW_HARD('probability table part files should be ready')
        call fopen(funit,binfname,access='STREAM',action='READ',status='OLD',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab; read_dense_part; file: '//binfname%to_char(), io_stat)
        read(funit,pos=1) file_header
        if( int(file_header(1)) /= self%nrefs ) THROW_HARD('nrefs mismatch in eul_prob_tab%read_dense_part')
        if( file_header(3) /= PROB_TAB_DENSE_FORMAT ) THROW_HARD('not a dense part file; eul_prob_tab%read_dense_part')
        nptcls_loc   = int(file_header(2))
        nrefined_loc = int(file_header(4))
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        read(funit,pos=addr,iostat=io_stat) pinds_loc
        call fileiochk('simple_eul_prob_tab; read_dense_part pinds; file: '//binfname%to_char(), io_stat)
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        if( nptcls_loc < 1 ) seed_nrots_loc = 0
        rows_addr = addr
        row_bytes = int(self%nrefs,int64) * int(storage_size(entry)/8,int64)
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        allocate(glob_of(nptcls_loc), source=0)
        do i_loc = 1,nptcls_loc
            pind = pinds_loc(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob < 1 ) cycle
            glob_of(i_loc) = i_glob
            !$omp atomic
            nparts_of(i_glob) = nparts_of(i_glob) + 1
            self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
            self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            read(funit,pos=rows_addr+int(i_loc-1,int64)*row_bytes,iostat=io_stat) self%rows(:,i_glob)
            call fileiochk('simple_eul_prob_tab; read_dense_part rows; file: '//binfname%to_char(), io_stat)
        enddo
        allocate(refined_loc(nrefined_loc))
        if( nrefined_loc > 0 )then
            read(funit,pos=rows_addr+int(nptcls_loc,int64)*row_bytes,iostat=io_stat) refined_loc
            call fileiochk('simple_eul_prob_tab; read_dense_part refined; file: '//binfname%to_char(), io_stat)
        endif
        call fclose(funit)
        n = 0
        do j = 1,nrefined_loc
            i_loc = refined_loc(j)%i
            if( i_loc < 1 .or. i_loc > nptcls_loc ) THROW_HARD('invalid particle of a refined entry')
            if( refined_loc(j)%ri < 1 .or. refined_loc(j)%ri > self%nrefs ) THROW_HARD('invalid reference of a refined entry')
            if( glob_of(i_loc) < 1 ) cycle
            n = n + 1
            refined_loc(n)   = refined_loc(j)
            refined_loc(n)%i = glob_of(i_loc)
        enddo
        refined = refined_loc(1:n)
        deallocate(refined_loc, pinds_loc, seed_shifts_loc, seed_has_sh_loc, pind2glob, glob_of)
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
        self%table_addr     = 1
        self%rows_addr      = 1
        self%table_nrefined = 0
        self%table_nfilled  = 0
        if( allocated(self%rows)            ) deallocate(self%rows)
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
