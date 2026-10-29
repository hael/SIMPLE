!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is geometrically sparse and searched independently per state.
module simple_eul_prob_tab_neigh
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_utils, only: build_pind_lookup, eulprob_dist_switch, eulprob_corr_switch, &
    &read_seed_shift_table, write_seed_shift_table, record_seed_shift, prob_candidate, &
    &prob_candidate_buffer, prob_candidate_store, prob_dense_entry, prob_refined_entry, posterior_temperature, &
    &draw_posterior_row, shift_prior_sigma2, prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN, draw_posterior_extras
use simple_decay_funs,         only: extremal_decay
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_defs,               only: PROB_INSERT_WMIN
use simple_polarft_calc,       only: SHERRSQ
use simple_ori,                only: ori
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

! Every direct stochastic search contributes this many active fine-reference
! candidates before its stochastic stopping rule may take effect.  When fewer
! active references exist, all of them are searched.
integer, parameter :: MIN_STOCH_FINE_REFS = 5

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    type(prob_candidate_store) :: candidate_store
    type(prob_candidate_buffer), allocatable :: candidate_buffers(:) !< per thread, flushed to the part file
    integer,        allocatable :: candidate_fill_counts(:)
    integer,        allocatable :: cand_counts(:)                    !< candidates written per particle (part file)
    integer(int64)              :: table_nnz     = 0_int64
    integer                     :: table_nchunks = 0
    logical                     :: l_direct_stoch_neigh = .false.
contains
    procedure :: new_worker     => new_neigh_worker
    procedure :: new_neigh
    procedure :: new_neigh_global
    procedure :: fill_tab_range => fill_tab_neigh_range
    procedure :: begin_write
    procedure :: write_rows    => write_rows_neigh
    procedure :: write_tab
    procedure :: flush_candidate_buffers
    procedure :: ref_assign    => ref_assign_neigh
    procedure :: kill          => kill_neigh
    procedure :: read_tabs_to_glob
    procedure, private :: read_sparse_tab_to_glob
end type eul_prob_tab_neigh

! Workspace for per-thread arrays for both the stochastic and subspace paths
type :: eval_ws
    integer, allocatable :: evaluated_ref_ids(:,:)     ! [nrefs,nthr]
    integer, allocatable :: best_evals(:,:)            ! [max_refine,nthr]
    integer, allocatable :: fullref_to_sparse_ref(:)   ! [nstates*nspace]
    real,    allocatable :: evaluated_ref_dists(:,:)   ! [nrefs,nthr]
    real,    allocatable :: state_eval_dists(:,:)      ! [nrefs,nthr]
    ! Subspace-path only:
    integer, allocatable :: sub_ref_counts(:,:)        ! [nsubs,nstates]
    integer, allocatable :: sub_ref_offsets(:,:)       ! [nsubs,nstates]
    integer, allocatable :: sub_ref_list(:)            ! [nrefs]
contains
    procedure, private :: alloc_eval_stoch_ws
    procedure, private :: alloc_eval_subspace_ws
    procedure :: init_eval_stoch_ws
    procedure :: init_eval_subspace_ws
    procedure :: dealloc_eval_ws
end type eval_ws

! Workspace for coarse subspaces used by the geom/state path
type :: coarse_search_ws
    real,    allocatable :: peak_subspace_dists(:,:,:)  ! [npeak,nstates,nthr]
    integer, allocatable :: peak_subspace_inds(:,:,:)   ! [npeak,nstates,nthr]
    integer, allocatable :: peak_subspace_count(:,:)    ! [nstates,nthr]
    integer, allocatable :: pooled_sub_inds(:,:,:)      ! [npooled,nstates,nthr]
    integer, allocatable :: pooled_sub_count(:,:)       ! [nstates,nthr]
contains
    procedure :: alloc_coarse_ws
    procedure :: dealloc_coarse_ws
end type coarse_search_ws

contains

    !> a part's table with its per-thread candidate buffers (an empty part needs no more)
    subroutine new_neigh_worker( self, params, build, pinds )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%kill
        call self%eul_prob_tab%new_worker(params,build,pinds)
        allocate(self%candidate_buffers(nthr_glob))
    end subroutine new_neigh_worker

    subroutine new_neigh( self, params, build, pinds )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_worker(params,build,pinds)
        call validate_neigh_mode(self)
    end subroutine new_neigh

    subroutine new_neigh_global( self, params, build, pinds )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%eul_prob_tab%new_compact_global(params,build,pinds)
        call validate_neigh_mode(self)
    end subroutine new_neigh_global

    subroutine validate_neigh_mode( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: nsubs
        select case(trim(self%p_ptr%prob_neigh_mode))
            case('shc','snhc')
                self%l_direct_stoch_neigh = .true.
            case('geom','state')
                self%l_direct_stoch_neigh = .false.
                if( .not. allocated(self%b_ptr%subspace_inds) )then
                    THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; missing subspace indices')
                endif
                nsubs = size(self%b_ptr%subspace_inds)
                if( nsubs < 1 )then
                    THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; empty subspace indices')
                endif
            case DEFAULT
                THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; unsupported prob_neigh_mode')
        end select
    end subroutine validate_neigh_mode

    ! Dispatches to the stochastic (shc/snhc) or subspace (geom/state) path
    subroutine fill_tab_neigh_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        ! Search paths
        select case(trim(self%p_ptr%prob_neigh_mode))
            case('shc','snhc')
                call fill_tab_stoch_range(self, i_first, i_last)
            case('geom','state')
                if( .not. allocated(self%b_ptr%subspace_inds) )&
                    &THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_subspace_range; missing subspace indices')
                call fill_tab_subspace_range(self, i_first, i_last)
            case DEFAULT
                THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh_range; unsupported prob_neigh_mode')
        end select
        if( self%table_is_open ) call self%flush_candidate_buffers
    end subroutine fill_tab_neigh_range

    ! STOCHASTIC NEIGHBORHOOD BRANCH

    subroutine fill_tab_stoch_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad)      :: grad_shsrch_obj(nthr_glob)
        type(eval_ws)               :: eval_work
        type(ran_tabu), allocatable :: direct_rts(:)
        integer,        allocatable :: direct_srch_order(:,:)
        real,           allocatable :: corrs_inpl(:,:)
        complex(sp),    allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        integer :: i, ithr, max_refs_to_refine, nfull_refs
        real    :: sigma2_prior  !< the shift prior's variance per axis (shift_prior=yes)
        integer :: iptcl, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_shc_neigh, l_snhc_neigh, l_seed_sh_first, l_bank(nthr_glob)
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob)
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_shc_neigh     = (trim(self%p_ptr%prob_neigh_mode) == 'shc')
        l_snhc_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'snhc')
        l_seed_sh_first = self%p_ptr%l_doshift .and. l_shc_neigh
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        l_bank    = .false.
        sigma2_prior = 0.
        call seed_rnd
        call alloc_seed_bank(self, seed_bank)
        nfull_refs = self%p_ptr%nstates * self%p_ptr%nspace
        allocate(corrs_inpl(nrots,nthr_glob))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_stoch_ws(self, max_refs_to_refine)
        allocate(direct_srch_order(nfull_refs,nthr_glob), direct_rts(nthr_glob))
        do ithr = 1,nthr_glob
            direct_rts(ithr) = ran_tabu(nfull_refs)
        enddo
        t_region = omp_get_wtime()
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new_alternating(self%b_ptr, lims, lims_init=lims_init, &
                    &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
            end do
            ! the shift prior's width, from the sampled particles' last shift increments (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
            ! dynamic: particles differ in cost (stopping rule, shift searches); each is processed
            ! whole by one thread with its own search object, so only its random draws change
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed)&
            !$omp proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_stoch(i, iptcl, ithr, .true., shift_seed)
                self%seed_shifts(:,i) = shift_seed(2:3)
                self%seed_has_sh(i)   = l_seed_sh_first
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed)&
            !$omp proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                shift_seed = 0.
                call process_particle_stoch(i, iptcl, ithr, .false., shift_seed)
            enddo
            !$omp end parallel do
        endif
        call record_neigh_fill_profile(self, grad_shsrch_obj, omp_get_wtime() - t_region, ts_seed, ts_score, ts_refine,&
            &nscored, nbank, i_to - i_from + 1)
        ! cleanup
        call eval_work%dealloc_eval_ws
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        deallocate(direct_srch_order,direct_rts)
        deallocate(corrs_inpl, seed_bank)

    contains

        ! One particle evaluation in the stochastic path
        subroutine process_particle_stoch(i_loc, iptcl_loc, ithr_loc, l_with_shift, shift_seed)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: shift_seed(3)
            type(ori) :: o_prev
            integer   :: prev_state_loc, prev_proj_loc, neval_loc
            real(dp)  :: t0
            call get_particle_context(self, iptcl_loc, o_prev, prev_state_loc, prev_proj_loc)
            t0 = omp_get_wtime()
            if( self%p_ptr%l_shift_prior .and. l_with_shift ) call grad_shsrch_obj(ithr_loc)%set_shift_prior(sigma2_prior,&
                &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz()))
            call estimate_shift_seed(self, grad_shsrch_obj, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc, o_prev,&
                &l_seed_sh_first, shift_seed)
            call memoize_seed_bank(self, iptcl_loc, shift_seed, seed_bank(:,:,ithr_loc), l_bank(ithr_loc), nbank(ithr_loc))
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            call evaluate_direct_stochastic_refs(i_loc, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc,&
                &shift_seed, l_with_shift, neval_loc)
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            nscored(ithr_loc)  = nscored(ithr_loc) + int(neval_loc, int64)
            if( l_with_shift )then
                t0 = omp_get_wtime()
                call refine_best_neighbors(self, eval_work, grad_shsrch_obj, l_seed_sh_first, max_refs_to_refine,&
                &i_loc, ithr_loc, iptcl_loc, shift_seed, neval_loc)
                ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            endif
            call o_prev%kill
        end subroutine process_particle_stoch

        ! Scores references drawn in random order
        subroutine evaluate_direct_stochastic_refs(i_loc, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc,&
            &shift_seed, l_with_shift, neval_loc)
            integer, intent(in)  :: i_loc, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc
            real,    intent(in)  :: shift_seed(3)
            logical, intent(in)  :: l_with_shift
            integer, intent(out) :: neval_loc
            integer :: full_ref_loc, prev_full_ref_loc, isample, ri_loc, irot_loc
            integer :: nrefs_bound, nmin_refs
            real    :: dist_loc, corr_loc, prev_corr_loc, neigh_frac
            logical :: l_greedy_first, l_prev_scored
            neval_loc     = 0
            l_prev_scored = .false.
            call direct_rts(ithr_loc)%ne_ran_iarr(direct_srch_order(:,ithr_loc))
            prev_full_ref_loc = 0
            if( prev_state_loc >= 1 .and. prev_state_loc <= self%p_ptr%nstates .and.&
                &prev_proj_loc >= 1 .and. prev_proj_loc <= self%p_ptr%nspace )then
                prev_full_ref_loc = (prev_state_loc - 1) * self%p_ptr%nspace + prev_proj_loc
                call put_last(prev_full_ref_loc, direct_srch_order(:,ithr_loc))
            endif
            ! Only force a dense first pass when the particle truly has no
            ! search result. Keying this off updatecnt makes iteration 2 behave
            ! like a full scan for first-time sampled particles.
            l_greedy_first = l_shc_neigh .and. (.not. self%b_ptr%spproj_field%has_been_searched(iptcl_loc))
            nrefs_bound    = nfull_refs
            if( self%nrefs < 1 ) THROW_HARD('direct stochastic neighborhood search has no active references')
            nmin_refs      = min(MIN_STOCH_FINE_REFS,self%nrefs)
            if( l_snhc_neigh )then
                neigh_frac  = extremal_decay(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
                nrefs_bound = max(2, min(nfull_refs, nint(real(nfull_refs) * (1. - neigh_frac))))
            endif
            prev_corr_loc = -huge(1.0)
            if( (l_shc_neigh .and. (.not. l_greedy_first)) )then
                call calc_previous_corr(ithr_loc, iptcl_loc, prev_full_ref_loc, prev_corr_loc)
            endif
            do isample = 1,nfull_refs
                ! SNHC may leave its stochastic neighbourhood only after the
                ! particle has the minimum number of active fine candidates.
                if( l_snhc_neigh .and. isample > nrefs_bound .and. neval_loc >= nmin_refs ) exit
                full_ref_loc = direct_srch_order(isample,ithr_loc)
                if( full_ref_loc < 1 .or. full_ref_loc > size(eval_work%fullref_to_sparse_ref) ) cycle
                ri_loc = eval_work%fullref_to_sparse_ref(full_ref_loc)
                if( ri_loc < 1 ) cycle
                call score_direct_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed, l_with_shift, dist_loc, corr_loc, irot_loc)
                call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                neval_loc = neval_loc + 1
                eval_work%evaluated_ref_ids(neval_loc,ithr_loc)   = ri_loc
                eval_work%evaluated_ref_dists(neval_loc,ithr_loc) = dist_loc
                if( full_ref_loc == prev_full_ref_loc ) l_prev_scored = .true.
                ! SHC acceptance is likewise deferred until the candidate
                ! floor has been established for this sampled particle.
                if( l_shc_neigh .and. (.not. l_greedy_first) .and. corr_loc > prev_corr_loc .and. &
                    &neval_loc >= nmin_refs ) exit
            enddo
            if( neval_loc < nmin_refs )then
                THROW_HARD('direct stochastic neighborhood search did not meet the minimum active-reference count')
            endif
            ! the current pose is a candidate: the posterior draw over the store is then an acceptance step over
            ! the proposals and the incumbent, not a forced move (the stop rule above leaves a set of references
            ! worse than the current pose plus the one that beat it; 2026-10-24)
            if( prev_full_ref_loc > 0 .and. .not. l_prev_scored )then
                ri_loc = eval_work%fullref_to_sparse_ref(prev_full_ref_loc)
                if( ri_loc > 0 )then
                    call score_direct_ref(prev_full_ref_loc, ithr_loc, iptcl_loc, shift_seed, l_with_shift,&
                        &dist_loc, corr_loc, irot_loc)
                    call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                    neval_loc = neval_loc + 1
                    eval_work%evaluated_ref_ids(neval_loc,ithr_loc)   = ri_loc
                    eval_work%evaluated_ref_dists(neval_loc,ithr_loc) = dist_loc
                endif
            endif
        end subroutine evaluate_direct_stochastic_refs

        ! Computes the maximum correlation of the previous orientation without shift
        subroutine calc_previous_corr(ithr_loc, iptcl_loc, prev_full_ref_loc, prev_corr_loc)
            integer, intent(in)  :: ithr_loc, iptcl_loc, prev_full_ref_loc
            real,    intent(out) :: prev_corr_loc
            prev_corr_loc = -huge(1.0)
            if( prev_full_ref_loc < 1 .or. prev_full_ref_loc > size(eval_work%fullref_to_sparse_ref) ) return
            if( eval_work%fullref_to_sparse_ref(prev_full_ref_loc) < 1 ) return
            call self%b_ptr%pftc%gen_objfun_vals(prev_full_ref_loc, iptcl_loc, [0.,0.], corrs_inpl(:,ithr_loc))
            prev_corr_loc = max(0., maxval(corrs_inpl(:,ithr_loc)))
        end subroutine calc_previous_corr

        ! Scores a single reference at the seed shift
        subroutine score_direct_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, dist_loc, corr_loc, irot_loc)
            integer, intent(in)  :: full_ref_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            real,    intent(out) :: dist_loc, corr_loc
            integer, intent(out) :: irot_loc
            real    :: sh_loc(2)
            sh_loc = 0.
            if( l_with_shift .and. l_seed_sh_first ) sh_loc = shift_seed_loc(2:3)
            call score_seeded_ref(self, full_ref_loc, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &dist_loc, corr_loc, irot_loc)
        end subroutine score_direct_ref

    end subroutine fill_tab_stoch_range

    ! NEIGHBORHOOD SUBSPACES BRANCH

    ! Identifies peak subspaces via coarse scoring, builds a pooled neighborhoo
    ! and evaluates all fine refs & refines the best
    subroutine fill_tab_subspace_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(coarse_search_ws) :: coarse_ws
        type(eval_ws)          :: eval_work
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        integer :: i, ithr, max_refs_to_refine, nsubs, npeak_target, npooled_capacity
        real    :: sigma2_prior  !< the shift prior's variance per axis (shift_prior=yes)
        integer :: iptcl, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_geom_neigh, l_state_neigh, l_seed_sh_first, l_bank(nthr_glob)
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob)
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_geom_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'geom')
        l_state_neigh   = (trim(self%p_ptr%prob_neigh_mode) == 'state')
        ! Keep shift-first seeding single-state only. In multi-state runs this
        ! can bias all states toward the same previous-state shift basin.
        l_seed_sh_first = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        l_bank    = .false.
        sigma2_prior = 0.
        call seed_rnd
        call alloc_seed_bank(self, seed_bank)
        nsubs        = size(self%b_ptr%subspace_inds)
        npeak_target = min(max(1, self%p_ptr%npeaks), nsubs)
        npooled_capacity = min(nsubs, max(1, npeak_target * max(1, self%nstates)))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_subspace_ws(self, nsubs, max_refs_to_refine)
        call coarse_ws%alloc_coarse_ws(npeak_target, npooled_capacity, self%p_ptr%nstates, nthr_glob)
        t_region = omp_get_wtime()
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new_alternating(self%b_ptr, lims, lims_init=lims_init, &
                    &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
            end do
            ! the shift prior's width, from the sampled particles' last shift increments (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
            ! dynamic: particles differ in cost (stopping rule, shift searches); each is processed
            ! whole by one thread with its own search object, so only its random draws change
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_subspace(i, iptcl, ithr, .true., shift_seed)
                self%seed_shifts(:,i) = shift_seed(2:3)
                self%seed_has_sh(i)   = l_seed_sh_first
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl      = self%pinds(i)
                ithr       = omp_get_thread_num() + 1
                shift_seed = 0.
                call process_particle_subspace(i, iptcl, ithr, .false., shift_seed)
            enddo
            !$omp end parallel do
        endif
        call record_neigh_fill_profile(self, grad_shsrch_obj, omp_get_wtime() - t_region, ts_seed, ts_score, ts_refine,&
            &nscored, nbank, i_to - i_from + 1)
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call coarse_ws%dealloc_coarse_ws
        call eval_work%dealloc_eval_ws
        deallocate(seed_bank)

    contains

        ! One particle evaluation
        subroutine process_particle_subspace(i_loc, iptcl_loc, ithr_loc, l_with_shift, shift_seed_loc)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: shift_seed_loc(3)
            type(ori) :: o_prev
            integer   :: prev_state, prev_proj, neval
            real(dp)  :: t0
            call get_particle_context(self, iptcl_loc, o_prev, prev_state, prev_proj)
            t0 = omp_get_wtime()
            if( self%p_ptr%l_shift_prior .and. l_with_shift ) call grad_shsrch_obj(ithr_loc)%set_shift_prior(sigma2_prior,&
                &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz()))
            call estimate_shift_seed(self, grad_shsrch_obj, ithr_loc, iptcl_loc, prev_state, prev_proj, o_prev,&
                &l_seed_sh_first, shift_seed_loc)
            call memoize_seed_bank(self, iptcl_loc, shift_seed_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &nbank(ithr_loc))
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            if( l_geom_neigh )then
                call build_geometric_neighborhood(ithr_loc, prev_proj)
            else
                if( l_state_neigh ) call find_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
                call build_pooled_neighborhood(ithr_loc, prev_proj)
            endif
            call evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval)
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            nscored(ithr_loc)  = nscored(ithr_loc) + int(neval, int64)
            if( l_with_shift )then
                t0 = omp_get_wtime()
                call refine_best_neighbors(self, eval_work, grad_shsrch_obj, l_seed_sh_first, max_refs_to_refine,&
                    &i_loc, ithr_loc, iptcl_loc, shift_seed_loc, neval)
                ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            endif
            call o_prev%kill
        end subroutine process_particle_subspace

        ! Scores all subspace representatives per active state independently and records the
        ! top-npeak_target subspace indices; prob_neigh_mode=state
        subroutine find_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            integer :: si_loc, istate_loc, isub_loc, full_ref_subspace, irot_loc, ri_loc, coarse_proj
            real    :: dist
            coarse_ws%peak_subspace_dists(:,:,ithr_loc) = huge(1.0)
            coarse_ws%peak_subspace_inds(:,:,ithr_loc)  = 0
            coarse_ws%peak_subspace_count(:,ithr_loc)   = 0
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                if( .not. self%state_exists(istate_loc) ) cycle
                do isub_loc = 1, nsubs
                    coarse_proj = self%b_ptr%subspace_inds(isub_loc)
                    full_ref_subspace = (istate_loc-1)*self%p_ptr%nspace + coarse_proj
                    call score_subspace_ref(full_ref_subspace, ithr_loc, iptcl_loc, shift_seed_loc,&
                        &l_with_shift, dist, irot_loc)
                    nscored(ithr_loc) = nscored(ithr_loc) + 1_int64
                    call consider_peak_subspace(isub_loc, istate_loc, ithr_loc, dist)
                    ri_loc = eval_work%fullref_to_sparse_ref(full_ref_subspace)
                    if( ri_loc > 0 )&
                        &call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist, irot_loc, 0., 0., .false.)
                enddo
            enddo
        end subroutine find_peak_subspaces

        ! Selects the subspace that contains the particle's previous best projection
        ! No coarse path; prob_neigh_mode=geom
        subroutine build_geometric_neighborhood(ithr_loc, prev_proj_loc)
            integer, intent(in) :: ithr_loc, prev_proj_loc
            integer :: si_loc, istate_loc, coarse_proj_loc, isub_loc
            coarse_proj_loc = max(1, min(self%p_ptr%nspace, prev_proj_loc))
            isub_loc = self%b_ptr%subspace_full2sub_map(coarse_proj_loc)
            if( isub_loc < 1 .or. isub_loc > nsubs ) isub_loc = 1
            coarse_ws%pooled_sub_count(:,ithr_loc) = 0
            do si_loc = 1,self%nstates
                istate_loc = self%ssinds(si_loc)
                coarse_ws%pooled_sub_inds(1,istate_loc,ithr_loc) = isub_loc
                coarse_ws%pooled_sub_count(istate_loc,ithr_loc)  = 1
            enddo
        end subroutine build_geometric_neighborhood

        ! Evaluates a coarse subspace representative; prob_neigh_mode=state
        subroutine score_subspace_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed_loc,&
            &l_with_shift, dist_loc, irot_loc)
            integer, intent(in)  :: full_ref_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            real    :: corr_loc, sh_loc(2)
            sh_loc = 0.
            if( l_with_shift ) sh_loc = shift_seed_loc(2:3)
            call score_seeded_ref(self, full_ref_loc, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &dist_loc, corr_loc, irot_loc)
        end subroutine score_subspace_ref

        ! Maintains a sorted top-npeak_target list of subspace indices per state; prob_neigh_mode=state
        subroutine consider_peak_subspace(isub_loc, istate_loc, ithr_loc, dist)
            integer, intent(in) :: isub_loc, istate_loc, ithr_loc
            real,    intent(in) :: dist
            integer :: nfound_loc, insert_loc, k_loc
            nfound_loc = coarse_ws%peak_subspace_count(istate_loc,ithr_loc)
            if( nfound_loc >= npeak_target )then
                if( dist >= coarse_ws%peak_subspace_dists(npeak_target,istate_loc,ithr_loc) ) return
                insert_loc = npeak_target
            else
                nfound_loc = nfound_loc + 1
                coarse_ws%peak_subspace_count(istate_loc,ithr_loc) = nfound_loc
                insert_loc = nfound_loc
            endif
            do k_loc = 1, nfound_loc
                if( dist < coarse_ws%peak_subspace_dists(k_loc,istate_loc,ithr_loc) )then
                    insert_loc = k_loc
                    exit
                endif
            enddo
            do k_loc = nfound_loc, insert_loc + 1, -1
                coarse_ws%peak_subspace_dists(k_loc,istate_loc,ithr_loc) =&
                    &coarse_ws%peak_subspace_dists(k_loc-1,istate_loc,ithr_loc)
                coarse_ws%peak_subspace_inds(k_loc,istate_loc,ithr_loc) =&
                    &coarse_ws%peak_subspace_inds(k_loc-1,istate_loc,ithr_loc)
            enddo
            coarse_ws%peak_subspace_dists(insert_loc,istate_loc,ithr_loc) = dist
            coarse_ws%peak_subspace_inds(insert_loc,istate_loc,ithr_loc)  = isub_loc
        end subroutine consider_peak_subspace

        ! Merges per-state peak subspace lists into one shared pooled set, always including
        ! the subspace of the previous best orientation. The same pooled directions are
        ! evaluated for every active state so state assignment compares like with like.
        subroutine build_pooled_neighborhood(ithr_loc, prev_proj_loc)
            integer, intent(in) :: ithr_loc, prev_proj_loc
            integer :: si_loc, istate_loc, src_state, npeak_found, ncopy
            integer :: prev_isub_loc, k_loc, pooled_count, candidate_isub
            logical :: state_has_prev
            coarse_ws%pooled_sub_count(:,ithr_loc) = 0
            pooled_count = 0
            prev_isub_loc = self%b_ptr%subspace_full2sub_map(max(1, min(self%p_ptr%nspace, prev_proj_loc)))
            if( prev_isub_loc < 1 .or. prev_isub_loc > nsubs ) prev_isub_loc = 1
            do si_loc = 1, self%nstates
                src_state = self%ssinds(si_loc)
                npeak_found = coarse_ws%peak_subspace_count(src_state,ithr_loc)
                state_has_prev = .false.
                do k_loc = 1, npeak_found
                    if( coarse_ws%peak_subspace_inds(k_loc,src_state,ithr_loc) == prev_isub_loc )then
                        state_has_prev = .true.
                        exit
                    endif
                enddo
                ncopy = npeak_found
                if( (.not. state_has_prev) .and. npeak_found >= npeak_target ) ncopy = max(0, npeak_target - 1)
                do k_loc = 1, ncopy
                    candidate_isub = coarse_ws%peak_subspace_inds(k_loc,src_state,ithr_loc)
                    call add_pooled_subspace(candidate_isub, ithr_loc, pooled_count)
                enddo
                if( .not. state_has_prev ) call add_pooled_subspace(prev_isub_loc, ithr_loc, pooled_count)
            enddo
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                if( pooled_count > 0 )then
                    coarse_ws%pooled_sub_inds(1:pooled_count,istate_loc,ithr_loc) =&
                        &coarse_ws%pooled_sub_inds(1:pooled_count,1,ithr_loc)
                    coarse_ws%pooled_sub_count(istate_loc,ithr_loc) = pooled_count
                endif
            enddo
        end subroutine build_pooled_neighborhood

        subroutine add_pooled_subspace(isub_loc, ithr_loc, pooled_count)
            integer, intent(in)    :: isub_loc, ithr_loc
            integer, intent(inout) :: pooled_count
            integer :: k_loc
            if( isub_loc < 1 .or. isub_loc > nsubs ) return
            do k_loc = 1, pooled_count
                if( coarse_ws%pooled_sub_inds(k_loc,1,ithr_loc) == isub_loc ) return
            enddo
            if( pooled_count >= size(coarse_ws%pooled_sub_inds,1) ) return
            pooled_count = pooled_count + 1
            coarse_ws%pooled_sub_inds(pooled_count,1,ithr_loc) = isub_loc
        end subroutine add_pooled_subspace

        ! Pass over fine references of the pooled neighborhoods
        subroutine evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval)
            integer, intent(in)  :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            integer, intent(out) :: neval
            integer :: jsub, kref_loc, nrefs_sub_loc, offset_loc
            integer :: si_loc, ri_loc, istate_loc, iproj_loc, irot_loc, iref_loc, isub_loc
            real    :: dist, corr_loc, sh_loc(2)
            neval  = 0
            sh_loc = 0.
            if( l_with_shift ) sh_loc = shift_seed_loc(2:3)
            do si_loc = 1,self%nstates
                istate_loc = self%ssinds(si_loc)
                do jsub = 1,coarse_ws%pooled_sub_count(istate_loc,ithr_loc)
                    isub_loc = coarse_ws%pooled_sub_inds(jsub, istate_loc, ithr_loc)
                    if( isub_loc < 1 .or. isub_loc > nsubs ) cycle
                    offset_loc    = eval_work%sub_ref_offsets(isub_loc,istate_loc)
                    nrefs_sub_loc = eval_work%sub_ref_counts(isub_loc,istate_loc)
                    do kref_loc = 1,nrefs_sub_loc
                        ri_loc = eval_work%sub_ref_list(offset_loc + kref_loc - 1)
                        if( ri_loc < 1 ) cycle
                        iproj_loc = self%ref_proj(ri_loc)
                        iref_loc  = self%ref_full(ri_loc)
                        call score_seeded_ref(self, iref_loc, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc),&
                            &l_bank(ithr_loc), dist, corr_loc, irot_loc)
                        call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist, irot_loc, 0., 0., .false.)
                        neval = neval + 1
                        eval_work%evaluated_ref_ids(neval,ithr_loc)   = ri_loc
                        eval_work%evaluated_ref_dists(neval,ithr_loc) = dist
                    enddo
                enddo
            enddo
        end subroutine evaluate_neighborhood

    end subroutine fill_tab_subspace_range

    ! ASSIGNMENT SUBROUTINES

    ! Sparse graph traversal on-the-fly for reference assignment
    !> ptcl -> (proj, state) assignment over the candidate store: a particle without a valid candidate is re-seeded
    !! at its previous orientation, then every particle draws one of its own candidates from its tempered posterior.
    !! Particles are assigned to references; no reference is forced to be populated (the reference-balanced frontier
    !! was deleted 2026-10-24: if every image shows one direction, every image is assigned to it)
    subroutine ref_assign_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: i, candidate_pos
        logical :: l_seed_fallback_shift
        l_seed_fallback_shift = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1 .and. &
            &trim(self%p_ptr%prob_neigh_mode) /= 'snhc'
        do i = 1, self%nptcls
            call seed_fallback_if_empty(i)
        enddo
        call assign_particles_posterior()

    contains

        !> ptcl -> (proj, state) by one draw per particle from its tempered posterior over its evaluated candidates
        !! The particle's state is the
        !! state of its closest valid candidate over all states (an active reference with an in-plane index; ties
        !! to the first in store order), as the dense table labels it; within that state p_c proportional to
        !! exp(-T_i (L_c - L_min)) over the valid candidates, T_i the particle's temperature from the fill (temps)
        !! scaled by prob_beta; the mean n_eff of a subsample of the rows at the draw's factor and at beta 1 is
        !! logged (prob_neff_readout). The drawn candidate is committed with its shift, the seed shift recorded, frac the
        !! evaluated fraction; no reference is forced to be populated. A particle without a valid candidate (a
        !! state emptied mid-run) keeps its orientation; one without a temperature (none expected) takes its
        !! argmin with a warning. The draws run in particle order on the master's random stream.
        subroutine assign_particles_posterior()
            type(prob_neff_readout) :: neff
            type(prob_candidate), allocatable :: cands(:)
            real,    allocatable :: w(:), dists_c(:), kw(:)
            integer, allocatable :: pos_c(:), kinds(:)
            logical, allocatable :: eligible(:)
            integer :: i_loc, n_loc, which, ibest, maxcand, nmode, ntemp_missing, nkept, istate_best, kmax, k, nk, nk_kept, r_loc
            real    :: temp, dist
            logical :: l_extras
            maxcand = 0
            do i_loc = 1,self%nptcls
                maxcand = max(maxcand, int(self%candidate_store%offsets(i_loc+1) - self%candidate_store%offsets(i_loc)))
            enddo
            allocate(w(max(1,maxcand)), dists_c(max(1,maxcand)), pos_c(max(1,maxcand)))
            ! the weighted insertion's extras (topk > 1): from the refined candidates of the state's row, in the
            ! modes whose candidates are shift-refined (state, geom); the shc/snhc neighbourhoods of the early
            ! stages insert the committed pose alone, their candidates being the hill climb's unrefined proposals
            kmax = max(1, self%p_ptr%topk)
            select case(trim(self%p_ptr%prob_neigh_mode))
                case('shc','snhc')
                    l_extras = .false.
                case DEFAULT
                    l_extras = kmax > 1
            end select
            call self%insert_tab%new(self%nptcls, self%pinds, merge(kmax, 1, l_extras))
            if( l_extras ) allocate(kinds(kmax), kw(kmax), cands(kmax), eligible(max(1,maxcand)))
            ! the n_eff readout: the mean effective number of candidates over a subsample of the rows the draws are from
            call neff%new(self%nptcls, PROB_NEFF_ROWS, PROB_NEFF_NMIN)
            do i_loc = 1,self%nptcls
                if( .not. neff%wanted(i_loc) ) cycle
                call gather_state_row(i_loc, dists_c, pos_c, n_loc, ibest, istate_best)
                if( istate_best == 0 ) cycle
                call neff%offer(n_loc, dists_c(1:n_loc), self%temps(i_loc))
            enddo
            call neff%report(self%p_ptr%prob_beta)
            call neff%kill
            nmode         = 0
            ntemp_missing = 0
            nkept         = 0
            do i_loc = 1,self%nptcls
                call gather_state_row(i_loc, dists_c, pos_c, n_loc, ibest, istate_best)
                if( istate_best == 0 )then
                    nkept = nkept + 1
                    cycle
                endif
                temp = self%p_ptr%prob_beta * self%temps(i_loc)
                if( .not. (temp > 0.) )then
                    ntemp_missing = ntemp_missing + 1
                    temp = huge(temp)
                endif
                call draw_posterior_row(n_loc, dists_c(1:n_loc), temp, w, which, dist)
                if( which < 1 ) which = ibest
                if( which == ibest ) nmode = nmode + 1
                call self%assign_candidate(i_loc, self%candidate_store%candidates(pos_c(which)))
                call record_seed_shift(self%assgn_map(i_loc), self%seed_shifts(:,i_loc), self%seed_has_sh(i_loc),&
                    &self%p_ptr%l_doshift)
                self%assgn_map(i_loc)%frac = search_frac(i_loc)
                ! the extras for the weighted insertion, without replacement from the refined candidates of the row
                if( l_extras )then
                    do r_loc = 1,n_loc
                        eligible(r_loc) = self%candidate_store%candidates(pos_c(r_loc))%has_sh
                    enddo
                    call draw_posterior_extras(n_loc, dists_c(1:n_loc), temp, eligible(1:n_loc), which, kmax, &
                        &PROB_INSERT_WMIN, w, nk_kept, kinds, kw)
                    nk = nk_kept - 1
                    do k = 1,nk
                        cands(k) = self%candidate_store%candidates(pos_c(kinds(k+1)))
                    enddo
                    call self%insert_tab%add(i_loc, kw(1), nk, cands(1:max(1,nk)), kw(2:max(2,nk+1)))
                endif
            enddo
            if( l_extras ) call self%insert_tab%report
            write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> PROB_TAB_ASSIGN: neighbourhood posterior draw at the particle temperature, ',&
                &self%nptcls - nkept, ' particles x ', self%nrefs, ' references; fraction drawn at the row mode ',&
                &real(nmode) / real(max(1, self%nptcls - nkept))
            if( nkept > 0 ) write(logfhandle,'(A,I0)') &
                &'WARNING: particles without a valid neighbourhood candidate, orientation kept: ', nkept
            if( ntemp_missing > 0 ) write(logfhandle,'(A,I0)') &
                &'WARNING: particles without a temperature in the neighbourhood part files, assigned to their argmin: ', ntemp_missing
            call flush(logfhandle)
            deallocate(w, dists_c, pos_c)
            if( allocated(kinds) ) deallocate(kinds, kw, cands, eligible)
        end subroutine assign_particles_posterior

        !> the particle's row: its state (the closest valid candidate over all states; 0 without one) and,
        !! into dists_row/pos_row, its valid candidates of that state with the position of the best
        subroutine gather_state_row( i_row, dists_row, pos_row, n_row, ibest_row, istate_row )
            integer, intent(in)    :: i_row
            real,    intent(inout) :: dists_row(:)
            integer, intent(inout) :: pos_row(:)
            integer, intent(out)   :: n_row, ibest_row, istate_row
            integer :: pos_loc, ri_loc
            real    :: best_loc
            istate_row = 0
            best_loc   = huge(1.0)
            do pos_loc = int(self%candidate_store%offsets(i_row)), int(self%candidate_store%offsets(i_row+1))-1
                ri_loc = valid_compact_ref(pos_loc)
                if( ri_loc < 1 ) cycle
                if( self%candidate_store%candidates(pos_loc)%dist < best_loc )then
                    best_loc   = self%candidate_store%candidates(pos_loc)%dist
                    istate_row = self%ref_state(ri_loc)
                endif
            enddo
            n_row     = 0
            ibest_row = 0
            if( istate_row == 0 ) return
            best_loc = huge(1.0)
            do pos_loc = int(self%candidate_store%offsets(i_row)), int(self%candidate_store%offsets(i_row+1))-1
                ri_loc = valid_compact_ref(pos_loc)
                if( ri_loc < 1 ) cycle
                if( self%ref_state(ri_loc) /= istate_row ) cycle
                n_row            = n_row + 1
                dists_row(n_row) = self%candidate_store%candidates(pos_loc)%dist
                pos_row(n_row)   = pos_loc
                if( dists_row(n_row) < best_loc )then
                    best_loc  = dists_row(n_row)
                    ibest_row = n_row
                endif
            enddo
        end subroutine gather_state_row

        !> the compact reference of a store entry when it is a valid candidate (an active reference with an
        !! in-plane index), 0 otherwise
        integer function valid_compact_ref( pos_loc ) result(ri_loc)
            integer, intent(in) :: pos_loc
            ri_loc = self%full_to_compact_ref(self%candidate_store%candidates(pos_loc)%iref)
            if( ri_loc < 1 .or. ri_loc > self%nrefs )then
                ri_loc = 0
                return
            endif
            if( self%candidate_store%candidates(pos_loc)%inpl <= 0 ) ri_loc = 0
        end function valid_compact_ref

        !> the particle holds at least one valid candidate
        logical function has_valid_candidate( iptcl_loc ) result( l_valid )
            integer, intent(in) :: iptcl_loc
            integer :: pos_loc
            l_valid = .false.
            do pos_loc = int(self%candidate_store%offsets(iptcl_loc)), int(self%candidate_store%offsets(iptcl_loc+1))-1
                if( valid_compact_ref(pos_loc) > 0 )then
                    l_valid = .true.
                    return
                endif
            enddo
        end function has_valid_candidate

        real function search_frac(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            if( self%l_direct_stoch_neigh .and. self%nrefs > 0 )then
                search_frac = 100.0 * real(particle_evaluated_count(iptcl_loc)) / real(self%nrefs)
            else
                search_frac = 100.0
            endif
        end function search_frac

        integer function particle_evaluated_count(iptcl_loc) result(nevaluated)
            integer, intent(in) :: iptcl_loc
            integer :: pos_loc
            nevaluated = 0
            do pos_loc = int(self%candidate_store%offsets(iptcl_loc)),&
                &int(self%candidate_store%offsets(iptcl_loc+1))-1
                if( self%candidate_store%candidates(pos_loc)%inpl > 0 ) nevaluated = nevaluated + 1
            enddo
        end function particle_evaluated_count

        integer function find_candidate_position(iptcl_loc, iref_loc) result(pos_found)
            integer, intent(in) :: iptcl_loc, iref_loc
            integer :: pos_loc, first_loc, last_loc
            pos_found = 0
            first_loc = int(self%candidate_store%offsets(iptcl_loc))
            last_loc  = int(self%candidate_store%offsets(iptcl_loc+1)) - 1
            do pos_loc = first_loc,last_loc
                if( self%full_to_compact_ref(self%candidate_store%candidates(pos_loc)%iref) == iref_loc )then
                    pos_found = pos_loc
                    return
                endif
            enddo
        end function find_candidate_position

        !> a particle without a valid candidate is re-seeded at its previous orientation (its state's reference
        !! at its closest direction), scored at the seed shift with the posterior in-plane draw of the fill
        subroutine seed_fallback_if_empty(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            type(ori) :: o_prev
            integer :: istate_loc, iproj_loc, irot_loc, fallback_state, fallback_proj, fallback_ref_full
            integer :: ri_loc, fallback_ref_loc
            real    :: sh_seed(2), dist
            if( has_valid_candidate(iptcl_loc) ) return
            call self%b_ptr%spproj_field%get_ori(self%pinds(iptcl_loc), o_prev)
            istate_loc = o_prev%get_state()
            if( istate_loc < 1 .or. istate_loc > self%p_ptr%nstates ) istate_loc = 1
            iproj_loc = self%b_ptr%eulspace%find_closest_proj(o_prev)
            iproj_loc = max(1, min(self%p_ptr%nspace, iproj_loc))
            fallback_ref_loc = 0
            do ri_loc = 1, self%nrefs
                if( self%ref_state(ri_loc) == istate_loc .and. self%ref_proj(ri_loc) == iproj_loc )then
                    fallback_ref_loc = ri_loc
                    exit
                endif
            enddo
            if( fallback_ref_loc == 0 )then
                do ri_loc = 1, self%nrefs
                    if( self%ref_state(ri_loc) == istate_loc )then
                        fallback_ref_loc = ri_loc
                        exit
                    endif
                enddo
            endif
            if( fallback_ref_loc == 0 )then
                ! previous state is no longer active (e.g. a state emptied mid-run);
                ! leave this particle unseeded so its orientation is kept unchanged
                call o_prev%kill
                return
            endif
            candidate_pos = find_candidate_position(iptcl_loc, fallback_ref_loc)
            if( candidate_pos < 1 ) candidate_pos = int(self%candidate_store%offsets(iptcl_loc))
            if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 )then
                fallback_state    = self%ref_state(fallback_ref_loc)
                fallback_proj     = self%ref_proj(fallback_ref_loc)
                fallback_ref_full = (fallback_state-1)*self%p_ptr%nspace + fallback_proj
                sh_seed = 0.
                if( l_seed_fallback_shift ) sh_seed = o_prev%get_2Dshift()
                call self%b_ptr%pftc%gen_posterior_val(fallback_ref_full, self%pinds(iptcl_loc), sh_seed,&
                    &self%b_ptr%pftc%ptcl_posterior_temp(self%pinds(iptcl_loc)), dist, irot_loc)
                if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                self%candidate_store%candidates(candidate_pos)%iref   = self%ref_full(fallback_ref_loc)
                self%candidate_store%candidates(candidate_pos)%dist   = dist
                self%candidate_store%candidates(candidate_pos)%inpl   = irot_loc
                self%candidate_store%candidates(candidate_pos)%x      = sh_seed(1)
                self%candidate_store%candidates(candidate_pos)%y      = sh_seed(2)
                self%candidate_store%candidates(candidate_pos)%has_sh = self%p_ptr%l_doshift
            endif
            call o_prev%kill
        end subroutine seed_fallback_if_empty

    end subroutine ref_assign_neigh

    ! Table & Assignment I/O

    ! The sparse part file: header (nrefs, nptcls, candidates, chunks), particle indices, candidate counts per
    ! particle (they size the master's store without a scan of the stream), seed table, then the candidate chunks
    subroutine begin_write( self, binfname )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: binfname
        integer        :: io_stat
        integer(int64) :: file_header(4), addr
        if( self%table_is_open ) THROW_HARD('probability candidate stream is already open')
        if( .not. allocated(self%candidate_buffers) ) THROW_HARD('candidate stream requested without worker buffers')
        file_header = [int(self%nrefs,int64),int(self%nptcls,int64),0_int64,0_int64]
        call fopen(self%table_unit,binfname,access='STREAM',action='WRITE',status='REPLACE',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab_neigh; begin_write; file: '//binfname%to_char(),io_stat)
        self%table_is_open = .true.
        allocate(self%cand_counts(self%nptcls), source=0)
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        write(self%table_unit,pos=addr) self%cand_counts
        addr = addr + sizeof(self%cand_counts)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        write(self%table_unit,pos=addr) self%temps
        addr = addr + sizeof(self%temps)
        self%table_addr    = addr
        self%table_nnz     = 0
        self%table_nchunks = 0
    end subroutine begin_write

    !> dense rows have no place in the sparse part file
    subroutine write_rows_neigh( self, i_first, rows, refined )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first
        type(prob_dense_entry),    intent(in)    :: rows(:,:)
        type(prob_refined_entry),  intent(in)    :: refined(:)
        THROW_HARD('dense rows cannot be written to a neighbourhood part file; write_rows_neigh')
    end subroutine write_rows_neigh

    subroutine flush_candidate_buffers( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ithr, nused, j, i_loc
        if( .not. self%table_is_open ) return
        do ithr = 1,size(self%candidate_buffers)
            nused = self%candidate_buffers(ithr)%nused
            if( nused < 1 ) cycle
            do j = 1,nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                self%cand_counts(i_loc) = self%cand_counts(i_loc) + 1
            enddo
            write(self%table_unit,pos=self%table_addr) nused
            self%table_addr = self%table_addr + sizeof(nused)
            write(self%table_unit,pos=self%table_addr) self%candidate_buffers(ithr)%particle_indices(1:nused)
            self%table_addr = self%table_addr + sizeof(self%candidate_buffers(ithr)%particle_indices(1:nused))
            write(self%table_unit,pos=self%table_addr) self%candidate_buffers(ithr)%candidates(1:nused)
            self%table_addr = self%table_addr + sizeof(self%candidate_buffers(ithr)%candidates(1:nused))
            self%table_nnz     = self%table_nnz + int(nused,int64)
            self%table_nchunks = self%table_nchunks + 1
            call self%candidate_buffers(ithr)%kill
        enddo
    end subroutine flush_candidate_buffers

    ! Completes the sparse part file: the header, the candidate counts, the seed table and the particle
    ! temperatures (posterior_temperature, one real per particle after the seed table; the posterior assignment
    ! draws at them)
    subroutine write_tab( self, binfname )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: binfname
        integer(int64) :: file_header(4), addr
        if( .not. self%table_is_open ) call self%begin_write(binfname)
        call self%flush_candidate_buffers
        ! a partition without particles writes a table without candidates
        if( self%table_nnz < 1 .and. self%nptcls > 0 ) THROW_HARD('eul_prob_tab_neigh%write_tab; empty candidate stream')
        file_header = [int(self%nrefs,int64),int(self%nptcls,int64),self%table_nnz,int(self%table_nchunks,int64)]
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        write(self%table_unit,pos=addr) self%cand_counts
        addr = addr + sizeof(self%cand_counts)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        write(self%table_unit,pos=addr) self%temps
        addr = addr + sizeof(self%temps)
        call fclose(self%table_unit)
        deallocate(self%cand_counts)
        self%table_unit    = -1
        self%table_is_open = .false.
        self%table_addr    = 1
    end subroutine write_tab

    ! Reads one partition candidate stream into the compact global store.
    subroutine read_sparse_tab_to_glob( self, binfname )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: binfname
        integer, parameter :: IO_CHUNK = 65536
        type(prob_candidate), allocatable :: candidates_loc(:)
        real,           allocatable :: seed_shifts_loc(:,:), temps_loc(:)
        logical,        allocatable :: seed_has_sh_loc(:)
        integer,        allocatable :: pinds_loc(:), particle_indices(:), pind2glob(:)
        integer :: funit, io_stat, nrefs_loc, nptcls_loc, nchunks
        integer :: i_loc, i_glob, ichunk, chunk_n, first, last, nread, j, pind, max_pind, seed_nrots_loc
        integer :: candidate_pos
        integer(int64) :: file_header(4), nnz, nnz_read, addr, chunk_indices_addr, chunk_candidates_addr
        if( file_exists(binfname) )then
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab_neigh; read_sparse_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD( 'sparse corr/rot files of partitions should be ready! ' )
        endif
        read(unit=funit,pos=1) file_header
        nrefs_loc  = int(file_header(1))
        nptcls_loc = int(file_header(2))
        nnz        = file_header(3)
        nchunks    = int(file_header(4))
        if( nrefs_loc .ne. self%nrefs ) THROW_HARD('nrefs mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        ! a partition without particles writes a table without candidates
        if( nnz < 1 .and. nptcls_loc > 0 ) THROW_HARD('empty sparse table in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc), temps_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        read(funit, pos=addr) pinds_loc
        addr = addr + 2 * sizeof(pinds_loc)   ! the particle indices and the candidate counts
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        if( nptcls_loc > 0 )then
            if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
            if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
            read(funit, pos=addr) temps_loc
        endif
        addr = addr + sizeof(temps_loc)
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        if( max_pind < 1 )then
            call fclose(funit)
            deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, temps_loc, pind2glob)
            return
        endif
        do i_loc = 1,nptcls_loc
            pind = pinds_loc(i_loc)
            i_glob = 0
            if( pind >= 1 .and. pind <= max_pind ) i_glob = pind2glob(pind)
            if( i_glob > 0 )then
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
                self%temps(i_glob)         = temps_loc(i_loc)
                self%candidate_store%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%candidate_store%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            endif
        enddo
        allocate(particle_indices(IO_CHUNK), candidates_loc(IO_CHUNK))
        nnz_read = 0
        do ichunk = 1,nchunks
            read(funit,pos=addr) chunk_n
            addr = addr + sizeof(chunk_n)
            if( chunk_n < 1 ) THROW_HARD('invalid empty candidate chunk')
            chunk_indices_addr    = addr
            chunk_candidates_addr = chunk_indices_addr + chunk_n * sizeof(chunk_n)
            do first = 1,chunk_n,IO_CHUNK
                last  = min(chunk_n, first + IO_CHUNK - 1)
                nread = last - first + 1
                read(funit,pos=chunk_indices_addr + (first-1)*sizeof(chunk_n)) particle_indices(1:nread)
                read(funit,pos=chunk_candidates_addr + (first-1)*sizeof(candidates_loc(1))) candidates_loc(1:nread)
                do j = 1,nread
                    i_loc = particle_indices(j)
                    if( i_loc < 1 .or. i_loc > nptcls_loc ) THROW_HARD('invalid particle index in sparse candidate stream')
                    pind = pinds_loc(i_loc)
                    i_glob = 0
                    if( pind >= 1 .and. pind <= max_pind ) i_glob = pind2glob(pind)
                    if( i_glob < 1 ) cycle
                    self%candidate_fill_counts(i_glob) = self%candidate_fill_counts(i_glob) + 1
                    candidate_pos = int(self%candidate_store%offsets(i_glob)) + self%candidate_fill_counts(i_glob) - 1
                    if( candidate_pos >= int(self%candidate_store%offsets(i_glob+1)) )&
                        &THROW_HARD('candidate-store overflow in read_sparse_tab_to_glob')
                    self%candidate_store%candidates(candidate_pos) = candidates_loc(j)
                enddo
            enddo
            addr = chunk_candidates_addr + chunk_n * sizeof(candidates_loc(1))
            nnz_read = nnz_read + int(chunk_n,int64)
        enddo
        if( nnz_read /= nnz ) THROW_HARD('candidate stream count mismatch')
        call fclose(funit)
        deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, temps_loc, particle_indices, candidates_loc, pind2glob)
    end subroutine read_sparse_tab_to_glob

    ! Aggregates nparts tables into the global table: the candidate counts every part writes after its
    ! particle indices size the store, then the parts' candidates are placed in file order
    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(string) :: fname
        integer, allocatable :: counts_glob(:), counts_loc(:), pinds_loc(:), pind2glob(:)
        integer :: ipart, funit, io_stat, nptcls_loc, i, pind, i_glob, max_pind
        integer(int64) :: header(4), addr
        allocate(counts_glob(self%nptcls), source=0)
        do ipart = 1,nparts
            fname = fbody//int2str_pad(ipart,numlen)//'.dat'
            call fopen(funit,fname,access='STREAM',action='READ',status='OLD',iostat=io_stat)
            call fileiochk('simple_eul_prob_tab_neigh; read_tabs_to_glob header; file: '//fname%to_char(), io_stat)
            read(funit,pos=1) header
            if( header(1) /= self%nrefs ) THROW_HARD('nrefs mismatch in read_tabs_to_glob')
            nptcls_loc = int(header(2))
            allocate(pinds_loc(nptcls_loc), counts_loc(nptcls_loc))
            addr = sizeof(header) + 1
            read(funit,pos=addr) pinds_loc
            addr = addr + sizeof(pinds_loc)
            read(funit,pos=addr) counts_loc
            call fclose(funit)
            if( sum(int(counts_loc,int64)) /= header(3) ) THROW_HARD('candidate counts disagree with the stream size')
            call build_pind_lookup(self%pinds,pinds_loc,pind2glob,max_pind)
            do i = 1,nptcls_loc
                pind = pinds_loc(i)
                if( pind < 1 .or. pind > max_pind ) cycle
                i_glob = pind2glob(pind)
                if( i_glob > 0 ) counts_glob(i_glob) = counts_glob(i_glob) + counts_loc(i)
            enddo
            deallocate(pinds_loc,counts_loc,pind2glob)
        enddo
        where( counts_glob < 1 ) counts_glob = 1
        call self%candidate_store%new_ragged(self%pinds,counts_glob)
        allocate(self%candidate_fill_counts(self%nptcls), source=0)
        deallocate(counts_glob)
        do ipart = 1, nparts
            fname = fbody//int2str_pad(ipart,numlen)//'.dat'
            call self%read_sparse_tab_to_glob(fname)
        enddo
        deallocate(self%candidate_fill_counts)
    end subroutine read_tabs_to_glob

    ! COMMON STOCHASTIC/SUBSPACE ROUTINES

    ! Reads the previous orientation metadata for one particle
    subroutine get_particle_context( self, iptcl, o_prev, prev_state, prev_proj )
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(in)    :: iptcl
        type(ori),                 intent(inout) :: o_prev
        integer,                   intent(out)   :: prev_state, prev_proj
        call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
        prev_state = o_prev%get_state()
        prev_proj  = self%b_ptr%eulspace%find_closest_proj(o_prev)
    end subroutine get_particle_context

    ! Records one particle vs. ref evaluation result in its thread-owned buffer.
    subroutine record_sparse_eval( self, i, ithr, ri, dist, irot, x, y, has_sh )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i, ithr, ri, irot
        real,                      intent(in)    :: dist, x, y
        logical,                   intent(in)    :: has_sh
        type(prob_candidate) :: candidate
        candidate%iref   = self%ref_full(ri)
        candidate%dist   = dist
        candidate%inpl   = irot
        candidate%x      = x
        candidate%y      = y
        candidate%has_sh = has_sh
        call self%candidate_buffers(ithr)%append_or_replace(i, candidate)
    end subroutine record_sparse_eval

    ! Minimizes shift at the previous best orientation
    !> the cost record of one neighbourhood fill (bench file): the wall time of the parallel region, thread-seconds of
    !! the seed, the score loop and the refinement, references scored, shifted particle memos built, and the search
    !! objects' counters
    subroutine record_neigh_fill_profile( self, grad_obj, wall, ts_seed, ts_score, ts_refine, nscored, nbank, nfilled )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad),    intent(inout) :: grad_obj(:)
        real(dp),                  intent(in)    :: wall, ts_seed(:), ts_score(:), ts_refine(:)
        integer(int64),            intent(in)    :: nscored(:), nbank(:)
        integer,                   intent(in)    :: nfilled
        integer(int64) :: obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr
        integer        :: ithr
        self%fill_profile%wall_region = self%fill_profile%wall_region + wall
        self%fill_profile%ts_seed     = self%fill_profile%ts_seed     + sum(ts_seed)
        self%fill_profile%ts_score    = self%fill_profile%ts_score    + sum(ts_score)
        self%fill_profile%ts_refine   = self%fill_profile%ts_refine   + sum(ts_refine)
        self%fill_profile%nthr        = nthr_glob
        self%fill_profile%nptcls      = self%fill_profile%nptcls  + int(nfilled, int64)
        self%fill_profile%nscored     = self%fill_profile%nscored + sum(nscored)
        self%fill_profile%bank_builds = self%fill_profile%bank_builds + sum(nbank)
        if( .not. self%p_ptr%l_doshift ) return
        do ithr = 1,size(grad_obj)
            call grad_obj(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
            self%fill_profile%objective_evals = self%fill_profile%objective_evals + obj_ithr
            self%fill_profile%gradient_evals  = self%fill_profile%gradient_evals  + grad_ithr
            self%fill_profile%alternations    = self%fill_profile%alternations    + alt_ithr
            self%fill_profile%grid_points     = self%fill_profile%grid_points     + grid_ithr
            self%fill_profile%minimize_calls  = self%fill_profile%minimize_calls  + min_ithr
        enddo
    end subroutine record_neigh_fill_profile

    subroutine estimate_shift_seed( self, grad_obj, ithr, iptcl, prev_state, prev_proj, o_prev, l_with_shift, shift_seed )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad),    intent(inout) :: grad_obj(:)
        integer,                   intent(in)    :: ithr, iptcl, prev_state, prev_proj
        type(ori),                 intent(inout) :: o_prev
        logical,                   intent(in)    :: l_with_shift
        real,                      intent(inout) :: shift_seed(3)
        integer :: irot, iref
        shift_seed = 0.
        if( .not. l_with_shift ) return
        if( prev_state >= 1 .and. prev_state <= self%p_ptr%nstates .and.&
            &prev_proj >= 1 .and. prev_proj <= self%p_ptr%nspace )then
            if( self%state_exists(prev_state) )then
                irot = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                iref = (prev_state-1)*self%p_ptr%nspace + prev_proj
                call grad_obj(ithr)%set_indices(iref, iptcl)
                shift_seed = grad_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                if( irot == 0 )then
                    shift_seed(2:3) = 0.
                else
                    ! one particle-frame shift from here on: the estimate converted at the rotation it ended at
                    shift_seed(2:3) = shift_ref2ptcl(shift_seed(2:3), self%b_ptr%pftc%get_rot(irot))
                endif
            endif
        endif
    end subroutine estimate_shift_seed

    !> the per-thread memo of the particle shifted by its particle-frame seed, which every reference of the particle is
    !! scored against (a seed below SHERRSQ scores against the plain particle memo)
    subroutine alloc_seed_bank( self, seed_bank )
        class(eul_prob_tab_neigh), intent(in)    :: self
        complex(sp), allocatable,  intent(inout) :: seed_bank(:,:,:)
        integer :: kfromto(2)
        kfromto = self%b_ptr%pftc%get_kfromto()
        allocate(seed_bank(self%b_ptr%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1, nthr_glob))
    end subroutine alloc_seed_bank

    !> shifts the particle once by its seed (l_bank: the bank holds it) and counts the memos built
    subroutine memoize_seed_bank( self, iptcl, shift_seed, bank, l_bank, nbank )
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(in)    :: iptcl
        real,                      intent(in)    :: shift_seed(3)
        complex(sp),               intent(inout) :: bank(:,:)
        logical,                   intent(out)   :: l_bank
        integer(int64),            intent(inout) :: nbank
        l_bank = sum(shift_seed(2:3)**2) > SHERRSQ
        if( .not. l_bank ) return
        call self%b_ptr%pftc%memoize_shifted_ptcl(iptcl, shift_seed(2:3), bank)
        nbank = nbank + 1_int64
    end subroutine memoize_seed_bank

    !> one reference scored against the particle shifted by shift_xy (particle frame): from the seed bank when it holds
    !! the shifted particle, with the angular Nyquist bin, otherwise by the shifted-reference kernel. The in-plane
    !! rotation is drawn from the row's tempered posterior at the particle's temperature over every grid rotation;
    !! corr is exp(-dist)
    subroutine score_seeded_ref( self, full_ref, iptcl, shift_xy, bank, l_bank, dist, corr, irot )
        class(eul_prob_tab_neigh), intent(in)  :: self
        integer,                   intent(in)  :: full_ref, iptcl
        real,                      intent(in)  :: shift_xy(2)
        complex(sp),               intent(in)  :: bank(:,:)
        logical,                   intent(in)  :: l_bank
        real,                      intent(out) :: dist, corr
        integer,                   intent(out) :: irot
        if( l_bank )then
            call self%b_ptr%pftc%gen_posterior_val_bank(full_ref, iptcl, bank,&
                &self%b_ptr%pftc%ptcl_posterior_temp(iptcl), dist, irot)
        else
            call self%b_ptr%pftc%gen_posterior_val(full_ref, iptcl, shift_xy,&
                &self%b_ptr%pftc%ptcl_posterior_temp(iptcl), dist, irot)
        endif
        corr = eulprob_corr_switch(dist, self%p_ptr%cc_objfun)
    end subroutine score_seeded_ref

    ! Optimizes shift & rotation on the top-npeaks_inpl references
    subroutine refine_best_neighbors( self, ew, grad_obj, l_seed_sh_first, max_refs_to_refine,&
        &i, ithr, iptcl, shift_seed, neval)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(eval_ws),             intent(inout) :: ew
        type(pftc_shsrch_grad),    intent(inout) :: grad_obj(:)
        logical,                   intent(in)    :: l_seed_sh_first
        integer,                   intent(in)    :: max_refs_to_refine, i, ithr, iptcl, neval
        real,                      intent(in)    :: shift_seed(3)
        integer :: nrefs_to_refine, nstate_eval, j, eval_slot, ri, istate, iproj, irot
        integer :: si, state_eval
        real    :: refined_shift(3), xy_start(2)
        ew%best_evals(:,ithr) = 0
        if( neval <= 0 ) return
        do si = 1, self%nstates
            state_eval = self%ssinds(si)
            ew%state_eval_dists(1:neval,ithr) = huge(1.0)
            nstate_eval = 0
            do j = 1, neval
                ri = ew%evaluated_ref_ids(j,ithr)
                if( ri < 1 ) cycle
                if( self%ref_state(ri) /= state_eval ) cycle
                nstate_eval = nstate_eval + 1
                ew%state_eval_dists(j,ithr) = ew%evaluated_ref_dists(j,ithr)
            enddo
            nrefs_to_refine = min(max_refs_to_refine, nstate_eval)
            if( nrefs_to_refine < 1 ) cycle
            ew%best_evals(1:nrefs_to_refine,ithr) = &
                &minnloc(ew%state_eval_dists(1:neval,ithr), nrefs_to_refine)
            do j = 1, nrefs_to_refine
                eval_slot = ew%best_evals(j,ithr)
                if( eval_slot < 1 ) cycle
                ri = ew%evaluated_ref_ids(eval_slot,ithr)
                if( ri < 1 ) cycle
                istate = self%ref_state(ri)
                iproj  = self%ref_proj(ri)
                call grad_obj(ithr)%set_indices(self%ref_full(ri), iptcl)
                irot = self%candidate_buffers(ithr)%get_inpl(i, self%ref_full(ri))
                if( l_seed_sh_first )then
                    ! the particle-frame seed as the shift of this reference at the rotation it was scored at, in the box
                    xy_start = shift_seed(2:3)
                    if( irot > 0 ) xy_start = shift_ptcl2ref(shift_seed(2:3), self%b_ptr%pftc%get_rot(irot))
                    xy_start = max(-self%p_ptr%trs, min(self%p_ptr%trs, xy_start))
                    refined_shift = grad_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=xy_start)
                else
                    refined_shift = grad_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                endif
                if( irot > 0 )then
                    call record_sparse_eval(self, i, ithr, ri, &
                        &eulprob_dist_switch(refined_shift(1), self%p_ptr%cc_objfun), &
                        &irot, refined_shift(2), refined_shift(3), .true.)
                endif
            enddo
        enddo
    end subroutine refine_best_neighbors

    ! DESTRUCTOR

    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ithr
        if( allocated(self%candidate_buffers) )then
            do ithr = 1,size(self%candidate_buffers)
                call self%candidate_buffers(ithr)%kill
            enddo
            deallocate(self%candidate_buffers)
        endif
        if( allocated(self%cand_counts) ) deallocate(self%cand_counts)
        self%table_nnz     = 0
        self%table_nchunks = 0
        call self%candidate_store%kill
        if( allocated(self%candidate_fill_counts) ) deallocate(self%candidate_fill_counts)
        self%l_direct_stoch_neigh = .false.
        call self%eul_prob_tab%kill
    end subroutine kill_neigh

    ! EVAL_WS bound procedures

    ! Allocates all stochastic-path arrays
    subroutine alloc_eval_stoch_ws( self, nrefs, nstates, nspace, max_refine, nthr )
        class(eval_ws), intent(inout) :: self
        integer,           intent(in)    :: nrefs, nstates, nspace, max_refine, nthr
        call self%dealloc_eval_ws
        allocate(self%best_evals(max(1,max_refine),nthr), self%evaluated_ref_ids(nrefs,nthr), source=0)
        allocate(self%fullref_to_sparse_ref(nstates*nspace), source=0)
        allocate(self%evaluated_ref_dists(nrefs,nthr), source=huge(1.0))
        allocate(self%state_eval_dists(nrefs,nthr),    source=huge(1.0))
    end subroutine alloc_eval_stoch_ws

    ! Allocates all subspace-path arrays
    subroutine alloc_eval_subspace_ws( self, nrefs, nstates, nspace, nsubs, max_refine, nthr )
        class(eval_ws), intent(inout) :: self
        integer,           intent(in)    :: nrefs, nstates, nspace, nsubs, max_refine, nthr
        call self%dealloc_eval_ws
        allocate(self%best_evals(max(1,max_refine),nthr), self%evaluated_ref_ids(nrefs,nthr), source=0)
        allocate(self%fullref_to_sparse_ref(nstates*nspace), source=0)
        allocate(self%sub_ref_counts(nsubs,nstates), self%sub_ref_offsets(nsubs,nstates),&
            &self%sub_ref_list(nrefs), source=0)
        allocate(self%evaluated_ref_dists(nrefs,nthr), source=huge(1.0))
        allocate(self%state_eval_dists(nrefs,nthr),    source=huge(1.0))
    end subroutine alloc_eval_subspace_ws

        ! Allocates and initialises eval_ws arrays required by the direct stochastic path.
    subroutine init_eval_stoch_ws( self, tabneigh, max_refs_to_refine )
        class(eval_ws),            intent(inout) :: self
        class(eul_prob_tab_neigh), intent(in)    :: tabneigh
        integer,                   intent(in)    :: max_refs_to_refine
        integer :: iref, iref_full
        call self%alloc_eval_stoch_ws(tabneigh%nrefs, tabneigh%p_ptr%nstates, tabneigh%p_ptr%nspace, max_refs_to_refine, nthr_glob)
        do iref = 1, tabneigh%nrefs
            iref_full = tabneigh%ref_full(iref)
            if( iref_full >= 1 .and. iref_full <= size(self%fullref_to_sparse_ref) )&
                &self%fullref_to_sparse_ref(iref_full) = iref
        enddo
    end subroutine init_eval_stoch_ws

    ! Allocates and initialises eval_ws arrays required by the subspace neighborhood path.
    subroutine init_eval_subspace_ws( self, tabneigh, nsubs, max_refs_to_refine )
        class(eval_ws),            intent(inout) :: self
        class(eul_prob_tab_neigh), intent(in)    :: tabneigh
        integer,                   intent(in)    :: nsubs, max_refs_to_refine
        integer :: i, ri, istate, isub, iref, iref_full
        call self%alloc_eval_subspace_ws(tabneigh%nrefs, tabneigh%p_ptr%nstates, tabneigh%p_ptr%nspace,&
            &nsubs, max_refs_to_refine, nthr_glob)
        do iref = 1, tabneigh%nrefs
            iref_full = tabneigh%ref_full(iref)
            if( iref_full >= 1 .and. iref_full <= size(self%fullref_to_sparse_ref) )&
                &self%fullref_to_sparse_ref(iref_full) = iref
        enddo
        do ri = 1,tabneigh%nrefs
            istate = tabneigh%ref_state(ri)
            isub   = tabneigh%b_ptr%subspace_full2sub_map(tabneigh%ref_proj(ri))
            if( isub < 1 .or. isub > nsubs ) cycle
            self%sub_ref_counts(isub,istate) = self%sub_ref_counts(isub,istate) + 1
        enddo
        i = 1
        do istate = 1,tabneigh%p_ptr%nstates
            do isub = 1,nsubs
                self%sub_ref_offsets(isub,istate) = i
                i = i + self%sub_ref_counts(isub,istate)
            enddo
        enddo
        self%sub_ref_counts = 0
        do ri = 1,tabneigh%nrefs
            istate = tabneigh%ref_state(ri)
            isub   = tabneigh%b_ptr%subspace_full2sub_map(tabneigh%ref_proj(ri))
            if( isub < 1 .or. isub > nsubs ) cycle
            i = self%sub_ref_offsets(isub,istate) + self%sub_ref_counts(isub,istate)
            self%sub_ref_list(i) = ri
            self%sub_ref_counts(isub,istate) = self%sub_ref_counts(isub,istate) + 1
        enddo
    end subroutine init_eval_subspace_ws

    ! Deallocates all allocatable fields of an eval_ws workspace
    subroutine dealloc_eval_ws( self )
        class(eval_ws), intent(inout) :: self
        if( allocated(self%evaluated_ref_ids)     ) deallocate(self%evaluated_ref_ids)
        if( allocated(self%best_evals)            ) deallocate(self%best_evals)
        if( allocated(self%fullref_to_sparse_ref) ) deallocate(self%fullref_to_sparse_ref)
        if( allocated(self%evaluated_ref_dists)   ) deallocate(self%evaluated_ref_dists)
        if( allocated(self%state_eval_dists)      ) deallocate(self%state_eval_dists)
        if( allocated(self%sub_ref_counts)        ) deallocate(self%sub_ref_counts)
        if( allocated(self%sub_ref_offsets)       ) deallocate(self%sub_ref_offsets)
        if( allocated(self%sub_ref_list)          ) deallocate(self%sub_ref_list)
    end subroutine dealloc_eval_ws

    ! COARSE_SEARCH_WS bound procedures

    ! Allocates all arrays of a coarse_search_ws workspace for the subspace path.
    subroutine alloc_coarse_ws( self, npeak_target, npooled_capacity, nstates, nthr )
        class(coarse_search_ws), intent(inout) :: self
        integer,                 intent(in)    :: npeak_target, npooled_capacity, nstates, nthr
        call self%dealloc_coarse_ws
        allocate(self%peak_subspace_dists(npeak_target,nstates,nthr), source=huge(1.0))
        allocate(self%peak_subspace_inds(npeak_target,nstates,nthr),  source=0)
        allocate(self%pooled_sub_inds(npooled_capacity,nstates,nthr), source=0)
        allocate(self%pooled_sub_count(nstates,nthr),                 source=0)
        allocate(self%peak_subspace_count(nstates,nthr),              source=0)
    end subroutine alloc_coarse_ws

    ! Safely deallocates all allocatable fields of a coarse_search_ws workspace.
    subroutine dealloc_coarse_ws( self )
        class(coarse_search_ws), intent(inout) :: self
        if( allocated(self%peak_subspace_count) ) deallocate(self%peak_subspace_count)
        if( allocated(self%pooled_sub_count)    ) deallocate(self%pooled_sub_count)
        if( allocated(self%pooled_sub_inds)     ) deallocate(self%pooled_sub_inds)
        if( allocated(self%peak_subspace_inds)  ) deallocate(self%peak_subspace_inds)
        if( allocated(self%peak_subspace_dists) ) deallocate(self%peak_subspace_dists)
    end subroutine dealloc_coarse_ws

end module simple_eul_prob_tab_neigh
