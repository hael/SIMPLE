!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is geometrically sparse and searched independently per state.
module simple_eul_prob_tab_neigh
use, intrinsic :: iso_fortran_env, only: int64
use omp_lib,                   only: omp_get_max_threads
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_utils, only: update_screen_best, eulprob_dist_switch, eulprob_corr_switch, &
    &record_seed_shift, prob_candidate, &
    &prob_candidate_buffer, prob_candidate_store, prob_dense_entry, prob_refined_entry, posterior_temperature, &
    &draw_posterior_row, prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN, top_posterior_extras, &
    &temperature_or_argmin, prob_rnd_scope, keyed_candidate_order, prob_readout, prob_readout_fname, &
    &write_prob_readout, report_prob_readout, PROB_READOUT_3D_NEIGH
use simple_rnd,                only: rnd_keyed, rnd_draw_key, RND_FAM_ASSIGN
use simple_decay_funs,         only: extremal_decay
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_eul_prob_tab_utils, only: shift_needs_disk, disk_radius_px, screen_pool_size, map_prior_weight, insert_pose_loss, &
    &shift_prior_sigma2
use simple_defs,               only: PROB_INSERT_WMIN
use simple_polarft_calc, only: SHERRSQ, SHIFTMAP_SCREEN, SHIFTMAP_WINDOW
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
    type(prob_candidate_buffer), allocatable :: candidate_buffers(:) !< per thread; the batch's store at the assignment
    integer,        allocatable :: candidate_fill_counts(:)
    logical                     :: l_direct_stoch_neigh = .false.
    integer                     :: nkept = 0   !< particles without a valid candidate, orientation kept (the assignment's count)
contains
    procedure :: new_neigh_pipeline
    procedure :: set_model_from_field => set_model_from_field_neigh
    procedure :: fill_tab_range => fill_tab_neigh_range
    procedure :: open_assignment  => open_assignment_neigh
    procedure :: assign_batch     => assign_batch_neigh
    procedure :: close_assignment => close_assignment_neigh
    procedure, private :: assign_rows_neigh
    procedure, private :: batch_store_from_buffers
    procedure :: kill          => kill_neigh
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

    !> the pipeline's neighbourhood table of a part (the shared-memory matcher): the assignment map, the thread
    !! buffers the fill keeps for the batch's own selection (assign_batch); no part file
    subroutine new_neigh_pipeline( self, params, build, pinds )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%eul_prob_tab%new_pipeline(params,build,pinds)
        allocate(self%candidate_buffers(nthr_glob))
        call validate_neigh_mode(self)
    end subroutine new_neigh_pipeline

    !> the neighbourhood fill has no refined-reference budget; the shift prior's width from the complete incoming field
    subroutine set_model_from_field_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        real :: sigma2_prior
        sigma2_prior = 0.
        if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
        call self%set_iteration_model(0, sigma2_prior)
    end subroutine set_model_from_field_neigh

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
        ! the polar shift map lives in the dense fill; this fill's increments are the window search's
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
    end subroutine fill_tab_neigh_range

    ! STOCHASTIC NEIGHBORHOOD BRANCH

    subroutine fill_tab_stoch_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad), target :: grad_shsrch_obj(nthr_glob)
        type(pftc_shsrch_grad), target :: grad_shsrch_disk(nthr_glob) !< the disk's (shift_sigma > 0)
        type(eval_ws)               :: eval_work
        integer,        allocatable :: direct_srch_order(:,:)
        real,           allocatable :: corrs_inpl(:,:)
        complex(sp),    allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        integer :: i, ithr, max_refs_to_refine, nfull_refs
        real    :: sigma2_prior  !< the shift prior's variance per axis (shift_prior=yes)
        integer :: iptcl, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_shc_neigh, l_snhc_neigh, l_seed_sh_first, l_bank(nthr_glob), l_map
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        real(dp)       :: prior_w
        real           :: step_win, step_disk, lims_disk(2,2)
        integer        :: kscreen
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob)
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_shc_neigh     = (trim(self%p_ptr%prob_neigh_mode) == 'shc')
        l_snhc_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'snhc')
        ! under the polar shift map's key the window map initialises every polish and there is no seed
        l_map           = self%p_ptr%l_centring
        l_seed_sh_first = self%p_ptr%l_doshift .and. l_shc_neigh .and. .not. l_map
        prior_w         = 0.d0
        kscreen         = 0
        step_win        = 0.
        step_disk       = 0.
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
        if( .not. self%l_model_set ) THROW_HARD('the iteration model is set by the caller before the fill (set_iteration_model); eul_prob_tab_neigh%fill_tab_range')
        ! the keyed draws of this fill: the candidate order per particle and the in-plane draw per scored reference
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .false.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
        call alloc_seed_bank(self, seed_bank)
        nfull_refs = self%p_ptr%nstates * self%p_ptr%nspace
        allocate(corrs_inpl(nrots,nthr_glob))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_stoch_ws(self, max_refs_to_refine)
        allocate(direct_srch_order(nfull_refs,nthr_glob))
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
            ! the shift prior's width of the iteration (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = self%sigma2_prior_model
            if( l_map )then
                call setup_shiftmap_neigh(self, prior_w, kscreen, step_win, step_disk)
                ! the disk's own search objects: the polish box clamps to the object's limits
                lims_disk(:,1) = -disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd, self%p_ptr%box)
                lims_disk(:,2) = -lims_disk(:,1)
                do ithr = 1,nthr_glob
                    call grad_shsrch_disk(ithr)%new_alternating(self%b_ptr, lims_disk, lims_init=lims_init, &
                        &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
                end do
            endif
            ! dynamic: particles differ in cost (stopping rule, shift searches); each is processed
            ! whole by one thread with its own search object; its draws are keyed on the particle
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
            &nscored, nbank, i_to - i_from + 1, grad_disk=grad_shsrch_disk)
        ! cleanup
        call eval_work%dealloc_eval_ws
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
            if( l_map ) call grad_shsrch_disk(ithr)%kill
        end do
        deallocate(direct_srch_order)
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
            logical   :: l_disk, l_prior
            real      :: shsig
            call get_particle_context(self, iptcl_loc, o_prev, prev_state_loc, prev_proj_loc)
            l_disk  = .false.
            l_prior = .false.
            if( l_map )then
                shsig   = self%b_ptr%spproj_field%get(iptcl_loc, 'shsig')
                l_disk  = shift_needs_disk(shsig, self%p_ptr%trs)
                l_prior = shsig < 0.
            endif
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
                &i_loc, ithr_loc, iptcl_loc, shift_seed, neval_loc, l_map, kscreen, step_win, step_disk, l_disk,&
                &l_prior, grad_shsrch_disk, prior_w)
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
            prev_full_ref_loc = 0
            if( prev_state_loc >= 1 .and. prev_state_loc <= self%p_ptr%nstates .and.&
                &prev_proj_loc >= 1 .and. prev_proj_loc <= self%p_ptr%nspace )then
                prev_full_ref_loc = (prev_state_loc - 1) * self%p_ptr%nspace + prev_proj_loc
            endif
            ! the candidate order: the previous reference last, the others keyed on the particle
            call keyed_candidate_order(self%rnd_scope, iptcl_loc, nfull_refs, prev_full_ref_loc, direct_srch_order(:,ithr_loc))
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
        type(pftc_shsrch_grad), target :: grad_shsrch_obj(nthr_glob)
        type(pftc_shsrch_grad), target :: grad_shsrch_disk(nthr_glob) !< the disk's (shift_sigma > 0)
        type(coarse_search_ws) :: coarse_ws
        type(eval_ws)          :: eval_work
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        integer :: i, ithr, max_refs_to_refine, nsubs, npeak_target, npooled_capacity
        real    :: sigma2_prior  !< the shift prior's variance per axis (shift_prior=yes)
        integer :: iptcl, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_geom_neigh, l_state_neigh, l_seed_sh_first, l_bank(nthr_glob), l_map
        real(dp)       :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        real(dp)       :: prior_w
        real           :: step_win, step_disk, lims_disk(2,2)
        integer        :: kscreen
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob)
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_geom_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'geom')
        l_state_neigh   = (trim(self%p_ptr%prob_neigh_mode) == 'state')
        ! Keep shift-first seeding single-state only. In multi-state runs this
        ! can bias all states toward the same previous-state shift basin. Under the polar shift map's key the
        ! window map initialises every polish and there is no seed
        l_map           = self%p_ptr%l_centring
        l_seed_sh_first = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1 .and. .not. l_map
        prior_w         = 0.d0
        kscreen         = 0
        step_win        = 0.
        step_disk       = 0.
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
        if( .not. self%l_model_set ) THROW_HARD('the iteration model is set by the caller before the fill (set_iteration_model); eul_prob_tab_neigh%fill_tab_subspace_range')
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .false.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
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
            ! the shift prior's width of the iteration (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = self%sigma2_prior_model
            if( l_map )then
                call setup_shiftmap_neigh(self, prior_w, kscreen, step_win, step_disk)
                ! the disk's own search objects: the polish box clamps to the object's limits
                lims_disk(:,1) = -disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd, self%p_ptr%box)
                lims_disk(:,2) = -lims_disk(:,1)
                do ithr = 1,nthr_glob
                    call grad_shsrch_disk(ithr)%new_alternating(self%b_ptr, lims_disk, lims_init=lims_init, &
                        &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
                end do
            endif
            ! dynamic: particles differ in cost (stopping rule, shift searches); each is processed
            ! whole by one thread with its own search object; its draws are keyed on the particle
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
            &nscored, nbank, i_to - i_from + 1, grad_disk=grad_shsrch_disk)
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
            if( l_map ) call grad_shsrch_disk(ithr)%kill
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
            logical   :: l_disk, l_prior
            real      :: shsig
            call get_particle_context(self, iptcl_loc, o_prev, prev_state, prev_proj)
            l_disk  = .false.
            l_prior = .false.
            if( l_map )then
                shsig   = self%b_ptr%spproj_field%get(iptcl_loc, 'shsig')
                l_disk  = shift_needs_disk(shsig, self%p_ptr%trs)
                l_prior = shsig < 0.
            endif
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
                    &i_loc, ithr_loc, iptcl_loc, shift_seed_loc, neval, l_map, kscreen, step_win, step_disk, l_disk,&
                    &l_prior, grad_shsrch_disk, prior_w)
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

    !> the assignment opened: the weighted insertion's table and its cone's geometry (topk > 1, the modes whose
    !! candidates are shift-refined), the n_eff readout over the rows of the global sample (nglob its size, rank0 the
    !! rank before this table's first row: a part's offset), the draws' scope on the table and the calculator (the
    !! fallback's in-plane draw), the counts close_assignment reports
    subroutine open_assignment_neigh( self, nglob, rank0 )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, optional,         intent(in)    :: nglob, rank0
        integer :: nglob_here, rank0_here
        nglob_here = self%nptcls
        rank0_here = 0
        if( present(nglob) ) nglob_here = nglob
        if( present(rank0) ) rank0_here = rank0
        if( self%l_assign_open ) THROW_HARD('the assignment is open; eul_prob_tab_neigh%open_assignment')
        ! the weighted insertion's extras (topk > 1): from the refined candidates of the state's row, in the
        ! modes whose candidates are shift-refined (state, geom); the shc/snhc neighbourhoods of the early
        ! stages insert the committed pose alone, their candidates being the hill climb's unrefined proposals
        self%kmax = max(1, self%p_ptr%topk)
        select case(trim(self%p_ptr%prob_neigh_mode))
            case('shc','snhc')
                self%l_extras = .false.
            case DEFAULT
                self%l_extras = self%kmax > 1
        end select
        ! the insertion cone's geometry: the direction space, the point group and the rotation grid of the table
        self%l_extras = self%l_extras .and. self%insert_cone%is_open()
        if( self%l_extras )then
            if( self%b_ptr%eulspace%get_noris() /= self%p_ptr%nspace ) THROW_HARD('the direction space is not built; the weighted insertion needs it; open_assignment_neigh')
            if( self%seed_nrots < 1 ) THROW_HARD('the table carries no rotation grid; the weighted insertion needs it; open_assignment_neigh')
            call self%insert_cone%set_space(self%b_ptr%eulspace, self%b_ptr%pgrpsyms, self%seed_nrots)
        endif
        call self%insert_tab%new(self%nptcls, self%pinds, merge(self%kmax, 1, self%l_extras))
        call self%neff%new(nglob_here, PROB_NEFF_ROWS, PROB_NEFF_NMIN, rank0=rank0_here)
        ! the keyed draws of the assignment: the fallback's in-plane draw on the calculator, the draw per particle
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .false.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
        self%nmode         = 0
        self%ntemp_missing = 0
        self%nkept         = 0
        self%l_assign_open = .true.
    end subroutine open_assignment_neigh

    !> the batch's selection: the candidates the fill kept in its thread buffers become the batch's store, in their
    !! append order, then the rows are assigned and the store released
    subroutine assign_batch_neigh( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        call self%batch_store_from_buffers(i_first, i_last)
        call self%assign_rows_neigh(i_first, i_last)
        call self%candidate_store%kill
    end subroutine assign_batch_neigh

    !> the thread buffers of the batch as the candidate store of the part: the particles outside the batch hold no
    !! candidate; a particle's candidates come in its thread's append order
    subroutine batch_store_from_buffers( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        integer, allocatable :: counts(:), filled(:)
        integer :: ithr, j, i_loc, pos
        allocate(counts(self%nptcls), filled(self%nptcls), source=0)
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                if( i_loc < i_first .or. i_loc > i_last ) THROW_HARD('a candidate of another batch in the buffers; eul_prob_tab_neigh%batch_store_from_buffers')
                counts(i_loc) = counts(i_loc) + 1
            enddo
        enddo
        call self%candidate_store%new_ragged(self%pinds, counts)
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                filled(i_loc) = filled(i_loc) + 1
                pos = int(self%candidate_store%offsets(i_loc)) + filled(i_loc) - 1
                self%candidate_store%candidates(pos) = self%candidate_buffers(ithr)%candidates(j)
            enddo
            call self%candidate_buffers(ithr)%kill
        enddo
        do i_loc = i_first,i_last
            self%candidate_store%seed_shifts(:,i_loc) = self%seed_shifts(:,i_loc)
            self%candidate_store%seed_has_sh(i_loc)   = self%seed_has_sh(i_loc)
        enddo
        deallocate(counts, filled)
    end subroutine batch_store_from_buffers

    !> the assignment closed: its readouts as a record, published by a distributed part or reported by the owner (the
    !! "no FSC yet" line only when no cone is open: the shc/snhc neighbourhood stages insert the committed pose alone)
    subroutine close_assignment_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(prob_readout) :: rec
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab_neigh%close_assignment')
        rec%route         = PROB_READOUT_3D_NEIGH
        rec%nptcls        = self%nptcls
        rec%nrefs         = self%nrefs
        rec%nmode         = self%nmode
        rec%ntemp_missing = self%ntemp_missing
        rec%nkept         = self%nkept
        rec%kmax          = self%kmax
        rec%l_extras      = self%l_extras
        if( .not. self%l_extras .and. self%insert_cone%is_open() ) rec%kmax = 1
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
    end subroutine close_assignment_neigh

    !> ptcl -> (proj, state) over the rows i_first..i_last of the batch's candidate store, by one draw per particle from its tempered posterior over its evaluated candidates. The
    !! particle's state is the state of its closest valid candidate over all states (an active reference with an
    !! in-plane index; ties to the first in store order), as the dense table labels it; within that state p_c
    !! proportional to exp(-T_i (L_c - L_min)) over the valid candidates, T_i the particle's temperature from the fill
    !! (temps) scaled by prob_beta; the mean n_eff of a subsample of the rows is logged at the close. The drawn
    !! candidate is committed with its shift, the seed shift recorded, frac the evaluated fraction; no reference is
    !! forced to be populated. A particle without a valid candidate (a state emptied mid-run) is first re-seeded at its
    !! previous orientation when its state is active, else keeps its orientation; one without a temperature (none
    !! expected) takes its argmin. The draws are keyed on the particle.
    subroutine assign_rows_neigh( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        integer :: i, candidate_pos
        logical :: l_seed_fallback_shift
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab_neigh%assign_rows_neigh')
        l_seed_fallback_shift = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1 .and. &
            &trim(self%p_ptr%prob_neigh_mode) /= 'snhc'
        do i = i_first,i_last
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
        !! argmin with a warning. The draws are keyed on the particle.
        subroutine assign_particles_posterior()
            real,    allocatable :: w(:), dists_c(:)
            integer, allocatable :: pos_c(:), which_pos(:)
            logical, allocatable :: l_refined(:)
            integer :: i_loc, n_loc, which, ibest, maxcand, istate_best
            real    :: temp, dist
            maxcand = 0
            do i_loc = i_first,i_last
                maxcand = max(maxcand, int(self%candidate_store%offsets(i_loc+1) - self%candidate_store%offsets(i_loc)))
            enddo
            allocate(w(max(1,maxcand)), dists_c(max(1,maxcand)), pos_c(max(1,maxcand)))
            if( self%l_extras )then
                allocate(l_refined(i_first:i_last), source=.false.)
                allocate(which_pos(i_first:i_last), source=0)
            endif
            do i_loc = i_first,i_last
                call gather_state_row(i_loc, dists_c, pos_c, n_loc, ibest, istate_best)
                if( istate_best == 0 )then
                    self%nkept = self%nkept + 1
                    cycle
                endif
                ! the n_eff readout takes the rows it wants in particle order
                if( self%neff%wanted(i_loc) ) call self%neff%offer(n_loc, dists_c(1:n_loc), self%temps(i_loc))
                temp = self%p_ptr%prob_beta * self%temps(i_loc)
                if( .not. (temp > 0.) )then
                    self%ntemp_missing = self%ntemp_missing + 1
                    temp = huge(temp)
                endif
                call draw_posterior_row(n_loc, dists_c(1:n_loc), temp, w, &
                    &rnd_keyed(rnd_draw_key(self%rnd_scope, self%pinds(i_loc), RND_FAM_ASSIGN)), which, dist)
                if( which < 1 ) which = ibest
                if( which == ibest ) self%nmode = self%nmode + 1
                call self%assign_candidate(i_loc, self%candidate_store%candidates(pos_c(which)))
                if( self%l_extras )then
                    ! the extras of the weighted insertion are weighed against the committed candidate's distance,
                    ! so the committed candidate must be a refined one too (an unrefined one's distance is at the
                    ! seed shift); its position in the row's gathered order is kept for the selection
                    l_refined(i_loc) = self%candidate_store%candidates(pos_c(which))%has_sh
                    which_pos(i_loc) = which
                endif
                call record_seed_shift(self%assgn_map(i_loc), self%seed_shifts(:,i_loc), self%seed_has_sh(i_loc),&
                    &self%p_ptr%l_doshift)
                self%assgn_map(i_loc)%frac = search_frac(i_loc)
            enddo
            if( self%l_extras ) call select_extras(self%kmax, l_refined, which_pos, maxcand)
            deallocate(w, dists_c, pos_c)
            if( allocated(l_refined) ) deallocate(l_refined, which_pos)
        end subroutine assign_particles_posterior

        !> the extras of the weighted insertion: for every particle whose committed candidate is refined, the refined
        !! candidates of its state's row whose pose lies within the insertion cone of the committed one (the rotation
        !! modulo the point group and the shift together moving the rim by less than the map resolves), the kmax - 1 of
        !! highest posterior at the particle's own temperature, with the weights that posterior renormalised over the
        !! kept set (top_posterior_extras): the draw decides the hypothesis, the extras only spread the particle over
        !! the poses the map resolves alike. Under the map the refined candidates are weighed by the dense table's
        !! insertion loss (insert_pose_loss: the map's prior on the pose's increment), not their data term alone.
        !! Deterministic given the committed poses, so the particles run in parallel after the draws, each thread
        !! gathering the row again on its own arrays; the table takes the result in particle order
        subroutine select_extras( kmax, l_refined, which_pos, maxcand )
            integer, intent(in) :: kmax, maxcand
            integer, intent(in) :: which_pos(i_first:i_last)
            logical, intent(in) :: l_refined(i_first:i_last)
            type(prob_candidate) :: cands(kmax)
            integer, allocatable :: kpos(:,:), kinds(:,:), nkeep(:), nelig(:), pos_c(:,:)
            real,    allocatable :: kw(:,:), dists_c(:,:)
            logical, allocatable :: eligible(:,:)
            integer  :: i_loc, n_loc, ibest, istate_best, r_loc, proj_r, k, ithr, nthr_team
            real(dp) :: prior_w_ins
            logical  :: l_prior_ins
            l_prior_ins = self%p_ptr%l_centring
            prior_w_ins = 0.d0
            if( l_prior_ins ) prior_w_ins = map_prior_weight(self%p_ptr%shift_sigma, self%p_ptr%smpd)
            nthr_team = omp_get_max_threads() ! the per-thread rows and masks are sized by the team the region runs
            allocate(kpos(kmax,i_first:i_last), nkeep(i_first:i_last), nelig(i_first:i_last), source=0)
            allocate(kw(kmax,i_first:i_last), source=0.)
            allocate(pos_c(max(1,maxcand),nthr_team), kinds(kmax,nthr_team), source=0)
            allocate(dists_c(max(1,maxcand),nthr_team), source=0.)
            allocate(eligible(max(1,maxcand),nthr_team), source=.false.)
            !$omp parallel do default(shared) proc_bind(close) schedule(static) &
            !$omp private(i_loc,n_loc,ibest,istate_best,r_loc,proj_r,k,ithr)
            do i_loc = i_first,i_last
                nkeep(i_loc) = 1
                kw(1,i_loc)  = 1.
                if( .not. l_refined(i_loc) ) cycle
                ithr = omp_get_thread_num() + 1
                call gather_state_row(i_loc, dists_c(:,ithr), pos_c(:,ithr), n_loc, ibest, istate_best)
                if( istate_best == 0 ) cycle
                do r_loc = 1,n_loc
                    eligible(r_loc,ithr) = self%candidate_store%candidates(pos_c(r_loc,ithr))%has_sh
                    if( .not. eligible(r_loc,ithr) ) cycle
                    if( l_prior_ins ) dists_c(r_loc,ithr) = insert_pose_loss(dists_c(r_loc,ithr),&
                        &real([self%candidate_store%candidates(pos_c(r_loc,ithr))%x,&
                        &self%candidate_store%candidates(pos_c(r_loc,ithr))%y], dp), prior_w_ins, self%temps(i_loc))
                    proj_r = modulo(self%candidate_store%candidates(pos_c(r_loc,ithr))%iref - 1, self%p_ptr%nspace) + 1
                    eligible(r_loc,ithr) = self%insert_cone%within(istate_best, self%assgn_map(i_loc)%iproj, &
                        &self%assgn_map(i_loc)%inpl, [self%assgn_map(i_loc)%x, self%assgn_map(i_loc)%y], proj_r, &
                        &self%candidate_store%candidates(pos_c(r_loc,ithr))%inpl, &
                        &[self%candidate_store%candidates(pos_c(r_loc,ithr))%x, self%candidate_store%candidates(pos_c(r_loc,ithr))%y])
                enddo
                call top_posterior_extras(n_loc, dists_c(1:n_loc,ithr), temperature_or_argmin(self%temps(i_loc)), &
                    &eligible(1:n_loc,ithr), which_pos(i_loc), kmax, PROB_INSERT_WMIN, nkeep(i_loc), kinds(:,ithr), &
                    &kw(:,i_loc), nelig(i_loc))
                ! the extras as store positions, which the serial pass reads without gathering the row again
                do k = 2,nkeep(i_loc)
                    kpos(k,i_loc) = pos_c(kinds(k,ithr),ithr)
                enddo
            enddo
            !$omp end parallel do
            do i_loc = i_first,i_last
                do k = 2,nkeep(i_loc)
                    cands(k-1) = self%candidate_store%candidates(kpos(k,i_loc))
                enddo
                call self%insert_tab%add(i_loc, kw(1,i_loc), nkeep(i_loc) - 1, cands(1:nkeep(i_loc)-1), &
                    &kw(2:nkeep(i_loc),i_loc), nelig(i_loc))
            enddo
            deallocate(kpos, kinds, nkeep, nelig, pos_c, kw, dists_c, eligible)
        end subroutine select_extras

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

    end subroutine assign_rows_neigh

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
    subroutine record_sparse_eval( self, i, ithr, ri, dist, irot, x, y, has_sh, sdev )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i, ithr, ri, irot
        real,                      intent(in)    :: dist, x, y
        logical,                   intent(in)    :: has_sh
        real,            optional, intent(in)    :: sdev
        type(prob_candidate) :: candidate
        candidate%iref   = self%ref_full(ri)
        candidate%dist   = dist
        candidate%inpl   = irot
        candidate%x      = x
        candidate%y      = y
        candidate%has_sh = has_sh
        if( present(sdev) ) candidate%sdev = sdev
        call self%candidate_buffers(ithr)%append_or_replace(i, candidate)
    end subroutine record_sparse_eval

    ! Minimizes shift at the previous best orientation
    !> the cost record of one neighbourhood fill (bench file): the wall time of the parallel region, thread-seconds of
    !! the seed, the score loop and the refinement, references scored, shifted particle memos built, and the search
    !! objects' counters
    subroutine record_neigh_fill_profile( self, grad_obj, wall, ts_seed, ts_score, ts_refine, nscored, nbank, nfilled,&
            &grad_disk )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad),    intent(inout) :: grad_obj(:)
        type(pftc_shsrch_grad), optional, intent(inout) :: grad_disk(:)
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
        if( .not. present(grad_disk) ) return
        if( .not. self%p_ptr%l_centring ) return
        do ithr = 1,size(grad_disk)
            call grad_disk(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
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
    !> the coarse maps of a fill (shift_sigma > 0): the disk of three sigma (disk_radius_px) and the window of trs,
    !! with each map's radial step, the box of the polish about a cell
    subroutine setup_shiftmap_neigh( self, prior_w, kscreen, step_win, step_disk )
        class(eul_prob_tab_neigh), intent(inout) :: self
        real(dp),                  intent(out)   :: prior_w
        integer,                   intent(out)   :: kscreen
        real,                      intent(out)   :: step_win, step_disk
        real(dp) :: sigma_px
        if( .not. (self%p_ptr%prob_beta > 0.) ) THROW_HARD('the shift map draws at prob_beta x T: prob_beta must be positive')
        sigma_px = real(self%p_ptr%shift_sigma, dp) / real(self%p_ptr%smpd, dp)
        call self%b_ptr%pftc%setup_shiftmaps(sigma_px, real(disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd,&
            &self%p_ptr%box), dp), real(self%p_ptr%trs, dp), prior_w, kscreen)
        step_win  = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 2, 1) -&
            &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 1, 1)))
        step_disk = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 2, 1) -&
            &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 1, 1)))
    end subroutine setup_shiftmap_neigh

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
    !> the shifts of the best candidates per state: the windowed refinement from the seed, or under the key (l_map)
    !! the polish in a box about a coarse map's optimum (kept if it beat the origin, or the cell if it beats the
    !! candidate's row): the window map at the best candidates; for a particle on the disk (l_disk) the screen of the
    !! state's candidates holding SCREEN_QUANTILE of its posterior at the origin with the disk map under the prior,
    !! the best of the screen polished in the disk's search object
    subroutine refine_best_neighbors( self, ew, grad_obj, l_seed_sh_first, max_refs_to_refine,&
        &i, ithr, iptcl, shift_seed, neval, l_map, kscreen, step_win, step_disk, l_disk, l_prior, grad_disk, prior_w)
        class(eul_prob_tab_neigh),      intent(inout) :: self
        type(eval_ws),                  intent(inout) :: ew
        type(pftc_shsrch_grad), target, intent(inout) :: grad_obj(:), grad_disk(:)
        logical,                        intent(in)    :: l_seed_sh_first
        integer,                        intent(in)    :: max_refs_to_refine, i, ithr, iptcl, neval
        real,                           intent(in)    :: shift_seed(3)
        logical,                        intent(in)    :: l_map, l_disk, l_prior
        integer,                        intent(in)    :: kscreen
        real,                           intent(in)    :: step_win, step_disk
        real(dp),                       intent(in)    :: prior_w
        type(pftc_shsrch_grad), pointer :: gso
        integer  :: nrefs_to_refine, nstate_eval, j, k, eval_slot, ri, istate, iproj, irot, irho, ialpha, imap, npool
        integer  :: si, state_eval, pool(neval), scr_irot(neval), loc(1)
        real     :: refined_shift(3), xy_start(2), loss_cell, dist, sdev, step, dists_sorted(neval)
        real     :: scr_xy(2,neval), scr_loss(neval), scr_sdev(neval), width_best(max(1,max_refs_to_refine))
        real(dp) :: temp_s, temp_i, prior_w_map
        logical  :: l_polished
        ew%best_evals(:,ithr) = 0
        if( neval <= 0 ) return
        gso         => grad_obj(ithr)
        imap        = SHIFTMAP_WINDOW
        prior_w_map = 0.d0
        step        = step_win
        if( l_map )then
            temp_s = real(self%p_ptr%prob_beta, dp) * real(posterior_temperature(&
                &self%b_ptr%pftc%wsqsum_ptcl_band(iptcl, kscreen), self%b_ptr%pftc%get_pftsz()), dp)
            temp_i = real(self%p_ptr%prob_beta, dp) * real(self%temps(i), dp)
            if( l_disk )then
                imap        = SHIFTMAP_SCREEN
                if( l_prior ) prior_w_map = prior_w
                step        = step_disk
                gso         => grad_disk(ithr)
            endif
        endif
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
            if( l_map .and. l_disk )then
                ! the screen over the state's candidates (their slots), then the best of the screen
                pool(1:nstate_eval)         = minnloc(ew%state_eval_dists(1:neval,ithr), nstate_eval)
                dists_sorted(1:nstate_eval) = ew%state_eval_dists(pool(1:nstate_eval),ithr)
                npool = min(nstate_eval, max(2 * nrefs_to_refine,&
                    &screen_pool_size(dists_sorted(1:nstate_eval), nstate_eval, temp_i)))
                scr_loss = huge(1.0)
                width_best = huge(1.)
                do k = 1, npool
                    eval_slot = pool(k)
                    if( eval_slot < 1 ) cycle
                    ri = ew%evaluated_ref_ids(eval_slot,ithr)
                    if( ri < 1 ) cycle
                    call self%b_ptr%pftc%gen_shiftmap(imap, self%ref_full(ri), iptcl)
                    call self%b_ptr%pftc%shiftmap_mode(imap, temp_s, prior_w_map, &
                        &irot, irho, ialpha, loss_cell, &
                        &temp_phys=temp_s / real(self%p_ptr%prob_beta, dp), sdev=scr_sdev(eval_slot), &
                        &width_cutoff=maxval(width_best(1:nrefs_to_refine)))
                    call update_screen_best(loss_cell, width_best(1:nrefs_to_refine))
                    scr_irot(eval_slot) = irot
                    scr_xy(:,eval_slot) = real(self%b_ptr%pftc%shiftmap_shift(imap, irho, ialpha))
                    scr_loss(eval_slot) = loss_cell
                enddo
                do j = 1, nrefs_to_refine
                    loc       = minloc(scr_loss)
                    eval_slot = loc(1)
                    if( scr_loss(eval_slot) >= huge(1.0)/2. ) exit
                    ri = ew%evaluated_ref_ids(eval_slot,ithr)
                    call gso%set_indices(self%ref_full(ri), iptcl)
                    irot = scr_irot(eval_slot)
                    call gso%minimize_in_cell(irot, scr_xy(:,eval_slot), step, refined_shift, l_polished)
                    dist = eulprob_dist_switch(refined_shift(1), self%p_ptr%cc_objfun)
                    if( l_polished .or. dist < ew%evaluated_ref_dists(eval_slot,ithr) )then
                        call record_sparse_eval(self, i, ithr, ri, dist, irot, refined_shift(2), refined_shift(3), .true.,&
                            &scr_sdev(eval_slot))
                    endif
                    scr_loss(eval_slot) = huge(1.0)
                enddo
                cycle
            endif
            ew%best_evals(1:nrefs_to_refine,ithr) = &
                &minnloc(ew%state_eval_dists(1:neval,ithr), nrefs_to_refine)
            do j = 1, nrefs_to_refine
                eval_slot = ew%best_evals(j,ithr)
                if( eval_slot < 1 ) cycle
                ri = ew%evaluated_ref_ids(eval_slot,ithr)
                if( ri < 1 ) cycle
                istate = self%ref_state(ri)
                iproj  = self%ref_proj(ri)
                call gso%set_indices(self%ref_full(ri), iptcl)
                irot = self%candidate_buffers(ithr)%get_inpl(i, self%ref_full(ri))
                if( l_map )then
                    call self%b_ptr%pftc%gen_shiftmap(SHIFTMAP_WINDOW, self%ref_full(ri), iptcl)
                    call self%b_ptr%pftc%shiftmap_mode(SHIFTMAP_WINDOW, temp_s, 0.d0, &
                        &irot, irho, ialpha, loss_cell, &
                        &temp_phys=temp_s / real(self%p_ptr%prob_beta, dp), sdev=sdev)
                    xy_start = real(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, irho, ialpha))
                    call gso%minimize_in_cell(irot, xy_start, step_win, refined_shift, l_polished)
                    dist = eulprob_dist_switch(refined_shift(1), self%p_ptr%cc_objfun)
                    if( l_polished .or. dist < ew%evaluated_ref_dists(eval_slot,ithr) )then
                        call record_sparse_eval(self, i, ithr, ri, dist, irot, refined_shift(2), refined_shift(3), .true., sdev)
                    endif
                    cycle
                endif
                if( l_seed_sh_first )then
                    ! the particle-frame seed as the shift of this reference at the rotation it was scored at, in the box
                    xy_start = shift_seed(2:3)
                    if( irot > 0 ) xy_start = shift_ptcl2ref(shift_seed(2:3), self%b_ptr%pftc%get_rot(irot))
                    xy_start = max(-self%p_ptr%trs, min(self%p_ptr%trs, xy_start))
                    refined_shift = gso%minimize(irot=irot, sh_rot=.true., xy_in=xy_start)
                else
                    refined_shift = gso%minimize(irot=irot, sh_rot=.true.)
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
        call self%candidate_store%kill
        if( allocated(self%candidate_fill_counts) ) deallocate(self%candidate_fill_counts)
        self%l_direct_stoch_neigh = .false.
        self%nkept = 0
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
