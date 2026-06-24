!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is geometrically sparse and searched independently per state.
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_utils, only: angle_sampling, build_pind_lookup, calc_athres, eulprob_dist_switch,&
    &materialize_seed_shift, read_seed_shift_table, write_seed_shift_table
use simple_decay_funs,        only: extremal_decay
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_ori,                only: ori
use simple_type_defs,          only: OBJFUN_EUCLID
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    integer, allocatable    :: eval_touched_refs(:,:)
    integer, allocatable    :: eval_touched_counts(:)
    integer                 :: eval_max_touched     = 0
    logical                 :: l_direct_stoch_neigh = .false.
contains
    procedure :: new_neigh
    procedure :: fill_tab      => fill_tab_neigh
    procedure :: fill_tab_range => fill_tab_neigh_range
    procedure :: ref_normalize => ref_normalize_neigh
    procedure :: ref_assign    => ref_assign_neigh
    procedure :: kill          => kill_neigh
    procedure :: write_tab     => write_tab_neigh
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

! Workspace for coarse subspaces used by the geom/state/sum path
type :: coarse_search_ws
    real,    allocatable :: peak_subspace_dists(:,:,:)  ! [npeak,nstates,nthr]
    integer, allocatable :: peak_subspace_inds(:,:,:)   ! [npeak,nstates,nthr]
    integer, allocatable :: peak_subspace_count(:,:)    ! [nstates,nthr]
    integer, allocatable :: pooled_sub_inds(:,:,:)      ! [npeak,nstates,nthr]
    integer, allocatable :: pooled_sub_count(:,:)       ! [nstates,nthr]
contains
    procedure :: alloc_coarse_ws
    procedure :: dealloc_coarse_ws
end type coarse_search_ws

contains

    subroutine new_neigh( self, params, build, pinds, empty_okay )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer :: nsubs
        logical :: l_empty_okay
        call self%kill
        select case(trim(params%prob_neigh_mode))
            case('shc','snhc')
                self%l_direct_stoch_neigh = .true.
                call self%eul_prob_tab%new(params, build, pinds, .true.)
            case('geom','state','sum')
                self%l_direct_stoch_neigh = .false.
                l_empty_okay = .false.
                if( present(empty_okay) ) l_empty_okay = empty_okay
                call self%eul_prob_tab%new(params, build, pinds, l_empty_okay)
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
        ! Sparse support now spans neighborhoods across all active states.
        self%eval_max_touched = max(1, self%nrefs)
        allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        allocate(self%eval_touched_counts(self%nptcls),                     source=0)
    end subroutine new_neigh

    ! Fills the full particle range by delegating to fill_tab_neigh_range(1, nptcls).
    subroutine fill_tab_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%fill_tab_range(1, self%nptcls)
    end subroutine fill_tab_neigh

    ! Dispatches to the stochastic (shc/snhc) or subspace (geom/state/sum) path
    subroutine fill_tab_neigh_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        if( self%eval_max_touched < 1 ) self%eval_max_touched = 1
        if( .not. allocated(self%eval_touched_refs) )&
            &allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        if( .not. allocated(self%eval_touched_counts) )&
            &allocate(self%eval_touched_counts(self%nptcls), source=0)
        ! Search paths
        select case(trim(self%p_ptr%prob_neigh_mode))
            case('shc','snhc')
                call fill_tab_stoch_range(self, i_first, i_last)
            case('geom','state','sum')
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
        type(pftc_shsrch_grad)      :: grad_shsrch_obj(nthr_glob)
        type(eval_ws)               :: eval_work
        type(ran_tabu), allocatable :: direct_rts(:)
        integer,        allocatable :: inds_sorted(:,:), direct_srch_order(:,:)
        real,           allocatable :: inpl_athres(:), dists_inpl(:,:), dists_inpl_sorted(:,:), corrs_inpl(:,:)
        integer :: i, istate, ithr, max_refs_to_refine, nfull_refs
        integer :: iptcl, si, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_prob_objfun, l_shc_neigh, l_snhc_neigh, l_seed_sh_first
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_prob_objfun   = (self%p_ptr%cc_objfun == OBJFUN_EUCLID)
        l_shc_neigh     = (trim(self%p_ptr%prob_neigh_mode) == 'shc')
        l_snhc_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'snhc')
        l_seed_sh_first = self%p_ptr%l_doshift .and. l_shc_neigh
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        nfull_refs = self%p_ptr%nstates * self%p_ptr%nspace
        call clear_sparse_eval_table_range(self, i_from, i_to)
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(nrots,nthr_glob), dists_inpl_sorted(nrots,nthr_glob),&
            &corrs_inpl(nrots,nthr_glob), inds_sorted(nrots,nthr_glob))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_stoch_ws(self, max_refs_to_refine)
        do si = 1, self%nstates
            istate = self%ssinds(si)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        allocate(direct_srch_order(nfull_refs,nthr_glob), direct_rts(nthr_glob))
        do ithr = 1,nthr_glob
            direct_rts(ithr) = ran_tabu(nfull_refs)
        enddo
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed)&
            !$omp proc_bind(close) schedule(static)
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
            !$omp proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                shift_seed = 0.
                call process_particle_stoch(i, iptcl, ithr, .false., shift_seed)
            enddo
            !$omp end parallel do
        endif
        ! cleanup
        call eval_work%dealloc_eval_ws
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        deallocate(direct_srch_order,direct_rts)
        deallocate(inds_sorted, dists_inpl_sorted, dists_inpl, corrs_inpl, inpl_athres)

    contains

        ! One particle evaluation in the stochastic path
        subroutine process_particle_stoch(i_loc, iptcl_loc, ithr_loc, l_with_shift, shift_seed)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: shift_seed(3)
            type(ori) :: o_prev
            integer   :: prev_state_loc, prev_proj_loc, neval_loc
            call get_particle_context(self, iptcl_loc, o_prev, prev_state_loc, prev_proj_loc)
            call estimate_shift_seed(self, grad_shsrch_obj, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc, o_prev,&
                &l_seed_sh_first, shift_seed)
            call evaluate_direct_stochastic_refs(i_loc, ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc,&
                &shift_seed, l_with_shift, neval_loc)
            if( l_with_shift )then
                call refine_best_neighbors(self, eval_work, grad_shsrch_obj, l_seed_sh_first, max_refs_to_refine,&
                &i_loc, ithr_loc, iptcl_loc, shift_seed, neval_loc)
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
            integer :: nrefs_bound, smpl_ninpl, nrots
            real    :: dist_loc, corr_loc, prev_corr_loc, neigh_frac, power_loc
            logical :: l_greedy_first
            neval_loc = 0
            call direct_rts(ithr_loc)%ne_ran_iarr(direct_srch_order(:,ithr_loc))
            prev_full_ref_loc = 0
            if( prev_state_loc >= 1 .and. prev_state_loc <= self%p_ptr%nstates .and.&
                &prev_proj_loc >= 1 .and. prev_proj_loc <= self%p_ptr%nspace )then
                prev_full_ref_loc = (prev_state_loc - 1) * self%p_ptr%nspace + prev_proj_loc
                call put_last(prev_full_ref_loc, direct_srch_order(:,ithr_loc))
            endif
            l_greedy_first = l_shc_neigh .and. self%b_ptr%spproj_field%is_first_update(self%p_ptr%which_iter, iptcl_loc)
            nrots          = self%b_ptr%pftc%get_nrots()
            nrefs_bound    = nfull_refs
            smpl_ninpl     = nrots
            power_loc      = POST_EXTR_POWER
            if( l_snhc_neigh )then
                neigh_frac  = extremal_decay(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
                nrefs_bound = max(2, min(nfull_refs, nint(real(nfull_refs) * (1. - neigh_frac))))
                smpl_ninpl  = neighfrac2nsmpl(neigh_frac, nrots)
                power_loc   = merge(EXTR_POWER, POST_EXTR_POWER, self%p_ptr%extr_iter <= self%p_ptr%extr_lim)
            endif
            prev_corr_loc = -huge(1.0)
            if( l_shc_neigh .and. (.not. l_greedy_first) )then
                call calc_previous_corr(ithr_loc, iptcl_loc, prev_full_ref_loc, prev_corr_loc)
            endif
            do isample = 1,nfull_refs
                if( l_snhc_neigh .and. isample > nrefs_bound ) exit
                full_ref_loc = direct_srch_order(isample,ithr_loc)
                if( full_ref_loc < 1 .or. full_ref_loc > size(eval_work%fullref_to_sparse_ref) ) cycle
                ri_loc = eval_work%fullref_to_sparse_ref(full_ref_loc)
                if( ri_loc < 1 ) cycle
                call score_direct_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed, l_with_shift,&
                    &l_snhc_neigh, power_loc, smpl_ninpl, dist_loc, corr_loc, irot_loc)
                call record_sparse_eval(self, i_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                neval_loc = neval_loc + 1
                eval_work%evaluated_ref_ids(neval_loc,ithr_loc)   = ri_loc
                eval_work%evaluated_ref_dists(neval_loc,ithr_loc) = dist_loc
                if( l_shc_neigh .and. (.not. l_greedy_first) .and. corr_loc > prev_corr_loc ) exit
            enddo
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

        ! Scores a single reference selecting
        subroutine score_direct_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift,&
            &l_snhc_style, power_loc, smpl_ninpl_loc, dist_loc, corr_loc, irot_loc)
            integer, intent(in)  :: full_ref_loc, ithr_loc, iptcl_loc, smpl_ninpl_loc
            real,    intent(in)  :: shift_seed_loc(3), power_loc
            logical, intent(in)  :: l_with_shift, l_snhc_style
            real,    intent(out) :: dist_loc, corr_loc
            integer, intent(out) :: irot_loc
            integer :: istate_loc, order_ind_loc
            real    :: sh_loc(2)
            sh_loc = 0.
            if( l_with_shift .and. l_seed_sh_first ) sh_loc = shift_seed_loc(2:3)
            if( l_snhc_style )then
                if( l_prob_objfun )then
                    call self%b_ptr%pftc%gen_prob_power_objfun_val(full_ref_loc, iptcl_loc, sh_loc, power_loc,&
                        &smpl_ninpl_loc, dist_loc, corr_loc, irot_loc, dists_inpl(:,ithr_loc), inds_sorted(:,ithr_loc))
                else
                    call self%b_ptr%pftc%gen_objfun_vals(full_ref_loc, iptcl_loc, sh_loc, corrs_inpl(:,ithr_loc))
                    call power_sampling(power_loc, self%b_ptr%pftc%get_nrots(), corrs_inpl(:,ithr_loc),&
                        &inds_sorted(:,ithr_loc), smpl_ninpl_loc, irot_loc, order_ind_loc, corr_loc)
                    if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                    dist_loc = eulprob_dist_switch(corr_loc, self%p_ptr%cc_objfun)
                endif
            else
                istate_loc = (full_ref_loc - 1) / self%p_ptr%nspace + 1
                if( l_prob_objfun )then
                    call self%b_ptr%pftc%gen_prob_objfun_val(full_ref_loc, iptcl_loc, sh_loc,&
                        &inpl_athres(istate_loc), self%p_ptr%prob_athres, dist_loc, irot_loc,&
                        &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
                    corr_loc = exp(-dist_loc)
                else
                    call self%b_ptr%pftc%gen_objfun_vals(full_ref_loc, iptcl_loc, sh_loc, corrs_inpl(:,ithr_loc))
                    dists_inpl(:,ithr_loc) = eulprob_dist_switch(corrs_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                    irot_loc = angle_sampling(dists_inpl(:,ithr_loc), dists_inpl_sorted(:,ithr_loc),&
                        &inds_sorted(:,ithr_loc), inpl_athres(istate_loc), self%p_ptr%prob_athres)
                    if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                    dist_loc = dists_inpl(irot_loc,ithr_loc)
                    corr_loc = corrs_inpl(irot_loc,ithr_loc)
                endif
            endif
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
        integer, allocatable   :: inds_sorted(:,:)
        real,    allocatable   :: inpl_athres(:), dists_inpl(:,:), dists_inpl_sorted(:,:)
        integer :: i, istate, ithr, max_refs_to_refine, nsubs, npeak_target
        integer :: iptcl, si, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_prob_objfun, l_geom_neigh, l_state_neigh, l_sum_neigh, l_seed_sh_first
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_prob_objfun   = (self%p_ptr%cc_objfun == OBJFUN_EUCLID)
        l_geom_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'geom')
        l_state_neigh   = (trim(self%p_ptr%prob_neigh_mode) == 'state')
        l_sum_neigh     = (trim(self%p_ptr%prob_neigh_mode) == 'sum')
        l_seed_sh_first = self%p_ptr%l_doshift
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        nsubs        = size(self%b_ptr%subspace_inds)
        npeak_target = min(max(1, self%p_ptr%npeaks), nsubs)
        call clear_sparse_eval_table_range(self, i_from, i_to)
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(nrots,nthr_glob), dists_inpl_sorted(nrots,nthr_glob), inds_sorted(nrots,nthr_glob))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_subspace_ws(self, nsubs, max_refs_to_refine)
        do si = 1, self%nstates
            istate = self%ssinds(si)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        call coarse_ws%alloc_coarse_ws(npeak_target, self%p_ptr%nstates, nthr_glob)
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_subspace(i, iptcl, ithr, .true., shift_seed)
                self%seed_shifts(:,i) = shift_seed(2:3)
                self%seed_has_sh(i)   = l_seed_sh_first
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl      = self%pinds(i)
                ithr       = omp_get_thread_num() + 1
                shift_seed = 0.
                call process_particle_subspace(i, iptcl, ithr, .false., shift_seed)
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call coarse_ws%dealloc_coarse_ws
        call eval_work%dealloc_eval_ws
        deallocate(inds_sorted, dists_inpl_sorted, dists_inpl, inpl_athres)

    contains

        ! One particle evaluation
        subroutine process_particle_subspace(i_loc, iptcl_loc, ithr_loc, l_with_shift, shift_seed_loc)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: shift_seed_loc(3)
            type(ori) :: o_prev
            integer   :: prev_state, prev_proj, neval
            call get_particle_context(self, iptcl_loc, o_prev, prev_state, prev_proj)
            call estimate_shift_seed(self, grad_shsrch_obj, ithr_loc, iptcl_loc, prev_state, prev_proj, o_prev,&
                &l_seed_sh_first, shift_seed_loc)
            if( l_geom_neigh )then
                call build_geometric_neighborhood(ithr_loc, prev_proj)
            else
                if( l_sum_neigh )   call find_sum_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
                if( l_state_neigh ) call find_peak_subspaces    (i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
                call build_pooled_neighborhood(ithr_loc, prev_proj)
            endif
            call evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval)
            if( l_with_shift )then
                call refine_best_neighbors(self, eval_work, grad_shsrch_obj, l_seed_sh_first, max_refs_to_refine,&
                    &i_loc, ithr_loc, iptcl_loc, shift_seed_loc, neval)
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
                    call consider_peak_subspace(isub_loc, istate_loc, ithr_loc, dist)
                    ri_loc = eval_work%fullref_to_sparse_ref(full_ref_subspace)
                    if( ri_loc > 0 )&
                        &call record_sparse_eval(self, i_loc, ri_loc, dist, irot_loc, 0., 0., .false.)
                enddo
            enddo
        end subroutine find_peak_subspaces

        ! Scores all subspace representatives with distance summed across all active states,
        ! keep the top-npeak_target and propagates to every state; prob_neigh_mode=sum
        subroutine find_sum_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            integer :: si_loc, istate_loc, isub_loc, full_ref_subspace, irot_loc, ri_loc, coarse_proj
            integer :: nvalid_loc, anchor_state, expected_states
            real    :: dist_loc, sum_dist_loc
            coarse_ws%peak_subspace_dists(:,:,ithr_loc) = huge(1.0)
            coarse_ws%peak_subspace_inds(:,:,ithr_loc)  = 0
            coarse_ws%peak_subspace_count(:,ithr_loc)   = 0
            anchor_state    = 0
            expected_states = 0
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                if( .not. self%state_exists(istate_loc) ) cycle
                expected_states = expected_states + 1
                if( anchor_state == 0 ) anchor_state = istate_loc
            enddo
            if( expected_states < 1 ) return
            do isub_loc = 1, nsubs
                coarse_proj = self%b_ptr%subspace_inds(isub_loc)
                sum_dist_loc = 0.
                nvalid_loc   = 0
                do si_loc = 1, self%nstates
                    istate_loc = self%ssinds(si_loc)
                    if( .not. self%state_exists(istate_loc) ) cycle
                    full_ref_subspace = (istate_loc-1)*self%p_ptr%nspace + coarse_proj
                    call score_subspace_ref(full_ref_subspace, ithr_loc, iptcl_loc, shift_seed_loc,&
                        &l_with_shift, dist_loc, irot_loc)
                    sum_dist_loc = sum_dist_loc + dist_loc
                    nvalid_loc = nvalid_loc + 1
                    ri_loc = eval_work%fullref_to_sparse_ref(full_ref_subspace)
                    if( ri_loc > 0 )&
                        &call record_sparse_eval(self, i_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                enddo
                if( nvalid_loc == expected_states )&
                    &call consider_peak_subspace(isub_loc, anchor_state, ithr_loc, sum_dist_loc)
            enddo
            call copy_anchor_peak_subspaces(anchor_state, ithr_loc)
        end subroutine find_sum_peak_subspaces

        ! Copies the peak subspace selection from the anchor state to all other active states
        ! prob_neigh_mode=sum
        subroutine copy_anchor_peak_subspaces(anchor_state, ithr_loc)
            integer, intent(in) :: anchor_state, ithr_loc
            integer :: si_loc, istate_loc, npeak_found_loc
            npeak_found_loc = coarse_ws%peak_subspace_count(anchor_state,ithr_loc)
            do si_loc = 1,self%nstates
                istate_loc = self%ssinds(si_loc)
                if( istate_loc == anchor_state ) cycle
                coarse_ws%peak_subspace_count(istate_loc,ithr_loc) = npeak_found_loc
                if( npeak_found_loc > 0 )then
                    coarse_ws%peak_subspace_dists(1:npeak_found_loc,istate_loc,ithr_loc) =&
                        &coarse_ws%peak_subspace_dists(1:npeak_found_loc,anchor_state,ithr_loc)
                    coarse_ws%peak_subspace_inds(1:npeak_found_loc,istate_loc,ithr_loc) =&
                        &coarse_ws%peak_subspace_inds(1:npeak_found_loc,anchor_state,ithr_loc)
                endif
            enddo
        end subroutine copy_anchor_peak_subspaces

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

        ! Evaluates a coarse subspace representative; prob_neigh_mode=state|sum
        subroutine score_subspace_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed_loc,&
            &l_with_shift, dist_loc, irot_loc)
            integer, intent(in)  :: full_ref_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            if( l_prob_objfun )then
                if( l_with_shift )then
                    call self%b_ptr%pftc%gen_best_objfun_val(full_ref_loc, iptcl_loc,&
                        &shift_seed_loc(2:3), dist_loc, irot_loc)
                else
                    call self%b_ptr%pftc%gen_best_objfun_val(full_ref_loc, iptcl_loc,&
                        &[0.,0.], dist_loc, irot_loc)
                endif
            else
                if( l_with_shift )then
                    call self%b_ptr%pftc%gen_objfun_vals(full_ref_loc, iptcl_loc,&
                        &shift_seed_loc(2:3), dists_inpl(:,ithr_loc))
                else
                    call self%b_ptr%pftc%gen_objfun_vals(full_ref_loc, iptcl_loc,&
                        &[0.,0.], dists_inpl(:,ithr_loc))
                endif
                dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                irot_loc = minloc(dists_inpl(:,ithr_loc), dim=1)
                dist_loc = dists_inpl(irot_loc,ithr_loc)
            endif
        end subroutine score_subspace_ref

        ! Maintains a sorted top-npeak_target list of subspace indices per state; prob_neigh_mode=state|sum
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

        ! Merges per-state peak subspace lists into the pooled set, always including the
        ! subspace of the previous best orientation
        subroutine build_pooled_neighborhood(ithr_loc, prev_proj_loc)
            integer, intent(in) :: ithr_loc, prev_proj_loc
            integer :: si_loc, istate_loc, npeak_found
            integer :: prev_isub_loc, k_loc, slot
            logical :: already_in_pool
            coarse_ws%pooled_sub_count(:,ithr_loc) = 0
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                npeak_found = coarse_ws%peak_subspace_count(istate_loc,ithr_loc)
                if( npeak_found > 0 )then
                    coarse_ws%pooled_sub_inds(1:npeak_found,istate_loc,ithr_loc) =&
                        &coarse_ws%peak_subspace_inds(1:npeak_found,istate_loc,ithr_loc)
                    coarse_ws%pooled_sub_count(istate_loc,ithr_loc) = npeak_found
                endif
                ! Subspace corresponding to the previous best orientation
                prev_isub_loc = self%b_ptr%subspace_full2sub_map(max(1, min(self%p_ptr%nspace, prev_proj_loc)))
                if( prev_isub_loc < 1 .or. prev_isub_loc > nsubs ) prev_isub_loc = 1
                ! Avoid duplicate evaluation
                already_in_pool = .false.
                do k_loc = 1, coarse_ws%pooled_sub_count(istate_loc,ithr_loc)
                    if( coarse_ws%pooled_sub_inds(k_loc,istate_loc,ithr_loc) == prev_isub_loc )then
                        already_in_pool = .true.
                        exit
                    endif
                enddo
                if( already_in_pool ) cycle
                ! Make the previous-orientation subspace part of the npeak_target peaks.
                ! Append or evict the worst (last) entry.
                if( coarse_ws%pooled_sub_count(istate_loc,ithr_loc) < npeak_target )then
                    slot = coarse_ws%pooled_sub_count(istate_loc,ithr_loc) + 1
                    coarse_ws%pooled_sub_count(istate_loc,ithr_loc) = slot
                else
                    slot = npeak_target
                endif
                coarse_ws%pooled_sub_inds(slot,istate_loc,ithr_loc) = prev_isub_loc
            enddo
        end subroutine build_pooled_neighborhood

        ! Pass over fine references of the pooled neighborhoods
        subroutine evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval)
            integer, intent(in)  :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            integer, intent(out) :: neval
            integer :: jsub, kref_loc, nrefs_sub_loc, offset_loc
            integer :: si_loc, ri_loc, istate_loc, iproj_loc, irot_loc, iref_loc, isub_loc
            real    :: dist
            neval = 0
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
                        iproj_loc = self%jinds(ri_loc)
                        iref_loc  = (istate_loc-1)*self%p_ptr%nspace + iproj_loc
                        if( l_prob_objfun )then
                            if( l_with_shift )then
                                call self%b_ptr%pftc%gen_prob_objfun_val(iref_loc, iptcl_loc, shift_seed_loc(2:3),&
                                    &inpl_athres(istate_loc), self%p_ptr%prob_athres, dist, irot_loc,&
                                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
                            else
                                call self%b_ptr%pftc%gen_prob_objfun_val(iref_loc, iptcl_loc, [0.,0.],&
                                    &inpl_athres(istate_loc), self%p_ptr%prob_athres, dist, irot_loc,&
                                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
                            endif
                        else
                            if( l_with_shift )then
                                call self%b_ptr%pftc%gen_objfun_vals(iref_loc, iptcl_loc,&
                                    &shift_seed_loc(2:3), dists_inpl(:,ithr_loc))
                            else
                                call self%b_ptr%pftc%gen_objfun_vals(iref_loc, iptcl_loc, [0.,0.], dists_inpl(:,ithr_loc))
                            endif
                            dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                            irot_loc = angle_sampling(dists_inpl(:,ithr_loc), dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc),&
                                &inpl_athres(istate_loc), self%p_ptr%prob_athres)
                            dist = dists_inpl(irot_loc,ithr_loc)
                        endif
                        call record_sparse_eval(self, i_loc, ri_loc, dist, irot_loc, 0., 0., .false.)
                        neval = neval + 1
                        eval_work%evaluated_ref_ids(neval,ithr_loc)   = ri_loc
                        eval_work%evaluated_ref_dists(neval,ithr_loc) = dist
                    enddo
                enddo
            enddo
        end subroutine evaluate_neighborhood

    end subroutine fill_tab_subspace_range

    ! ASSIGNMENT SUBROUTINES

    ! Normalises sparse loc_tab distances
    subroutine ref_normalize_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, iref, k
        logical :: any_eval
        ! normalize only over evaluated refs (inpl > 0)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,k,iref,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = 0.
            do k = 1,self%eval_touched_counts(i)
                iref = self%eval_touched_refs(k,i)
                if( iref < 1 .or. iref > self%nrefs ) cycle
                if( self%loc_tab(iref,i)%inpl > 0 ) sum_dist_all = sum_dist_all + self%loc_tab(iref,i)%dist
            enddo
            if( sum_dist_all < TINY )then
                do k = 1,self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl > 0 ) self%loc_tab(iref,i)%dist = 0.
                enddo
            else
                do k = 1,self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl > 0 ) self%loc_tab(iref,i)%dist = self%loc_tab(iref,i)%dist / sum_dist_all
                enddo
            endif
        enddo
        !$omp end parallel do
        min_dist = huge(1.0)
        max_dist = -huge(1.0)
        any_eval = .false.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,k,iref)&
        !$omp reduction(min:min_dist) reduction(max:max_dist) reduction(.or.:any_eval)
        do i = 1,self%nptcls
            do k = 1,self%eval_touched_counts(i)
                iref = self%eval_touched_refs(k,i)
                if( iref < 1 .or. iref > self%nrefs ) cycle
                if( self%loc_tab(iref,i)%inpl <= 0 ) cycle
                min_dist = min(min_dist, self%loc_tab(iref,i)%dist)
                max_dist = max(max_dist, self%loc_tab(iref,i)%dist)
                any_eval = .true.
            enddo
        enddo
        !$omp end parallel do
        if( .not. any_eval ) return
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab_neigh normalize')
            ! dense-style stochastic tie break over evaluated entries only.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,k,iref)
            do i = 1,self%nptcls
                do k = 1,self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl > 0 ) self%loc_tab(iref,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,k,iref)
            do i = 1,self%nptcls
                do k = 1,self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl > 0 )&
                    &self%loc_tab(iref,i)%dist = (self%loc_tab(iref,i)%dist - min_dist) / (max_dist - min_dist)
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine ref_normalize_neigh

    ! Sparse graph traversal on-the-fly for reference assignment
    subroutine ref_assign_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        ! Sparse global assignment:
        ! 1) ensure every particle has at least one evaluated candidate
        ! 2) normalize sparse scores
        ! 3) build ref->particle and particle->ref adjacency
        ! 4) repeatedly assign the best currently available particle per active ref
        ! 5) fall back to best evaluated sparse ref for leftovers

        type :: assign_graph_ws
            integer, allocatable :: ref_counts(:), ref_offsets(:), ref_fill(:), ref_list(:), ref_pos(:), active_refs(:)
            integer, allocatable :: ptcl_counts(:), ptcl_offsets(:), ptcl_fill(:), ptcl_refs(:)
            real,    allocatable :: ref_dists(:)
        end type assign_graph_ws

        type :: assign_frontier_ws
            integer, allocatable :: inds_sorted(:), order(:), work_ptcl(:), sel_refs(:), sel_pos(:)
            real,    allocatable :: iref_dist(:), dists_sorted(:), work_d(:), sel_dists(:), sel_dists_sorted(:)
            logical, allocatable :: ptcl_avail(:)
        end type assign_frontier_ws

        type(assign_graph_ws)    :: graph
        type(assign_frontier_ws) :: frontier
        real,    allocatable :: dists_inpl(:), dists_inpl_sorted(:)
        integer, allocatable :: inds_sorted(:)
        integer   :: i, iref, assigned_iref, assigned_ptcl, istate, si, fallback_ref
        integer   :: k, idx, nactive, total, m, start, maxref, nleft, assigned_idx, nsel, pos, last_ref
        integer   :: greedy_state(self%nptcls)
        real      :: projs_athres, state_projs_athres(self%p_ptr%nstates)
        real      :: huge_val
        logical   :: l_prob_objfun, final_assigned(self%nptcls), l_filter_greedy_state, l_seed_fallback_shift
        huge_val = huge(1.0)
        l_prob_objfun = (self%p_ptr%cc_objfun == OBJFUN_EUCLID)
        l_seed_fallback_shift = self%p_ptr%l_doshift .and. trim(self%p_ptr%prob_neigh_mode) /= 'snhc'
        allocate(dists_inpl(self%b_ptr%pftc%get_nrots()), dists_inpl_sorted(self%b_ptr%pftc%get_nrots()), inds_sorted(self%b_ptr%pftc%get_nrots()))
        do i = 1, self%nptcls
            ! Seed particle with no evaluated candidate with previous orientation
            call seed_fallback_if_empty(i)
        enddo
        call self%ref_normalize()
        call build_sparse_assignment_graph()
        call assign_particles_globally()
        call assign_remaining_particles_from_best_touched_ref()
        deallocate(inds_sorted, dists_inpl_sorted, dists_inpl)
    contains

        subroutine build_sparse_assignment_graph()
            allocate(graph%ref_counts(self%nrefs), graph%ref_offsets(self%nrefs+1), graph%ref_fill(self%nrefs), graph%ref_pos(self%nrefs),&
            &frontier%iref_dist(self%nrefs), frontier%dists_sorted(self%nrefs), frontier%inds_sorted(self%nrefs), frontier%ptcl_avail(self%nptcls),&
            &graph%ptcl_counts(self%nptcls), graph%ptcl_offsets(self%nptcls+1), graph%ptcl_fill(self%nptcls))
            graph%ref_counts = 0
            graph%ptcl_counts = 0
            do i = 1, self%nptcls
                do k = 1, self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs       ) cycle
                    if( self%loc_tab(iref,i)%inpl <= 0        ) cycle
                    graph%ref_counts(iref) = graph%ref_counts(iref) + 1
                    graph%ptcl_counts(i)   = graph%ptcl_counts(i) + 1
                enddo
            enddo
            nactive = count(graph%ref_counts > 0)
            allocate(graph%active_refs(max(1,nactive)), source=0)
            if( nactive > 0 ) graph%active_refs(1:nactive) = pack((/(iref, iref=1,self%nrefs)/), graph%ref_counts > 0)
            graph%ref_offsets(1) = 1
            do iref = 1, self%nrefs
                graph%ref_offsets(iref+1) = graph%ref_offsets(iref) + graph%ref_counts(iref)
            enddo
            total = graph%ref_offsets(self%nrefs+1) - 1
            allocate(graph%ref_list(max(1,total)), graph%ref_dists(max(1,total)))
            graph%ptcl_offsets(1) = 1
            do i = 1, self%nptcls
                graph%ptcl_offsets(i+1) = graph%ptcl_offsets(i) + graph%ptcl_counts(i)
            enddo
            allocate(graph%ptcl_refs(max(1,graph%ptcl_offsets(self%nptcls+1)-1)))
            graph%ref_fill = graph%ref_offsets(1:self%nrefs)
            graph%ptcl_fill = graph%ptcl_offsets(1:self%nptcls)
            do i = 1, self%nptcls
                do k = 1, self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl <= 0   ) cycle
                    idx = graph%ref_fill(iref)
                    graph%ref_list(idx)  = i
                    graph%ref_dists(idx) = self%loc_tab(iref,i)%dist
                    graph%ref_fill(iref) = idx + 1
                    idx = graph%ptcl_fill(i)
                    graph%ptcl_refs(idx) = iref
                    graph%ptcl_fill(i)   = idx + 1
                enddo
            enddo
            maxref = max(1, maxval(graph%ref_counts))
            allocate(frontier%work_d(maxref), frontier%order(maxref), frontier%work_ptcl(maxref))
            allocate(frontier%sel_refs(max(1,nactive)), frontier%sel_pos(self%nrefs), frontier%sel_dists(max(1,nactive)),&
            &frontier%sel_dists_sorted(max(1,nactive)))
            frontier%sel_pos = 0
            do idx = 1, nactive
                iref = graph%active_refs(idx)
                m    = graph%ref_counts(iref)
                if( m <= 1 ) cycle
                start = graph%ref_offsets(iref)
                frontier%work_d(1:m)    = graph%ref_dists(start:start+m-1)
                frontier%work_ptcl(1:m) = graph%ref_list(start:start+m-1)
                frontier%order(1:m)     = (/(k,k=1,m)/)
                call hpsort(frontier%work_d(1:m), frontier%order(1:m))
                do k = 1, m
                    graph%ref_dists(start+k-1) = frontier%work_d(k)
                    graph%ref_list(start+k-1)  = frontier%work_ptcl(frontier%order(k))
                enddo
            enddo
        end subroutine build_sparse_assignment_graph

        subroutine assign_particles_globally()
            projs_athres = 0.
            state_projs_athres = 0.
            do si = 1, self%nstates
                istate = self%ssinds(si)
                state_projs_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist',&
                    &self%p_ptr%prob_athres, state=istate)
                projs_athres = max(projs_athres, state_projs_athres(istate))
            enddo
            final_assigned = .false.
            if( self%nstates > 1 )then
                call assign_greedy_state_labels()
                do si = 1,self%nstates
                    call assign_particles_for_state(self%ssinds(si))
                enddo
            else
                call assign_particles_single_state()
            endif
            frontier%ptcl_avail = .not. final_assigned
        end subroutine assign_particles_globally

        subroutine init_frontier()
            graph%ref_pos       = 1
            frontier%iref_dist  = huge_val
            frontier%sel_pos    = 0
            nsel                = 0
            do idx = 1, nactive
                iref = graph%active_refs(idx)
                call advance_ref_head(iref)
                call sync_frontier_ref(iref)
            enddo
        end subroutine init_frontier

        subroutine assign_greedy_state_labels()
            greedy_state = 0
            l_filter_greedy_state = .false.
            graph%ref_pos       = 1
            frontier%iref_dist  = huge_val
            frontier%ptcl_avail = .true.
            nleft = self%nptcls
            call init_frontier()
            do while( nleft > 0 )
                if( nsel == 0 ) exit
                assigned_idx  = minloc(frontier%sel_dists(1:nsel), dim=1)
                assigned_iref = frontier%sel_refs(assigned_idx)
                assigned_ptcl = graph%ref_list(graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1)
                greedy_state(assigned_ptcl) = self%sinds(assigned_iref)
                frontier%ptcl_avail(assigned_ptcl) = .false.
                nleft = nleft - 1
                call update_frontier_after_assignment(assigned_ptcl)
            enddo
            do i = 1,self%nptcls
                if( greedy_state(i) == 0 )then
                    fallback_ref = pick_best_evaluated_ref(i)
                    if( fallback_ref > 0 ) greedy_state(i) = self%sinds(fallback_ref)
                endif
            enddo
        end subroutine assign_greedy_state_labels

        subroutine assign_particles_single_state()
            l_filter_greedy_state = .false.
            frontier%ptcl_avail = .true.
            nleft = self%nptcls
            call init_frontier()
            do while( nleft > 0 )
                if( nsel == 0 ) exit
                assigned_idx = angle_sampling(frontier%sel_dists(1:nsel), frontier%sel_dists_sorted(1:nsel),&
                    &frontier%inds_sorted(1:nsel), projs_athres, self%p_ptr%prob_athres)
                call commit_selected_assignment()
            enddo
        end subroutine assign_particles_single_state

        subroutine assign_particles_for_state( state_filter )
            integer, intent(in) :: state_filter
            l_filter_greedy_state = .true.
            frontier%ptcl_avail = (greedy_state == state_filter) .and. (.not. final_assigned)
            nleft = count(frontier%ptcl_avail)
            call init_frontier()
            do while( nleft > 0 )
                if( nsel == 0 ) exit
                assigned_idx = angle_sampling(frontier%sel_dists(1:nsel), frontier%sel_dists_sorted(1:nsel),&
                    &frontier%inds_sorted(1:nsel), state_projs_athres(state_filter), self%p_ptr%prob_athres)
                call commit_selected_assignment()
            enddo
        end subroutine assign_particles_for_state

        subroutine commit_selected_assignment()
            assigned_iref = frontier%sel_refs(assigned_idx)
            assigned_ptcl = graph%ref_list(graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1)
            frontier%ptcl_avail(assigned_ptcl) = .false.
            final_assigned(assigned_ptcl) = .true.
            nleft = nleft - 1
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
            self%assgn_map(assigned_ptcl)%frac = search_frac(assigned_ptcl)
            call update_frontier_after_assignment(assigned_ptcl)
        end subroutine commit_selected_assignment

        subroutine update_frontier_after_assignment( iptcl_assigned )
            integer, intent(in) :: iptcl_assigned
            do idx = graph%ptcl_offsets(iptcl_assigned), graph%ptcl_offsets(iptcl_assigned+1)-1
                iref = graph%ptcl_refs(idx)
                m    = graph%ref_counts(iref)
                if( graph%ref_pos(iref) > m ) cycle
                start = graph%ref_offsets(iref)
                if( graph%ref_list(start + graph%ref_pos(iref) - 1) /= iptcl_assigned ) cycle
                call advance_ref_head(iref)
                call sync_frontier_ref(iref)
            enddo
        end subroutine update_frontier_after_assignment

        subroutine assign_remaining_particles_from_best_touched_ref()
            do i = 1, self%nptcls
                if( .not. frontier%ptcl_avail(i) ) cycle
                if( self%nstates > 1 .and. greedy_state(i) > 0 )then
                    fallback_ref = pick_best_evaluated_ref(i, greedy_state(i))
                else
                    fallback_ref = pick_best_evaluated_ref(i)
                endif
                if( fallback_ref == 0 )then
                    call seed_fallback_if_empty(i)
                    if( self%nstates > 1 .and. greedy_state(i) > 0 )then
                        fallback_ref = pick_best_evaluated_ref(i, greedy_state(i))
                    else
                        fallback_ref = pick_best_evaluated_ref(i)
                    endif
                endif
                if( fallback_ref == 0 ) fallback_ref = pick_best_evaluated_ref(i)
                if( fallback_ref == 0 ) fallback_ref = 1
                self%assgn_map(i) = self%loc_tab(fallback_ref,i)
                call materialize_seed_shift(self%assgn_map(i), self%seed_shifts(:,i),&
                    &self%seed_has_sh(i), self%p_ptr%l_doshift, self%seed_nrots)
                self%assgn_map(i)%frac = search_frac(i)
                final_assigned(i) = .true.
            enddo
        end subroutine assign_remaining_particles_from_best_touched_ref

        real function search_frac(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            if( self%l_direct_stoch_neigh .and. self%nrefs > 0 )then
                search_frac = 100.0 * real(self%eval_touched_counts(iptcl_loc)) / real(self%nrefs)
            else
                search_frac = 100.0
            endif
        end function search_frac

        integer function pick_best_evaluated_ref(iptcl_loc, state_filter) result(ri_best)
            integer, intent(in) :: iptcl_loc
            integer, intent(in), optional :: state_filter
            integer :: kt, ri_loc
            real    :: best_dist
            ri_best = 0
            best_dist = huge(1.0)
            do kt = 1, self%eval_touched_counts(iptcl_loc)
                ri_loc = self%eval_touched_refs(kt,iptcl_loc)
                if( ri_loc < 1 .or. ri_loc > self%nrefs ) cycle
                if( self%loc_tab(ri_loc,iptcl_loc)%inpl <= 0 ) cycle
                if( present(state_filter) )then
                    if( self%sinds(ri_loc) /= state_filter ) cycle
                endif
                if( self%loc_tab(ri_loc,iptcl_loc)%dist < best_dist )then
                    best_dist = self%loc_tab(ri_loc,iptcl_loc)%dist
                    ri_best   = ri_loc
                endif
            enddo
        end function pick_best_evaluated_ref

        subroutine seed_fallback_if_empty(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            type(ori) :: o_prev
            integer :: istate_loc, iproj_loc, irot_loc, fallback_state, fallback_proj, fallback_ref_full
            integer :: ri_loc, fallback_ref_loc
            real    :: inpl_athres_state, sh_seed(2), dist
            if( pick_best_evaluated_ref(iptcl_loc) > 0 ) return
            call self%b_ptr%spproj_field%get_ori(self%pinds(iptcl_loc), o_prev)
            istate_loc = o_prev%get_state()
            if( istate_loc < 1 .or. istate_loc > self%p_ptr%nstates ) istate_loc = 1
            iproj_loc = self%b_ptr%eulspace%find_closest_proj(o_prev)
            iproj_loc = max(1, min(self%p_ptr%nspace, iproj_loc))
            fallback_ref_loc = 0
            do ri_loc = 1, self%nrefs
                if( self%sinds(ri_loc) == istate_loc .and. self%jinds(ri_loc) == iproj_loc )then
                    fallback_ref_loc = ri_loc
                    exit
                endif
            enddo
            if( fallback_ref_loc == 0 )then
                do ri_loc = 1, self%nrefs
                    if( self%sinds(ri_loc) == istate_loc )then
                        fallback_ref_loc = ri_loc
                        exit
                    endif
                enddo
            endif
            if( fallback_ref_loc == 0 ) fallback_ref_loc = 1
            if( self%loc_tab(fallback_ref_loc,iptcl_loc)%inpl <= 0 )then
                fallback_state    = self%sinds(fallback_ref_loc)
                fallback_proj     = self%jinds(fallback_ref_loc)
                fallback_ref_full = (fallback_state-1)*self%p_ptr%nspace + fallback_proj
                sh_seed = 0.
                if( l_seed_fallback_shift ) sh_seed = o_prev%get_2Dshift()
                if( l_prob_objfun )then
                    if( self%p_ptr%l_prob_inpl )then
                        inpl_athres_state = calc_athres(self%b_ptr%spproj_field, 'dist_inpl',&
                            &self%p_ptr%prob_athres, state=fallback_state)
                        call self%b_ptr%pftc%gen_prob_objfun_val(fallback_ref_full, self%pinds(iptcl_loc), sh_seed,&
                            &inpl_athres_state, self%p_ptr%prob_athres, dist, irot_loc, dists_inpl_sorted, inds_sorted)
                    else
                        call self%b_ptr%pftc%gen_best_objfun_val(fallback_ref_full, self%pinds(iptcl_loc), sh_seed,&
                            &dist, irot_loc)
                    endif
                else
                    call self%b_ptr%pftc%gen_objfun_vals(fallback_ref_full, self%pinds(iptcl_loc), sh_seed, dists_inpl)
                    dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                    irot_loc = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                    if( self%p_ptr%l_prob_inpl )then
                        inpl_athres_state = calc_athres(self%b_ptr%spproj_field, 'dist_inpl',&
                            &self%p_ptr%prob_athres, state=fallback_state)
                        irot_loc = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted,&
                            &inpl_athres_state, self%p_ptr%prob_athres)
                    else
                        irot_loc = minloc(dists_inpl, dim=1)
                    endif
                    if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                    dist = dists_inpl(irot_loc)
                endif
                if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                self%loc_tab(fallback_ref_loc,iptcl_loc)%dist   = dist
                self%loc_tab(fallback_ref_loc,iptcl_loc)%inpl   = irot_loc
                self%loc_tab(fallback_ref_loc,iptcl_loc)%x      = sh_seed(1)
                self%loc_tab(fallback_ref_loc,iptcl_loc)%y      = sh_seed(2)
                self%loc_tab(fallback_ref_loc,iptcl_loc)%has_sh = self%p_ptr%l_doshift
            endif
            call mark_ref_touched(self, iptcl_loc, fallback_ref_loc)
            call o_prev%kill
        end subroutine seed_fallback_if_empty

        subroutine advance_ref_head(iref)
            integer, intent(in) :: iref
            integer :: mloc, sloc, cand
            mloc = graph%ref_counts(iref)
            if( graph%ref_pos(iref) > mloc )then
                frontier%iref_dist(iref) = huge_val
                return
            endif
            sloc = graph%ref_offsets(iref)
            do while( graph%ref_pos(iref) <= mloc )
                cand = graph%ref_list(sloc + graph%ref_pos(iref) - 1)
                if( frontier%ptcl_avail(cand) )then
                    if( .not. l_filter_greedy_state )then
                        frontier%iref_dist(iref) = graph%ref_dists(sloc + graph%ref_pos(iref) - 1)
                        return
                    endif
                    if( greedy_state(cand) == self%sinds(iref) )then
                        frontier%iref_dist(iref) = graph%ref_dists(sloc + graph%ref_pos(iref) - 1)
                        return
                    endif
                endif
                graph%ref_pos(iref) = graph%ref_pos(iref) + 1
            enddo
            frontier%iref_dist(iref) = huge_val
        end subroutine advance_ref_head

        subroutine sync_frontier_ref(iref)
            integer, intent(in) :: iref
            if( frontier%sel_pos(iref) > 0 )then
                if( frontier%iref_dist(iref) >= huge_val )then
                    pos = frontier%sel_pos(iref)
                    if( pos < nsel )then
                        last_ref                   = frontier%sel_refs(nsel)
                        frontier%sel_refs(pos)     = last_ref
                        frontier%sel_dists(pos)    = frontier%sel_dists(nsel)
                        frontier%sel_pos(last_ref) = pos
                    endif
                    frontier%sel_pos(iref) = 0
                    nsel = nsel - 1
                else
                    frontier%sel_dists(frontier%sel_pos(iref)) = frontier%iref_dist(iref)
                endif
            else
                if( frontier%iref_dist(iref) < huge_val )then
                    nsel                     = nsel + 1
                    frontier%sel_refs(nsel)  = iref
                    frontier%sel_dists(nsel) = frontier%iref_dist(iref)
                    frontier%sel_pos(iref)   = nsel
                endif
            endif
        end subroutine sync_frontier_ref

    end subroutine ref_assign_neigh

    ! Table & Assignment I/O

    ! Serialises the sparse evaluation table to a binary file
    subroutine write_tab_neigh( self, binfname )
        class(eul_prob_tab_neigh), intent(in) :: self
        class(string),             intent(in) :: binfname
        type(ptcl_ref), allocatable :: sparse_tab(:)
        integer,        allocatable :: sparse_refs(:), sparse_counts(:)
        integer :: funit, addr, io_stat, file_header(3), i, k, ri, nnz, pos
        nnz = 0
        allocate(sparse_counts(self%nptcls), source=0)
        do i = 1,self%nptcls
            do k = 1,self%eval_touched_counts(i)
                ri = self%eval_touched_refs(k,i)
                if( ri < 1 .or. ri > self%nrefs ) cycle
                if( self%loc_tab(ri,i)%inpl <= 0 ) cycle
                sparse_counts(i) = sparse_counts(i) + 1
                nnz = nnz + 1
            enddo
        enddo
        if( nnz < 1 ) THROW_HARD('eul_prob_tab_neigh%write_tab_neigh; empty sparse table')
        allocate(sparse_refs(nnz), sparse_tab(nnz))
        pos = 0
        do i = 1,self%nptcls
            do k = 1,self%eval_touched_counts(i)
                ri = self%eval_touched_refs(k,i)
                if( ri < 1 .or. ri > self%nrefs ) cycle
                if( self%loc_tab(ri,i)%inpl <= 0 ) cycle
                pos = pos + 1
                sparse_refs(pos) = ri
                sparse_tab(pos)  = self%loc_tab(ri,i)
            enddo
        enddo
        file_header(1) = self%nrefs
        file_header(2) = self%nptcls
        file_header(3) = nnz
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit, pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(funit, addr, self%seed_nrots, self%seed_shifts, self%seed_has_sh)
        write(funit, pos=addr) sparse_counts
        addr = addr + sizeof(sparse_counts)
        write(funit, pos=addr) sparse_refs
        addr = addr + sizeof(sparse_refs)
        write(funit, pos=addr) sparse_tab
        call fclose(funit)
        deallocate(sparse_tab, sparse_refs, sparse_counts)
    end subroutine write_tab_neigh

    ! Reads one partition sparse table and merges it into self%loc_tab
    subroutine read_sparse_tab_to_glob( self, binfname )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: binfname
        type(ptcl_ref), allocatable :: sparse_tab(:)
        real,           allocatable :: seed_shifts_loc(:,:)
        logical,        allocatable :: seed_has_sh_loc(:)
        integer,        allocatable :: pinds_loc(:), sparse_counts(:), sparse_refs(:), pind2glob(:)
        integer :: funit, addr, io_stat, file_header(3), nrefs_loc, nptcls_loc, nnz
        integer :: i_loc, i_glob, k, pos, ri, pind, max_pind, seed_nrots_loc
        if( file_exists(binfname) )then
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab_neigh; read_sparse_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD( 'sparse corr/rot files of partitions should be ready! ' )
        endif
        read(unit=funit,pos=1) file_header
        nrefs_loc  = file_header(1)
        nptcls_loc = file_header(2)
        nnz        = file_header(3)
        if( nrefs_loc .ne. self%nrefs ) THROW_HARD('nrefs mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        if( nnz < 1 ) THROW_HARD('empty sparse table in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        allocate(sparse_counts(nptcls_loc), sparse_refs(nnz), sparse_tab(nnz))
        addr = sizeof(file_header) + 1
        read(funit, pos=addr) pinds_loc
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        read(funit, pos=addr) sparse_counts
        addr = addr + sizeof(sparse_counts)
        read(funit, pos=addr) sparse_refs
        addr = addr + sizeof(sparse_refs)
        read(funit, pos=addr) sparse_tab
        call fclose(funit)
        if( any(sparse_counts < 0) .or. sum(sparse_counts) /= nnz )&
            &THROW_HARD('sparse table count mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, sparse_counts, sparse_refs, sparse_tab, pind2glob)
            return
        endif
        pos = 0
        do i_loc = 1,nptcls_loc
            pind = pinds_loc(i_loc)
            i_glob = 0
            if( pind >= 1 .and. pind <= max_pind ) i_glob = pind2glob(pind)
            if( i_glob > 0 )then
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            endif
            do k = 1,sparse_counts(i_loc)
                pos = pos + 1
                if( i_glob < 1 ) cycle
                ri = sparse_refs(pos)
                if( ri < 1 .or. ri > self%nrefs ) cycle
                self%loc_tab(ri,i_glob) = sparse_tab(pos)
                if( allocated(self%eval_touched_counts) .and. allocated(self%eval_touched_refs) )then
                    if( self%eval_touched_counts(i_glob) >= size(self%eval_touched_refs,1) )&
                        &THROW_HARD('eval_touched overflow in eul_prob_tab_neigh%read_sparse_tab_to_glob')
                    self%eval_touched_counts(i_glob) = self%eval_touched_counts(i_glob) + 1
                    self%eval_touched_refs(self%eval_touched_counts(i_glob),i_glob) = ri
                endif
            enddo
        enddo
        if( pos /= nnz ) THROW_HARD('sparse table count mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, sparse_counts, sparse_refs, sparse_tab, pind2glob)
    end subroutine read_sparse_tab_to_glob

    ! Aggregates nparts tables into the global table
    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(string) :: fname
        integer :: ipart
        if( .not. allocated(self%eval_touched_counts) ) allocate(self%eval_touched_counts(self%nptcls), source=0)
        if( .not. allocated(self%eval_touched_refs)   ) allocate(self%eval_touched_refs(self%nrefs,self%nptcls), source=0)
        if( allocated(self%eval_touched_counts) ) self%eval_touched_counts = 0
        if( allocated(self%eval_touched_refs)   ) self%eval_touched_refs   = 0
        do ipart = 1, nparts
            fname = fbody//int2str_pad(ipart,numlen)//'.dat'
            call self%read_sparse_tab_to_glob(fname)
        enddo
        self%eval_max_touched = size(self%eval_touched_refs,1)
    end subroutine read_tabs_to_glob

    ! COMMON STOCHASTIC/SUBSPACE ROUTINES

    ! Resets entries for the particle range
    subroutine clear_sparse_eval_table_range( self, i_from, i_to )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_from, i_to
        integer :: i, k, iref
        !$omp parallel do default(shared) private(i,k,iref) proc_bind(close) schedule(static)
        do i = i_from, i_to
            do k = 1, self%eval_touched_counts(i)
                iref = self%eval_touched_refs(k,i)
                if( iref > 0 )then
                    self%loc_tab(iref,i)%dist   = huge(1.0)
                    self%loc_tab(iref,i)%inpl   = 0
                    self%loc_tab(iref,i)%x      = 0.
                    self%loc_tab(iref,i)%y      = 0.
                    self%loc_tab(iref,i)%has_sh = .false.
                endif
            enddo
            self%eval_touched_counts(i) = 0
        enddo
        !$omp end parallel do
    end subroutine clear_sparse_eval_table_range

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

    ! Records one particle vs. ref evaluation result
    subroutine record_sparse_eval( self, i, ri, dist, irot, x, y, has_sh )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i, ri, irot
        real,                      intent(in)    :: dist, x, y
        logical,                   intent(in)    :: has_sh
        self%loc_tab(ri,i)%dist   = dist
        self%loc_tab(ri,i)%inpl   = irot
        self%loc_tab(ri,i)%x      = x
        self%loc_tab(ri,i)%y      = y
        self%loc_tab(ri,i)%has_sh = has_sh
        call mark_ref_touched(self, i, ri)
    end subroutine record_sparse_eval

    ! Minimizes shift at the previous best orientation
    subroutine estimate_shift_seed( self, grad_obj, ithr, iptcl, prev_state, prev_proj, o_prev, l_with_shift, shift_seed )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad),    intent(inout) :: grad_obj(:)
        integer,                   intent(in)    :: ithr, iptcl, prev_state, prev_proj
        type(ori),                 intent(inout) :: o_prev
        logical,                   intent(in)    :: l_with_shift
        real,                      intent(inout) :: shift_seed(3)
        integer :: irot, iref_start
        shift_seed = 0.
        if( .not. l_with_shift ) return
        irot = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
        if( prev_state >= 1 .and. prev_state <= self%p_ptr%nstates .and.&
            &prev_proj >= 1 .and. prev_proj <= self%p_ptr%nspace )then
            if( self%state_exists(prev_state) )then
                iref_start = (prev_state-1)*self%p_ptr%nspace
                call grad_obj(ithr)%set_indices(iref_start + prev_proj, iptcl)
                shift_seed = grad_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                if( irot == 0 ) shift_seed(2:3) = 0.
            endif
        endif
    end subroutine estimate_shift_seed

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
        real    :: refined_shift(3)
        ew%best_evals(:,ithr) = 0
        if( neval <= 0 ) return
        do si = 1, self%nstates
            state_eval = self%ssinds(si)
            ew%state_eval_dists(1:neval,ithr) = huge(1.0)
            nstate_eval = 0
            do j = 1, neval
                ri = ew%evaluated_ref_ids(j,ithr)
                if( ri < 1 ) cycle
                if( self%sinds(ri) /= state_eval ) cycle
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
                istate = self%sinds(ri)
                iproj  = self%jinds(ri)
                call grad_obj(ithr)%set_indices((istate-1)*self%p_ptr%nspace + iproj, iptcl)
                irot = self%loc_tab(ri,i)%inpl
                if( l_seed_sh_first )then
                    refined_shift = grad_obj(ithr)%minimize(irot=irot, sh_rot=.true.,&
                        &xy_in=shift_seed(2:3))
                else
                    refined_shift = grad_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                endif
                if( irot > 0 )then
                    call record_sparse_eval(self, i, ri, &
                        &eulprob_dist_switch(refined_shift(1), self%p_ptr%cc_objfun), &
                        &irot, refined_shift(2), refined_shift(3), .true.)
                endif
            enddo
        enddo
    end subroutine refine_best_neighbors

    ! Appends refernce index to the particle list of evaluated references
    subroutine mark_ref_touched( self, i, ri )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i, ri
        integer :: kt
        do kt = 1, self%eval_touched_counts(i)
            if( self%eval_touched_refs(kt,i) == ri ) return
        enddo
        if( self%eval_touched_counts(i) >= self%eval_max_touched )then
            THROW_HARD('simple_eul_prob_tab_neigh::mark_ref_touched; eval_touched overflow')
        endif
        self%eval_touched_counts(i) = self%eval_touched_counts(i) + 1
        self%eval_touched_refs(self%eval_touched_counts(i),i) = ri
    end subroutine mark_ref_touched

    ! DESTRUCTOR
    
    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        if( allocated(self%eval_touched_refs)   ) deallocate(self%eval_touched_refs)
        if( allocated(self%eval_touched_counts) ) deallocate(self%eval_touched_counts)
        self%eval_max_touched     = 0
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
            iref_full = (tabneigh%sinds(iref)-1)*tabneigh%p_ptr%nspace + tabneigh%jinds(iref)
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
            iref_full = (tabneigh%sinds(iref)-1)*tabneigh%p_ptr%nspace + tabneigh%jinds(iref)
            if( iref_full >= 1 .and. iref_full <= size(self%fullref_to_sparse_ref) )&
                &self%fullref_to_sparse_ref(iref_full) = iref
        enddo
        do ri = 1,tabneigh%nrefs
            istate = tabneigh%sinds(ri)
            isub   = tabneigh%b_ptr%subspace_full2sub_map(tabneigh%jinds(ri))
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
            istate = tabneigh%sinds(ri)
            isub   = tabneigh%b_ptr%subspace_full2sub_map(tabneigh%jinds(ri))
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
    subroutine alloc_coarse_ws( self, npeak_target, nstates, nthr )
        class(coarse_search_ws), intent(inout) :: self
        integer,                 intent(in)    :: npeak_target, nstates, nthr
        call self%dealloc_coarse_ws
        allocate(self%peak_subspace_dists(npeak_target,nstates,nthr), source=huge(1.0))
        allocate(self%peak_subspace_inds(npeak_target,nstates,nthr),  source=0)
        allocate(self%pooled_sub_inds(npeak_target,nstates,nthr),     source=0)
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
