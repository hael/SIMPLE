!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is geometrically sparse and searched independently per state.
module simple_eul_prob_tab_neigh
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_eul_prob_tab_utils
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab, POSTERIOR3D_FDR_Q
use simple_prob_posterior3D,    only: posterior3d_reader, posterior3d_writer, POSTERIOR3D_FNAME
use simple_decay_funs,         only: extremal_decay
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_ori,                only: ori
use simple_oris,               only: oris
use simple_syslib,              only: simple_rename
use simple_segmentation,        only: detect_peak_thres_fdr
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

integer, parameter :: POSTERIOR_EXPLORE_NNEIGH = 3

! Persistent sparse-posterior artifact state for one worker object.  The
! builder's full-to-subspace map is materialized as compact per-subspace
! lists so posterior rows can expand to fine-grid indices without rescanning
! the complete target angular grid for every particle.
type :: posterior3d_ws
    type(posterior3d_reader), allocatable :: readers(:)
    integer, allocatable :: pind_lookup(:)
    logical, allocatable :: posterior_available(:)
    integer, allocatable :: sub_offsets(:), sub_full_inds(:)
    integer, allocatable :: source_to_target_sub(:)
    type(oris) :: target_subspace
    integer :: source_nspace = 0
    integer :: coverage_count = 0
    logical :: initialized = .false.
    logical :: valid = .false.
    logical :: warning_emitted = .false.
contains
    procedure :: kill => kill_posterior3d_ws
end type posterior3d_ws

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    type(prob_candidate_store) :: candidate_store
    integer, allocatable       :: candidate_fill_counts(:)
    type(posterior3d_ws)        :: posterior_ws
    logical                    :: l_direct_stoch_neigh = .false.
contains
    procedure :: new_neigh
    procedure :: new_neigh_global
    procedure :: fill_tab       => fill_tab_neigh
    procedure :: fill_tab_range => fill_tab_neigh_range
    procedure :: ref_assign     => ref_assign_neigh
    procedure :: write_posterior3D_candidates
    procedure :: kill           => kill_neigh
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

    subroutine init_posterior3d_ws( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        logical :: ok
        integer :: i, ithr, max_artifact_pind, nsub, source_nspace, isrc, isub, target_proj
        type(ori) :: o_tgt
        if( self%posterior_ws%initialized ) return
        self%posterior_ws%initialized = .true.
        allocate(self%posterior_ws%readers(nthr_glob))
        self%posterior_ws%valid = .true.
        do ithr = 1, nthr_glob
            call self%posterior_ws%readers(ithr)%open(POSTERIOR3D_FNAME, ok)
            self%posterior_ws%valid = self%posterior_ws%valid .and. ok
        enddo
        if( self%posterior_ws%valid )then
            ! Rows are sparse and direct-access records are keyed by particle
            ! id.  The next stage may sample a different particle subset.
            self%posterior_ws%valid = self%posterior_ws%readers(1)%nptcls > 0 .and.&
                &self%posterior_ws%readers(1)%nstates == self%p_ptr%nstates .and.&
                &self%posterior_ws%readers(1)%source_nspace > 0 .and.&
                &self%posterior_ws%readers(1)%source_nspace <= self%p_ptr%nspace .and.&
                &self%posterior_ws%readers(1)%source_nspace_sub > 0 .and.&
                &self%posterior_ws%readers(1)%source_nspace_sub <= self%posterior_ws%readers(1)%source_nspace .and.&
                &self%p_ptr%nspace > 1 .and.&
                &trim(self%posterior_ws%readers(1)%pgrp) == trim(self%p_ptr%pgrp) .and.&
                    &allocated(self%b_ptr%subspace_inds) .and. allocated(self%b_ptr%subspace_full2sub_map) .and.&
                    &allocated(self%b_ptr%subspace_neighbors)
        endif
        if( self%posterior_ws%valid )then
            nsub = size(self%b_ptr%subspace_inds)
            self%posterior_ws%valid = self%posterior_ws%readers(1)%source_nspace > 0 .and.&
                &self%posterior_ws%readers(1)%source_nspace_sub > 0 .and.&
                &self%posterior_ws%readers(1)%source_nspace_sub <= self%posterior_ws%readers(1)%source_nspace .and.&
                &size(self%b_ptr%subspace_full2sub_map) == self%p_ptr%nspace .and.&
                &minval(self%b_ptr%subspace_full2sub_map) >= 1 .and.&
                &maxval(self%b_ptr%subspace_full2sub_map) <= nsub .and.&
                &size(self%b_ptr%subspace_neighbors,2) == nsub .and.&
                &size(self%b_ptr%subspace_neighbors,1) >= min(POSTERIOR_EXPLORE_NNEIGH+1,nsub)
            if( self%posterior_ws%valid )then
                call build_posterior_subspace_lists(self)
                source_nspace = self%posterior_ws%readers(1)%source_nspace
                allocate(self%posterior_ws%source_to_target_sub(source_nspace), source=0)
                if( source_nspace == self%p_ptr%nspace )then
                    ! Same full angular grid: use the target builder map directly.
                    do isrc = 1, source_nspace
                        self%posterior_ws%source_to_target_sub(isrc) = self%b_ptr%subspace_full2sub_map(isrc)
                    enddo
                else
                    ! Changed grid: retain target coarse representatives. Source
                    ! Euler coordinates are mapped lazily and cached below.
                    call self%posterior_ws%target_subspace%new(nsub, is_ptcl=.false.)
                    do isub = 1, nsub
                        target_proj = self%b_ptr%subspace_inds(isub)
                        call self%b_ptr%eulspace%get_ori(target_proj, o_tgt)
                        call self%posterior_ws%target_subspace%set_ori(isub, o_tgt)
                    enddo
                    call o_tgt%kill
                endif
            endif
        endif
        if( self%posterior_ws%valid )then
            do ithr = 2, nthr_glob
                self%posterior_ws%valid = self%posterior_ws%valid .and.&
                    &self%posterior_ws%readers(ithr)%source_nspace == self%posterior_ws%readers(1)%source_nspace .and.&
                    &self%posterior_ws%readers(ithr)%source_nspace_sub == self%posterior_ws%readers(1)%source_nspace_sub
            enddo
        endif
        if( .not. self%posterior_ws%valid )then
            do ithr = 1, nthr_glob
                call self%posterior_ws%readers(ithr)%kill
            enddo
            return
        endif
        self%posterior_ws%source_nspace = self%posterior_ws%readers(1)%source_nspace
        call build_pind_lookup(self%posterior_ws%readers(1)%pinds, self%pinds,&
            &self%posterior_ws%pind_lookup, max_artifact_pind)
        allocate(self%posterior_ws%posterior_available(self%nptcls), source=.false.)
        do i = 1, self%nptcls
            if( self%pinds(i) >= 1 .and. self%pinds(i) <= max_artifact_pind )then
                self%posterior_ws%posterior_available(i) = self%posterior_ws%pind_lookup(self%pinds(i)) > 0
            endif
        enddo
        self%posterior_ws%coverage_count = count(self%posterior_ws%posterior_available)
        if( self%posterior_ws%coverage_count == 0 ) self%posterior_ws%valid = .false.
    end subroutine init_posterior3d_ws

    subroutine build_posterior_subspace_lists( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: counts(:), cursor(:)
        integer :: nsub, nfull, i, isub
        nsub  = size(self%b_ptr%subspace_inds)
        nfull = self%p_ptr%nspace
        allocate(counts(nsub), cursor(nsub), self%posterior_ws%sub_offsets(nsub+1),&
            &self%posterior_ws%sub_full_inds(nfull))
        counts = 0
        do i = 1, nfull
            counts(self%b_ptr%subspace_full2sub_map(i)) = counts(self%b_ptr%subspace_full2sub_map(i)) + 1
        enddo
        self%posterior_ws%sub_offsets(1) = 1
        do isub = 1, nsub
            self%posterior_ws%sub_offsets(isub+1) = self%posterior_ws%sub_offsets(isub) + counts(isub)
        enddo
        cursor = self%posterior_ws%sub_offsets(1:nsub)
        do i = 1, nfull
            isub = self%b_ptr%subspace_full2sub_map(i)
            self%posterior_ws%sub_full_inds(cursor(isub)) = i
            cursor(isub) = cursor(isub) + 1
        enddo
        deallocate(counts, cursor)
    end subroutine build_posterior_subspace_lists

    ! Maps one source full-grid projection to the target coarse subspace.  For
    ! changed grids the stored source Euler coordinates are authoritative; the
    ! result is cached by source projection identity so the expensive angular
    ! lookup is performed at most once per source direction.
    subroutine ensure_source_to_target_sub( self, source_proj, source_euls )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, intent(in) :: source_proj
        real, intent(in) :: source_euls(3)
        type(ori) :: o_src
        integer :: cached_sub
        if( source_proj < 1 .or. source_proj > size(self%posterior_ws%source_to_target_sub) ) return
        if( .not. all(ieee_is_finite(source_euls)) ) return
        !$omp atomic read
        cached_sub = self%posterior_ws%source_to_target_sub(source_proj)
        if( cached_sub > 0 ) return
        !$omp critical(posterior_source_to_target_map)
        if( self%posterior_ws%source_to_target_sub(source_proj) == 0 )then
            call o_src%new(is_ptcl=.false.)
            call o_src%set_euler(source_euls)
            cached_sub = self%b_ptr%pgrpsyms%find_closest_proj(self%posterior_ws%target_subspace, o_src)
            !$omp atomic write
            self%posterior_ws%source_to_target_sub(source_proj) = cached_sub
            call o_src%kill
        endif
        !$omp end critical(posterior_source_to_target_map)
    end subroutine ensure_source_to_target_sub

    subroutine new_neigh( self, params, build, pinds )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%eul_prob_tab%new_worker(params,build,pinds)
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
            case('posterior')
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

    ! Fills the full particle range by delegating to fill_tab_neigh_range(1, nptcls).
    subroutine fill_tab_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%fill_tab_range(1, self%nptcls)
    end subroutine fill_tab_neigh

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
            case('posterior')
                call fill_tab_posterior_range(self, i_first, i_last)
            case DEFAULT
                THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh_range; unsupported prob_neigh_mode')
        end select
        if( self%table_is_open ) call self%flush_candidate_buffers
    end subroutine fill_tab_neigh_range

    ! SPARSE PROFILED POSTERIOR BRANCH
    !
    ! The source artifact is a sparse posterior on the preceding dense grid.
    ! The source posterior is indexed on the builder's coarse subspace.  The
    ! existing full-to-subspace map expands each retained coarse direction to
    ! the fine-grid directions assigned to that subspace cell.  This keeps the
    ! number of target candidates proportional to the resolution ratio rather
    ! than to an angular-radius union over the full target grid.
    subroutine fill_tab_posterior_range( self, i_first, i_last )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, intent(in) :: i_first, i_last
        integer, allocatable :: mapped_stamp(:,:), map_generation(:)
        logical, allocatable :: posterior_usable(:), recovery_mask(:)
        integer, allocatable :: inds_sorted(:,:)
        real, allocatable :: dists_inpl(:,:)
        real :: shift(2), dist, corr
        integer :: i, j, k, iproj, ri, full_ref, irot, ithr, istate, nrots, ineigh
        integer :: i_from, i_to, ncovered, nrange, source_proj, source_sub, coarse_proj
        logical :: okrow
        character(len=8) :: fallback_mode
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        ! Preserve the workflow's ordinary neighborhood policy for the
        ! conditional probe.  refine3D_multi/input-orientation workflows use
        ! the geometric path; single/independent ab-initio workflows use state.
        fallback_mode = 'state'
        select case(trim(self%p_ptr%multivol_mode))
            case('docked','geom','input_oris_refine','input_oris_fixed')
                fallback_mode = 'geom'
        end select
        call init_posterior3d_ws(self)
        if( .not. self%posterior_ws%valid )then
            if( .not. self%posterior_ws%warning_emitted )then
                THROW_WARN('prob_neigh_mode=posterior: support artifact missing or incompatible; falling back to bounded prob_neigh_mode='//trim(fallback_mode)//' search')
                self%posterior_ws%warning_emitted = .true.
            endif
            call fill_tab_subspace_range(self, i_first, i_last, fallback_mode)
            return
        endif
        nrange = i_to - i_from + 1
        allocate(posterior_usable(self%nptcls), source=self%posterior_ws%posterior_available)
        allocate(recovery_mask(self%nptcls), source=.false.)
        ncovered = count(posterior_usable(i_from:i_to))
        if( ncovered == 0 )then
            call fill_tab_subspace_range(self, i_first, i_last, fallback_mode)
            deallocate(posterior_usable, recovery_mask)
            return
        endif
        nrots = self%b_ptr%pftc%get_nrots()
        ! Every posterior worker writes the deferred-shift metadata, even when
        ! no particle in this partition needs the bounded ordinary-search
        ! recovery.  Set this before the posterior-only path can write a table;
        ! otherwise partitions that happened to enter the recovery path publish
        ! a nonzero seed grid while the others publish the default zero.
        self%seed_nrots = nrots
        allocate(mapped_stamp(self%nrefs,nthr_glob), map_generation(nthr_glob),&
            &dists_inpl(nrots,nthr_glob), inds_sorted(nrots,nthr_glob))
        mapped_stamp = 0
        map_generation = 0
        dists_inpl = 0.
        inds_sorted = 0
        !$omp parallel do default(shared) private(i,ithr,j,k,iproj,ri,full_ref,irot,istate,shift,dist,corr,okrow,source_proj,source_sub,coarse_proj,ineigh)&
        !$omp proc_bind(close) schedule(static)
        do i = i_from, i_to
            if( .not. posterior_usable(i) ) cycle
            ithr = omp_get_thread_num() + 1
            map_generation(ithr) = map_generation(ithr) + 1
            call self%posterior_ws%readers(ithr)%read_row(self%pinds(i), okrow)
            if( .not. okrow )then
                posterior_usable(i) = .false.
                cycle
            endif
            recovery_mask(i) = posterior_row_needs_recovery(self, self%pinds(i), self%posterior_ws%readers(ithr))
            do k = 1, self%posterior_ws%readers(ithr)%nsel
                istate = self%posterior_ws%readers(ithr)%state(k)
                if( istate < 1 .or. istate > self%p_ptr%nstates ) cycle
                if( .not. self%state_exists(istate) ) cycle
                source_proj = self%posterior_ws%readers(ithr)%proj(k)
                if( source_proj < 1 .or. source_proj > size(self%posterior_ws%source_to_target_sub) ) cycle
                call ensure_source_to_target_sub(self, source_proj, self%posterior_ws%readers(ithr)%euls(:,k))
                source_sub = self%posterior_ws%source_to_target_sub(source_proj)
                if( source_sub < 1 .or. source_sub > size(self%posterior_ws%sub_offsets)-1 ) cycle
                do ineigh = 1, min(POSTERIOR_EXPLORE_NNEIGH+1, size(self%b_ptr%subspace_neighbors,1))
                    coarse_proj = self%b_ptr%subspace_neighbors(ineigh,source_sub)
                    if( coarse_proj < 1 .or. coarse_proj > size(self%posterior_ws%sub_offsets)-1 ) cycle
                    do j = self%posterior_ws%sub_offsets(coarse_proj), self%posterior_ws%sub_offsets(coarse_proj+1)-1
                        iproj = self%posterior_ws%sub_full_inds(j)
                        full_ref = (istate-1)*self%p_ptr%nspace + iproj
                        ri = self%full_to_compact_ref(full_ref)
                        if( ri < 1 .or. mapped_stamp(ri,ithr) == map_generation(ithr) ) cycle
                        mapped_stamp(ri,ithr) = map_generation(ithr)
                        shift = 0.
                        if( self%posterior_ws%readers(ithr)%has_sh(k) /= 0 .and. self%p_ptr%l_doshift )&
                            &shift = [self%posterior_ws%readers(ithr)%x(k), self%posterior_ws%readers(ithr)%y(k)]
                        call self%b_ptr%pftc%gen_likelihood_val(full_ref, self%pinds(i), shift,&
                            &nrots, dist, corr, irot, dists_inpl(:,ithr), inds_sorted(:,ithr))
                        if( irot < 1 .or. irot > nrots ) irot = 1
                        call record_sparse_eval(self, i, ithr, ri, dist, irot, shift(1), shift(2),&
                            &(self%posterior_ws%readers(ithr)%has_sh(k) /= 0 .and. self%p_ptr%l_doshift))
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        ! A malformed row may only become visible during the direct read.  Recount
        ! after the parallel section so those particles enter the bounded fallback.
        ncovered = count(posterior_usable(i_from:i_to))
        if( ncovered < nrange )then
            ! Particle sampling may change between producer and consumer.
            ! Particles without a valid posterior row or distance estimate use
            ! the bounded workflow-specific fallback.
            call fill_tab_subspace_range(self, i_first, i_last, fallback_mode,&
                &active_mask=.not.posterior_usable)
        endif
        if( any(recovery_mask(i_from:i_to)) )then
            write(logfhandle,'(A,I0,A,A)') '>>> POSTERIOR BOUNDED RECOVERY PROBES: ',&
                &count(recovery_mask(i_from:i_to)), ' particles using prob_neigh_mode=', trim(fallback_mode)
            call fill_tab_subspace_range(self, i_first, i_last, fallback_mode, active_mask=recovery_mask)
        endif
        deallocate(mapped_stamp, map_generation, posterior_usable, recovery_mask, dists_inpl, inds_sorted)
    end subroutine fill_tab_posterior_range

    ! A posterior row is considered in need of an ordinary recovery probe only
    ! when it is empty or its coarse support no longer contains the particle's
    ! current projection cell (including the local three-neighbor ring).
    ! The ordinary state/geom path is used only for these particles; it is not a
    ! global fallback for a valid posterior artifact.
    logical function posterior_row_needs_recovery( self, pind, reader ) result(l_recover)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, intent(in) :: pind
        type(posterior3d_reader), intent(in) :: reader
        type(ori) :: o_prev
        integer :: prev_state, prev_proj, prev_sub, k, ineigh, source_proj, source_sub, nkeep
        logical :: l_found
        ! Low retained mass requests bounded posterior exploration, but does
        ! not by itself launch the expensive ordinary neighborhood search.
        l_recover = reader%nsel == 0
        call get_particle_context(self, pind, o_prev, prev_state, prev_proj)
        prev_proj = max(1, min(self%p_ptr%nspace, prev_proj))
        prev_sub = self%b_ptr%subspace_full2sub_map(prev_proj)
        nkeep = min(POSTERIOR_EXPLORE_NNEIGH+1, size(self%b_ptr%subspace_neighbors,1))
        l_found = .false.
        do k = 1, reader%nsel
            if( prev_state >= 1 .and. prev_state <= self%p_ptr%nstates )then
                if( reader%state(k) /= prev_state ) cycle
            endif
            source_proj = reader%proj(k)
            if( source_proj < 1 .or. source_proj > size(self%posterior_ws%source_to_target_sub) ) cycle
            call ensure_source_to_target_sub(self, source_proj, reader%euls(:,k))
            source_sub = self%posterior_ws%source_to_target_sub(source_proj)
            if( source_sub < 1 .or. source_sub > size(self%b_ptr%subspace_neighbors,2) ) cycle
            do ineigh = 1, nkeep
                if( self%b_ptr%subspace_neighbors(ineigh,source_sub) == prev_sub )then
                    l_found = .true.
                    exit
                endif
            enddo
            if( l_found ) exit
        enddo
        l_recover = l_recover .or. .not. l_found
        call o_prev%kill
    end function posterior_row_needs_recovery

    ! Rebuilds the sparse posterior from the merged current candidate evidence.
    ! This is called only by the global prob_neigh owner, after all distributed
    ! candidate tables have been merged.  The source rows are therefore based
    ! on the current reference and current calibrated distances, not on a
    ! cumulative product of evidence from previous iterations.
    subroutine write_posterior3D_candidates( self, binfname )
        class(eul_prob_tab_neigh), intent(in) :: self
        character(len=*), intent(in) :: binfname
        type(posterior3d_writer) :: writer
        integer, allocatable :: work_pos(:), seen_stamp(:), seen_index(:)
        integer, allocatable :: out_state(:), out_proj(:), out_inpl(:), out_has_sh(:)
        real, allocatable :: work_dist(:), work_weight(:), out_dist(:), out_weight(:), out_x(:), out_y(:), out_euls(:,:)
        character(len=STDLEN) :: tmpfname, tmpmeta, metafname
        real :: dmin, norm, euls(3), dist_thres
        integer :: i, pos, pick, ncand, k, k_this, kmax, iref, iproj, istate, nsub, npeaks_detected
        if( .not. allocated(self%candidate_store%offsets) .or. self%nptcls < 1 .or. self%nrefs < 1 )then
            THROW_WARN('write_posterior3D_candidates requires merged candidate evidence; no artifact written')
            return
        endif
        nsub = 0
        if( allocated(self%b_ptr%subspace_inds) ) nsub = size(self%b_ptr%subspace_inds)
        if( nsub < 1 )then
            THROW_WARN('write_posterior3D_candidates requires a coarse subspace map; no artifact written')
            return
        endif
        kmax = min(self%nrefs, max(3, min(128, max(1, 8*self%p_ptr%npeaks_inpl))))
        allocate(work_pos(self%nrefs), seen_stamp(self%nrefs), seen_index(self%nrefs),&
            &work_dist(self%nrefs), work_weight(self%nrefs))
        seen_stamp = 0
        allocate(out_state(kmax), out_proj(kmax), out_inpl(kmax), out_has_sh(kmax), out_dist(kmax), out_weight(kmax),&
            &out_x(kmax), out_y(kmax), out_euls(3,kmax))
        tmpfname = trim(binfname)//'.tmp'
        tmpmeta  = trim(tmpfname)//'.meta'
        metafname = trim(binfname)//'.meta'
        call writer%open(trim(tmpfname), self%nptcls, self%p_ptr%nstates, self%p_ptr%nspace, nsub, kmax,&
            &self%p_ptr%pgrp, self%pinds)
        do i = 1, self%nptcls
            ncand = 0
            do pos = int(self%candidate_store%offsets(i)), int(self%candidate_store%offsets(i+1))-1
                if( self%candidate_store%candidates(pos)%inpl < 1 ) cycle
                if( .not. ieee_is_finite(self%candidate_store%candidates(pos)%dist) ) cycle
                if( self%candidate_store%candidates(pos)%iref < 1 .or.&
                    &self%candidate_store%candidates(pos)%iref > self%p_ptr%nstates*self%p_ptr%nspace ) cycle
                iref = self%full_to_compact_ref(self%candidate_store%candidates(pos)%iref)
                if( iref < 1 .or. iref > self%nrefs ) cycle
                if( seen_stamp(iref) == i )then
                    pick = seen_index(iref)
                    if( self%candidate_store%candidates(pos)%dist < work_dist(pick) )then
                        work_pos(pick)  = pos
                        work_dist(pick) = self%candidate_store%candidates(pos)%dist
                    endif
                else
                    ncand = ncand + 1
                    seen_stamp(iref) = i
                    seen_index(iref) = ncand
                    work_pos(ncand)  = pos
                    work_dist(ncand) = self%candidate_store%candidates(pos)%dist
                endif
            enddo
            if( ncand < 1 )then
                call writer%write_row(self%pinds(i), 0, 0., out_state, out_proj, out_inpl, out_has_sh,&
                    &out_dist, out_weight, out_x, out_y, out_euls)
                cycle
            endif
            dmin = minval(work_dist(1:ncand))
            work_weight(1:ncand) = exp(-(work_dist(1:ncand)-dmin))
            norm = sum(work_weight(1:ncand))
            if( norm > 0. ) work_weight(1:ncand) = work_weight(1:ncand) / norm
            call detect_peak_thres_fdr(ncand, work_dist(1:ncand), POSTERIOR3D_FDR_Q, 1, min(kmax,ncand),&
                &dist_thres, npeaks_detected, lower_tail=.true.)
            k_this = min(kmax, ncand, max(1, npeaks_detected))
            out_state = 0
            out_proj = 0
            out_inpl = 0
            out_has_sh = 0
            out_dist = 0.
            out_weight = 0.
            out_x = 0.
            out_y = 0.
            out_euls = 0.
            do k = 1, k_this
                pick = minloc(work_dist(1:ncand), dim=1)
                pos = work_pos(pick)
                iref = self%full_to_compact_ref(self%candidate_store%candidates(pos)%iref)
                if( iref < 1 .or. iref > self%nrefs ) THROW_HARD('posterior refresh contains an invalid reference')
                istate = self%ref_state(iref)
                iproj  = self%ref_proj(iref)
                out_state(k) = istate
                out_proj(k) = iproj
                out_inpl(k) = self%candidate_store%candidates(pos)%inpl
                out_has_sh(k) = merge(1, 0, self%candidate_store%candidates(pos)%has_sh)
                out_dist(k) = self%candidate_store%candidates(pos)%dist
                out_weight(k) = work_weight(pick)
                out_x(k) = self%candidate_store%candidates(pos)%x
                out_y(k) = self%candidate_store%candidates(pos)%y
                euls = self%b_ptr%eulspace%get_euler(iproj)
                out_euls(:,k) = euls
                work_dist(pick) = huge(1.0)
            enddo
            call writer%write_row(self%pinds(i), k_this, sum(out_weight(1:k_this)), out_state, out_proj, out_inpl,&
                &out_has_sh, out_dist, out_weight, out_x, out_y, out_euls)
        enddo
        call writer%close
        call writer%kill
        call simple_rename(trim(tmpfname), trim(binfname), overwrite=.true.)
        call simple_rename(trim(tmpmeta), trim(metafname), overwrite=.true.)
        write(logfhandle,'(A,I0,A,I0)') '>>> POSTERIOR REFRESH PUBLISHED: particles=', self%nptcls,&
            &', maximum support width=', kmax
        deallocate(work_pos, seen_stamp, seen_index, work_dist, work_weight, out_state, out_proj, out_inpl, out_has_sh, out_dist, out_weight,&
            &out_x, out_y, out_euls)
    end subroutine write_posterior3D_candidates

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
        logical :: l_shc_neigh, l_snhc_neigh, l_seed_sh_first
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        l_shc_neigh     = (trim(self%p_ptr%prob_neigh_mode) == 'shc')
        l_snhc_neigh    = (trim(self%p_ptr%prob_neigh_mode) == 'snhc')
        l_seed_sh_first = self%p_ptr%l_doshift .and. l_shc_neigh
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        nfull_refs = self%p_ptr%nstates * self%p_ptr%nspace
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
            real    :: dist_loc, corr_loc, prev_corr_loc, neigh_frac
            logical :: l_greedy_first
            neval_loc = 0
            call direct_rts(ithr_loc)%ne_ran_iarr(direct_srch_order(:,ithr_loc))
            prev_full_ref_loc = 0
            if( prev_state_loc >= 1 .and. prev_state_loc <= self%p_ptr%nstates .and.&
                &prev_proj_loc >= 1 .and. prev_proj_loc <= self%p_ptr%nspace )then
                prev_full_ref_loc = (prev_state_loc - 1) * self%p_ptr%nspace + prev_proj_loc
                call put_last(prev_full_ref_loc, direct_srch_order(:,ithr_loc))
            endif
            ! Only force a dense first pass when the particle truly has no previous
            ! search result. Keying this off updatecnt makes iteration 2 behave
            ! like a full scan for first-time sampled particles.
            l_greedy_first = l_shc_neigh .and. (.not. self%b_ptr%spproj_field%has_been_searched(iptcl_loc))
            nrots          = self%b_ptr%pftc%get_nrots()
            nrefs_bound    = nfull_refs
            smpl_ninpl     = nrots
            if( l_snhc_neigh )then
                neigh_frac  = extremal_decay(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
                nrefs_bound = max(2, min(nfull_refs, nint(real(nfull_refs) * (1. - neigh_frac))))
                smpl_ninpl  = neighfrac2nsmpl(neigh_frac, nrots)
            endif
            prev_corr_loc = -huge(1.0)
            if( (l_shc_neigh .and. (.not. l_greedy_first)) )then
                call calc_previous_corr(ithr_loc, iptcl_loc, prev_full_ref_loc, prev_corr_loc)
            endif
            do isample = 1,nfull_refs
                if( l_snhc_neigh .and. isample > nrefs_bound ) exit
                full_ref_loc = direct_srch_order(isample,ithr_loc)
                if( full_ref_loc < 1 .or. full_ref_loc > size(eval_work%fullref_to_sparse_ref) ) cycle
                ri_loc = eval_work%fullref_to_sparse_ref(full_ref_loc)
                if( ri_loc < 1 ) cycle
                call score_direct_ref(full_ref_loc, ithr_loc, iptcl_loc, shift_seed, l_with_shift,&
                    &l_snhc_neigh, smpl_ninpl, dist_loc, corr_loc, irot_loc)
                call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
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
            &l_snhc_style, smpl_ninpl_loc, dist_loc, corr_loc, irot_loc)
            integer, intent(in)  :: full_ref_loc, ithr_loc, iptcl_loc, smpl_ninpl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift, l_snhc_style
            real,    intent(out) :: dist_loc, corr_loc
            integer, intent(out) :: irot_loc
            integer :: istate_loc, nsample_loc
            real    :: sh_loc(2)
            sh_loc = 0.
            if( l_with_shift .and. l_seed_sh_first ) sh_loc = shift_seed_loc(2:3)
            istate_loc = (full_ref_loc - 1) / self%p_ptr%nspace + 1
            nsample_loc = merge(smpl_ninpl_loc, inpl_likelihood_nsample(istate_loc), l_snhc_style)
            call self%b_ptr%pftc%gen_likelihood_val(full_ref_loc, iptcl_loc, sh_loc,&
                &nsample_loc, dist_loc, corr_loc, irot_loc, dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
        end subroutine score_direct_ref

        integer function inpl_likelihood_nsample( istate_loc ) result(nsample)
            integer, intent(in) :: istate_loc
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, inpl_athres(istate_loc))
            nsample = min(nrots, max(1, int(athres_ub * real(nrots) / 180.)))
        end function inpl_likelihood_nsample

    end subroutine fill_tab_stoch_range

    ! NEIGHBORHOOD SUBSPACES BRANCH

    ! Identifies peak subspaces via coarse scoring, builds a pooled neighborhoo
    ! and evaluates all fine refs & refines the best
    subroutine fill_tab_subspace_range( self, i_first, i_last, search_mode, active_mask )
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: i_first, i_last
        character(len=*), optional, intent(in)   :: search_mode
        logical, optional, intent(in)            :: active_mask(:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(coarse_search_ws) :: coarse_ws
        type(eval_ws)          :: eval_work
        integer, allocatable   :: inds_sorted(:,:)
        real,    allocatable   :: inpl_athres(:), dists_inpl(:,:), dists_inpl_sorted(:,:)
        integer :: i, istate, ithr, max_refs_to_refine, nsubs, npeak_target, npooled_capacity
        integer :: iptcl, si, i_from, i_to, nrots
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_geom_neigh, l_state_neigh, l_seed_sh_first
        character(len=STDLEN) :: active_search_mode
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        active_search_mode = trim(self%p_ptr%prob_neigh_mode)
        if( present(search_mode) ) active_search_mode = trim(search_mode)
        l_geom_neigh    = (trim(active_search_mode) == 'geom')
        l_state_neigh   = (trim(active_search_mode) == 'state')
        ! Keep shift-first seeding single-state only. In multi-state runs this
        ! can bias all states toward the same previous-state shift basin.
        l_seed_sh_first = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        if( present(active_mask) )then
            if( size(active_mask) < i_to ) THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_subspace_range; active mask is too small')
        endif
        call seed_rnd
        nsubs        = size(self%b_ptr%subspace_inds)
        npeak_target = min(max(1, self%p_ptr%npeaks), nsubs)
        npooled_capacity = min(nsubs, max(1, npeak_target * max(1, self%nstates)))
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(nrots,nthr_glob), dists_inpl_sorted(nrots,nthr_glob), inds_sorted(nrots,nthr_glob))
        max_refs_to_refine = max(1, self%p_ptr%npeaks_inpl)
        call eval_work%init_eval_subspace_ws(self, nsubs, max_refs_to_refine)
        do si = 1, self%nstates
            istate = self%ssinds(si)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        call coarse_ws%alloc_coarse_ws(npeak_target, npooled_capacity, self%p_ptr%nstates, nthr_glob)
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
                if( present(active_mask) )then
                    if( .not. active_mask(i) ) cycle
                endif
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
                if( present(active_mask) )then
                    if( .not. active_mask(i) ) cycle
                endif
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
            integer :: istate_loc
            real    :: corr_loc
            istate_loc = (full_ref_loc - 1) / self%p_ptr%nspace + 1
            if( l_with_shift )then
                call self%b_ptr%pftc%gen_likelihood_val(full_ref_loc, iptcl_loc,&
                    &shift_seed_loc(2:3), inpl_likelihood_nsample(istate_loc), dist_loc, corr_loc, irot_loc,&
                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
            else
                call self%b_ptr%pftc%gen_likelihood_val(full_ref_loc, iptcl_loc,&
                    &[0.,0.], inpl_likelihood_nsample(istate_loc), dist_loc, corr_loc, irot_loc,&
                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
            endif
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
            real    :: dist, corr_loc
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
                        iproj_loc = self%ref_proj(ri_loc)
                        iref_loc  = self%ref_full(ri_loc)
                        if( l_with_shift )then
                            call self%b_ptr%pftc%gen_likelihood_val(iref_loc, iptcl_loc,&
                                &shift_seed_loc(2:3), inpl_likelihood_nsample(istate_loc), dist, corr_loc, irot_loc,&
                                &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
                        else
                            call self%b_ptr%pftc%gen_likelihood_val(iref_loc, iptcl_loc,&
                                &[0.,0.], inpl_likelihood_nsample(istate_loc), dist, corr_loc, irot_loc,&
                                &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
                        endif
                        call record_sparse_eval(self, i_loc, ithr_loc, ri_loc, dist, irot_loc, 0., 0., .false.)
                        neval = neval + 1
                        eval_work%evaluated_ref_ids(neval,ithr_loc)   = ri_loc
                        eval_work%evaluated_ref_dists(neval,ithr_loc) = dist
                    enddo
                enddo
            enddo
        end subroutine evaluate_neighborhood

        integer function inpl_likelihood_nsample( istate_loc ) result(nsample)
            integer, intent(in) :: istate_loc
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, inpl_athres(istate_loc))
            nsample = min(nrots, max(1, int(athres_ub * real(nrots) / 180.)))
        end function inpl_likelihood_nsample

    end subroutine fill_tab_subspace_range

    ! ASSIGNMENT SUBROUTINES

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
            integer, allocatable :: ref_candidate_pos(:)
            real,    allocatable :: ref_dists(:)
        end type assign_graph_ws

        type :: assign_frontier_ws
            integer, allocatable :: inds_sorted(:), order(:), work_ptcl(:), work_cpos(:), sel_refs(:), sel_pos(:)
            real,    allocatable :: iref_dist(:), dists_sorted(:), work_d(:), sel_dists(:), sel_dists_sorted(:)
            logical, allocatable :: ptcl_avail(:)
        end type assign_frontier_ws

        type(assign_graph_ws)    :: graph
        type(assign_frontier_ws) :: frontier
        real,    allocatable :: dists_inpl_sorted(:)
        integer, allocatable :: inds_sorted(:)
        integer   :: i, iref, assigned_iref, assigned_ptcl, istate, si, fallback_ref
        integer   :: k, idx, nactive, total, m, start, maxref, nleft, assigned_idx, nsel, pos, last_ref
        integer   :: candidate_pos, first_candidate, last_candidate
        integer   :: greedy_state(self%nptcls)
        real      :: projs_athres, state_projs_athres(self%p_ptr%nstates)
        real      :: huge_val, dist_tmp, corr_tmp
        logical   :: final_assigned(self%nptcls), l_filter_greedy_state, l_seed_fallback_shift, l_report_npeaks
        huge_val = huge(1.0)
        l_report_npeaks = trim(self%p_ptr%refine) == 'prob_neigh' .and.&
            &trim(self%p_ptr%prob_neigh_mode) == 'posterior'
        l_seed_fallback_shift = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1 .and. &
            &trim(self%p_ptr%prob_neigh_mode) /= 'snhc'
        allocate(dists_inpl_sorted(self%b_ptr%pftc%get_nrots()), inds_sorted(self%b_ptr%pftc%get_nrots()))
        do i = 1, self%nptcls
            ! Seed particle with no evaluated candidate with previous orientation
            call seed_fallback_if_empty(i)
        enddo
        call build_sparse_assignment_graph()
        call assign_particles_globally()
        call assign_remaining_particles_from_best_touched_ref()
        deallocate(inds_sorted, dists_inpl_sorted)
    contains

        subroutine build_sparse_assignment_graph()
            allocate(graph%ref_counts(self%nrefs), graph%ref_offsets(self%nrefs+1), graph%ref_fill(self%nrefs), graph%ref_pos(self%nrefs),&
            &frontier%iref_dist(self%nrefs), frontier%dists_sorted(self%nrefs), frontier%inds_sorted(self%nrefs), frontier%ptcl_avail(self%nptcls),&
            &graph%ptcl_counts(self%nptcls), graph%ptcl_offsets(self%nptcls+1), graph%ptcl_fill(self%nptcls))
            graph%ref_counts = 0
            graph%ptcl_counts = 0
            do i = 1, self%nptcls
                first_candidate = int(self%candidate_store%offsets(i))
                last_candidate  = int(self%candidate_store%offsets(i+1)) - 1
                do candidate_pos = first_candidate,last_candidate
                    iref = self%full_to_compact_ref(self%candidate_store%candidates(candidate_pos)%iref)
                    if( iref < 1 .or. iref > self%nrefs       ) cycle
                    if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 ) cycle
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
            allocate(graph%ref_list(max(1,total)), graph%ref_candidate_pos(max(1,total)), graph%ref_dists(max(1,total)))
            graph%ptcl_offsets(1) = 1
            do i = 1, self%nptcls
                graph%ptcl_offsets(i+1) = graph%ptcl_offsets(i) + graph%ptcl_counts(i)
            enddo
            allocate(graph%ptcl_refs(max(1,graph%ptcl_offsets(self%nptcls+1)-1)))
            graph%ref_fill = graph%ref_offsets(1:self%nrefs)
            graph%ptcl_fill = graph%ptcl_offsets(1:self%nptcls)
            do i = 1, self%nptcls
                first_candidate = int(self%candidate_store%offsets(i))
                last_candidate  = int(self%candidate_store%offsets(i+1)) - 1
                do candidate_pos = first_candidate,last_candidate
                    iref = self%full_to_compact_ref(self%candidate_store%candidates(candidate_pos)%iref)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 ) cycle
                    idx = graph%ref_fill(iref)
                    graph%ref_list(idx)          = i
                    graph%ref_candidate_pos(idx) = candidate_pos
                    graph%ref_dists(idx)         = self%candidate_store%candidates(candidate_pos)%dist
                    graph%ref_fill(iref) = idx + 1
                    idx = graph%ptcl_fill(i)
                    graph%ptcl_refs(idx) = iref
                    graph%ptcl_fill(i)   = idx + 1
                enddo
            enddo
            maxref = max(1, maxval(graph%ref_counts))
            allocate(frontier%work_d(maxref), frontier%order(maxref), frontier%work_ptcl(maxref), frontier%work_cpos(maxref))
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
                frontier%work_cpos(1:m) = graph%ref_candidate_pos(start:start+m-1)
                frontier%order(1:m)     = (/(k,k=1,m)/)
                call hpsort(frontier%work_d(1:m), frontier%order(1:m))
                do k = 1, m
                    graph%ref_dists(start+k-1) = frontier%work_d(k)
                    graph%ref_list(start+k-1)  = frontier%work_ptcl(frontier%order(k))
                    graph%ref_candidate_pos(start+k-1) = frontier%work_cpos(frontier%order(k))
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
                ! Multi-state state labelling is deterministic (argmin distance) for both weighting
                ! schemes; probabilistic exploration is confined to the within-state projection
                ! assignment (assign_particles_for_state). Only refine=prob_state samples the state.
                assigned_idx = minloc(frontier%sel_dists(1:nsel), dim=1)
                assigned_iref = frontier%sel_refs(assigned_idx)
                assigned_ptcl = graph%ref_list(graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1)
                greedy_state(assigned_ptcl) = self%ref_state(assigned_iref)
                frontier%ptcl_avail(assigned_ptcl) = .false.
                nleft = nleft - 1
                call update_frontier_after_assignment(assigned_ptcl)
            enddo
            do i = 1,self%nptcls
                if( greedy_state(i) == 0 )then
                    fallback_ref = pick_best_evaluated_ref(i)
                    if( fallback_ref > 0 ) greedy_state(i) = self%ref_state(fallback_ref)
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
                call sample_likelihood_dist(nsel, selected_frontier_dist, likelihood_nsample(nsel, projs_athres),&
                    &dist_tmp, corr_tmp, assigned_idx, frontier%sel_dists_sorted, frontier%inds_sorted)
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
                call sample_likelihood_dist(nsel, selected_frontier_dist, likelihood_nsample(nsel, state_projs_athres(state_filter)),&
                    &dist_tmp, corr_tmp, assigned_idx, frontier%sel_dists_sorted, frontier%inds_sorted)
                call commit_selected_assignment()
            enddo
        end subroutine assign_particles_for_state

        subroutine commit_selected_assignment()
            assigned_iref = frontier%sel_refs(assigned_idx)
            idx = graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1
            assigned_ptcl = graph%ref_list(idx)
            frontier%ptcl_avail(assigned_ptcl) = .false.
            final_assigned(assigned_ptcl) = .true.
            nleft = nleft - 1
            call self%assign_candidate(assigned_ptcl, self%candidate_store%candidates(graph%ref_candidate_pos(idx)))
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
            self%assgn_map(assigned_ptcl)%frac = search_frac(assigned_ptcl)
            ! Report the distinct valid references used for this posterior-guided neighborhood.
            if( l_report_npeaks ) self%assgn_map(assigned_ptcl)%npeaks = particle_evaluated_count(assigned_ptcl)
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
                    fallback_ref = pick_best_evaluated_ref_for_state(i,greedy_state(i))
                else
                    fallback_ref = pick_best_evaluated_ref(i)
                endif
                if( fallback_ref == 0 )then
                    call seed_fallback_if_empty(i)
                    if( self%nstates > 1 .and. greedy_state(i) > 0 )then
                        fallback_ref = pick_best_evaluated_ref_for_state(i,greedy_state(i))
                    else
                        fallback_ref = pick_best_evaluated_ref(i)
                    endif
                endif
                if( fallback_ref == 0 ) fallback_ref = pick_best_evaluated_ref(i)
                if( fallback_ref == 0 ) fallback_ref = 1
                candidate_pos = find_candidate_position(i, fallback_ref)
                if( candidate_pos < 1 ) THROW_HARD('missing fallback candidate in sparse assignment')
                call self%assign_candidate(i, self%candidate_store%candidates(candidate_pos))
                call materialize_seed_shift(self%assgn_map(i), self%seed_shifts(:,i),&
                    &self%seed_has_sh(i), self%p_ptr%l_doshift, self%seed_nrots)
                self%assgn_map(i)%frac = search_frac(i)
                if( l_report_npeaks ) self%assgn_map(i)%npeaks = particle_evaluated_count(i)
                final_assigned(i) = .true.
            enddo
        end subroutine assign_remaining_particles_from_best_touched_ref

        real function search_frac(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            if( self%l_direct_stoch_neigh .and. self%nrefs > 0 )then
                search_frac = 100.0 * real(particle_evaluated_count(iptcl_loc)) / real(self%nrefs)
            else
                search_frac = 100.0
            endif
        end function search_frac

        integer function pick_best_evaluated_ref(iptcl_loc) result(ri_best)
            integer, intent(in) :: iptcl_loc
            ri_best = pick_best_evaluated_ref_filtered(iptcl_loc,0)
        end function pick_best_evaluated_ref

        integer function pick_best_evaluated_ref_for_state(iptcl_loc,state_filter) result(ri_best)
            integer, intent(in) :: iptcl_loc, state_filter
            ri_best = pick_best_evaluated_ref_filtered(iptcl_loc,state_filter)
        end function pick_best_evaluated_ref_for_state

        integer function pick_best_evaluated_ref_filtered(iptcl_loc,state_filter) result(ri_best)
            integer, intent(in) :: iptcl_loc
            integer, intent(in) :: state_filter
            integer :: pos_loc, ri_loc, first_loc, last_loc
            real    :: best_dist
            ri_best = 0
            best_dist = huge(1.0)
            first_loc = int(self%candidate_store%offsets(iptcl_loc))
            last_loc  = int(self%candidate_store%offsets(iptcl_loc+1)) - 1
            do pos_loc = first_loc,last_loc
                ri_loc = self%full_to_compact_ref(self%candidate_store%candidates(pos_loc)%iref)
                if( ri_loc < 1 .or. ri_loc > self%nrefs ) cycle
                if( self%candidate_store%candidates(pos_loc)%inpl <= 0 ) cycle
                if( state_filter > 0 .and. self%ref_state(ri_loc) /= state_filter ) cycle
                if( self%candidate_store%candidates(pos_loc)%dist < best_dist )then
                    best_dist = self%candidate_store%candidates(pos_loc)%dist
                    ri_best   = ri_loc
                endif
            enddo
        end function pick_best_evaluated_ref_filtered

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

        integer function particle_evaluated_count(iptcl_loc) result(nevaluated)
            integer, intent(in) :: iptcl_loc
            integer :: pos_loc
            nevaluated = 0
            do pos_loc = int(self%candidate_store%offsets(iptcl_loc)),&
                &int(self%candidate_store%offsets(iptcl_loc+1))-1
                if( self%candidate_store%candidates(pos_loc)%inpl > 0 ) nevaluated = nevaluated + 1
            enddo
        end function particle_evaluated_count

        subroutine seed_fallback_if_empty(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            type(ori) :: o_prev
            integer :: istate_loc, iproj_loc, irot_loc, fallback_state, fallback_proj, fallback_ref_full
            integer :: ri_loc, fallback_ref_loc
            real    :: inpl_athres_state, sh_seed(2), dist, corr_loc
            if( pick_best_evaluated_ref(iptcl_loc) > 0 ) return
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
            if( fallback_ref_loc == 0 ) fallback_ref_loc = 1
            candidate_pos = find_candidate_position(iptcl_loc, fallback_ref_loc)
            if( candidate_pos < 1 ) candidate_pos = int(self%candidate_store%offsets(iptcl_loc))
            if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 )then
                fallback_state    = self%ref_state(fallback_ref_loc)
                fallback_proj     = self%ref_proj(fallback_ref_loc)
                fallback_ref_full = (fallback_state-1)*self%p_ptr%nspace + fallback_proj
                sh_seed = 0.
                if( l_seed_fallback_shift ) sh_seed = o_prev%get_2Dshift()
                if( self%p_ptr%l_prob_inpl )then
                    inpl_athres_state = calc_athres(self%b_ptr%spproj_field, 'dist_inpl',&
                        &self%p_ptr%prob_athres, state=fallback_state)
                    call self%b_ptr%pftc%gen_likelihood_val(fallback_ref_full, self%pinds(iptcl_loc),&
                        &sh_seed, likelihood_nsample(self%b_ptr%pftc%get_nrots(), inpl_athres_state),&
                        &dist, corr_loc, irot_loc, dists_inpl_sorted, inds_sorted)
                else
                    call self%b_ptr%pftc%gen_likelihood_val(fallback_ref_full, self%pinds(iptcl_loc),&
                        &sh_seed, self%b_ptr%pftc%get_nrots(), dist, corr_loc, irot_loc, dists_inpl_sorted, inds_sorted)
                endif
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
                    if( greedy_state(cand) == self%ref_state(iref) )then
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

        integer function likelihood_nsample( n, athres_ub_in ) result(nsample)
            integer, intent(in) :: n
            real,    intent(in) :: athres_ub_in
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, athres_ub_in)
            nsample = min(n, max(1, int(athres_ub * real(n) / 180.)))
        end function likelihood_nsample

        real function selected_frontier_dist( isel ) result(dist)
            integer, intent(in) :: isel
            dist = frontier%sel_dists(isel)
        end function selected_frontier_dist

    end subroutine ref_assign_neigh

    ! Table & Assignment I/O

    ! Reads one partition candidate stream into the compact global store.
    subroutine read_sparse_tab_to_glob( self, binfname )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: binfname
        integer, parameter :: IO_CHUNK = 65536
        type(prob_candidate), allocatable :: candidates_loc(:)
        real,           allocatable :: seed_shifts_loc(:,:)
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
        if( nnz < 1 ) THROW_HARD('empty sparse table in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        read(funit, pos=addr) pinds_loc
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab_neigh%read_sparse_tab_to_glob')
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        if( max_pind < 1 )then
            call fclose(funit)
            deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, pind2glob)
            return
        endif
        do i_loc = 1,nptcls_loc
            pind = pinds_loc(i_loc)
            i_glob = 0
            if( pind >= 1 .and. pind <= max_pind ) i_glob = pind2glob(pind)
            if( i_glob > 0 )then
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
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
        deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, particle_indices, candidates_loc, pind2glob)
    end subroutine read_sparse_tab_to_glob

    ! Aggregates nparts tables into the global table
    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(string) :: fname
        integer, parameter :: IO_CHUNK = 65536
        integer, allocatable :: counts_glob(:), particle_indices(:), pinds_loc(:), pind2glob(:)
        real,    allocatable :: seed_shifts_loc(:,:)
        logical, allocatable :: seed_has_sh_loc(:)
        type(prob_candidate) :: candidate_dummy
        integer :: ipart, funit, io_stat, nptcls_loc, nchunks, ichunk, chunk_n
        integer :: i, j, first, last, nread
        integer(int64) :: header(4), nnz, nnz_read, addr, chunk_indices_addr, chunk_candidates_addr
        integer :: pind, i_glob, max_pind, seed_nrots_loc
        allocate(counts_glob(self%nptcls), source=0)
        do ipart = 1,nparts
                fname = fbody//int2str_pad(ipart,numlen)//'.dat'
                call fopen(funit,fname,access='STREAM',action='READ',status='OLD',iostat=io_stat)
                call fileiochk('simple_eul_prob_tab_neigh; read_tabs_to_glob header; file: '//fname%to_char(), io_stat)
                read(funit,pos=1) header
                if( header(1) /= self%nrefs ) THROW_HARD('nrefs mismatch in read_tabs_to_glob')
                nptcls_loc = int(header(2))
                nnz = header(3)
                nchunks = int(header(4))
                allocate(pinds_loc(nptcls_loc))
                allocate(seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
                addr = sizeof(header) + 1
                read(funit,pos=addr) pinds_loc
                addr = addr + sizeof(pinds_loc)
                call read_seed_shift_table(funit,addr,seed_nrots_loc,seed_shifts_loc,seed_has_sh_loc)
                allocate(particle_indices(IO_CHUNK))
                call build_pind_lookup(self%pinds,pinds_loc,pind2glob,max_pind)
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
                        do j = 1,nread
                            i = particle_indices(j)
                            if( i < 1 .or. i > nptcls_loc ) THROW_HARD('invalid particle index in sparse candidate stream')
                            pind = pinds_loc(i)
                            if( pind < 1 .or. pind > max_pind ) cycle
                            i_glob = pind2glob(pind)
                            if( i_glob > 0 ) counts_glob(i_glob) = counts_glob(i_glob) + 1
                        enddo
                    enddo
                    addr = chunk_candidates_addr + chunk_n * sizeof(candidate_dummy)
                    nnz_read = nnz_read + int(chunk_n,int64)
                enddo
                if( nnz_read /= nnz ) THROW_HARD('candidate stream count mismatch')
                call fclose(funit)
                deallocate(pinds_loc,particle_indices,seed_shifts_loc,seed_has_sh_loc,pind2glob)
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
                    refined_shift = grad_obj(ithr)%minimize(irot=irot, sh_rot=.true.,&
                        &xy_in=shift_seed(2:3))
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

    subroutine kill_posterior3d_ws( self )
        class(posterior3d_ws), intent(inout) :: self
        integer :: ithr
        if( allocated(self%readers) )then
            do ithr = 1, size(self%readers)
                call self%readers(ithr)%kill
            enddo
            deallocate(self%readers)
        endif
        if( allocated(self%pind_lookup) ) deallocate(self%pind_lookup)
        if( allocated(self%posterior_available) ) deallocate(self%posterior_available)
        if( allocated(self%sub_offsets) ) deallocate(self%sub_offsets)
        if( allocated(self%sub_full_inds) ) deallocate(self%sub_full_inds)
        if( allocated(self%source_to_target_sub) ) deallocate(self%source_to_target_sub)
        call self%target_subspace%kill
        self%source_nspace = 0
        self%coverage_count = 0
        self%initialized = .false.
        self%valid = .false.
        self%warning_emitted = .false.
    end subroutine kill_posterior3d_ws

    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%candidate_store%kill
        if( allocated(self%candidate_fill_counts) ) deallocate(self%candidate_fill_counts)
        call self%posterior_ws%kill
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
