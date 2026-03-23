! Refactored direct-CSR version of simple_eul_prob_tab_neigh.
!
! Main changes:
!   * neighborhood builders now emit sparse rows directly
!   * no global logical mask(:,:) is built
!   * non-ptree neighborhood pooling uses disjoint subspace blocks
!   * shift refinement runs internally on already-scored edges
!
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_eul_prob_tab,       only: angle_sampling, calc_num2sample, calc_athres, eulprob_dist_switch
use simple_eulspace_neigh_map, only: eulspace_neigh_map
implicit none
private
#include "simple_local_flags.inc"

public :: eul_prob_tab_neigh

type :: eul_prob_tab_neigh
    integer,        allocatable :: pinds(:)
    type(ptcl_ref), allocatable :: assgn_map(:)
    logical,        allocatable :: proj_exists(:,:)
    logical,        allocatable :: state_exists(:)
    integer,        allocatable :: active_state_indices(:)
    integer,        allocatable :: ref_proj_indices(:)
    integer,        allocatable :: ref_state_indices(:)
    integer,        allocatable :: ref_index_map(:,:)
    integer,        allocatable :: proj_active_state_count(:)
    integer,        allocatable :: ptcl_off(:)
    integer,        allocatable :: edge_ref_index(:)
    integer,        allocatable :: edge_ptcl(:)
    type(ptcl_ref), allocatable :: edge_val(:)
    integer,        allocatable :: ref_edge_offsets(:)
    integer,        allocatable :: ref_edge_indices(:)
    class(builder),     pointer :: b_ptr => null()
    class(parameters),  pointer :: p_ptr => null()
    integer                     :: nptcls = 0
    integer                     :: nstates = 0
    integer                     :: nrefs   = 0
    integer                     :: nedges      = 0
    integer                     :: maxdeg_ptcl = 0
contains
    procedure          :: build_sparse_neigh_graph
    procedure, private :: init_common
    procedure, private :: build_ref_lists_and_map
    procedure, private :: check_subspace_prereqs
    procedure, private :: init_edge_default
    procedure, private :: materialize_row_to_csr
    procedure, private :: finalize_rows_to_csr_geom
    procedure, private :: finalize_rows_to_csr_subspace
    generic,   private :: finalize_rows_to_csr => finalize_rows_to_csr_geom, finalize_rows_to_csr_subspace
    procedure, private :: build_graph_from_prev_geom
    procedure, private :: build_graph_from_subspace_peaks
    procedure, private :: build_graph_from_subspace_peaks_one_state
    procedure, private :: build_graph_from_subspace_peaks_impl
    procedure, private :: build_ref_adjacency
    procedure, private :: refine_shift_edges
    procedure, private :: ref_normalize_sparse
    procedure, private :: sort_ref_lists_by_dist
    procedure, private :: sort_eidx_by_dist
    procedure, private :: sift_down
    procedure          :: ref_assign => ref_assign_sparse
    procedure          :: write_tab
    procedure          :: read_tabs_to_glob
    procedure          :: write_assignment
    procedure          :: read_assignment
    procedure          :: kill
end type eul_prob_tab_neigh

type :: proj_sel_tmp
    integer, allocatable :: proj_idx(:)
    real                 :: cxy_shift(2) = [0., 0.]
    logical              :: shifts_valid = .false.
end type proj_sel_tmp

type :: geom_thread_buf
    real,           allocatable :: sub_dists(:)
    real,           allocatable :: score_buf(:)
    real,           allocatable :: sorted_dist_scratch(:)
    integer,        allocatable :: sorted_idx_scratch(:)
    integer,        allocatable :: peak_sub_idxs(:)
    integer,        allocatable :: proj_sel(:)
end type geom_thread_buf

type :: subspace_thread_buf
    real,           allocatable :: score_buf(:)
    real,           allocatable :: sorted_dist_scratch(:)
    integer,        allocatable :: sorted_idx_scratch(:)
    integer,        allocatable :: proj_sel(:)
end type subspace_thread_buf

! Benchmark timings for route-internal steps (set by builder implementations).
real(timer_int_kind) :: rt_route_setup = 0.
real(timer_int_kind) :: rt_route_shift_init = 0.
real(timer_int_kind) :: rt_route_particle_loop = 0.
real(timer_int_kind) :: rt_route_finalize = 0.

contains

    ! Initialize object state and optionally build the sparse neighborhood graph.
    subroutine build_sparse_neigh_graph(self, params, build, pinds, neigh_type, empty_okay, state, build_sparse_graph)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        character(len=*),          intent(in)    :: neigh_type
        logical, optional,         intent(in)    :: empty_okay
        integer, optional,         intent(in)    :: state
        logical, optional,         intent(in)    :: build_sparse_graph
        type(string)            :: benchfname
        logical                 :: do_build_sparse_graph
        integer(timer_int_kind) :: t_route_tot, t_build_ref_adjacency, t_refine_shift_edges
        integer                 :: fnr
        real(timer_int_kind)    :: rt_route_build
        real(timer_int_kind)    :: rt_build_ref_adjacency, rt_refine_shift_edges, rt_route_tot
        do_build_sparse_graph = .true.
        if (present(build_sparse_graph)) do_build_sparse_graph = build_sparse_graph
        rt_route_build         = 0.
        rt_route_setup         = 0.
        rt_route_shift_init    = 0.
        rt_route_particle_loop = 0.
        rt_route_finalize      = 0.
        rt_build_ref_adjacency = 0.
        rt_refine_shift_edges  = 0.
        rt_route_tot           = 0.
        call self%init_common(params, build, pinds, empty_okay)
        call self%build_ref_lists_and_map
        if (do_build_sparse_graph) then
            call self%check_subspace_prereqs
            if (L_BENCH_GLOB) t_route_tot = tic()
            select case(trim(neigh_type))
                case('geom')
                    call self%build_graph_from_prev_geom
                case('subspace_srch')
                    if (present(state)) then
                        call self%build_graph_from_subspace_peaks_impl(state)
                    else
                        call self%build_graph_from_subspace_peaks_impl()
                    endif
                case default
                    THROW_HARD('build_sparse_neigh_graph: unsupported neigh_type='//trim(neigh_type)//'; expected geom or subspace_srch')
            end select
            if (L_BENCH_GLOB) then
                rt_route_build = toc(t_route_tot)
                t_build_ref_adjacency = tic()
            endif
            call self%build_ref_adjacency
            if (L_BENCH_GLOB) then
                rt_build_ref_adjacency = toc(t_build_ref_adjacency)
                t_refine_shift_edges = tic()
            endif
            call self%refine_shift_edges
            if (L_BENCH_GLOB) then
                rt_refine_shift_edges = toc(t_refine_shift_edges)
                rt_route_tot = rt_route_build + rt_build_ref_adjacency + rt_refine_shift_edges
                if (rt_route_tot > TINY) then
                    benchfname = 'PROB_TAB_NEIGH_'//trim(neigh_type)//'_BENCH_PART'//int2str_pad(self%p_ptr%part, self%p_ptr%numlen)//'.txt'
                    call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                    write(fnr,'(a)') '*** TIMINGS (s) ***'
                    write(fnr,'(a,1x,f9.3)') trim(neigh_type)//' route total      : ', rt_route_build
                    write(fnr,'(a,1x,f9.3)') '  setup                   : ', rt_route_setup
                    write(fnr,'(a,1x,f9.3)') '  shift init              : ', rt_route_shift_init
                    write(fnr,'(a,1x,f9.3)') '  particle loop           : ', rt_route_particle_loop
                    write(fnr,'(a,1x,f9.3)') '  finalize                : ', rt_route_finalize
                    write(fnr,'(a,1x,f9.3)') 'build_ref_adjacency       : ', rt_build_ref_adjacency
                    write(fnr,'(a,1x,f9.3)') 'refine_shift_edges        : ', rt_refine_shift_edges
                    write(fnr,'(a,1x,f9.3)') 'route total               : ', rt_route_tot
                    write(fnr,'(a)') ''
                    write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                    write(fnr,'(a,1x,f9.2)') trim(neigh_type)//' route total      : ', (rt_route_build/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') '  setup                   : ', (rt_route_setup/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') '  shift init              : ', (rt_route_shift_init/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') '  particle loop           : ', (rt_route_particle_loop/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') '  finalize                : ', (rt_route_finalize/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') 'build_ref_adjacency       : ', (rt_build_ref_adjacency/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') 'refine_shift_edges        : ', (rt_refine_shift_edges/rt_route_tot) * 100.
                    write(fnr,'(a,1x,f9.2)') '% accounted for           : ', &
                        ((rt_route_setup+rt_route_shift_init+rt_route_particle_loop+rt_route_finalize+ &
                          rt_build_ref_adjacency+rt_refine_shift_edges)/rt_route_tot) * 100.
                    call fclose(fnr)
                endif
            endif
        endif
    end subroutine build_sparse_neigh_graph

    ! Initialize common per-particle and per-reference metadata used by all workflows.
    subroutine init_common(self, params, build, pinds, empty_okay)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer, parameter :: MIN_POP = 5
        logical :: l_empty
        integer :: istate
        call self%kill
        self%p_ptr => params
        self%b_ptr => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        self%assgn_map(:)%pind   = self%pinds(:)
        self%assgn_map(:)%istate = 0
        self%assgn_map(:)%iproj  = 0
        self%assgn_map(:)%inpl   = 0
        self%assgn_map(:)%dist   = huge(1.)
        self%assgn_map(:)%x      = 0.
        self%assgn_map(:)%y      = 0.
        self%assgn_map(:)%has_sh = .false.
        l_empty = (trim(params%empty3Dcavgs) .eq. 'yes')
        if (present(empty_okay)) l_empty = empty_okay
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        if (self%nstates == 0) then
            THROW_HARD('No valid states available after state existence filtering')
        endif
        if (l_empty) then
            allocate(self%proj_exists(self%p_ptr%nspace, self%p_ptr%nstates), source=.true.)
        else
            self%proj_exists = self%b_ptr%spproj_field%projs_exist(self%p_ptr%nstates, self%p_ptr%nspace, thres=MIN_POP)
        endif
        self%nrefs = 0
        do istate = 1, self%p_ptr%nstates
            if (.not. self%state_exists(istate)) cycle
            self%nrefs = self%nrefs + count(self%proj_exists(:, istate))
        enddo
        if (self%nrefs == 0) then
            THROW_HARD('No valid references available after state/projection existence filtering')
        endif
    end subroutine init_common

    ! Build compressed reference lists and (iproj,istate)->ref index lookup map.
    subroutine build_ref_lists_and_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: state_ref_counts(:), state_ref_offsets(:)
        integer :: cnt, iproj, istate, local_ref_idx, state_i, state_list_idx
        if (self%nstates <= 0) then
            THROW_HARD('build_ref_lists_and_map: no active states available')
        endif
        if (self%nrefs <= 0) then
            THROW_HARD('build_ref_lists_and_map: no active references available')
        endif
        allocate(self%active_state_indices(self%nstates), self%ref_proj_indices(self%nrefs), self%ref_state_indices(self%nrefs))
        allocate(self%ref_index_map(self%p_ptr%nspace, self%p_ptr%nstates), self%proj_active_state_count(self%p_ptr%nspace), source=0)
        state_list_idx = 0
        do istate = 1, self%p_ptr%nstates
            if (.not. self%state_exists(istate)) cycle
            state_list_idx = state_list_idx + 1
            self%active_state_indices(state_list_idx) = istate
        enddo
        allocate(state_ref_counts(self%nstates), state_ref_offsets(self%nstates), source=0)
        !$omp parallel do default(shared) private(state_i,istate,iproj,cnt) &
        !$omp proc_bind(close) schedule(static)
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            cnt = 0
            do iproj = 1, self%p_ptr%nspace
                if (self%proj_exists(iproj, istate)) cnt = cnt + 1
            enddo
            state_ref_counts(state_i) = cnt
        enddo
        !$omp end parallel do
        state_ref_offsets(1) = 1
        do state_i = 2, self%nstates
            state_ref_offsets(state_i) = state_ref_offsets(state_i-1) + state_ref_counts(state_i-1)
        enddo
        if (state_ref_offsets(self%nstates) + state_ref_counts(self%nstates) - 1 /= self%nrefs) then
            THROW_HARD('build_ref_lists_and_map: reference count mismatch while building sparse mapping')
        endif
        !$omp parallel default(shared) private(state_i,istate,iproj,local_ref_idx,cnt) proc_bind(close) 
        !$omp do schedule(static)
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            local_ref_idx = state_ref_offsets(state_i)
            do iproj = 1, self%p_ptr%nspace
                if (.not. self%proj_exists(iproj, istate)) cycle
                self%ref_proj_indices(local_ref_idx) = iproj
                self%ref_state_indices(local_ref_idx) = istate
                self%ref_index_map(iproj, istate) = local_ref_idx
                local_ref_idx = local_ref_idx + 1
            enddo
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do iproj = 1, self%p_ptr%nspace
            cnt = 0
            do state_i = 1, self%nstates
                istate = self%active_state_indices(state_i)
                if (self%proj_exists(iproj, istate)) cnt = cnt + 1
            enddo
            self%proj_active_state_count(iproj) = cnt
        enddo
        !$omp end do nowait
        !$omp end parallel
        deallocate(state_ref_counts, state_ref_offsets)
    end subroutine build_ref_lists_and_map

    ! Validate required subspace mappings before any neighborhood graph construction.
    subroutine check_subspace_prereqs(self)
        class(eul_prob_tab_neigh), intent(in) :: self
        if (.not. allocated(self%b_ptr%subspace_inds)) then
            THROW_HARD('check_subspace_prereqs: subspace_inds not allocated; enable l_neigh and set nspace_sub')
        endif
        if (size(self%b_ptr%subspace_inds) /= self%p_ptr%nspace_sub) then
            THROW_HARD('check_subspace_prereqs: size(subspace_inds) must equal nspace_sub')
        endif
        if (.not. allocated(self%b_ptr%subspace_full2sub_map)) then
            THROW_HARD('check_subspace_prereqs: subspace_full2sub_map not allocated')
        endif
        if (size(self%b_ptr%subspace_full2sub_map) /= self%p_ptr%nspace) then
            THROW_HARD('check_subspace_prereqs: size(subspace_full2sub_map) must equal nspace')
        endif
    end subroutine check_subspace_prereqs

    ! Fill a ptcl_ref with default values for one particle/reference edge.
    subroutine init_edge_default(self, iptcl, ref_idx, v)
        class(eul_prob_tab_neigh), intent(in)  :: self
        integer,                   intent(in)  :: iptcl, ref_idx
        type(ptcl_ref),            intent(out) :: v
        v = ptcl_ref( &
            pind   = iptcl, &
            istate = self%ref_state_indices(ref_idx), &
            iproj  = self%ref_proj_indices(ref_idx), &
            inpl   = 0, &
            dist   = huge(1.), &
            x      = 0., &
            y      = 0., &
            has_sh = .false. )
    end subroutine init_edge_default

    ! Score one particle's selected projections and write directly into CSR edge arrays.
    subroutine materialize_row_to_csr(self, iptcl, ptcl_local_idx, proj_idx, nproj, edge_first, cxy_shift, shifts_valid, &
                                      score_buf, sorted_dist_scratch, sorted_idx_scratch, inpl_angle_thres)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: iptcl, ptcl_local_idx, nproj, edge_first
        integer,                   intent(in)    :: proj_idx(:)
        real,                      intent(in)    :: cxy_shift(2)
        logical,                   intent(in)    :: shifts_valid
        real,                      intent(inout) :: score_buf(:), sorted_dist_scratch(:)
        integer,                   intent(inout) :: sorted_idx_scratch(:)
        real,                      intent(in)    :: inpl_angle_thres(:)
        integer :: i, j, iproj, istate, ref_idx, irot, iref_full, e, nrow_expect
        real    :: rotmat_loc(2,2), rot_xy_loc(2)
        nrow_expect = 0
        do i = 1, nproj
            iproj = proj_idx(i)
            nrow_expect = nrow_expect + self%proj_active_state_count(iproj)
        enddo
        if (nrow_expect <= 0) then
            THROW_HARD('materialize_row_to_csr: zero-row materialization')
        endif
        e = edge_first - 1
        do i = 1, nproj
            iproj = proj_idx(i)
            if (self%proj_active_state_count(iproj) <= 0) cycle
            do j = 1, self%nstates
                istate = self%active_state_indices(j)
                ref_idx = self%ref_index_map(iproj, istate)
                if (ref_idx <= 0) cycle
                e = e + 1
                self%edge_ref_index(e) = ref_idx
                self%edge_ptcl(e)      = ptcl_local_idx
                iref_full = (istate - 1) * self%p_ptr%nspace + iproj
                call self%b_ptr%pftc%gen_objfun_vals(iref_full, iptcl, cxy_shift, score_buf)
                score_buf = eulprob_dist_switch(score_buf, self%p_ptr%cc_objfun)
                irot = angle_sampling(score_buf, sorted_dist_scratch, sorted_idx_scratch, inpl_angle_thres(istate), self%p_ptr%prob_athres)
                call self%init_edge_default(iptcl, ref_idx, self%edge_val(e))
                self%edge_val(e)%dist = score_buf(irot)
                self%edge_val(e)%inpl = irot
                if (shifts_valid) then
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat_loc)
                    rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                    self%edge_val(e)%x = rot_xy_loc(1)
                    self%edge_val(e)%y = rot_xy_loc(2)
                    self%edge_val(e)%has_sh = .true.
                endif
            enddo
        enddo
        if (e /= edge_first + nrow_expect - 1) then
            THROW_HARD('materialize_row_to_csr: edge fill mismatch')
        endif
    end subroutine materialize_row_to_csr

    ! Finalize staged row selections into CSR storage for the geom route.
    subroutine finalize_rows_to_csr_geom(self, sel_rows, row_nedges, tbuf, inpl_angle_thres)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(proj_sel_tmp),        intent(in)    :: sel_rows(:)
        integer,                   intent(in)    :: row_nedges(:)
        type(geom_thread_buf),     intent(inout) :: tbuf(:)
        real,                      intent(in)    :: inpl_angle_thres(:)
        integer :: i, iptcl, ithr, nproj
        if (allocated(self%ptcl_off))       deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index)) deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))      deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))       deallocate(self%edge_val)
        allocate(self%ptcl_off(self%nptcls + 1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            if (row_nedges(i) <= 0) then
                THROW_HARD('finalize_rows_to_csr_geom: particle has zero valid neighbors')
            endif
            self%ptcl_off(i+1) = self%ptcl_off(i) + row_nedges(i)
        enddo
        self%nedges = self%ptcl_off(self%nptcls + 1) - 1
        self%maxdeg_ptcl = maxval(row_nedges)
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        !$omp parallel do default(shared) private(i,iptcl,ithr,nproj) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            nproj = size(sel_rows(i)%proj_idx)
            call self%materialize_row_to_csr(iptcl, i, sel_rows(i)%proj_idx, nproj, self%ptcl_off(i), sel_rows(i)%cxy_shift, &
                                             sel_rows(i)%shifts_valid, tbuf(ithr)%score_buf, tbuf(ithr)%sorted_dist_scratch, &
                                             tbuf(ithr)%sorted_idx_scratch, inpl_angle_thres)
        enddo
        !$omp end parallel do
    end subroutine finalize_rows_to_csr_geom

    ! Finalize staged row selections into CSR storage for the subspace route.
    subroutine finalize_rows_to_csr_subspace(self, sel_rows, row_nedges, tbuf, inpl_angle_thres)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(proj_sel_tmp),        intent(in)    :: sel_rows(:)
        integer,                   intent(in)    :: row_nedges(:)
        type(subspace_thread_buf), intent(inout) :: tbuf(:)
        real,                      intent(in)    :: inpl_angle_thres(:)
        integer :: i, iptcl, ithr, nproj
        if (allocated(self%ptcl_off))       deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index)) deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))      deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))       deallocate(self%edge_val)
        allocate(self%ptcl_off(self%nptcls + 1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            if (row_nedges(i) <= 0) then
                THROW_HARD('finalize_rows_to_csr_subspace: particle has zero valid neighbors')
            endif
            self%ptcl_off(i+1) = self%ptcl_off(i) + row_nedges(i)
        enddo
        self%nedges = self%ptcl_off(self%nptcls + 1) - 1
        self%maxdeg_ptcl = maxval(row_nedges)
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        !$omp parallel do default(shared) private(i,iptcl,ithr,nproj) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            nproj = size(sel_rows(i)%proj_idx)
            call self%materialize_row_to_csr(iptcl, i, sel_rows(i)%proj_idx, nproj, self%ptcl_off(i), sel_rows(i)%cxy_shift, &
                                             sel_rows(i)%shifts_valid, tbuf(ithr)%score_buf, tbuf(ithr)%sorted_dist_scratch, &
                                             tbuf(ithr)%sorted_idx_scratch, inpl_angle_thres)
        enddo
        !$omp end parallel do
    end subroutine finalize_rows_to_csr_subspace

    ! Release per-particle staged projection selections.
    subroutine cleanup_staged_rows(sel_rows)
        type(proj_sel_tmp), allocatable, intent(inout) :: sel_rows(:)
        integer :: i
        if (.not. allocated(sel_rows)) return
        do i = 1, size(sel_rows)
            if (allocated(sel_rows(i)%proj_idx)) deallocate(sel_rows(i)%proj_idx)
        enddo
        deallocate(sel_rows)
    end subroutine cleanup_staged_rows

    ! Release temporary buffers used by the geom builder.
    subroutine cleanup_staged_buffers_geom(row_nedges, inpl_angle_thres, tbuf)
        integer,               allocatable, intent(inout) :: row_nedges(:)
        real,                  allocatable, intent(inout) :: inpl_angle_thres(:)
        type(geom_thread_buf), allocatable, intent(inout) :: tbuf(:)
        integer :: ithr
        if (allocated(tbuf)) then
            do ithr = 1, size(tbuf)
                if (allocated(tbuf(ithr)%sub_dists))            deallocate(tbuf(ithr)%sub_dists)
                if (allocated(tbuf(ithr)%peak_sub_idxs))        deallocate(tbuf(ithr)%peak_sub_idxs)
                if (allocated(tbuf(ithr)%score_buf))            deallocate(tbuf(ithr)%score_buf)
                if (allocated(tbuf(ithr)%sorted_dist_scratch))  deallocate(tbuf(ithr)%sorted_dist_scratch)
                if (allocated(tbuf(ithr)%sorted_idx_scratch))   deallocate(tbuf(ithr)%sorted_idx_scratch)
                if (allocated(tbuf(ithr)%proj_sel))             deallocate(tbuf(ithr)%proj_sel)
            enddo
            deallocate(tbuf)
        endif
        if (allocated(row_nedges))       deallocate(row_nedges)
        if (allocated(inpl_angle_thres)) deallocate(inpl_angle_thres)
    end subroutine cleanup_staged_buffers_geom

    ! Release temporary buffers used by the subspace builder.
    subroutine cleanup_staged_buffers_subspace(inpl_dists, coarse_best_dist, peak_sub_idxs, row_nedges, inpl_angle_thres, tbuf)
        real,                     allocatable, intent(inout) :: inpl_dists(:,:), coarse_best_dist(:,:)
        integer,                  allocatable, intent(inout) :: peak_sub_idxs(:,:), row_nedges(:)
        real,                     allocatable, intent(inout) :: inpl_angle_thres(:)
        type(subspace_thread_buf), allocatable, intent(inout) :: tbuf(:)
        integer :: ithr
        if (allocated(inpl_dists))       deallocate(inpl_dists)
        if (allocated(coarse_best_dist)) deallocate(coarse_best_dist)
        if (allocated(peak_sub_idxs))    deallocate(peak_sub_idxs)
        if (allocated(tbuf)) then
            do ithr = 1, size(tbuf)
                if (allocated(tbuf(ithr)%score_buf))            deallocate(tbuf(ithr)%score_buf)
                if (allocated(tbuf(ithr)%sorted_dist_scratch))  deallocate(tbuf(ithr)%sorted_dist_scratch)
                if (allocated(tbuf(ithr)%sorted_idx_scratch))   deallocate(tbuf(ithr)%sorted_idx_scratch)
                if (allocated(tbuf(ithr)%proj_sel))             deallocate(tbuf(ithr)%proj_sel)
            enddo
            deallocate(tbuf)
        endif
        if (allocated(row_nedges))       deallocate(row_nedges)
        if (allocated(inpl_angle_thres)) deallocate(inpl_angle_thres)
    end subroutine cleanup_staged_buffers_subspace

    ! Build graph from geometric neighbors around previous particle orientations.
    subroutine build_graph_from_prev_geom(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        real, allocatable             :: inpl_angle_thres(:)
        integer, allocatable          :: neigh_proj(:), row_nedges(:)
        type(proj_sel_tmp), allocatable :: sel_rows(:)
        type(geom_thread_buf), allocatable :: tbuf(:)
        type(eulspace_neigh_map)      :: neigh_map
        type(ori)                     :: o_prev, o_sub, osym
        integer :: i, iptcl, iproj_full, isub, istate, ithr, j, npeak_use, nproj
        integer :: nspace_sub, nrots, prev_iproj, prev_sub, state_i
        real    :: dtmp, inplrotdist
        nspace_sub = self%p_ptr%nspace_sub
        npeak_use  = max(1, min(self%p_ptr%npeaks, nspace_sub))
        nrots      = self%b_ptr%pftc%get_nrots()
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(sel_rows(self%nptcls))
        allocate(row_nedges(self%nptcls), source=0)
        allocate(inpl_angle_thres(self%p_ptr%nstates), source=0.)
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            inpl_angle_thres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        allocate(tbuf(nthr_glob))
        do ithr = 1, nthr_glob
            allocate(tbuf(ithr)%sub_dists(nspace_sub), tbuf(ithr)%peak_sub_idxs(npeak_use), &
                     tbuf(ithr)%score_buf(nrots), tbuf(ithr)%sorted_dist_scratch(nrots))
            allocate(tbuf(ithr)%sorted_idx_scratch(nrots), tbuf(ithr)%proj_sel(self%p_ptr%nspace), source=0)
        enddo
        !$omp parallel do default(shared) &
        !$omp private(i,iptcl,ithr,o_prev,isub,iproj_full,dtmp,inplrotdist,o_sub,osym, &
        !$omp& prev_iproj,prev_sub,j,nproj,neigh_proj) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            tbuf(ithr)%sub_dists = huge(1.)
            call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
            do isub = 1, nspace_sub
                iproj_full = self%b_ptr%subspace_inds(isub)
                call self%b_ptr%eulspace%get_ori(iproj_full, o_sub)
                call self%b_ptr%pgrpsyms%sym_dists(o_prev, o_sub, osym, dtmp, inplrotdist)
                tbuf(ithr)%sub_dists(isub) = dtmp
                call o_sub%kill
                call osym%kill
            enddo
            tbuf(ithr)%peak_sub_idxs = minnloc(tbuf(ithr)%sub_dists, npeak_use)
            nproj = 0
            do j = 1, npeak_use
                call neigh_map%get_neighbors_list(tbuf(ithr)%peak_sub_idxs(j), neigh_proj)
                if (allocated(neigh_proj)) then
                    call append_proj_block(tbuf(ithr)%proj_sel, nproj, neigh_proj)
                    deallocate(neigh_proj)
                endif
            enddo
            prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
            prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
            if (.not. any(tbuf(ithr)%peak_sub_idxs(1:npeak_use) == prev_sub)) then
                call neigh_map%get_neighbors_list(prev_sub, neigh_proj)
                if (allocated(neigh_proj)) then
                    call append_proj_block(tbuf(ithr)%proj_sel, nproj, neigh_proj)
                    deallocate(neigh_proj)
                endif
            endif
            if (nproj <= 0) then
                call add_all_active_projs(self, tbuf(ithr)%proj_sel, nproj)
            else if (.not. any(self%proj_active_state_count(tbuf(ithr)%proj_sel(1:nproj)) > 0)) then
                call add_all_active_projs(self, tbuf(ithr)%proj_sel, nproj)
            endif
            if (allocated(sel_rows(i)%proj_idx)) deallocate(sel_rows(i)%proj_idx)
            allocate(sel_rows(i)%proj_idx(nproj), source=tbuf(ithr)%proj_sel(1:nproj))
            sel_rows(i)%cxy_shift = [0., 0.]
            sel_rows(i)%shifts_valid = .false.
            row_nedges(i) = 0
            do j = 1, nproj
                row_nedges(i) = row_nedges(i) + self%proj_active_state_count(sel_rows(i)%proj_idx(j))
            enddo
            call o_prev%kill
        enddo
        !$omp end parallel do
        call neigh_map%kill
        call self%finalize_rows_to_csr(sel_rows, row_nedges, tbuf, inpl_angle_thres)
        call cleanup_staged_rows(sel_rows)
        call cleanup_staged_buffers_geom(row_nedges, inpl_angle_thres, tbuf)
    end subroutine build_graph_from_prev_geom

    ! Multi-state wrapper for subspace-peak driven graph construction.
    subroutine build_graph_from_subspace_peaks(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%build_graph_from_subspace_peaks_impl()
    end subroutine build_graph_from_subspace_peaks

    ! Fixed-state wrapper for subspace-peak driven graph construction.
    subroutine build_graph_from_subspace_peaks_one_state(self, state)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: state
        if (state < 1 .or. state > self%p_ptr%nstates) then
            THROW_HARD('build_graph_from_subspace_peaks_one_state: state index out of range')
        endif
        if (.not. self%state_exists(state)) then
            THROW_HARD('build_graph_from_subspace_peaks_one_state: requested state does not exist')
        endif
        call self%build_graph_from_subspace_peaks_impl(state)
    end subroutine build_graph_from_subspace_peaks_one_state

    ! Build graph from best coarse subspace peaks, with optional fixed-state mode.
    subroutine build_graph_from_subspace_peaks_impl(self, state_opt)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, optional,         intent(in)    :: state_opt
        real,                    allocatable :: coarse_best_dist(:,:), inpl_angle_thres(:), inpl_dists(:,:)
        integer,                 allocatable :: neigh_proj(:), peak_sub_idxs(:,:), row_nedges(:)
        type(proj_sel_tmp), allocatable :: sel_rows(:)
        type(subspace_thread_buf), allocatable :: tbuf(:)
        type(eulspace_neigh_map)      :: neigh_map
        type(pftc_shsrch_grad)        :: grad_shsrch_obj(nthr_glob)
        type(ori)                     :: o_prev
        logical :: do_shift_first, has_state_opt, prev_included
        integer :: i, iref_prev, iref_start, iproj_full, iptcl, irot, isub, istate, ithr
        integer :: j, npeak_use, nproj, nspace_sub, nvalid_sub, nrots, prev_iproj, prev_sub
        integer :: state_fixed, state_i
        integer(timer_int_kind) :: t_setup, t_shift_init, t_particle_loop, t_finalize
        real    :: cxy(3), cxy_shift(2), huge_dist, lims(2,2), lims_init(2,2)
        real(timer_int_kind)    :: rt_setup, rt_shift_init, rt_particle_loop, rt_finalize
        rt_setup = 0.
        rt_shift_init = 0.
        rt_particle_loop = 0.
        rt_finalize = 0.
        if (L_BENCH_GLOB) t_setup = tic()
        call self%check_subspace_prereqs()
        nspace_sub     = self%p_ptr%nspace_sub
        npeak_use      = max(1, min(self%p_ptr%npeaks, nspace_sub))
        nrots          = self%b_ptr%pftc%get_nrots()
        huge_dist      = huge(1.)
        do_shift_first = self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift
        has_state_opt  = present(state_opt)
        if (has_state_opt) then
            state_fixed = state_opt
        else
            state_fixed = 1
        endif
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(sel_rows(self%nptcls))
        allocate(row_nedges(self%nptcls), source=0)
        allocate(inpl_dists(nrots, nthr_glob), coarse_best_dist(nspace_sub, nthr_glob), peak_sub_idxs(npeak_use, nthr_glob), tbuf(nthr_glob))
        allocate(inpl_angle_thres(self%p_ptr%nstates), source=0.)
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            inpl_angle_thres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        do ithr = 1, nthr_glob
            allocate(tbuf(ithr)%score_buf(nrots), tbuf(ithr)%sorted_dist_scratch(nrots))
            allocate(tbuf(ithr)%sorted_idx_scratch(nrots), tbuf(ithr)%proj_sel(self%p_ptr%nspace), source=0)
        enddo
        if (L_BENCH_GLOB) rt_setup = toc(t_setup)
        if (do_shift_first) then
            if (L_BENCH_GLOB) t_shift_init = tic()
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
            if (L_BENCH_GLOB) rt_shift_init = toc(t_shift_init)
        endif
        if (L_BENCH_GLOB) t_particle_loop = tic()
        !$omp parallel do default(shared) &
        !$omp private(i,ithr,iptcl,o_prev,istate,isub,iproj_full,irot,iref_start, &
        !$omp& iref_prev,cxy,cxy_shift,state_i,nvalid_sub,prev_sub,prev_iproj,j,nproj,neigh_proj,prev_included) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
            cxy_shift = [0., 0.]
            prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
            prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
            if (has_state_opt) then
                istate = state_fixed
                if (do_shift_first) then
                    if (self%proj_exists(prev_iproj, istate)) then
                        iref_start = (istate - 1) * self%p_ptr%nspace
                        iref_prev  = iref_start + prev_iproj
                        irot = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                        call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        cxy_shift = 0.
                        if (irot /= 0) cxy_shift = cxy(2:3)
                    endif
                endif
            else
                istate = o_prev%get_state()
                if (do_shift_first) then
                    if (istate >= 1 .and. istate <= self%p_ptr%nstates .and. &
                        self%proj_exists(prev_iproj, istate)) then
                        iref_start = (istate - 1) * self%p_ptr%nspace
                        iref_prev  = iref_start + prev_iproj
                        irot = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                        call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        cxy_shift = 0.
                        if (irot /= 0) cxy_shift = cxy(2:3)
                    endif
                endif
            endif
            coarse_best_dist(:,ithr) = huge_dist
            if (has_state_opt) then
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    if (.not. self%proj_exists(iproj_full, state_fixed)) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((state_fixed - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                    inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), self%p_ptr%cc_objfun)
                    irot = minloc(inpl_dists(:,ithr), dim=1)
                    coarse_best_dist(isub,ithr) = inpl_dists(irot,ithr)
                enddo
            else
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    do state_i = 1, self%nstates
                        istate = self%active_state_indices(state_i)
                        if (.not. self%proj_exists(iproj_full, istate)) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                        inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), self%p_ptr%cc_objfun)
                        irot = minloc(inpl_dists(:,ithr), dim=1)
                        coarse_best_dist(isub,ithr) = min(coarse_best_dist(isub,ithr), inpl_dists(irot,ithr))
                    enddo
                enddo
            endif
            nproj = 0
            prev_included = .false.
            nvalid_sub = count(coarse_best_dist(:,ithr) < huge_dist / 2.)
            if (nvalid_sub > 0) then
                if (nvalid_sub <= npeak_use) then
                    do isub = 1, nspace_sub
                        if (coarse_best_dist(isub,ithr) >= huge_dist / 2.) cycle
                        call neigh_map%get_neighbors_list(isub, neigh_proj)
                        if (allocated(neigh_proj)) then
                            call append_proj_block(tbuf(ithr)%proj_sel, nproj, neigh_proj)
                            deallocate(neigh_proj)
                        endif
                        if (isub == prev_sub) prev_included = .true.
                    enddo
                else
                    peak_sub_idxs(:,ithr) = minnloc(coarse_best_dist(:,ithr), npeak_use)
                    do j = 1, npeak_use
                        call neigh_map%get_neighbors_list(peak_sub_idxs(j,ithr), neigh_proj)
                        if (allocated(neigh_proj)) then
                            call append_proj_block(tbuf(ithr)%proj_sel, nproj, neigh_proj)
                            deallocate(neigh_proj)
                        endif
                    enddo
                    prev_included = any(peak_sub_idxs(:,ithr) == prev_sub)
                endif
            endif
            if (.not. prev_included) then
                call neigh_map%get_neighbors_list(prev_sub, neigh_proj)
                if (allocated(neigh_proj)) then
                    call append_proj_block(tbuf(ithr)%proj_sel, nproj, neigh_proj)
                    deallocate(neigh_proj)
                endif
            endif
            if (nproj <= 0) then
                call add_all_active_projs(self, tbuf(ithr)%proj_sel, nproj)
            else if (.not. any(self%proj_active_state_count(tbuf(ithr)%proj_sel(1:nproj)) > 0)) then
                call add_all_active_projs(self, tbuf(ithr)%proj_sel, nproj)
            endif
            if (allocated(sel_rows(i)%proj_idx)) deallocate(sel_rows(i)%proj_idx)
            allocate(sel_rows(i)%proj_idx(nproj), source=tbuf(ithr)%proj_sel(1:nproj))
            sel_rows(i)%cxy_shift = cxy_shift
            sel_rows(i)%shifts_valid = do_shift_first
            row_nedges(i) = 0
            do j = 1, nproj
                row_nedges(i) = row_nedges(i) + self%proj_active_state_count(sel_rows(i)%proj_idx(j))
            enddo
            call o_prev%kill
        enddo
        !$omp end parallel do
        if (L_BENCH_GLOB) then
            rt_particle_loop = toc(t_particle_loop)
            t_finalize = tic()
        endif
        if (do_shift_first) then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call neigh_map%kill
        call self%finalize_rows_to_csr(sel_rows, row_nedges, tbuf, inpl_angle_thres)
        call cleanup_staged_rows(sel_rows)
        call cleanup_staged_buffers_subspace(inpl_dists, coarse_best_dist, peak_sub_idxs, row_nedges, inpl_angle_thres, tbuf)
        if (L_BENCH_GLOB) rt_finalize = toc(t_finalize)
        rt_route_setup         = rt_setup
        rt_route_shift_init    = rt_shift_init
        rt_route_particle_loop = rt_particle_loop
        rt_route_finalize      = rt_finalize
    end subroutine build_graph_from_subspace_peaks_impl

    ! Build reference-major adjacency over edge indices for fast per-reference scans.
    subroutine build_ref_adjacency(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: edge_count_by_ref(:)
        integer, allocatable :: thread_ref_counts(:,:), thread_ref_offsets(:,:)
        integer :: edge_idx, ithr, ref_idx, slot
        allocate(edge_count_by_ref(self%nrefs), source=0)
        allocate(thread_ref_counts(self%nrefs, nthr_glob), source=0)
        !$omp parallel default(shared) private(edge_idx,ref_idx,ithr) proc_bind(close)
        !$omp do schedule(static)
        do edge_idx = 1, self%nedges
            ithr = omp_get_thread_num() + 1
            ref_idx = self%edge_ref_index(edge_idx)
            thread_ref_counts(ref_idx, ithr) = thread_ref_counts(ref_idx, ithr) + 1
        enddo
        !$omp end do
        !$omp do schedule(static)
        do ref_idx = 1, self%nrefs
            do ithr = 1, nthr_glob
                edge_count_by_ref(ref_idx) = edge_count_by_ref(ref_idx) + thread_ref_counts(ref_idx, ithr)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        if (allocated(self%ref_edge_offsets)) deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices)) deallocate(self%ref_edge_indices)
        allocate(self%ref_edge_offsets(self%nrefs+1), self%ref_edge_indices(self%nedges))
        self%ref_edge_offsets(1) = 1
        do ref_idx = 1, self%nrefs
            self%ref_edge_offsets(ref_idx+1) = self%ref_edge_offsets(ref_idx) + edge_count_by_ref(ref_idx)
        enddo
        allocate(thread_ref_offsets(self%nrefs, nthr_glob), source=0)
        !$omp parallel default(shared) private(edge_idx,ref_idx,slot,ithr) proc_bind(close)
        !$omp do schedule(static)
        do ref_idx = 1, self%nrefs
            slot = self%ref_edge_offsets(ref_idx)
            do ithr = 1, nthr_glob
                thread_ref_offsets(ref_idx, ithr) = slot
                slot = slot + thread_ref_counts(ref_idx, ithr)
            enddo
        enddo
        !$omp end do
        !$omp do schedule(static)
        do edge_idx = 1, self%nedges
            ithr = omp_get_thread_num() + 1
            ref_idx = self%edge_ref_index(edge_idx)
            slot = thread_ref_offsets(ref_idx, ithr)
            thread_ref_offsets(ref_idx, ithr) = slot + 1
            self%ref_edge_indices(slot) = edge_idx
        enddo
        !$omp end do
        !$omp end parallel
        deallocate(edge_count_by_ref, thread_ref_counts, thread_ref_offsets)
    end subroutine build_ref_adjacency

    ! Refine shifts on already-scored sparse edges.
    subroutine refine_shift_edges(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        real,    allocatable   :: candidate_dist(:,:)
        integer, allocatable   :: candidate_edge(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        logical :: do_shift_refine
        integer :: e, i, iref_full, iref_start, iproj, iptcl, irot, istate, ithr
        integer :: n_candidates, n_refine, n_refs_to_refine, n_samples, nrots
        integer :: refine_rank, ri, state_i
        real    :: cxy_prob(3), lims(2,2), lims_init(2,2)
        do_shift_refine = self%p_ptr%l_prob_sh .and. self%p_ptr%l_doshift
        if (.not. do_shift_refine) return
        if (any(self%edge_val(:)%dist >= huge(1.) / 2.)) then
            THROW_HARD('refine_shift_edges: unresolved sparse edges present; builders must score all edges before shift refinement')
        endif
        call seed_rnd
        nrots = self%b_ptr%pftc%get_nrots()
        n_refs_to_refine = 0
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n_samples, self%p_ptr%prob_athres, state=istate)
            n_refs_to_refine = max(n_refs_to_refine, n_samples)
        enddo
        allocate(candidate_dist(self%maxdeg_ptcl, nthr_glob), candidate_edge(self%maxdeg_ptcl, nthr_glob))
        lims(:,1)      = -self%p_ptr%trs
        lims(:,2)      =  self%p_ptr%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                           maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
        enddo
        !$omp parallel do default(shared) &
        !$omp private(i,iptcl,ithr,n_candidates,e,ri,istate,iproj,irot,n_refine,refine_rank,cxy_prob,iref_full,iref_start) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            n_candidates = self%ptcl_off(i+1) - self%ptcl_off(i)
            do e = 1, n_candidates
                candidate_dist(e,ithr) = self%edge_val(self%ptcl_off(i) + e - 1)%dist
                candidate_edge(e,ithr) = self%ptcl_off(i) + e - 1
            enddo
            n_refine = min(n_refs_to_refine, n_candidates)
            if (n_refine > 0) then
                call hpsort(candidate_dist(1:n_candidates,ithr), candidate_edge(1:n_candidates,ithr))
                do refine_rank = 1, n_refine
                    e      = candidate_edge(refine_rank,ithr)
                    ri     = self%edge_ref_index(e)
                    istate = self%ref_state_indices(ri)
                    iproj  = self%ref_proj_indices(ri)
                    iref_start = (istate - 1) * self%p_ptr%nspace
                    iref_full  = iref_start + iproj
                    call grad_shsrch_obj(ithr)%set_indices(iref_full, iptcl)
                    irot = self%edge_val(e)%inpl
                    cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if (irot > 0) then
                        self%edge_val(e)%inpl   = irot
                        self%edge_val(e)%dist   = eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun)
                        self%edge_val(e)%x      = cxy_prob(2)
                        self%edge_val(e)%y      = cxy_prob(3)
                        self%edge_val(e)%has_sh = .true.
                    endif
                enddo
            endif
        enddo
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        enddo
        deallocate(candidate_dist, candidate_edge)
    end subroutine refine_shift_edges

    ! Normalize per-particle edge distances and map them to a robust global range.
    subroutine ref_normalize_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ptcl_local_idx, edge_idx, lo, hi
        real    :: sum_dist, min_dist, max_dist
        !$omp parallel do default(shared) private(ptcl_local_idx,lo,hi,sum_dist) &
        !$omp proc_bind(close) schedule(static)
        do ptcl_local_idx = 1, self%nptcls
            lo = self%ptcl_off(ptcl_local_idx)
            hi = self%ptcl_off(ptcl_local_idx+1)-1
            sum_dist = sum(self%edge_val(lo:hi)%dist)
            if (sum_dist < TINY) then
                self%edge_val(lo:hi)%dist = 0.
            else
                self%edge_val(lo:hi)%dist = self%edge_val(lo:hi)%dist / sum_dist
            endif
        enddo
        !$omp end parallel do
        min_dist = minval(self%edge_val(:)%dist)
        max_dist = maxval(self%edge_val(:)%dist)
        if ((max_dist - min_dist) < TINY) then
            !$omp parallel do default(shared) private(edge_idx) &
            !$omp proc_bind(close) schedule(static)
            do edge_idx = 1, self%nedges
                self%edge_val(edge_idx)%dist = ran3()
            enddo
            !$omp end parallel do
        else
            self%edge_val(:)%dist = (self%edge_val(:)%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine ref_normalize_sparse

    ! Sort each reference's candidate edge list by ascending distance.
    subroutine sort_ref_lists_by_dist(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ref_idx, ref_first_ptr, ref_last_ptr, n_ref_edges
        !$omp parallel do default(shared) &
        !$omp private(ref_idx,ref_first_ptr,ref_last_ptr,n_ref_edges) &
        !$omp proc_bind(close) schedule(static)
        do ref_idx = 1, self%nrefs
            ref_first_ptr = self%ref_edge_offsets(ref_idx)
            ref_last_ptr  = self%ref_edge_offsets(ref_idx+1) - 1
            n_ref_edges   = ref_last_ptr - ref_first_ptr + 1
            if (n_ref_edges > 1) call self%sort_eidx_by_dist(self%ref_edge_indices(ref_first_ptr:ref_last_ptr))
        enddo
        !$omp end parallel do
    end subroutine sort_ref_lists_by_dist

    ! In-place heap sort of edge indices using edge distance as key.
    subroutine sort_eidx_by_dist(self, edge_idx_list)
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(inout) :: edge_idx_list(:)
        integer :: n_edges, root_idx, end_idx, tmp_edge_idx
        n_edges = size(edge_idx_list)
        if (n_edges <= 1) return
        do root_idx = n_edges/2, 1, -1
            call self%sift_down(edge_idx_list, root_idx, n_edges)
        enddo
        do end_idx = n_edges, 2, -1
            tmp_edge_idx           = edge_idx_list(1)
            edge_idx_list(1)       = edge_idx_list(end_idx)
            edge_idx_list(end_idx) = tmp_edge_idx
            call self%sift_down(edge_idx_list, 1, end_idx-1)
        enddo
    end subroutine sort_eidx_by_dist

    ! Heap helper: restore max-heap property in edge index list.
    subroutine sift_down(self, edge_idx_list, root_idx, last_idx)
        class(eul_prob_tab_neigh), intent(in) :: self
        integer, intent(inout) :: edge_idx_list(:)
        integer, intent(in) :: root_idx, last_idx
        integer :: scan_idx, left_child_idx, swap_idx, tmp_edge_idx
        real :: swap_dist, left_dist, right_dist
        scan_idx = root_idx
        do
            left_child_idx = 2 * scan_idx
            if (left_child_idx > last_idx) exit
            swap_idx  = scan_idx
            swap_dist = self%edge_val(edge_idx_list(swap_idx))%dist
            left_dist = self%edge_val(edge_idx_list(left_child_idx))%dist
            if (swap_dist < left_dist) then
                swap_idx  = left_child_idx
                swap_dist = left_dist
            endif
            if (left_child_idx + 1 <= last_idx) then
                right_dist = self%edge_val(edge_idx_list(left_child_idx+1))%dist
                if (swap_dist < right_dist) swap_idx = left_child_idx + 1
            endif
            if (swap_idx == scan_idx) exit
            tmp_edge_idx = edge_idx_list(scan_idx)
            edge_idx_list(scan_idx) = edge_idx_list(swap_idx)
            edge_idx_list(swap_idx) = tmp_edge_idx
            scan_idx = swap_idx
        enddo
    end subroutine sift_down

    ! Greedy global assignment from references to still-unassigned particles.
    subroutine ref_assign_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        logical, allocatable :: particle_available(:)
        integer, allocatable :: ref_frontier_edge(:), active_pos_by_ref(:), active_refs(:)
        real,    allocatable :: best_dist_by_ref(:), active_best_dist(:), sorted_dist_scratch(:)
        integer, allocatable :: sorted_idx_scratch(:)
        real    :: projs_athres
        integer :: edge_idx, istate, k, n_active_refs, n_assigned, neighbor_edge
        integer :: ptcl_local_idx, ref_idx, sampled_active_pos
        call self%ref_normalize_sparse()
        call self%sort_ref_lists_by_dist()
        allocate(particle_available(self%nptcls), source=.true.)
        allocate(ref_frontier_edge(self%nrefs), best_dist_by_ref(self%nrefs), active_pos_by_ref(self%nrefs), &
                 active_refs(self%nrefs), active_best_dist(self%nrefs), sorted_dist_scratch(self%nrefs), &
                 sorted_idx_scratch(self%nrefs))
        active_pos_by_ref = 0
        do ref_idx = 1, self%nrefs
            ref_frontier_edge(ref_idx) = self%ref_edge_offsets(ref_idx)
            best_dist_by_ref(ref_idx)  = huge(best_dist_by_ref(ref_idx))
        enddo
        n_active_refs = 0
        do ref_idx = 1, self%nrefs
            call update_ref_frontier(ref_idx)
        enddo
        projs_athres = 0.
        do k = 1, self%nstates
            istate = self%active_state_indices(k)
            projs_athres = max(projs_athres, calc_athres(self%b_ptr%spproj_field, 'dist', self%p_ptr%prob_athres, state=istate))
        enddo
        n_assigned = 0
        do while (n_assigned < self%nptcls)
            if (n_active_refs <= 0) then
                THROW_HARD('neighborhood too restrictive: no active references but particles remain unassigned')
            endif
            sampled_active_pos = angle_sampling(active_best_dist(1:n_active_refs), sorted_dist_scratch(1:n_active_refs), &
                                               sorted_idx_scratch(1:n_active_refs), projs_athres, self%p_ptr%prob_athres)
            ref_idx  = active_refs(sampled_active_pos)
            edge_idx = self%ref_edge_indices(ref_frontier_edge(ref_idx))
            ptcl_local_idx = self%edge_ptcl(edge_idx)
            if (particle_available(ptcl_local_idx)) then
                particle_available(ptcl_local_idx) = .false.
                self%assgn_map(ptcl_local_idx) = self%edge_val(edge_idx)
                n_assigned = n_assigned + 1
            endif
            do neighbor_edge = self%ptcl_off(ptcl_local_idx), self%ptcl_off(ptcl_local_idx+1)-1
                call update_ref_frontier(self%edge_ref_index(neighbor_edge))
            enddo
        enddo
        deallocate(particle_available, ref_frontier_edge, best_dist_by_ref, active_pos_by_ref, active_refs, active_best_dist, sorted_dist_scratch, sorted_idx_scratch)

    contains

        subroutine update_ref_frontier(ref_idx)
            integer, intent(in) :: ref_idx
            integer :: edge_ptr, last_edge_ptr, edge_idx, ptcl_local_idx, active_pos, last_active_pos, swapped_ref_idx
            edge_ptr = ref_frontier_edge(ref_idx)
            last_edge_ptr = self%ref_edge_offsets(ref_idx+1) - 1
            do while (edge_ptr <= last_edge_ptr)
                edge_idx = self%ref_edge_indices(edge_ptr)
                ptcl_local_idx = self%edge_ptcl(edge_idx)
                if (particle_available(ptcl_local_idx)) exit
                edge_ptr = edge_ptr + 1
            enddo
            ref_frontier_edge(ref_idx) = edge_ptr
            if (edge_ptr > last_edge_ptr) then
                best_dist_by_ref(ref_idx) = huge(best_dist_by_ref(ref_idx))
                active_pos = active_pos_by_ref(ref_idx)
                if (active_pos /= 0) then
                    last_active_pos = n_active_refs
                    swapped_ref_idx = active_refs(last_active_pos)
                    active_refs(active_pos) = swapped_ref_idx
                    active_best_dist(active_pos) = active_best_dist(last_active_pos)
                    active_pos_by_ref(swapped_ref_idx) = active_pos
                    active_pos_by_ref(ref_idx) = 0
                    n_active_refs = n_active_refs - 1
                endif
            else
                best_dist_by_ref(ref_idx) = self%edge_val(self%ref_edge_indices(edge_ptr))%dist
                active_pos = active_pos_by_ref(ref_idx)
                if (active_pos == 0) then
                    n_active_refs = n_active_refs + 1
                    active_pos = n_active_refs
                    active_refs(active_pos) = ref_idx
                    active_pos_by_ref(ref_idx) = active_pos
                endif
                active_best_dist(active_pos) = best_dist_by_ref(ref_idx)
            endif
        end subroutine update_ref_frontier

    end subroutine ref_assign_sparse

    ! Write sparse table payload for this partition/object to a binary stream file.
    subroutine write_tab(self, binfname)
        class(eul_prob_tab_neigh), intent(in) :: self
        class(string),             intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(3)
        file_header(1) = self%nrefs
        file_header(2) = self%nptcls
        file_header(3) = self%nedges
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1) file_header
        addr = size(file_header) * sizeof(file_header(1)) + 1
        write(unit=funit, pos=addr) self%pinds
        addr = addr + self%nptcls * sizeof(self%pinds(1))
        write(unit=funit, pos=addr) self%ptcl_off
        addr = addr + (self%nptcls + 1) * sizeof(self%ptcl_off(1))
        write(unit=funit, pos=addr) self%edge_val
        call fclose(funit)
    end subroutine write_tab

    ! Merge partitioned sparse tables into this global object.
    subroutine read_tabs_to_glob(self, fbody, nparts, numlen)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(ptcl_ref), allocatable :: edge_val_loc(:), stage_edge_val(:), tmp_stage_edge_val(:)
        integer,         allocatable :: degree_by_ptcl(:), fill_ptr(:), pind_to_local(:), pinds_loc(:), ptcl_off_loc(:)
        integer,         allocatable :: stage_global_i(:), tmp_stage_global_i(:)
        type(string) :: binfname
        integer      :: addr, file_header(3), funit, global_e, global_i, header_bytes, io_stat
        integer      :: ipart, local_e, local_i, needed_cap, nedges_loc, new_cap, nptcls_loc
        integer      :: nrefs_loc, ref_idx, row_deg, stage_cap, stage_nused
        if (nparts < 1) then
            THROW_HARD('read_tabs_to_glob: nparts must be >= 1')
        endif
        if (allocated(self%ptcl_off))         deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index))   deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))        deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))         deallocate(self%edge_val)
        if (allocated(self%ref_edge_offsets)) deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices)) deallocate(self%ref_edge_indices)
        allocate(degree_by_ptcl(self%nptcls), pind_to_local(self%p_ptr%nptcls), source=0)
        stage_nused = 0
        stage_cap   = 0
        do global_i = 1, self%nptcls
            if (self%pinds(global_i) < 1 .or. self%pinds(global_i) > self%p_ptr%nptcls) then
                THROW_HARD('read_tabs_to_glob: global particle index out of range')
            endif
            if (pind_to_local(self%pinds(global_i)) /= 0) then
                THROW_HARD('read_tabs_to_glob: duplicate particle index in global sampled set')
            endif
            pind_to_local(self%pinds(global_i)) = global_i
        enddo
        header_bytes = size(file_header) * sizeof(file_header(1))
        do ipart = 1, nparts
            binfname = fbody//int2str_pad(ipart, numlen)//'.dat'
            if (.not. file_exists(binfname)) then
                THROW_HARD('file '//binfname%to_char()//' does not exists!')
            else
                call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            endif
            call fileiochk('simple_eul_prob_tab_neigh; read_tabs_to_glob; file: '//binfname%to_char(), io_stat)
            read(unit=funit, pos=1) file_header
            nrefs_loc  = file_header(1)
            nptcls_loc = file_header(2)
            nedges_loc = file_header(3)
            if (nrefs_loc /= self%nrefs) then
                THROW_HARD('read_tabs_to_glob: nrefs mismatch between sparse partition file and global object')
            endif
            allocate(pinds_loc(nptcls_loc), ptcl_off_loc(nptcls_loc + 1), edge_val_loc(nedges_loc))
            addr = header_bytes + 1
            read(unit=funit, pos=addr) pinds_loc
            addr = addr + nptcls_loc * sizeof(pinds_loc(1))
            read(unit=funit, pos=addr) ptcl_off_loc
            addr = addr + (nptcls_loc + 1) * sizeof(ptcl_off_loc(1))
            read(unit=funit, pos=addr) edge_val_loc
            call fclose(funit)
            if (ptcl_off_loc(1) /= 1) then
                THROW_HARD('read_tabs_to_glob: invalid sparse partition offsets')
            endif
            if (ptcl_off_loc(nptcls_loc + 1) - 1 /= nedges_loc) then
                THROW_HARD('read_tabs_to_glob: sparse partition header/offset mismatch')
            endif
            do local_i = 1, nptcls_loc
                if (pinds_loc(local_i) < 1 .or. pinds_loc(local_i) > self%p_ptr%nptcls) then
                    THROW_HARD('read_tabs_to_glob: partition file contains particle index out of range')
                endif
                global_i = pind_to_local(pinds_loc(local_i))
                if (global_i <= 0) then
                    THROW_HARD('read_tabs_to_glob: partition file contains particle outside the sampled set')
                endif
                row_deg = ptcl_off_loc(local_i + 1) - ptcl_off_loc(local_i)
                if (row_deg <= 0) then
                    THROW_HARD('read_tabs_to_glob: sparse partition row has no edges')
                endif
                degree_by_ptcl(global_i) = degree_by_ptcl(global_i) + row_deg
                needed_cap = stage_nused + row_deg
                if (needed_cap > stage_cap) then
                    new_cap = max(needed_cap, max(1, stage_cap * 2))
                    allocate(tmp_stage_global_i(new_cap), tmp_stage_edge_val(new_cap))
                    if (stage_nused > 0) then
                        tmp_stage_global_i(1:stage_nused) = stage_global_i(1:stage_nused)
                        tmp_stage_edge_val(1:stage_nused) = stage_edge_val(1:stage_nused)
                    endif
                    call move_alloc(tmp_stage_global_i, stage_global_i)
                    call move_alloc(tmp_stage_edge_val, stage_edge_val)
                    stage_cap = new_cap
                endif
                stage_global_i(stage_nused + 1:stage_nused + row_deg) = global_i
                stage_edge_val(stage_nused + 1:stage_nused + row_deg) = &
                    edge_val_loc(ptcl_off_loc(local_i):ptcl_off_loc(local_i + 1) - 1)
                stage_nused = stage_nused + row_deg
            enddo
            deallocate(pinds_loc, ptcl_off_loc, edge_val_loc)
        enddo
        if (any(degree_by_ptcl <= 0)) then
            THROW_HARD('read_tabs_to_glob: missing sparse neighborhood rows for one or more particles')
        endif
        allocate(self%ptcl_off(self%nptcls + 1))
        self%ptcl_off(1) = 1
        do global_i = 1, self%nptcls
            self%ptcl_off(global_i + 1) = self%ptcl_off(global_i) + degree_by_ptcl(global_i)
        enddo
        self%nedges = self%ptcl_off(self%nptcls + 1) - 1
        self%maxdeg_ptcl = maxval(degree_by_ptcl)
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        if (stage_nused /= self%nedges) then
            THROW_HARD('read_tabs_to_glob: staged edge count does not match reconstructed global size')
        endif
        allocate(fill_ptr(self%nptcls), source=self%ptcl_off(1:self%nptcls))
        do global_e = 1, stage_nused
            global_i = stage_global_i(global_e)
            ref_idx = self%ref_index_map(stage_edge_val(global_e)%iproj, stage_edge_val(global_e)%istate)
            if (ref_idx <= 0) then
                THROW_HARD('read_tabs_to_glob: sparse partition edge refers to an inactive reference')
            endif
            local_e = fill_ptr(global_i)
            self%edge_ref_index(local_e) = ref_idx
            self%edge_ptcl(local_e)      = global_i
            self%edge_val(local_e)       = stage_edge_val(global_e)
            fill_ptr(global_i) = local_e + 1
        enddo
        if (any(fill_ptr /= self%ptcl_off(2:self%nptcls + 1))) then
            THROW_HARD('read_tabs_to_glob: failed to reconstruct the sparse global table consistently')
        endif
        call self%build_ref_adjacency()
        deallocate(degree_by_ptcl, fill_ptr, pind_to_local)
        if (allocated(stage_global_i)) deallocate(stage_global_i)
        if (allocated(stage_edge_val)) deallocate(stage_edge_val)
    end subroutine read_tabs_to_glob

    ! Write current particle assignments to a binary stream file.
    subroutine write_assignment(self, binfname)
        class(eul_prob_tab_neigh), intent(in) :: self
        class(string),             intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1)          self%nptcls
        write(unit=funit, pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! Read global assignment file and scatter into this local sampled set.
    subroutine read_assignment(self, binfname)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in) :: binfname
        type(ptcl_ref), allocatable :: assgn_glob(:)
        integer,        allocatable :: pind_to_local(:)
        integer :: funit, global_idx, headsz, io_stat, local_idx, nptcls_glob, pind
        headsz = sizeof(nptcls_glob)
        if (.not. file_exists(binfname)) then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        endif
        call fileiochk('simple_eul_prob_tab_neigh; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit, pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit, pos=headsz + 1) assgn_glob
        call fclose(funit)
        allocate(pind_to_local(self%p_ptr%nptcls), source=0)
        do local_idx = 1, self%nptcls
            pind = self%assgn_map(local_idx)%pind
            if (pind < 1 .or. pind > self%p_ptr%nptcls) cycle
            pind_to_local(pind) = local_idx
        enddo
        do global_idx = 1, nptcls_glob
            pind = assgn_glob(global_idx)%pind
            if (pind < 1 .or. pind > self%p_ptr%nptcls) cycle
            local_idx = pind_to_local(pind)
            if (local_idx > 0) self%assgn_map(local_idx) = assgn_glob(global_idx)
        enddo
        deallocate(pind_to_local, assgn_glob)
    end subroutine read_assignment

    ! Release all allocated storage and reset object to default state.
    subroutine kill(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        if (allocated(self%pinds))                   deallocate(self%pinds)
        if (allocated(self%assgn_map))               deallocate(self%assgn_map)
        if (allocated(self%proj_exists))             deallocate(self%proj_exists)
        if (allocated(self%state_exists))            deallocate(self%state_exists)
        if (allocated(self%active_state_indices))    deallocate(self%active_state_indices)
        if (allocated(self%ref_proj_indices))        deallocate(self%ref_proj_indices)
        if (allocated(self%ref_state_indices))       deallocate(self%ref_state_indices)
        if (allocated(self%ref_index_map))           deallocate(self%ref_index_map)
        if (allocated(self%proj_active_state_count)) deallocate(self%proj_active_state_count)
        if (allocated(self%ptcl_off))                deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index))          deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))               deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))                deallocate(self%edge_val)
        if (allocated(self%ref_edge_offsets))        deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices))        deallocate(self%ref_edge_indices)
        self%nptcls      = 0
        self%nstates     = 0
        self%nrefs       = 0
        self%nedges      = 0
        self%maxdeg_ptcl = 0
        self%b_ptr       => null()
        self%p_ptr       => null()
    end subroutine kill

    ! Append a projection block as-is (subspace neighborhoods are assumed disjoint).
    subroutine append_proj_block(proj_sel, nused, vals)
        integer, intent(inout) :: proj_sel(:)
        integer, intent(inout) :: nused
        integer, intent(in)    :: vals(:)
        integer :: nadd
        nadd = size(vals)
        if (nadd <= 0) return
        if (nused + nadd > size(proj_sel)) then
            THROW_HARD('append_proj_block: projection buffer overflow')
        endif
        proj_sel(nused+1:nused+nadd) = vals
        nused = nused + nadd
    end subroutine append_proj_block

    ! Fill projection buffer with all active projections (no dedup needed).
    subroutine add_all_active_projs(self, proj_sel, nused)
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(inout) :: proj_sel(:)
        integer,                   intent(inout) :: nused
        integer :: iproj
        nused = 0
        do iproj = 1, self%p_ptr%nspace
            if (self%proj_active_state_count(iproj) <= 0) cycle
            nused = nused + 1
            if (nused > size(proj_sel)) then
                THROW_HARD('add_all_active_projs: projection buffer overflow')
            endif
            proj_sel(nused) = iproj
        enddo
    end subroutine add_all_active_projs

end module simple_eul_prob_tab_neigh