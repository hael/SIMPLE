!@descr: builds sparse neighborhood graph directly in projection-major order
!
!   * neighborhood builders emit sparse rows directly
!   * after any build_graph_from_*, every edge in edge_val must have scored dist/inpl
!   * ptree neighborhood construction caches already-evaluated refs directly
!   * ref_assign assumes all edges are already scored. 
!
! 2do: after testing the ptree mode, make a ptree_srch_geom neigh_type that overcomes the need for initial coarse subpace evaluation
!
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,               only: builder
use simple_pftc_shsrch_grad,      only: pftc_shsrch_grad
use simple_eul_prob_tab,          only: angle_sampling, calc_num2sample, calc_athres, eulprob_dist_switch
use simple_eulspace_neigh_map,    only: eulspace_neigh_map
use simple_strategy3D_tree_utils, only: select_peak_trees, trace_tree_prob
implicit none
private
#include "simple_local_flags.inc"

public :: eul_prob_tab_neigh

type :: eul_prob_tab_neigh
    class(builder),     pointer :: b_ptr => null()
    class(parameters),  pointer :: p_ptr => null()
    integer                     :: nptcls = 0
    integer,        allocatable :: pinds(:)
    type(ptcl_ref), allocatable :: assgn_map(:)
    logical,        allocatable :: proj_exists(:,:)
    logical,        allocatable :: state_exists(:)
    integer                     :: nstates = 0
    integer                     :: nrefs   = 0
    integer,        allocatable :: active_state_indices(:)
    integer,        allocatable :: ref_proj_indices(:)
    integer,        allocatable :: ref_state_indices(:)
    integer,        allocatable :: ref_index_map(:,:)
    integer,        allocatable :: proj_active_state_count(:)
    integer,        allocatable :: proj_ref_offsets(:)
    integer,        allocatable :: proj_ref_indices(:)
    integer,        allocatable :: active_proj_indices(:)
    integer                     :: nactive_proj = 0
    integer                     :: nedges       = 0
    integer                     :: maxdeg_ptcl  = 0
    integer,        allocatable :: ptcl_off(:)
    integer,        allocatable :: edge_ref_index(:)
    integer,        allocatable :: edge_ptcl(:)
    type(ptcl_ref), allocatable :: edge_val(:)
    integer,        allocatable :: ref_edge_offsets(:)
    integer,        allocatable :: ref_edge_indices(:)
contains
    procedure          :: build_sparse_neigh_graph
    procedure, private :: init_common
    procedure, private :: build_ref_lists_and_map
    procedure, private :: check_subspace_prereqs
    procedure, private :: init_edge_default
    procedure, private :: score_and_materialize_row_from_proj_list
    procedure, private :: flatten_sparse_rows
    procedure, private :: build_proj_list_from_subspace_selection
    procedure, private :: fill_proj_set_from_active
    procedure, private :: get_ptcl_prev_geom
    procedure, private :: apply_shift_seed
    procedure, private :: build_graph_from_prev_geom
    procedure, private :: build_graph_from_subspace_peaks
    procedure, private :: build_graph_from_subspace_peaks_one_state
    procedure, private :: build_graph_from_subspace_peaks_impl
    procedure, private :: build_graph_from_ptree_srch
    procedure, private :: build_graph_from_ptree_srch_one_state
    procedure, private :: build_graph_from_ptree_srch_impl
    procedure, private :: build_ref_adjacency
    procedure, private :: refine_shift_sparse_edges
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

type :: sparse_row_tmp
    integer,        allocatable :: ref_idx(:)
    type(ptcl_ref), allocatable :: val(:)
end type sparse_row_tmp

type :: proj_set_buf
    integer, allocatable :: list(:)
    integer, allocatable :: stamp(:)
    integer              :: nused = 0
    integer              :: gen   = 0
end type proj_set_buf

type :: row_build_thread_buf
    type(proj_set_buf)          :: proj_set
    integer,        allocatable :: cached_ref_idx(:)
    type(ptcl_ref), allocatable :: cached_edge_val(:)
    integer                     :: ncached        = 0
end type row_build_thread_buf

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
        logical :: do_build_sparse_graph
        do_build_sparse_graph = .true.
        if (present(build_sparse_graph)) do_build_sparse_graph = build_sparse_graph
        call self%init_common(params, build, pinds, empty_okay)
        call self%build_ref_lists_and_map
        if (do_build_sparse_graph) then
            call self%check_subspace_prereqs
            select case(trim(neigh_type))
                case('geom')
                    call self%build_graph_from_prev_geom
                case('subspace_srch')
                    if (present(state)) then
                        call self%build_graph_from_subspace_peaks_one_state(state)
                    else
                        call self%build_graph_from_subspace_peaks
                    endif
                case('ptree_srch')
                    if (present(state)) then
                        call self%build_graph_from_ptree_srch_one_state(state)
                    else
                        call self%build_graph_from_ptree_srch
                    endif
                case default
                    THROW_HARD('build_sparse_neigh_graph: unsupported neigh_type='//trim(neigh_type)//'; expected geom, subspace_srch or ptree_srch')
            end select
            call self%build_ref_adjacency
            call self%refine_shift_sparse_edges
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
        call self%kill
        self%p_ptr => params
        self%b_ptr => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        !$omp workshare
        self%assgn_map(:)%pind   = self%pinds(:)
        self%assgn_map(:)%istate = 0
        self%assgn_map(:)%iproj  = 0
        self%assgn_map(:)%inpl   = 0
        self%assgn_map(:)%dist   = huge(1.)
        self%assgn_map(:)%x      = 0.
        self%assgn_map(:)%y      = 0.
        self%assgn_map(:)%has_sh = .false.
        !$omp end workshare
        l_empty = (trim(params%empty3Dcavgs) .eq. 'yes')
        if (present(empty_okay)) l_empty = empty_okay
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        if (l_empty) then
            allocate(self%proj_exists(self%p_ptr%nspace, self%p_ptr%nstates), source=.true.)
        else
            self%proj_exists = self%b_ptr%spproj_field%projs_exist(self%p_ptr%nstates, self%p_ptr%nspace, thres=MIN_POP)
        endif
        self%nrefs = count(self%proj_exists .eqv. .true.)
        if (self%nrefs == 0) then
            THROW_HARD('No valid references available after state/projection existence filtering')
        endif
    end subroutine init_common

    ! Build reference lists plus projection-major CSR and active projection list.
    subroutine build_ref_lists_and_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: write_ptr(:)
        integer :: state_list_idx, ref_idx, istate, iproj, pos
        allocate(self%active_state_indices(self%nstates), self%ref_proj_indices(self%nrefs), self%ref_state_indices(self%nrefs))
        allocate(self%ref_index_map(self%p_ptr%nspace, self%p_ptr%nstates), self%proj_active_state_count(self%p_ptr%nspace), source=0)
        state_list_idx = 0
        ref_idx = 0
        do istate = 1, self%p_ptr%nstates
            if (.not. self%state_exists(istate)) cycle
            state_list_idx = state_list_idx + 1
            self%active_state_indices(state_list_idx) = istate
            do iproj = 1, self%p_ptr%nspace
                if (.not. self%proj_exists(iproj, istate)) cycle
                ref_idx = ref_idx + 1
                self%ref_proj_indices(ref_idx) = iproj
                self%ref_state_indices(ref_idx) = istate
                self%ref_index_map(iproj, istate) = ref_idx
                self%proj_active_state_count(iproj) = self%proj_active_state_count(iproj) + 1
            enddo
        enddo
        allocate(self%proj_ref_offsets(self%p_ptr%nspace + 1))
        self%proj_ref_offsets(1) = 1
        do iproj = 1, self%p_ptr%nspace
            self%proj_ref_offsets(iproj + 1) = self%proj_ref_offsets(iproj) + self%proj_active_state_count(iproj)
        enddo
        if (self%proj_ref_offsets(self%p_ptr%nspace + 1) - 1 /= self%nrefs) then
            THROW_HARD('build_ref_lists_and_map: projection-major CSR size mismatch')
        endif
        allocate(self%proj_ref_indices(self%nrefs), write_ptr(self%p_ptr%nspace))
        write_ptr = self%proj_ref_offsets(1:self%p_ptr%nspace)
        do ref_idx = 1, self%nrefs
            iproj = self%ref_proj_indices(ref_idx)
            pos = write_ptr(iproj)
            self%proj_ref_indices(pos) = ref_idx
            write_ptr(iproj) = pos + 1
        enddo
        if (any(write_ptr /= self%proj_ref_offsets(2:self%p_ptr%nspace + 1))) then
            THROW_HARD('build_ref_lists_and_map: failed to build projection-major CSR')
        endif
        deallocate(write_ptr)
        self%nactive_proj = count(self%proj_active_state_count > 0)
        allocate(self%active_proj_indices(self%nactive_proj))
        pos = 0
        do iproj = 1, self%p_ptr%nspace
            if (self%proj_active_state_count(iproj) <= 0) cycle
            pos = pos + 1
            self%active_proj_indices(pos) = iproj
        enddo
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

    ! Build one sparse row in final projection-major order, scoring cache misses on demand.
    subroutine score_and_materialize_row_from_proj_list(self, iptcl, proj_idx, nproj, cxy_shift, shifts_valid, score_buf, row, &
                                                        cached_ref_idx, cached_edge_val, ncached)
        class(eul_prob_tab_neigh),   intent(in)    :: self
        integer,                     intent(in)    :: iptcl, nproj
        integer,                     intent(in)    :: proj_idx(:)
        real,                        intent(in)    :: cxy_shift(2)
        logical,                     intent(in)    :: shifts_valid
        real,                        intent(inout) :: score_buf(:)
        type(sparse_row_tmp),        intent(inout) :: row
        integer,        allocatable, intent(inout) :: cached_ref_idx(:)
        type(ptcl_ref), allocatable, intent(inout) :: cached_edge_val(:)
        integer,                     intent(inout) :: ncached
        integer :: i, nrow, iproj, istate, ref_idx, slot, irot, iref_full, p, state_i
        type(ptcl_ref) :: v
        real :: rotmat_loc(2,2), rot_xy_loc(2)
        real :: sorted_dist_scratch(size(score_buf)), inpl_angle_thres(self%p_ptr%nstates)
        integer :: sorted_idx_scratch(size(score_buf))
        if (allocated(row%ref_idx)) deallocate(row%ref_idx)
        if (allocated(row%val))     deallocate(row%val)
        nrow = 0
        do i = 1, nproj
            iproj = proj_idx(i)
            if (iproj < 1 .or. iproj > self%p_ptr%nspace) cycle
            nrow = nrow + (self%proj_ref_offsets(iproj + 1) - self%proj_ref_offsets(iproj))
        enddo
        if (nrow <= 0) then
            THROW_HARD('score_and_materialize_row_from_proj_list: zero-row materialization')
        endif
        allocate(row%ref_idx(nrow), row%val(nrow))
        inpl_angle_thres = 0.
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            inpl_angle_thres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        nrow = 0
        do i = 1, nproj
            iproj = proj_idx(i)
            if (iproj < 1 .or. iproj > self%p_ptr%nspace) cycle
            do p = self%proj_ref_offsets(iproj), self%proj_ref_offsets(iproj + 1) - 1
                ref_idx = self%proj_ref_indices(p)
                nrow = nrow + 1
                row%ref_idx(nrow) = ref_idx
                slot = 0
                if (ncached > 0) slot = find_int_buf(cached_ref_idx, ref_idx, ncached)
                if (slot > 0) then
                    row%val(nrow) = cached_edge_val(slot)
                else
                    istate = self%ref_state_indices(ref_idx)
                    iref_full = (istate - 1) * self%p_ptr%nspace + iproj
                    call self%b_ptr%pftc%gen_objfun_vals(iref_full, iptcl, cxy_shift, score_buf)
                    score_buf = eulprob_dist_switch(score_buf, self%p_ptr%cc_objfun)
                    irot = angle_sampling(score_buf, sorted_dist_scratch, sorted_idx_scratch, inpl_angle_thres(istate), self%p_ptr%prob_athres)
                    call self%init_edge_default(iptcl, ref_idx, v)
                    v%dist = score_buf(irot)
                    v%inpl = irot
                    if (shifts_valid) then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat_loc)
                        rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                        v%x = rot_xy_loc(1)
                        v%y = rot_xy_loc(2)
                        v%has_sh = .true.
                    endif
                    row%val(nrow) = v
                    call append_or_improve_ref(cached_ref_idx, cached_edge_val, ncached, ref_idx, v)
                endif
            enddo
        enddo
    end subroutine score_and_materialize_row_from_proj_list

    ! Flatten per-particle temporary rows into particle-major CSR edge arrays.
    subroutine flatten_sparse_rows(self, rows)
        class(eul_prob_tab_neigh),         intent(inout) :: self
        type(sparse_row_tmp), allocatable, intent(inout) :: rows(:)
        integer, allocatable :: nrow_by_ptcl(:)
        integer :: i, j, e
        if (allocated(self%ptcl_off))       deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index)) deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))      deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))       deallocate(self%edge_val)
        allocate(nrow_by_ptcl(self%nptcls), source=0)
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            if (allocated(rows(i)%ref_idx)) nrow_by_ptcl(i) = size(rows(i)%ref_idx)
        enddo
        !$omp end parallel do
        if (any(nrow_by_ptcl <= 0)) then
            THROW_HARD('flatten_sparse_rows: particle has zero valid neighbors')
        endif
        allocate(self%ptcl_off(self%nptcls + 1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            self%ptcl_off(i+1) = self%ptcl_off(i) + nrow_by_ptcl(i)
        enddo
        self%nedges = self%ptcl_off(self%nptcls + 1) - 1
        self%maxdeg_ptcl = maxval(self%ptcl_off(2:self%nptcls+1) - self%ptcl_off(1:self%nptcls))
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        !$omp parallel do default(shared) private(i,j,e) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            e = self%ptcl_off(i)
            do j = 1, size(rows(i)%ref_idx)
                self%edge_ref_index(e) = rows(i)%ref_idx(j)
                self%edge_ptcl(e)      = i
                self%edge_val(e)       = rows(i)%val(j)
                e = e + 1
            enddo
        enddo
        !$omp end parallel do
        do i = 1, size(rows)
            if (allocated(rows(i)%ref_idx)) deallocate(rows(i)%ref_idx)
            if (allocated(rows(i)%val))     deallocate(rows(i)%val)
        enddo
        deallocate(nrow_by_ptcl)
        deallocate(rows)
    end subroutine flatten_sparse_rows

    ! Build a projection list from selected subspace representations.
    ! Gathers neighborhood lists for selected subspace reps, unions them,
    ! unions with previous neighborhood, fallback to active_proj_indices if needed.
    ! proj_set%list/stamp must be pre-allocated to at least nspace by the caller.
    ! proj_set%gen is advanced here and proj_set%nused is produced.
    subroutine build_proj_list_from_subspace_selection(self, neigh_map, selected_subs, nsel, prev_sub, proj_set)
        class(eul_prob_tab_neigh),  intent(in)    :: self
        type(eulspace_neigh_map),   intent(in)    :: neigh_map
        integer,                    intent(in)    :: selected_subs(:), nsel, prev_sub
        type(proj_set_buf),         intent(inout) :: proj_set
        integer, allocatable :: neigh_proj(:)
        integer :: j, p
        proj_set%gen   = proj_set%gen + 1
        proj_set%nused = 0
        ! Gather neighborhoods for selected subspaces
        do j = 1, nsel
            call neigh_map%get_neighbors_list(selected_subs(j), neigh_proj)
            if (allocated(neigh_proj)) then
                do p = 1, size(neigh_proj)
                    if (proj_set%stamp(neigh_proj(p)) /= proj_set%gen) then
                        proj_set%stamp(neigh_proj(p)) = proj_set%gen
                        proj_set%nused = proj_set%nused + 1
                        proj_set%list(proj_set%nused) = neigh_proj(p)
                    endif
                enddo
                deallocate(neigh_proj)
            endif
        enddo
        ! Add neighborhoods of previous projection
        call neigh_map%get_neighbors_list(prev_sub, neigh_proj)
        if (allocated(neigh_proj)) then
            do p = 1, size(neigh_proj)
                if (proj_set%stamp(neigh_proj(p)) /= proj_set%gen) then
                    proj_set%stamp(neigh_proj(p)) = proj_set%gen
                    proj_set%nused = proj_set%nused + 1
                    proj_set%list(proj_set%nused) = neigh_proj(p)
                endif
            enddo
            deallocate(neigh_proj)
        endif
        ! Fallback to all active projections if no valid neighborhoods found
        if (proj_set%nused <= 0 .or. &
            &.not. any(self%proj_ref_offsets(proj_set%list(1:proj_set%nused) + 1) > self%proj_ref_offsets(proj_set%list(1:proj_set%nused)))) then
            call self%fill_proj_set_from_active(proj_set)
        endif
    end subroutine build_proj_list_from_subspace_selection

    ! Copy the full active-projection list into proj_set (fallback helper).
    subroutine fill_proj_set_from_active(self, proj_set)
        class(eul_prob_tab_neigh), intent(in)    :: self
        type(proj_set_buf),        intent(inout) :: proj_set
        proj_set%nused = self%nactive_proj
        proj_set%list(1:proj_set%nused) = self%active_proj_indices(1:self%nactive_proj)
    end subroutine fill_proj_set_from_active

    ! Extract previous particle geometry: orientation, closest projection, subspace index, and state.
    ! o_prev is returned live; caller is responsible for o_prev%kill.
    subroutine get_ptcl_prev_geom(self, iptcl, state_fixed, has_state_opt, &
                                   o_prev, prev_iproj, prev_sub, istate)
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(in)    :: iptcl, state_fixed
        logical,                   intent(in)    :: has_state_opt
        type(ori),                 intent(out)   :: o_prev
        integer,                   intent(out)   :: prev_iproj, prev_sub, istate
        call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
        prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
        prev_sub   = self%b_ptr%subspace_full2sub_map(prev_iproj)
        if (has_state_opt) then
            istate = state_fixed
        else
            istate = o_prev%get_state()
        endif
    end subroutine get_ptcl_prev_geom

    ! Seed cxy_shift from a coarse in-plane+shift search at the previous projection.
    ! Guard conditions (istate bounds, proj_exists) are checked internally; cxy_shift is
    ! set to [0.,0.] when the search is skipped or returns no valid rotation.
    ! grad_shsrch must already be initialized by the caller.
    subroutine apply_shift_seed(self, iptcl, prev_iproj, istate, o_prev, grad_shsrch, cxy_shift)
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer,                   intent(in)    :: iptcl, prev_iproj, istate
        type(ori),                 intent(in)    :: o_prev
        type(pftc_shsrch_grad),    intent(inout) :: grad_shsrch
        real,                      intent(out)   :: cxy_shift(2)
        integer :: iref_prev, irot
        real    :: cxy(3)
        cxy_shift = [0., 0.]
        if (istate < 1 .or. istate > self%p_ptr%nstates) return
        if (.not. self%proj_exists(prev_iproj, istate)) return
        iref_prev = (istate - 1) * self%p_ptr%nspace + prev_iproj
        irot      = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
        call grad_shsrch%set_indices(iref_prev, iptcl)
        cxy = grad_shsrch%minimize(irot=irot, sh_rot=.false.)
        if (irot /= 0) cxy_shift = cxy(2:3)
    end subroutine apply_shift_seed

    ! Build graph from geometric neighbors around previous particle orientations.
    subroutine build_graph_from_prev_geom(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(sparse_row_tmp), allocatable :: rows(:)
        type(eulspace_neigh_map)          :: neigh_map
        type(row_build_thread_buf), allocatable :: tbuf(:)
        type(proj_set_buf),    allocatable :: proj_buf(:)
        type(ori), allocatable :: sub_oris(:)        ! precomputed subspace orientations
        type(ori) :: o_prev, osym
        integer   :: nspace, nspace_sub, npeak_use, nrots
        integer   :: i, ithr, iptcl, isub, prev_iproj, prev_sub
        real      :: dtmp, inplrotdist
        real      :: cxy_shift(2)
        real,    allocatable :: score_buf(:,:)
        real,    allocatable :: sub_dists_thr(:,:)       ! thread-local distance scratch
        integer, allocatable :: peak_sub_idxs_thr(:,:)   ! thread-local peak-index scratch
        nspace     = self%p_ptr%nspace
        nspace_sub = self%p_ptr%nspace_sub
        npeak_use  = max(1, min(self%p_ptr%npeaks, nspace_sub))
        nrots      = self%b_ptr%pftc%get_nrots()
        ! Precompute subspace representative orientations; shared read-only inside OMP region
        allocate(sub_oris(nspace_sub))
        do isub = 1, nspace_sub
            call self%b_ptr%eulspace%get_ori(self%b_ptr%subspace_inds(isub), sub_oris(isub))
        enddo
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(rows(self%nptcls), score_buf(nrots, nthr_glob), tbuf(nthr_glob), proj_buf(nthr_glob), &
                 sub_dists_thr(nspace_sub, nthr_glob), peak_sub_idxs_thr(npeak_use, nthr_glob))
        do ithr = 1, nthr_glob
            allocate(proj_buf(ithr)%list(nspace), proj_buf(ithr)%stamp(nspace), source=0)
        enddo
        !$omp parallel do default(shared) &
        !$omp private(i,ithr,iptcl,o_prev,isub,dtmp,inplrotdist,osym, &
        !$omp& prev_iproj,prev_sub,cxy_shift) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
            do isub = 1, nspace_sub
                call self%b_ptr%pgrpsyms%sym_dists(o_prev, sub_oris(isub), osym, dtmp, inplrotdist)
                sub_dists_thr(isub, ithr) = dtmp
                call osym%kill
            enddo
            peak_sub_idxs_thr(1:npeak_use, ithr) = minnloc(sub_dists_thr(1:nspace_sub, ithr), npeak_use)
            prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
            prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
            call self%build_proj_list_from_subspace_selection(neigh_map, peak_sub_idxs_thr(1:npeak_use, ithr), npeak_use, prev_sub, proj_buf(ithr))
            tbuf(ithr)%ncached = 0
            cxy_shift = [0., 0.]
            call self%score_and_materialize_row_from_proj_list(iptcl, proj_buf(ithr)%list(1:proj_buf(ithr)%nused), proj_buf(ithr)%nused, cxy_shift, .false., &
                                                                score_buf(:,ithr), rows(i), tbuf(ithr)%cached_ref_idx, tbuf(ithr)%cached_edge_val, tbuf(ithr)%ncached)
            call o_prev%kill
        enddo
        !$omp end parallel do
        call neigh_map%kill
        call self%flatten_sparse_rows(rows)
        do ithr = 1, nthr_glob
            if (allocated(tbuf(ithr)%cached_ref_idx)  ) deallocate(tbuf(ithr)%cached_ref_idx)
            if (allocated(tbuf(ithr)%cached_edge_val) ) deallocate(tbuf(ithr)%cached_edge_val)
            if (allocated(proj_buf(ithr)%list)        ) deallocate(proj_buf(ithr)%list)
            if (allocated(proj_buf(ithr)%stamp)       ) deallocate(proj_buf(ithr)%stamp)
        enddo
        do isub = 1, nspace_sub
            call sub_oris(isub)%kill
        enddo
        if (allocated(tbuf)              ) deallocate(tbuf)
        if (allocated(proj_buf)          ) deallocate(proj_buf)
        if (allocated(score_buf)         ) deallocate(score_buf)
        if (allocated(sub_oris)          ) deallocate(sub_oris)
        if (allocated(sub_dists_thr)     ) deallocate(sub_dists_thr)
        if (allocated(peak_sub_idxs_thr) ) deallocate(peak_sub_idxs_thr)
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
        type(sparse_row_tmp), allocatable :: rows(:)
        real,                 allocatable :: inpl_dists(:,:), coarse_best_dist(:,:)
        integer,              allocatable :: peak_sub_idxs(:,:)
        integer,              allocatable :: valid_subs_thr(:,:)
        type(row_build_thread_buf), allocatable :: tbuf(:)
        type(proj_set_buf),   allocatable :: proj_buf(:)
        type(eulspace_neigh_map) :: neigh_map
        type(pftc_shsrch_grad)   :: grad_shsrch_obj(nthr_glob)
        type(ori)                :: o_prev
        real    :: lims(2,2), lims_init(2,2), cxy_shift(2), huge_dist
        logical :: do_shift_first, has_state_opt
        integer :: state_fixed, nspace, nrots, nspace_sub, npeak_use
        integer :: i, isub, ithr, iptcl, istate, iproj_full, irot
        integer :: state_i, nvalid_sub, prev_sub, prev_iproj
        call self%check_subspace_prereqs()
        nspace         = self%p_ptr%nspace
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
        allocate(rows(self%nptcls), inpl_dists(nrots, nthr_glob), coarse_best_dist(nspace_sub, nthr_glob), peak_sub_idxs(npeak_use, nthr_glob), &
             valid_subs_thr(nspace_sub, nthr_glob), tbuf(nthr_glob), proj_buf(nthr_glob))
        do ithr = 1, nthr_glob
            allocate(proj_buf(ithr)%list(nspace), proj_buf(ithr)%stamp(nspace), source=0)
        enddo
        if (do_shift_first) then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
        endif
        !$omp parallel do default(shared) &
        !$omp private(i,ithr,iptcl,o_prev,istate,isub,iproj_full,irot, &
        !$omp& cxy_shift,state_i,nvalid_sub,prev_sub,prev_iproj) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            block
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call self%get_ptcl_prev_geom(iptcl, state_fixed, has_state_opt, o_prev, prev_iproj, prev_sub, istate)
                cxy_shift = [0., 0.]
                if (do_shift_first) call self%apply_shift_seed(iptcl, prev_iproj, istate, o_prev, grad_shsrch_obj(ithr), cxy_shift)
                call o_prev%kill
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
                ! Determine which subspaces to use for projection-gathering
                nvalid_sub = count(dist_is_scored(coarse_best_dist(:,ithr)))
                if (nvalid_sub <= npeak_use .and. nvalid_sub > 0) then
                    ! Use all valid subspaces
                    valid_subs_thr(1:nvalid_sub, ithr) = pack([(isub, isub=1,nspace_sub)], &
                                                               dist_is_scored(coarse_best_dist(1:nspace_sub,ithr)))
                    call self%build_proj_list_from_subspace_selection(neigh_map, valid_subs_thr(1:nvalid_sub, ithr), nvalid_sub, prev_sub, proj_buf(ithr))
                else if (nvalid_sub > npeak_use) then
                    ! Use top npeak_use subspaces by score
                    peak_sub_idxs(:,ithr) = minnloc(coarse_best_dist(:,ithr), npeak_use)
                    call self%build_proj_list_from_subspace_selection(neigh_map, peak_sub_idxs(1:npeak_use,ithr), npeak_use, prev_sub, proj_buf(ithr))
                else
                    ! No valid subspaces, fallback to previous neighborhood only
                    call self%build_proj_list_from_subspace_selection(neigh_map, [(0)], 0, prev_sub, proj_buf(ithr))
                endif
                tbuf(ithr)%ncached = 0
                call self%score_and_materialize_row_from_proj_list(iptcl, proj_buf(ithr)%list(1:proj_buf(ithr)%nused), proj_buf(ithr)%nused, cxy_shift, do_shift_first, &
                                                                    inpl_dists(:,ithr), rows(i), tbuf(ithr)%cached_ref_idx, tbuf(ithr)%cached_edge_val, tbuf(ithr)%ncached)
            end block
        enddo
        !$omp end parallel do
        if (do_shift_first) then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call neigh_map%kill
        call self%flatten_sparse_rows(rows)
        if (allocated(inpl_dists))       deallocate(inpl_dists)
        if (allocated(coarse_best_dist)) deallocate(coarse_best_dist)
        if (allocated(peak_sub_idxs))    deallocate(peak_sub_idxs)
        if (allocated(valid_subs_thr))   deallocate(valid_subs_thr)
        do ithr = 1, nthr_glob
            if (allocated(tbuf(ithr)%cached_ref_idx)  ) deallocate(tbuf(ithr)%cached_ref_idx)
            if (allocated(tbuf(ithr)%cached_edge_val) ) deallocate(tbuf(ithr)%cached_edge_val)
            if (allocated(proj_buf(ithr)%list)        ) deallocate(proj_buf(ithr)%list)
            if (allocated(proj_buf(ithr)%stamp)       ) deallocate(proj_buf(ithr)%stamp)
        enddo
        if (allocated(tbuf))     deallocate(tbuf)
        if (allocated(proj_buf)) deallocate(proj_buf)
    end subroutine build_graph_from_subspace_peaks_impl

    ! Multi-state wrapper for tree-guided graph construction.
    subroutine build_graph_from_ptree_srch(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%build_graph_from_ptree_srch_impl()
    end subroutine build_graph_from_ptree_srch

    ! Fixed-state wrapper for tree-guided graph construction.
    subroutine build_graph_from_ptree_srch_one_state(self, state)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: state
        if (state < 1 .or. state > self%p_ptr%nstates) then
            THROW_HARD('build_graph_from_ptree_srch_one_state: state index out of range')
        endif
        if (.not. self%state_exists(state)) then
            THROW_HARD('build_graph_from_ptree_srch_one_state: requested state does not exist')
        endif
        call self%build_graph_from_ptree_srch_impl(state)
    end subroutine build_graph_from_ptree_srch_one_state

    ! Build graph from probabilistic tree descent and cache evaluated edge values.
    subroutine build_graph_from_ptree_srch_impl(self, state_opt)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, optional,         intent(in)    :: state_opt
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori)              :: o_prev
        type(sparse_row_tmp),   allocatable :: rows(:)
        integer,                allocatable :: peak_trees(:,:), inds_sorted(:,:)
        real,                   allocatable :: inpl_corrs(:,:), tree_best_corrs(:,:), peak_tree_corrs(:,:), dists_sorted(:,:), inpl_angle_thres(:)
        type(row_build_thread_buf), allocatable :: tbuf(:)
        real    :: lims(2,2), lims_init(2,2), cxy_shift(2)
        logical :: do_shift_first, has_state_opt
        integer :: nrots, nspace_sub, npeak_use, ntrees, state_fixed
        integer :: i, isub, ithr, iptcl, istate, iproj_full, itree, state_i
        integer :: npeak_trees, ipeak, prev_iproj, prev_sub
        real    :: corr_best
        call self%check_subspace_prereqs()
        nspace_sub     = self%p_ptr%nspace_sub
        ntrees         = self%b_ptr%block_tree%get_n_trees()
        if (ntrees <= 0) then
            THROW_HARD('build_graph_from_ptree_srch_impl: block_tree has no trees')
        endif
        npeak_use      = max(1, min(self%p_ptr%npeaks, ntrees))
        nrots          = self%b_ptr%pftc%get_nrots()
        do_shift_first = self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift
        has_state_opt  = present(state_opt)
        if (has_state_opt) then
            state_fixed = state_opt
        else
            state_fixed = 1
        endif
        allocate(rows(self%nptcls), inpl_corrs(nrots, nthr_glob), tree_best_corrs(ntrees, nthr_glob), &
                &peak_trees(npeak_use, nthr_glob), peak_tree_corrs(npeak_use, nthr_glob), &
                &dists_sorted(nrots, nthr_glob), inds_sorted(nrots, nthr_glob), tbuf(nthr_glob))
        do ithr = 1, nthr_glob
            allocate(tbuf(ithr)%proj_set%list(self%p_ptr%nspace))
            allocate(tbuf(ithr)%proj_set%stamp(self%p_ptr%nspace), source=0)
        enddo
        allocate(inpl_angle_thres(self%p_ptr%nstates), source=0.)
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            inpl_angle_thres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if (do_shift_first) then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
        endif
        !$omp parallel do default(shared) &
        !$omp private(i,iptcl,ithr,cxy_shift,o_prev,istate, &
        !$omp& isub,iproj_full,itree,corr_best,state_i,npeak_trees,ipeak, &
        !$omp& prev_iproj,prev_sub) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            cxy_shift = [0., 0.]
            tbuf(ithr)%proj_set%nused   = 0
            tbuf(ithr)%ncached          = 0
            tbuf(ithr)%proj_set%gen     = tbuf(ithr)%proj_set%gen + 1
            if (do_shift_first) then
                call self%get_ptcl_prev_geom(iptcl, state_fixed, has_state_opt, o_prev, prev_iproj, prev_sub, istate)
                call self%apply_shift_seed(iptcl, prev_iproj, istate, o_prev, grad_shsrch_obj(ithr), cxy_shift)
                call o_prev%kill
            endif
            tree_best_corrs(:,ithr) = -huge(1.0)
            if (has_state_opt) then
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    itree = self%b_ptr%subspace_full2sub_map(iproj_full)
                    if (.not. self%proj_exists(iproj_full, state_fixed)) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((state_fixed - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                    corr_best = maxval(inpl_corrs(:,ithr))
                    tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                enddo
            else
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    itree = self%b_ptr%subspace_full2sub_map(iproj_full)
                    do state_i = 1, self%nstates
                        istate = self%active_state_indices(state_i)
                        if (.not. self%proj_exists(iproj_full, istate)) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                        corr_best = maxval(inpl_corrs(:,ithr))
                        tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                    enddo
                enddo
            endif
            peak_trees(:,ithr)      = 0
            peak_tree_corrs(:,ithr) = -huge(1.0)
            call select_peak_trees(tree_best_corrs(:,ithr), peak_trees(:,ithr), peak_tree_corrs(:,ithr), npeak_trees)
            do ipeak = 1, npeak_trees
                call trace_tree_prob(self%b_ptr%block_tree, peak_trees(ipeak,ithr), iptcl, ithr, score_ptree_ref)
            enddo
            if (tbuf(ithr)%proj_set%nused <= 0) then
                call self%fill_proj_set_from_active(tbuf(ithr)%proj_set)
            endif
            call self%score_and_materialize_row_from_proj_list(iptcl, tbuf(ithr)%proj_set%list, tbuf(ithr)%proj_set%nused, cxy_shift, do_shift_first, &
                                                                inpl_corrs(:,ithr), rows(i), tbuf(ithr)%cached_ref_idx, tbuf(ithr)%cached_edge_val, tbuf(ithr)%ncached)
        enddo
        !$omp end parallel do
        if (do_shift_first) then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call self%flatten_sparse_rows(rows)
        do ithr = 1, nthr_glob
            if (allocated(tbuf(ithr)%proj_set%list)  ) deallocate(tbuf(ithr)%proj_set%list)
            if (allocated(tbuf(ithr)%proj_set%stamp) ) deallocate(tbuf(ithr)%proj_set%stamp)
            if (allocated(tbuf(ithr)%cached_ref_idx) ) deallocate(tbuf(ithr)%cached_ref_idx)
            if (allocated(tbuf(ithr)%cached_edge_val)) deallocate(tbuf(ithr)%cached_edge_val)
        enddo
        if (allocated(tbuf))             deallocate(tbuf)
        if (allocated(inpl_corrs))       deallocate(inpl_corrs)
        if (allocated(tree_best_corrs))  deallocate(tree_best_corrs)
        if (allocated(peak_trees))       deallocate(peak_trees)
        if (allocated(peak_tree_corrs))  deallocate(peak_tree_corrs)
        if (allocated(dists_sorted))     deallocate(dists_sorted)
        if (allocated(inds_sorted))      deallocate(inds_sorted)
        if (allocated(inpl_angle_thres)) deallocate(inpl_angle_thres)

    contains

        subroutine score_ptree_ref(iproj_eval, iptcl_eval, ithr_eval, best_corr)
            integer, intent(in)  :: iproj_eval
            integer, intent(in)  :: iptcl_eval
            integer, intent(in)  :: ithr_eval
            real,    intent(out) :: best_corr
            integer :: istate_loc, state_j, irot_loc, ref_idx_loc
            real    :: dist_obj(nrots), rotmat_loc(2,2), rot_xy_loc(2)
            type(ptcl_ref) :: v
            logical :: have_valid
            best_corr  = -huge(1.0)
            have_valid = .false.
            if (iproj_eval < 1 .or. iproj_eval > self%p_ptr%nspace) return
            if (has_state_opt) then
                if (.not. self%proj_exists(iproj_eval, state_fixed)) return
                call self%b_ptr%pftc%gen_objfun_vals((state_fixed - 1) * self%p_ptr%nspace + iproj_eval, iptcl_eval, cxy_shift, inpl_corrs(:,ithr_eval))
                best_corr = maxval(inpl_corrs(:,ithr_eval))
                dist_obj = eulprob_dist_switch(inpl_corrs(:,ithr_eval), self%p_ptr%cc_objfun)
                irot_loc = angle_sampling(dist_obj, dists_sorted(:,ithr_eval), inds_sorted(:,ithr_eval), inpl_angle_thres(state_fixed), self%p_ptr%prob_athres)
                ref_idx_loc = self%ref_index_map(iproj_eval, state_fixed)
                if (ref_idx_loc > 0) then
                    call self%init_edge_default(iptcl_eval, ref_idx_loc, v)
                    v%dist = dist_obj(irot_loc)
                    v%inpl = irot_loc
                    if (do_shift_first) then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                        rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                        v%x = rot_xy_loc(1)
                        v%y = rot_xy_loc(2)
                        v%has_sh = .true.
                    endif
                    call append_or_improve_ref(tbuf(ithr_eval)%cached_ref_idx, tbuf(ithr_eval)%cached_edge_val, tbuf(ithr_eval)%ncached, ref_idx_loc, v)
                    have_valid = .true.
                endif
            else
                do state_j = 1, self%nstates
                    istate_loc = self%active_state_indices(state_j)
                    if (.not. self%proj_exists(iproj_eval, istate_loc)) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((istate_loc - 1) * self%p_ptr%nspace + iproj_eval, iptcl_eval, cxy_shift, inpl_corrs(:,ithr_eval))
                    best_corr = max(best_corr, maxval(inpl_corrs(:,ithr_eval)))
                    dist_obj = eulprob_dist_switch(inpl_corrs(:,ithr_eval), self%p_ptr%cc_objfun)
                    irot_loc = angle_sampling(dist_obj, dists_sorted(:,ithr_eval), inds_sorted(:,ithr_eval), inpl_angle_thres(istate_loc), self%p_ptr%prob_athres)
                    ref_idx_loc = self%ref_index_map(iproj_eval, istate_loc)
                    if (ref_idx_loc <= 0) cycle
                    call self%init_edge_default(iptcl_eval, ref_idx_loc, v)
                    v%dist = dist_obj(irot_loc)
                    v%inpl = irot_loc
                    if (do_shift_first) then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                        rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                        v%x = rot_xy_loc(1)
                        v%y = rot_xy_loc(2)
                        v%has_sh = .true.
                    endif
                    call append_or_improve_ref(tbuf(ithr_eval)%cached_ref_idx, tbuf(ithr_eval)%cached_edge_val, tbuf(ithr_eval)%ncached, ref_idx_loc, v)
                    have_valid = .true.
                enddo
            endif
            if (have_valid) then
                if (tbuf(ithr_eval)%proj_set%stamp(iproj_eval) /= tbuf(ithr_eval)%proj_set%gen) then
                    tbuf(ithr_eval)%proj_set%stamp(iproj_eval) = tbuf(ithr_eval)%proj_set%gen
                    tbuf(ithr_eval)%proj_set%nused             = tbuf(ithr_eval)%proj_set%nused + 1
                    tbuf(ithr_eval)%proj_set%list(tbuf(ithr_eval)%proj_set%nused) = iproj_eval
                endif
            endif
        end subroutine score_ptree_ref

    end subroutine build_graph_from_ptree_srch_impl

    ! Build reference-major adjacency over edge indices for fast per-reference scans.
    subroutine build_ref_adjacency(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: edge_count_by_ref(:), edge_count_by_ref_thr(:,:), ref_write_ptr(:)
        integer :: ref_idx, edge_idx, ithr
        allocate(edge_count_by_ref(self%nrefs), source=0)
        allocate(edge_count_by_ref_thr(self%nrefs, nthr_glob), source=0)
        !$omp parallel default(shared) private(edge_idx,ref_idx,ithr) proc_bind(close)
        !$omp do schedule(static)
        do edge_idx = 1, self%nedges
            ithr = omp_get_thread_num() + 1
            ref_idx = self%edge_ref_index(edge_idx)
            edge_count_by_ref_thr(ref_idx, ithr) = edge_count_by_ref_thr(ref_idx, ithr) + 1
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do ref_idx = 1, self%nrefs
            do ithr = 1, nthr_glob
                edge_count_by_ref(ref_idx) = edge_count_by_ref(ref_idx) + edge_count_by_ref_thr(ref_idx, ithr)
            enddo
        enddo
        !$omp end do nowait
        !$omp end parallel
        deallocate(edge_count_by_ref_thr)
        if (allocated(self%ref_edge_offsets)) deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices)) deallocate(self%ref_edge_indices)
        allocate(self%ref_edge_offsets(self%nrefs+1), self%ref_edge_indices(self%nedges), ref_write_ptr(self%nrefs))
        self%ref_edge_offsets(1) = 1
        do ref_idx = 1, self%nrefs
            self%ref_edge_offsets(ref_idx+1) = self%ref_edge_offsets(ref_idx) + edge_count_by_ref(ref_idx)
        enddo
        ref_write_ptr = self%ref_edge_offsets(1:self%nrefs)
        do edge_idx = 1, self%nedges
            ref_idx = self%edge_ref_index(edge_idx)
            self%ref_edge_indices(ref_write_ptr(ref_idx)) = edge_idx
            ref_write_ptr(ref_idx) = ref_write_ptr(ref_idx) + 1
        enddo
        deallocate(edge_count_by_ref, ref_write_ptr)
    end subroutine build_ref_adjacency

    ! Shift-refine scored sparse edges only; unresolved edges are treated as builder invariant violations.
    subroutine refine_shift_sparse_edges(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real,    allocatable   :: candidate_dist(:,:)
        integer, allocatable   :: candidate_edge(:,:)
        integer :: i, e, ithr, iptcl, ri, istate, iproj, irot, iref_full, iref_start
        integer :: lo, hi, n_candidates, n_refine, n_refs_to_refine, n_samples, state_i, refine_rank
        real    :: lims(2,2), lims_init(2,2), cxy_prob(3)
        logical :: do_shift_refine
        do_shift_refine = self%p_ptr%l_prob_sh .and. self%p_ptr%l_doshift
        if (.not. do_shift_refine) return
        call seed_rnd
        if (any(.not. dist_is_scored(self%edge_val(:)%dist))) then
            THROW_HARD('refine_shift_sparse_edges: invariant violated - builders emitted unresolved sparse edges')
        endif
        allocate(candidate_dist(self%maxdeg_ptcl, nthr_glob), candidate_edge(self%maxdeg_ptcl, nthr_glob))
        ! Compute number of references to refine (for shift optimization)
        n_refs_to_refine = 0
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n_samples, self%p_ptr%prob_athres, state=istate)
            n_refs_to_refine = max(n_refs_to_refine, n_samples)
        enddo
        lims(:,1)      = -self%p_ptr%trs
        lims(:,2)      =  self%p_ptr%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                           maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
        enddo
        ! Collect scored edges and refine top-ranked ones.
        !$omp parallel do default(shared) &
        !$omp private(i,iptcl,ithr,istate,iproj,lo,hi,n_candidates,e,ri,n_refine,refine_rank, &
        !$omp& cxy_prob,iref_full,iref_start,irot) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            lo = self%ptcl_off(i)
            hi = self%ptcl_off(i+1)-1
            n_candidates = hi - lo + 1
            if (n_candidates > 0) then
                do e = 1, n_candidates
                    candidate_dist(e,ithr) = self%edge_val(lo + e - 1)%dist
                    candidate_edge(e,ithr) = lo + e - 1
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
            endif
        enddo
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        enddo
        deallocate(candidate_dist, candidate_edge)
    end subroutine refine_shift_sparse_edges

    ! Normalize per-particle edge distances and map them to a robust global range.
    subroutine ref_normalize_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ptcl_local_idx, edge_idx, lo, hi
        real    :: sum_dist, min_dist, max_dist
        if (any(.not. dist_is_scored(self%edge_val(:)%dist))) then
            THROW_HARD('Unscored sparse edges remain before normalization/assignment')
        endif
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
        integer :: ref_idx, sampled_active_pos, edge_idx, ptcl_local_idx, neighbor_edge
        integer :: n_active_refs, n_assigned
        real    :: projs_athres
        integer :: k, istate
        if (any(.not. dist_is_scored(self%edge_val(:)%dist))) then
            THROW_HARD('Unscored sparse edges remain before normalization/assignment')
        endif
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
        deallocate(particle_available, ref_frontier_edge, best_dist_by_ref, active_pos_by_ref)
        deallocate(active_refs, active_best_dist, sorted_dist_scratch, sorted_idx_scratch)

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
        integer      :: funit, addr, io_stat, file_header(3), header_bytes
        integer      :: ipart, nrefs_loc, nptcls_loc, nedges_loc, local_i, global_i, local_e, global_e, ref_idx
        integer      :: row_deg, stage_nused, stage_cap, needed_cap, new_cap
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
                if (degree_by_ptcl(global_i) /= 0) then
                    THROW_HARD('read_tabs_to_glob: partition invariant violated - particle appears in multiple partition files')
                endif
                row_deg = ptcl_off_loc(local_i + 1) - ptcl_off_loc(local_i)
                if (row_deg <= 0) then
                    THROW_HARD('read_tabs_to_glob: sparse partition row has no edges')
                endif
                degree_by_ptcl(global_i) = row_deg
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
            local_e                      = fill_ptr(global_i)
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
        integer :: funit, io_stat, nptcls_glob, headsz, local_idx, global_idx, pind
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
        if (allocated(self%proj_ref_offsets))        deallocate(self%proj_ref_offsets)
        if (allocated(self%proj_ref_indices))        deallocate(self%proj_ref_indices)
        if (allocated(self%active_proj_indices))     deallocate(self%active_proj_indices)
        if (allocated(self%ptcl_off))                deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index))          deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))               deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))                deallocate(self%edge_val)
        if (allocated(self%ref_edge_offsets))        deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices))        deallocate(self%ref_edge_indices)
        self%nptcls      = 0
        self%nstates     = 0
        self%nrefs       = 0
        self%nactive_proj = 0
        self%nedges      = 0
        self%maxdeg_ptcl = 0
        self%b_ptr       => null()
        self%p_ptr       => null()
    end subroutine kill

    ! Private helper functions for managing dynamic integer buffers used in reference edge lists

    ! Returns .true. if a distance value has been scored (i.e. was not left at the sentinel huge(1.) value).
    elemental logical function dist_is_scored(d)
        real, intent(in) :: d
        dist_is_scored = (d < huge(d) / 2.)
    end function dist_is_scored

    ! Return the 1-based position of val within buf(1:nused), or 0 if not present.
    integer function find_int_buf(buf, val, nused) result(pos)
        integer,           intent(in) :: buf(:)
        integer,           intent(in) :: val
        integer, optional, intent(in) :: nused
        integer :: i, nscan
        pos = 0
        nscan = size(buf)
        if (present(nused)) nscan = min(max(nused, 0), nscan)
        do i = 1, nscan
            if (buf(i) == val) then
                pos = i
                return
            endif
        enddo
    end function find_int_buf

    ! Insert or update one preferred ref entry, keeping the smaller dist on duplicates.
    subroutine append_or_improve_ref(cached_ref_idx, cached_edge_val, nused, ref_idx, v)
        integer,        allocatable, intent(inout) :: cached_ref_idx(:)
        type(ptcl_ref), allocatable, intent(inout) :: cached_edge_val(:)
        integer,                     intent(inout) :: nused
        integer,                     intent(in)    :: ref_idx
        type(ptcl_ref),              intent(in)    :: v
        integer,        allocatable :: tmp_ref(:)
        type(ptcl_ref), allocatable :: tmp_val(:)
        integer :: pos, used_count, cap, new_cap
        if (allocated(cached_ref_idx) .neqv. allocated(cached_edge_val)) then
            THROW_HARD('append_or_improve_ref: cached_ref_idx/cached_edge_val allocation mismatch')
        endif
        used_count = 0
        if (allocated(cached_ref_idx)) used_count = min(max(nused, 0), min(size(cached_ref_idx), size(cached_edge_val)))
        pos = 0
        if (used_count > 0) pos = find_int_buf(cached_ref_idx, ref_idx, used_count)
        if (pos > 0) then
            if (v%dist < cached_edge_val(pos)%dist) cached_edge_val(pos) = v
            return
        endif
        cap = 0
        if (allocated(cached_ref_idx)) cap = min(size(cached_ref_idx), size(cached_edge_val))
        if (used_count < cap) then
            cached_ref_idx(used_count + 1) = ref_idx
            cached_edge_val(used_count + 1) = v
            nused = used_count + 1
            return
        endif
        new_cap = max(1, max(cap * 2, used_count + 1))
        allocate(tmp_ref(new_cap), tmp_val(new_cap))
        if (used_count > 0) then
            tmp_ref(1:used_count) = cached_ref_idx(1:used_count)
            tmp_val(1:used_count) = cached_edge_val(1:used_count)
        endif
        tmp_ref(used_count + 1) = ref_idx
        tmp_val(used_count + 1) = v
        call move_alloc(tmp_ref, cached_ref_idx)
        call move_alloc(tmp_val, cached_edge_val)
        nused = used_count + 1
    end subroutine append_or_improve_ref

end module simple_eul_prob_tab_neigh
