! Refactored direct-CSR version of simple_eul_prob_tab_neigh.
!
! Main changes:
!   * neighborhood builders now emit sparse rows directly
!   * no global logical mask(:,:) is built
!   * ptree neighborhood construction caches already-evaluated refs directly
!   * fill_tab_sparse skips refs whose edge_val(:)%dist is already populated
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

type :: sparse_row_tmp
    integer,        allocatable :: ref_idx(:)
    type(ptcl_ref), allocatable :: val(:)
end type sparse_row_tmp

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
    integer                     :: nedges      = 0
    integer                     :: maxdeg_ptcl = 0
    integer,        allocatable :: ptcl_off(:)
    integer,        allocatable :: edge_ref_index(:)
    integer,        allocatable :: edge_ptcl(:)
    type(ptcl_ref), allocatable :: edge_val(:)
    integer,        allocatable :: ref_edge_offsets(:)
    integer,        allocatable :: ref_edge_indices(:)
contains
    procedure          :: new
    procedure, private :: init_common
    procedure, private :: build_ref_lists_and_map
    procedure, private :: check_subspace_prereqs
    procedure, private :: init_edge_default
    procedure, private :: collect_all_active_projs
    procedure, private :: materialize_row_from_proj_list
    procedure, private :: flatten_sparse_rows
    procedure, private :: build_graph_from_prev_geom
    procedure, private :: build_graph_from_subspace_peaks
    procedure, private :: build_graph_from_subspace_peaks_one_state
    procedure, private :: build_graph_from_subspace_peaks_impl
    procedure, private :: build_graph_from_ptree_srch
    procedure, private :: build_graph_from_ptree_srch_one_state
    procedure, private :: build_graph_from_ptree_srch_impl
    procedure, private :: build_ref_adjacency
    procedure          :: fill_tab   => fill_tab_sparse
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

contains

    ! Initialize object state and optionally build the sparse neighborhood graph.
    subroutine new(self, params, build, pinds, neigh_type, empty_okay, state, build_sparse_graph)
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
                    THROW_HARD('new: unsupported neigh_type='//trim(neigh_type)//'; expected geom, subspace_srch or ptree_srch')
            end select
            call self%build_ref_adjacency
        endif
    end subroutine new

    ! Initialize common per-particle and per-reference metadata used by all workflows.
    subroutine init_common(self, params, build, pinds, empty_okay)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer, parameter :: MIN_POP = 5
        integer :: ptcl_local_idx, iptcl
        logical :: l_empty
        real    :: x
        call self%kill
        self%p_ptr => params
        self%b_ptr => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        !$omp parallel do default(shared) private(ptcl_local_idx,iptcl) proc_bind(close) schedule(static)
        do ptcl_local_idx = 1, self%nptcls
            iptcl = self%pinds(ptcl_local_idx)
            self%assgn_map(ptcl_local_idx)%pind   = iptcl
            self%assgn_map(ptcl_local_idx)%istate = 0
            self%assgn_map(ptcl_local_idx)%iproj  = 0
            self%assgn_map(ptcl_local_idx)%inpl   = 0
            self%assgn_map(ptcl_local_idx)%dist   = huge(x)
            self%assgn_map(ptcl_local_idx)%x      = 0.
            self%assgn_map(ptcl_local_idx)%y      = 0.
            self%assgn_map(ptcl_local_idx)%has_sh = .false.
        enddo
        !$omp end parallel do
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

    ! Build compressed reference lists and (iproj,istate)->ref index lookup map.
    subroutine build_ref_lists_and_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: state_list_idx, ref_idx, istate, iproj
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
        real :: x
        v%pind   = iptcl
        v%istate = self%ref_state_indices(ref_idx)
        v%iproj  = self%ref_proj_indices(ref_idx)
        v%inpl   = 0
        v%dist   = huge(x)
        v%x      = 0.
        v%y      = 0.
        v%has_sh = .false.
    end subroutine init_edge_default

    ! Collect all projections that have at least one active state.
    subroutine collect_all_active_projs(self, proj_idx, nproj)
        class(eul_prob_tab_neigh), intent(in)    :: self
        integer, allocatable,      intent(inout) :: proj_idx(:)
        integer,                   intent(out)   :: nproj
        integer :: iproj
        if (allocated(proj_idx)) deallocate(proj_idx)
        nproj = count(self%proj_active_state_count > 0)
        if (nproj <= 0) then
            THROW_HARD('collect_all_active_projs: no active projections available')
        endif
        allocate(proj_idx(nproj))
        nproj = 0
        do iproj = 1, self%p_ptr%nspace
            if (self%proj_active_state_count(iproj) <= 0) cycle
            nproj = nproj + 1
            proj_idx(nproj) = iproj
        enddo
    end subroutine collect_all_active_projs

    ! Expand a projection list into one sparse row over all active states.
    subroutine materialize_row_from_proj_list(self, iptcl, proj_idx, nproj, row, pref_ref, pref_val, npref)
        class(eul_prob_tab_neigh),             intent(in)    :: self
        integer,                               intent(in)    :: iptcl, nproj
        integer,                               intent(in)    :: proj_idx(:)
        type(sparse_row_tmp),                  intent(inout) :: row
        integer,        allocatable, optional, intent(in)    :: pref_ref(:)
        type(ptcl_ref), allocatable, optional, intent(in)    :: pref_val(:)
        integer,         optional, intent(in)    :: npref
        integer :: i, j, nrow, iproj, istate, ref_idx, slot, npref_use
        if (allocated(row%ref_idx)) deallocate(row%ref_idx)
        if (allocated(row%val))     deallocate(row%val)
        nrow = 0
        do i = 1, nproj
            iproj = proj_idx(i)
            if (iproj < 1 .or. iproj > self%p_ptr%nspace) cycle
            nrow = nrow + self%proj_active_state_count(iproj)
        enddo
        if (nrow <= 0) then
            THROW_HARD('materialize_row_from_proj_list: zero-row materialization')
        endif
        allocate(row%ref_idx(nrow), row%val(nrow))
        nrow = 0
        npref_use = 0
        if (present(npref)) npref_use = npref
        do i = 1, nproj
            iproj = proj_idx(i)
            if (iproj < 1 .or. iproj > self%p_ptr%nspace) cycle
            if (self%proj_active_state_count(iproj) <= 0) cycle
            do j = 1, self%nstates
                istate = self%active_state_indices(j)
                ref_idx = self%ref_index_map(iproj, istate)
                if (ref_idx <= 0) cycle
                nrow = nrow + 1
                row%ref_idx(nrow) = ref_idx
                slot = 0
                if (present(pref_ref) .and. present(pref_val)) then
                    if (allocated(pref_ref) .and. allocated(pref_val)) then
                        if (npref_use > 0) slot = find_int_buf(pref_ref, npref_use, ref_idx)
                    endif
                endif
                if (slot > 0) then
                    row%val(nrow) = pref_val(slot)
                else
                    call self%init_edge_default(iptcl, ref_idx, row%val(nrow))
                endif
            enddo
        enddo
    end subroutine materialize_row_from_proj_list

    ! Flatten per-particle temporary rows into particle-major CSR edge arrays.
    subroutine flatten_sparse_rows(self, rows)
        class(eul_prob_tab_neigh),         intent(inout) :: self
        type(sparse_row_tmp), allocatable, intent(inout) :: rows(:)
        integer :: i, j, e, nrow
        if (allocated(self%ptcl_off))       deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index)) deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))      deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))       deallocate(self%edge_val)
        allocate(self%ptcl_off(self%nptcls + 1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            nrow = 0
            if (allocated(rows(i)%ref_idx)) nrow = size(rows(i)%ref_idx)
            if (nrow <= 0) then
                THROW_HARD('flatten_sparse_rows: particle has zero valid neighbors')
            endif
            self%ptcl_off(i+1) = self%ptcl_off(i) + nrow
        enddo
        self%nedges = self%ptcl_off(self%nptcls + 1) - 1
        self%maxdeg_ptcl = maxval(self%ptcl_off(2:self%nptcls+1) - self%ptcl_off(1:self%nptcls))
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        do i = 1, self%nptcls
            e = self%ptcl_off(i)
            do j = 1, size(rows(i)%ref_idx)
                self%edge_ref_index(e) = rows(i)%ref_idx(j)
                self%edge_ptcl(e)      = i
                self%edge_val(e)       = rows(i)%val(j)
                e = e + 1
            enddo
        enddo
        do i = 1, size(rows)
            if (allocated(rows(i)%ref_idx)) deallocate(rows(i)%ref_idx)
            if (allocated(rows(i)%val))     deallocate(rows(i)%val)
        enddo
        deallocate(rows)
    end subroutine flatten_sparse_rows

    ! Build graph from geometric neighbors around previous particle orientations.
    subroutine build_graph_from_prev_geom(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(sparse_row_tmp), allocatable :: rows(:)
        type(eulspace_neigh_map)          :: neigh_map
        type(ori) :: o_prev, o_sub, osym
        integer   :: nspace_sub, npeak_use
        integer   :: i, iptcl, isub, iproj_full, prev_iproj, prev_sub
        integer   :: j, nproj
        real      :: dtmp, inplrotdist
        nspace_sub = self%p_ptr%nspace_sub
        npeak_use  = max(1, min(self%p_ptr%npeaks, nspace_sub))
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(rows(self%nptcls))
        !$omp parallel do default(shared) private(i,iptcl,o_prev,isub,iproj_full,dtmp,inplrotdist,o_sub,osym,prev_iproj,prev_sub,j,nproj) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            block
                integer, allocatable :: peak_sub_idxs(:), proj_sel(:), neigh_proj(:)
                real,    allocatable :: sub_dists(:)
                iptcl = self%pinds(i)
                allocate(sub_dists(nspace_sub), peak_sub_idxs(npeak_use))
                sub_dists = huge(1.)
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    if (iproj_full < 1 .or. iproj_full > self%p_ptr%nspace) then
                        THROW_HARD('build_graph_from_prev_geom: representative projection index out of range')
                    endif
                    call self%b_ptr%eulspace%get_ori(iproj_full, o_sub)
                    call self%b_ptr%pgrpsyms%sym_dists(o_prev, o_sub, osym, dtmp, inplrotdist)
                    sub_dists(isub) = dtmp
                    call o_sub%kill
                    call osym%kill
                enddo
                peak_sub_idxs = minnloc(sub_dists, npeak_use)
                nproj = 0
                do j = 1, npeak_use
                    call neigh_map%get_neighbors_list(peak_sub_idxs(j), neigh_proj)
                    if (allocated(neigh_proj)) then
                        call append_unique_int_list(proj_sel, nproj, neigh_proj)
                        deallocate(neigh_proj)
                    endif
                enddo
                prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
                if (prev_iproj >= 1 .and. prev_iproj <= self%p_ptr%nspace) then
                    prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
                    if (prev_sub >= 1 .and. prev_sub <= nspace_sub) then
                        call neigh_map%get_neighbors_list(prev_sub, neigh_proj)
                        if (allocated(neigh_proj)) then
                            call append_unique_int_list(proj_sel, nproj, neigh_proj)
                            deallocate(neigh_proj)
                        endif
                    endif
                endif
                if (nproj <= 0) then
                    call self%collect_all_active_projs(proj_sel, nproj)
                else if (.not. any(self%proj_active_state_count(proj_sel(1:nproj)) > 0)) then
                    call self%collect_all_active_projs(proj_sel, nproj)
                endif
                call self%materialize_row_from_proj_list(iptcl, proj_sel, nproj, rows(i))
                call o_prev%kill
            end block
        enddo
        !$omp end parallel do
        call neigh_map%kill
        call self%flatten_sparse_rows(rows)
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
        type(eulspace_neigh_map) :: neigh_map
        type(pftc_shsrch_grad)   :: grad_shsrch_obj(nthr_glob)
        type(ori)                :: o_prev
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_shift(2), huge_dist
        logical :: do_shift_first, has_state_opt
        integer :: state_fixed, nrots, nspace_sub, npeak_use
        integer :: i, isub, ithr, iptcl, istate, iproj_full, irot, iref_start, iref_prev
        integer :: state_i, nvalid_sub, prev_sub, prev_iproj, j, nproj
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
        allocate(rows(self%nptcls), inpl_dists(nrots, nthr_glob), coarse_best_dist(nspace_sub, nthr_glob), peak_sub_idxs(npeak_use, nthr_glob))
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
        !$omp parallel do default(shared) private(i,ithr,iptcl,o_prev,istate,isub,iproj_full,irot,iref_start,iref_prev,cxy,cxy_shift,state_i,nvalid_sub,prev_sub,prev_iproj,j,nproj) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            block
                integer, allocatable :: proj_sel(:), neigh_proj(:)
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                cxy_shift = [0., 0.]
                prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
                prev_sub = 0
                if (prev_iproj >= 1 .and. prev_iproj <= self%p_ptr%nspace) then
                    prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
                endif
                if (has_state_opt) then
                    istate = state_fixed
                    if (do_shift_first) then
                        if (prev_iproj >= 1 .and. prev_iproj <= self%p_ptr%nspace) then
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
                    endif
                else
                    istate = o_prev%get_state()
                    if (do_shift_first) then
                        if (istate >= 1 .and. istate <= self%p_ptr%nstates) then
                            if (self%state_exists(istate)) then
                                if (prev_iproj >= 1 .and. prev_iproj <= self%p_ptr%nspace) then
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
                            endif
                        endif
                    endif
                endif
                coarse_best_dist(:,ithr) = huge_dist
                if (has_state_opt) then
                    do isub = 1, nspace_sub
                        iproj_full = self%b_ptr%subspace_inds(isub)
                        if (iproj_full < 1 .or. iproj_full > self%p_ptr%nspace) then
                            THROW_HARD('build_graph_from_subspace_peaks_impl: representative projection index out of range')
                        endif
                        if (.not. self%proj_exists(iproj_full, state_fixed)) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((state_fixed - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                        inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), self%p_ptr%cc_objfun)
                        irot = minloc(inpl_dists(:,ithr), dim=1)
                        coarse_best_dist(isub,ithr) = inpl_dists(irot,ithr)
                    enddo
                else
                    do isub = 1, nspace_sub
                        iproj_full = self%b_ptr%subspace_inds(isub)
                        if (iproj_full < 1 .or. iproj_full > self%p_ptr%nspace) then
                            THROW_HARD('build_graph_from_subspace_peaks_impl: representative projection index out of range')
                        endif
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
                nvalid_sub = count(coarse_best_dist(:,ithr) < huge_dist / 2.)
                if (nvalid_sub > 0) then
                    if (nvalid_sub <= npeak_use) then
                        do isub = 1, nspace_sub
                            if (coarse_best_dist(isub,ithr) >= huge_dist / 2.) cycle
                            call neigh_map%get_neighbors_list(isub, neigh_proj)
                            if (allocated(neigh_proj)) then
                                call append_unique_int_list(proj_sel, nproj, neigh_proj)
                                deallocate(neigh_proj)
                            endif
                        enddo
                    else
                        peak_sub_idxs(:,ithr) = minnloc(coarse_best_dist(:,ithr), npeak_use)
                        do j = 1, npeak_use
                            call neigh_map%get_neighbors_list(peak_sub_idxs(j,ithr), neigh_proj)
                            if (allocated(neigh_proj)) then
                                call append_unique_int_list(proj_sel, nproj, neigh_proj)
                                deallocate(neigh_proj)
                            endif
                        enddo
                    endif
                endif
                if (prev_sub >= 1 .and. prev_sub <= nspace_sub) then
                    call neigh_map%get_neighbors_list(prev_sub, neigh_proj)
                    if (allocated(neigh_proj)) then
                        call append_unique_int_list(proj_sel, nproj, neigh_proj)
                        deallocate(neigh_proj)
                    endif
                endif
                if (nproj <= 0) then
                    call self%collect_all_active_projs(proj_sel, nproj)
                else if (.not. any(self%proj_active_state_count(proj_sel(1:nproj)) > 0)) then
                    call self%collect_all_active_projs(proj_sel, nproj)
                endif
                call self%materialize_row_from_proj_list(iptcl, proj_sel, nproj, rows(i))
                call o_prev%kill
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
        type(sparse_row_tmp), allocatable :: rows(:)
        integer,              allocatable :: peak_trees(:,:), proj_sel(:), pref_ref(:)
        real,                 allocatable :: inpl_corrs(:,:), tree_best_corrs(:,:), peak_tree_corrs(:,:), dists_sorted(:,:), inpl_angle_thres(:)
        type(ptcl_ref),       allocatable :: pref_val(:)
        integer,              allocatable :: inds_sorted(:,:)
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_shift(2)
        logical :: do_shift_first
        integer :: nrots, nspace_sub, npeak_use, ntrees
        integer :: i, isub, ithr, iptcl, istate, iproj_full, irot, iref_start, iref_prev, itree, state_i
        integer :: npeak_trees, ipeak, nproj, npref, prev_iproj
        real    :: corr_best
        call self%check_subspace_prereqs()
        nspace_sub     = self%p_ptr%nspace_sub
        ntrees         = self%b_ptr%block_tree%get_n_trees()
        if (ntrees <= 0) then
            THROW_HARD('build_graph_from_ptree_srch_impl: block_tree has no trees')
        endif
        npeak_use      = max(1, min(self%p_ptr%npeaks, ntrees))
        nrots          = self%b_ptr%pftc%get_nrots()

        print *, '************************ nrots = ', nrots

        do_shift_first = self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift
        allocate(rows(self%nptcls))
        allocate(inpl_corrs(nrots, nthr_glob), tree_best_corrs(ntrees, nthr_glob), &
                 peak_trees(npeak_use, nthr_glob), peak_tree_corrs(npeak_use, nthr_glob), &
                 dists_sorted(nrots, nthr_glob), inds_sorted(nrots, nthr_glob))
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
        !$omp parallel do default(shared) private(i,iptcl,ithr,cxy_shift,o_prev,istate,iref_start,iref_prev,irot,cxy,isub,iproj_full,itree,corr_best,state_i,npeak_trees,ipeak,nproj,npref,prev_iproj,proj_sel,pref_ref,pref_val) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            cxy_shift = [0., 0.]
            nproj = 0
            npref = 0
            if (allocated(proj_sel)) deallocate(proj_sel)
            if (allocated(pref_ref)) deallocate(pref_ref)
            if (allocated(pref_val)) deallocate(pref_val)
            if (do_shift_first) then
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
                istate     = o_prev%get_state()
                if (present(state_opt)) istate = state_opt
                if (istate >= 1 .and. istate <= self%p_ptr%nstates) then
                    if (self%state_exists(istate) .and. prev_iproj > 0) then
                        if (self%proj_exists(prev_iproj, istate)) then
                            iref_start = (istate - 1) * self%p_ptr%nspace
                            iref_prev  = iref_start + prev_iproj
                            irot       = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                            call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                            cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                            cxy_shift = 0.
                            if (irot /= 0) cxy_shift = cxy(2:3)
                        endif
                    endif
                endif
                call o_prev%kill
            endif
            tree_best_corrs(:,ithr) = -huge(1.0)
            do isub = 1, nspace_sub
                iproj_full = self%b_ptr%subspace_inds(isub)
                if (iproj_full < 1 .or. iproj_full > self%p_ptr%nspace) then
                    THROW_HARD('build_graph_from_ptree_srch_impl: representative projection index out of range')
                endif
                itree = self%b_ptr%subspace_full2sub_map(iproj_full)
                if (itree < 1 .or. itree > ntrees) then
                    THROW_HARD('build_graph_from_ptree_srch_impl: tree index out of range')
                endif
                if (present(state_opt)) then
                    if (.not. self%proj_exists(iproj_full, state_opt)) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((state_opt - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                    corr_best = maxval(inpl_corrs(:,ithr))
                    tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                else
                    do state_i = 1, self%nstates
                        istate = self%active_state_indices(state_i)
                        if (.not. self%proj_exists(iproj_full, istate)) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                        corr_best = maxval(inpl_corrs(:,ithr))
                        tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                    enddo
                endif
            enddo
            peak_trees(:,ithr)      = 0
            peak_tree_corrs(:,ithr) = -huge(1.0)
            call select_peak_trees(tree_best_corrs(:,ithr), peak_trees(:,ithr), peak_tree_corrs(:,ithr), npeak_trees)
            do ipeak = 1, npeak_trees
                call trace_tree_prob(self%b_ptr%block_tree, peak_trees(ipeak,ithr), ithr, score_ptree_ref)
            enddo
            if (nproj <= 0) call self%collect_all_active_projs(proj_sel, nproj)
            call self%materialize_row_from_proj_list(iptcl, proj_sel, nproj, rows(i), pref_ref, pref_val, npref)
            if (allocated(proj_sel)) deallocate(proj_sel)
            if (allocated(pref_ref)) deallocate(pref_ref)
            if (allocated(pref_val)) deallocate(pref_val)
        enddo
        !$omp end parallel do
        if (do_shift_first) then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call self%flatten_sparse_rows(rows)
        if (allocated(inpl_corrs))       deallocate(inpl_corrs)
        if (allocated(tree_best_corrs))  deallocate(tree_best_corrs)
        if (allocated(peak_trees))       deallocate(peak_trees)
        if (allocated(peak_tree_corrs))  deallocate(peak_tree_corrs)
        if (allocated(dists_sorted))     deallocate(dists_sorted)
        if (allocated(inds_sorted))      deallocate(inds_sorted)
        if (allocated(inpl_angle_thres)) deallocate(inpl_angle_thres)

    contains

        subroutine score_ptree_ref(iproj_eval, ithr_eval, best_corr)
            integer, intent(in)  :: iproj_eval
            integer, intent(in)  :: ithr_eval
            real,    intent(out) :: best_corr
            integer :: istate_loc, state_j, irot_loc, ref_idx_loc
            real    :: dist_obj(nrots), rotmat_loc(2,2), rot_xy_loc(2)
            type(ptcl_ref) :: v
            logical :: have_valid
            best_corr  = -huge(1.0)
            have_valid = .false.
            if (iproj_eval < 1 .or. iproj_eval > self%p_ptr%nspace) return
            if (present(state_opt)) then
                if (.not. self%proj_exists(iproj_eval, state_opt)) return
                call self%b_ptr%pftc%gen_objfun_vals((state_opt - 1) * self%p_ptr%nspace + iproj_eval, iptcl, cxy_shift, inpl_corrs(:,ithr_eval))
                best_corr = maxval(inpl_corrs(:,ithr_eval))
                dist_obj = eulprob_dist_switch(inpl_corrs(:,ithr_eval), self%p_ptr%cc_objfun)
                irot_loc = angle_sampling(dist_obj, dists_sorted(:,ithr_eval), inds_sorted(:,ithr_eval), inpl_angle_thres(state_opt), self%p_ptr%prob_athres)
                ref_idx_loc = self%ref_index_map(iproj_eval, state_opt)
                if (ref_idx_loc > 0) then
                    call self%init_edge_default(iptcl, ref_idx_loc, v)
                    v%dist = dist_obj(irot_loc)
                    v%inpl = irot_loc
                    if (do_shift_first) then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                        rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                        v%x = rot_xy_loc(1)
                        v%y = rot_xy_loc(2)
                        v%has_sh = .true.
                    endif
                    call append_or_improve_ref(pref_ref, pref_val, npref, ref_idx_loc, v)
                    have_valid = .true.
                endif
            else
                do state_j = 1, self%nstates
                    istate_loc = self%active_state_indices(state_j)
                    if (.not. self%proj_exists(iproj_eval, istate_loc)) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((istate_loc - 1) * self%p_ptr%nspace + iproj_eval, iptcl, cxy_shift, inpl_corrs(:,ithr_eval))
                    best_corr = max(best_corr, maxval(inpl_corrs(:,ithr_eval)))
                    dist_obj = eulprob_dist_switch(inpl_corrs(:,ithr_eval), self%p_ptr%cc_objfun)
                    irot_loc = angle_sampling(dist_obj, dists_sorted(:,ithr_eval), inds_sorted(:,ithr_eval), inpl_angle_thres(istate_loc), self%p_ptr%prob_athres)
                    ref_idx_loc = self%ref_index_map(iproj_eval, istate_loc)
                    if (ref_idx_loc <= 0) cycle
                    call self%init_edge_default(iptcl, ref_idx_loc, v)
                    v%dist = dist_obj(irot_loc)
                    v%inpl = irot_loc
                    if (do_shift_first) then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                        rot_xy_loc = matmul(cxy_shift, rotmat_loc)
                        v%x = rot_xy_loc(1)
                        v%y = rot_xy_loc(2)
                        v%has_sh = .true.
                    endif
                    call append_or_improve_ref(pref_ref, pref_val, npref, ref_idx_loc, v)
                    have_valid = .true.
                enddo
            endif

            if (have_valid) call append_unique_int(proj_sel, nproj, iproj_eval)
        end subroutine score_ptree_ref

    end subroutine build_graph_from_ptree_srch_impl

    ! Build reference-major adjacency over edge indices for fast per-reference scans.
    subroutine build_ref_adjacency(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: edge_count_by_ref(:), ref_write_ptr(:)
        integer :: ref_idx, edge_idx
        allocate(edge_count_by_ref(self%nrefs), source=0)
        do edge_idx = 1, self%nedges
            edge_count_by_ref(self%edge_ref_index(edge_idx)) = edge_count_by_ref(self%edge_ref_index(edge_idx)) + 1
        enddo
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

    ! Evaluate and optionally refine all unresolved sparse edges.
    subroutine fill_tab_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori)              :: o_prev
        real,    allocatable   :: inpl_angle_thres(:)
        real,    allocatable   :: dists_inpl(:,:), dists_inpl_sorted(:,:), candidate_dist(:,:)
        integer, allocatable   :: inds_sorted(:,:), candidate_edge(:,:)
        integer :: nrots, i, e, ithr, iptcl, ri, istate, iproj, irot, iref_full, iref_start
        integer :: n_candidates, candidate_i, n_refine, state_i, refine_rank
        integer :: n_refs_to_refine, n_samples
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rotmat(2,2), rot_xy(2), huge_dist
        logical :: edge_ready
        call seed_rnd
        nrots = self%b_ptr%pftc%get_nrots()
        huge_dist = huge(1.)
        allocate(inpl_angle_thres(self%p_ptr%nstates), source=0.)
        n_refs_to_refine = 0
        do state_i = 1, self%nstates
            istate = self%active_state_indices(state_i)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n_samples, self%p_ptr%prob_athres, state=istate)
            n_refs_to_refine = max(n_refs_to_refine, n_samples)
            inpl_angle_thres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        allocate(dists_inpl(nrots, nthr_glob), dists_inpl_sorted(nrots, nthr_glob), inds_sorted(nrots, nthr_glob))
        allocate(candidate_dist(self%maxdeg_ptcl, nthr_glob), candidate_edge(self%maxdeg_ptcl, nthr_glob))
        if (self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift) then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,istate,iproj,irot,iref_start,cxy,rotmat,rot_xy,n_candidates,e,ri,candidate_i,n_refine,refine_rank,cxy_prob,iref_full,edge_ready) &
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                istate = o_prev%get_state()
                iproj  = self%b_ptr%eulspace%find_closest_proj(o_prev)
                irot   = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                if (istate >= 1 .and. istate <= self%p_ptr%nstates .and. iproj >= 1 .and. iproj <= self%p_ptr%nspace) then
                    if (self%ref_index_map(iproj,istate) > 0) then
                        iref_start = (istate - 1) * self%p_ptr%nspace
                        call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        if (irot == 0) cxy(2:3) = 0.
                    else
                        cxy(2:3) = 0.
                    endif
                else
                    cxy(2:3) = 0.
                endif
                n_candidates = self%ptcl_off(i+1) - self%ptcl_off(i)
                candidate_i = 0
                do e = self%ptcl_off(i), self%ptcl_off(i+1)-1
                    edge_ready = (self%edge_val(e)%dist < huge_dist / 2.)
                    if (.not. edge_ready) then
                        ri     = self%edge_ref_index(e)
                        istate = self%ref_state_indices(ri)
                        iproj  = self%ref_proj_indices(ri)
                        iref_start = (istate - 1) * self%p_ptr%nspace
                        iref_full  = iref_start + iproj
                        call self%b_ptr%pftc%gen_objfun_vals(iref_full, iptcl, cxy(2:3), dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                        irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), &
                                              inpl_angle_thres(istate), self%p_ptr%prob_athres)
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                        rot_xy = matmul(cxy(2:3), rotmat)
                        self%edge_val(e)%dist   = dists_inpl(irot,ithr)
                        self%edge_val(e)%inpl   = irot
                        self%edge_val(e)%x      = rot_xy(1)
                        self%edge_val(e)%y      = rot_xy(2)
                        self%edge_val(e)%has_sh = .true.
                    endif
                    candidate_i = candidate_i + 1
                    candidate_dist(candidate_i,ithr) = self%edge_val(e)%dist
                    candidate_edge(candidate_i,ithr) = e
                enddo
                call o_prev%kill
                if (self%p_ptr%l_prob_sh) then
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
                            cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=cxy(2:3))
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
        else
            if (self%p_ptr%l_prob_sh .and. self%p_ptr%l_doshift) then
                lims(:,1)      = -self%p_ptr%trs
                lims(:,2)      =  self%p_ptr%trs
                lims_init(:,1) = -SHC_INPL_TRSHWDTH
                lims_init(:,2) =  SHC_INPL_TRSHWDTH
                do ithr = 1, nthr_glob
                    call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                                   maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
                enddo
                !$omp parallel do default(shared) private(i,iptcl,ithr,n_candidates,e,ri,istate,iproj,irot,candidate_i,n_refine,refine_rank,cxy,iref_full,iref_start,edge_ready) &
                !$omp proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    n_candidates = self%ptcl_off(i+1) - self%ptcl_off(i)
                    candidate_i = 0
                    do e = self%ptcl_off(i), self%ptcl_off(i+1)-1
                        edge_ready = (self%edge_val(e)%dist < huge_dist / 2.)
                        if (.not. edge_ready) then
                            ri     = self%edge_ref_index(e)
                            istate = self%ref_state_indices(ri)
                            iproj  = self%ref_proj_indices(ri)
                            iref_start = (istate - 1) * self%p_ptr%nspace
                            iref_full  = iref_start + iproj
                            call self%b_ptr%pftc%gen_objfun_vals(iref_full, iptcl, [0.,0.], dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), &
                                                  inpl_angle_thres(istate), self%p_ptr%prob_athres)
                            self%edge_val(e)%dist = dists_inpl(irot,ithr)
                            self%edge_val(e)%inpl = irot
                        endif
                        candidate_i = candidate_i + 1
                        candidate_dist(candidate_i,ithr) = self%edge_val(e)%dist
                        candidate_edge(candidate_i,ithr) = e
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
                            cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                            if (irot > 0) then
                                self%edge_val(e)%inpl = irot
                                self%edge_val(e)%dist = eulprob_dist_switch(cxy(1), self%p_ptr%cc_objfun)
                                self%edge_val(e)%x    = cxy(2)
                                self%edge_val(e)%y    = cxy(3)
                            endif
                            self%edge_val(e)%has_sh = .true.
                        enddo
                    endif
                enddo
                !$omp end parallel do
            else
                !$omp parallel do default(shared) private(i,iptcl,ithr,e,ri,istate,iproj,irot,iref_full,iref_start,edge_ready) &
                !$omp proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    do e = self%ptcl_off(i), self%ptcl_off(i+1)-1
                        edge_ready = (self%edge_val(e)%dist < huge_dist / 2.)
                        if (.not. edge_ready) then
                            ri     = self%edge_ref_index(e)
                            istate = self%ref_state_indices(ri)
                            iproj  = self%ref_proj_indices(ri)
                            iref_start = (istate - 1) * self%p_ptr%nspace
                            iref_full  = iref_start + iproj
                            call self%b_ptr%pftc%gen_objfun_vals(iref_full, iptcl, [0.,0.], dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), &
                                                  inpl_angle_thres(istate), self%p_ptr%prob_athres)
                            self%edge_val(e)%dist = dists_inpl(irot,ithr)
                            self%edge_val(e)%inpl = irot
                        endif
                    enddo
                enddo
                !$omp end parallel do
            endif
        endif
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        enddo
        deallocate(inpl_angle_thres, dists_inpl, dists_inpl_sorted, inds_sorted, candidate_dist, candidate_edge)
    end subroutine fill_tab_sparse

    ! Normalize per-particle edge distances and map them to a robust global range.
    subroutine ref_normalize_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ptcl_local_idx, edge_idx
        real    :: sum_dist, min_dist, max_dist
        !$omp parallel do default(shared) private(ptcl_local_idx,edge_idx,sum_dist) proc_bind(close) schedule(static)
        do ptcl_local_idx = 1, self%nptcls
            sum_dist = 0.
            do edge_idx = self%ptcl_off(ptcl_local_idx), self%ptcl_off(ptcl_local_idx+1)-1
                sum_dist = sum_dist + self%edge_val(edge_idx)%dist
            enddo
            if (sum_dist < TINY) then
                do edge_idx = self%ptcl_off(ptcl_local_idx), self%ptcl_off(ptcl_local_idx+1)-1
                    self%edge_val(edge_idx)%dist = 0.
                enddo
            else
                do edge_idx = self%ptcl_off(ptcl_local_idx), self%ptcl_off(ptcl_local_idx+1)-1
                    self%edge_val(edge_idx)%dist = self%edge_val(edge_idx)%dist / sum_dist
                enddo
            endif
        enddo
        !$omp end parallel do
        min_dist = minval(self%edge_val(:)%dist)
        max_dist = maxval(self%edge_val(:)%dist)
        if ((max_dist - min_dist) < TINY) then
            !$omp parallel do default(shared) private(edge_idx) proc_bind(close) schedule(static)
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
        !$omp parallel do default(shared) private(ref_idx,ref_first_ptr,ref_last_ptr,n_ref_edges) proc_bind(close) schedule(static)
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
         type(ptcl_ref), allocatable :: edge_val_loc(:)
        integer,         allocatable :: degree_by_ptcl(:), fill_ptr(:), pind_to_local(:), pinds_loc(:), ptcl_off_loc(:)
        type(string) :: binfname
        integer      :: funit, addr, io_stat, file_header(3), header_bytes
        integer      :: ipart, nrefs_loc, nptcls_loc, nedges_loc, local_i, global_i, local_e, global_e, ref_idx
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
            allocate(pinds_loc(nptcls_loc), ptcl_off_loc(nptcls_loc + 1))
            addr = header_bytes + 1
            read(unit=funit, pos=addr) pinds_loc
            addr = addr + nptcls_loc * sizeof(pinds_loc(1))
            read(unit=funit, pos=addr) ptcl_off_loc
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
                    THROW_HARD('read_tabs_to_glob: duplicate particle across sparse partition files')
                endif
                degree_by_ptcl(global_i) = ptcl_off_loc(local_i + 1) - ptcl_off_loc(local_i)
                if (degree_by_ptcl(global_i) <= 0) then
                    THROW_HARD('read_tabs_to_glob: sparse partition row has no edges')
                endif
            enddo
            deallocate(pinds_loc, ptcl_off_loc)
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
        allocate(fill_ptr(self%nptcls), source=self%ptcl_off(1:self%nptcls))
        do ipart = 1, nparts
            binfname = fbody//int2str_pad(ipart, numlen)//'.dat'
            if (.not. file_exists(binfname)) then
                THROW_HARD('file '//binfname%to_char()//' does not exists!')
            else
                call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            endif
            call fileiochk('simple_eul_prob_tab_neigh; read_tabs_to_glob; file: '//binfname%to_char(), io_stat)
            read(unit=funit, pos=1) file_header
            nptcls_loc = file_header(2)
            nedges_loc = file_header(3)
            allocate(pinds_loc(nptcls_loc), ptcl_off_loc(nptcls_loc + 1), edge_val_loc(nedges_loc))
            addr = header_bytes + 1
            read(unit=funit, pos=addr) pinds_loc
            addr = addr + nptcls_loc * sizeof(pinds_loc(1))
            read(unit=funit, pos=addr) ptcl_off_loc
            addr = addr + (nptcls_loc + 1) * sizeof(ptcl_off_loc(1))
            read(unit=funit, pos=addr) edge_val_loc
            call fclose(funit)
            do local_i = 1, nptcls_loc
                global_i = pind_to_local(pinds_loc(local_i))
                global_e = fill_ptr(global_i)
                do local_e = ptcl_off_loc(local_i), ptcl_off_loc(local_i + 1) - 1
                    if (edge_val_loc(local_e)%iproj < 1 .or. edge_val_loc(local_e)%iproj > self%p_ptr%nspace) then
                        THROW_HARD('read_tabs_to_glob: sparse partition edge has iproj out of range')
                    endif
                    if (edge_val_loc(local_e)%istate < 1 .or. edge_val_loc(local_e)%istate > self%p_ptr%nstates) then
                        THROW_HARD('read_tabs_to_glob: sparse partition edge has istate out of range')
                    endif
                    ref_idx = self%ref_index_map(edge_val_loc(local_e)%iproj, edge_val_loc(local_e)%istate)
                    if (ref_idx <= 0) then
                        THROW_HARD('read_tabs_to_glob: sparse partition edge refers to an inactive reference')
                    endif
                    self%edge_ref_index(global_e) = ref_idx
                    self%edge_ptcl(global_e)      = global_i
                    self%edge_val(global_e)       = edge_val_loc(local_e)
                    global_e = global_e + 1
                enddo
                fill_ptr(global_i) = global_e
            enddo
            deallocate(pinds_loc, ptcl_off_loc, edge_val_loc)
        enddo
        if (any(fill_ptr /= self%ptcl_off(2:self%nptcls + 1))) then
            THROW_HARD('read_tabs_to_glob: failed to reconstruct the sparse global table consistently')
        endif
        call self%build_ref_adjacency()
        deallocate(degree_by_ptcl, fill_ptr, pind_to_local)
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

    ! Helper functions for managing dynamic integer buffers used in reference edge lists

    ! Return the 1-based position of val within buf(1:n), or 0 if not present.
    integer function find_int_buf(buf, n, val) result(pos)
        integer, intent(in) :: buf(:)
        integer, intent(in) :: n, val
        integer :: i
        pos = 0
        do i = 1, n
            if (buf(i) == val) then
                pos = i
                return
            endif
        enddo
    end function find_int_buf

    ! Append a value to a dynamic integer buffer only if it is not already present.
    subroutine append_unique_int(buf, n, val)
        integer, allocatable, intent(inout) :: buf(:)
        integer,              intent(inout) :: n
        integer,              intent(in)    :: val
        integer, allocatable :: tmp(:)
        integer :: newcap
        if (allocated(buf)) then
            if (n > 0) then
                if (find_int_buf(buf, n, val) > 0) return
            endif
            if (n == size(buf)) then
                newcap = max(8, 2 * size(buf))
                allocate(tmp(newcap))
                if (n > 0) tmp(1:n) = buf(1:n)
                call move_alloc(tmp, buf)
            endif
        else
            allocate(buf(8))
        endif
        n = n + 1
        buf(n) = val
    end subroutine append_unique_int

    ! Append all values from vals into buf while preserving uniqueness.
    subroutine append_unique_int_list(buf, n, vals)
        integer, allocatable, intent(inout) :: buf(:)
        integer,              intent(inout) :: n
        integer,              intent(in)    :: vals(:)
        integer :: i
        do i = 1, size(vals)
            call append_unique_int(buf, n, vals(i))
        enddo
    end subroutine append_unique_int_list

    ! Insert or update one preferred ref entry, keeping the smaller dist on duplicates.
    subroutine append_or_improve_ref(pref_ref, pref_val, npref, ref_idx, v)
        integer,        allocatable, intent(inout) :: pref_ref(:)
        type(ptcl_ref), allocatable, intent(inout) :: pref_val(:)
        integer,                     intent(inout) :: npref
        integer,                     intent(in)    :: ref_idx
        type(ptcl_ref),              intent(in)    :: v
        integer, allocatable :: tmp_ref(:)
        type(ptcl_ref), allocatable :: tmp_val(:)
        integer :: pos, newcap
        if (allocated(pref_ref)) then
            pos = 0
            if (npref > 0) pos = find_int_buf(pref_ref, npref, ref_idx)
            if (pos > 0) then
                if (v%dist < pref_val(pos)%dist) pref_val(pos) = v
                return
            endif
            if (npref == size(pref_ref)) then
                newcap = max(8, 2 * size(pref_ref))
                allocate(tmp_ref(newcap), tmp_val(newcap))
                if (npref > 0) then
                    tmp_ref(1:npref) = pref_ref(1:npref)
                    tmp_val(1:npref) = pref_val(1:npref)
                endif
                call move_alloc(tmp_ref, pref_ref)
                call move_alloc(tmp_val, pref_val)
            endif
        else
            allocate(pref_ref(8), pref_val(8))
        endif
        npref = npref + 1
        pref_ref(npref) = ref_idx
        pref_val(npref) = v
    end subroutine append_or_improve_ref

end module simple_eul_prob_tab_neigh
