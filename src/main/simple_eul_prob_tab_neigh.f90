!> Neighborhood sparse probabilistic ref table + globally coupled assignment
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
    class(builder),    pointer :: b_ptr => null()
    class(parameters), pointer :: p_ptr => null()
    ! particles in this object
    integer                     :: nptcls = 0
    integer,        allocatable :: pinds(:)                 ! global ptcl indices for processing
    type(ptcl_ref), allocatable :: assgn_map(:)             ! assignment map (size nptcls)
    ! existence maps + compressed reference list
    logical,        allocatable :: proj_exists(:,:)         ! (nspace, nstates_total)
    logical,        allocatable :: state_exists(:)          ! (nstates_total)
    integer                     :: nstates = 0              ! count(state_exists)
    integer                     :: nrefs   = 0              ! count(proj_exists)
    integer,        allocatable :: active_state_indices(:)  ! (nstates) active state IDs
    integer,        allocatable :: ref_proj_indices(:)      ! (nrefs) projection index for each compressed reference
    integer,        allocatable :: ref_state_indices(:)     ! (nrefs) state index for each compressed reference
    ! map from (projection,state) -> compressed reference index in 1..nrefs, 0 if missing
    integer,        allocatable :: ref_index_map(:,:)       ! (nspace, nstates_total)
    integer,        allocatable :: proj_active_state_count(:) ! (nspace), number of active states available for each projection
    ! sparse candidate graph: particle-major CSR edges (CSR=Compressed Sparse Row)
    integer                     :: nedges      = 0
    integer                     :: maxdeg_ptcl = 0
    integer,        allocatable :: ptcl_off(:)              ! (nptcls+1), 1-based offsets into edge arrays
    integer,        allocatable :: edge_ref_index(:)        ! (nedges), compressed reference index per edge
    integer,        allocatable :: edge_ptcl(:)             ! (nedges), local particle index (1..nptcls)
    type(ptcl_ref), allocatable :: edge_val(:)              ! (nedges), stores dist/inpl/x/y/...
    ! reference-major adjacency: CSR of edge indices grouped by reference
    integer,        allocatable :: ref_edge_offsets(:)      ! (nrefs+1)
    integer,        allocatable :: ref_edge_indices(:)      ! (nedges), edge indices sorted by reference for fast access to all particles assigned to a given reference
contains
    procedure          :: new
    procedure          :: fill_tab       => fill_tab_sparse
    procedure          :: ref_assign     => ref_assign_sparse
    procedure          :: write_tab
    procedure          :: read_tabs_to_glob
    procedure          :: write_assignment
    procedure          :: read_assignment
    procedure          :: kill
    procedure, private :: init_common
    procedure, private :: build_neigh_mask_from_subspace_peaks
    procedure, private :: build_neigh_mask_from_subspace_peaks_one_state
    procedure, private :: build_neigh_mask_from_subspace_peaks_impl
    procedure, private :: build_neigh_mask_from_ptree_srch
    procedure, private :: build_neigh_mask_from_ptree_srch_one_state
    procedure, private :: build_neigh_mask_from_ptree_srch_impl
    procedure, private :: build_neigh_mask_from_prev_geom
    procedure, private :: check_subspace_prereqs
    procedure, private :: build_ref_lists_and_map
    procedure, private :: build_sparse_from_mask
    procedure, private :: build_ref_adjacency
    procedure, private :: ref_normalize_sparse
    procedure, private :: sort_ref_lists_by_dist
    procedure, private :: sort_eidx_by_dist
    procedure, private :: sift_down
end type eul_prob_tab_neigh

contains

    !===========================================================
    ! Constructor for the global driver object.
    ! Initializes the global particle/reference maps and builds
    ! sparse neighborhoods
    !===========================================================
    subroutine new(self, params, build, pinds, neigh_type, empty_okay, state, build_sparse_graph )
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
        if( present(build_sparse_graph) ) do_build_sparse_graph = build_sparse_graph
        call self%init_common(params, build, pinds, empty_okay)
        call self%build_ref_lists_and_map
        if( do_build_sparse_graph )then
            call self%check_subspace_prereqs
            select case(trim(neigh_type))
                case('geom')
                    call self%build_neigh_mask_from_prev_geom
                case('subspace_srch')
                    if( present(state) ) then
                        call self%build_neigh_mask_from_subspace_peaks_one_state(state)
                    else
                        call self%build_neigh_mask_from_subspace_peaks
                    endif
                case('ptree_srch')
                    if( present(state) ) then
                        call self%build_neigh_mask_from_ptree_srch_one_state(state)
                    else
                        call self%build_neigh_mask_from_ptree_srch
                    endif
                case default
                    THROW_HARD('new: unsupported neigh_type='//trim(neigh_type)//'; expected geom, subspace_srch or ptree_srch')
            end select
            call self%build_ref_adjacency
        endif
    end subroutine new

    !===========================================================
    ! Common constructor logic for both local-part and driver
    ! objects. Builds the global reference/state maps, but does
    ! not yet assemble any sparse neighborhood graph.
    !===========================================================
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
        ! init global reference/state existence maps and compressed reference lists/maps
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

    subroutine build_ref_lists_and_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: state_list_idx, ref_idx, istate, iproj
        allocate(self%active_state_indices(self%nstates), self%ref_proj_indices(self%nrefs), self%ref_state_indices(self%nrefs))
        allocate(self%ref_index_map(self%p_ptr%nspace, self%p_ptr%nstates), source=0)
        allocate(self%proj_active_state_count(self%p_ptr%nspace), source=0)
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

    !===========================================================
    ! Validate that the subspace data structures required by all
    ! build_neigh_mask_* routines are present and consistent.
    !===========================================================
    subroutine check_subspace_prereqs(self)
        class(eul_prob_tab_neigh), intent(in) :: self
        if( .not. allocated(self%b_ptr%subspace_inds) )then
            THROW_HARD('check_subspace_prereqs: subspace_inds not allocated; enable l_neigh and set nspace_sub')
        endif
        if( size(self%b_ptr%subspace_inds) /= self%p_ptr%nspace_sub )then
            THROW_HARD('check_subspace_prereqs: size(subspace_inds) must equal nspace_sub')
        endif
        if( .not. allocated(self%b_ptr%subspace_full2sub_map) )then
            THROW_HARD('check_subspace_prereqs: subspace_full2sub_map not allocated')
        endif
        if( size(self%b_ptr%subspace_full2sub_map) /= self%p_ptr%nspace )then
            THROW_HARD('check_subspace_prereqs: size(subspace_full2sub_map) must equal nspace')
        endif
    end subroutine check_subspace_prereqs

    !===========================================================
    ! Build per-particle neighborhood mask from previous
    ! assigned orientations only (no objective evaluation,
    ! no gradient-based shift search).
    !
    ! For each particle:
    !   1) get previous assigned orientation
    !   2) compute symmetry-aware distance to each coarse representative
    !   3) take npeak_use nearest coarse representatives
    !   4) pool their neighborhoods into mask(:,i)
    !===========================================================
    subroutine build_neigh_mask_from_prev_geom(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(eulspace_neigh_map) :: neigh_map
        type(ori) :: o_prev
        logical, allocatable :: mask(:,:)
        integer :: nspace_sub, npeak_use, i, iptcl, isub, iproj_full
        real :: dtmp, inplrotdist
        type(ori) :: o_sub, osym
        nspace_sub = self%p_ptr%nspace_sub
        npeak_use  = max(1, min(self%p_ptr%npeaks, nspace_sub))
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(mask(self%p_ptr%nspace, self%nptcls), source=.false.)
        !$omp parallel do default(shared) private(i,iptcl,o_prev,isub,iproj_full,dtmp,inplrotdist,o_sub,osym) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            block
                logical :: neigh_mask(self%p_ptr%nspace)
                integer :: peak_sub_idxs(npeak_use)
                integer :: prev_iproj, prev_sub
                logical :: prev_mask(self%p_ptr%nspace), have_prev_mask
                real    :: sub_dists(nspace_sub)
                iptcl = self%pinds(i)
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                mask(:,i) = .false.
                sub_dists = huge(1.)
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    if( iproj_full < 1 .or. iproj_full > self%p_ptr%nspace )then
                        THROW_HARD('build_neigh_mask_from_prev_geom: representative projection index out of range')
                    endif
                    call self%b_ptr%eulspace%get_ori(iproj_full, o_sub)
                    call self%b_ptr%pgrpsyms%sym_dists(o_prev, o_sub, osym, dtmp, inplrotdist)
                    sub_dists(isub) = dtmp
                    call o_sub%kill
                    call osym%kill
                enddo
                peak_sub_idxs = minnloc(sub_dists, npeak_use)
                call neigh_map%get_neighbors_mask_pooled(peak_sub_idxs, neigh_mask)
                mask(:,i) = neigh_mask
                ! explicit previous-neighborhood OR logic: derive previous projection from orientation,
                ! map to subspace, and OR in full neighborhood mask for true continuity
                prev_mask      = .false.
                have_prev_mask = .false.
                prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
                if( prev_iproj >= 1 .and. prev_iproj <= self%p_ptr%nspace )then
                    prev_sub = self%b_ptr%subspace_full2sub_map(prev_iproj)
                    if( prev_sub >= 1 .and. prev_sub <= nspace_sub )then
                        call neigh_map%get_neighbors_mask(prev_sub, prev_mask)
                        have_prev_mask = .true.
                    endif
                endif
                if( have_prev_mask ) mask(:,i) = mask(:,i) .or. prev_mask
                if( .not. any(mask(:,i) .and. (self%proj_active_state_count > 0)) )then
                    mask(:,i) = (self%proj_active_state_count > 0)
                endif
                call o_prev%kill
            end block
        enddo
        !$omp end parallel do
        call neigh_map%kill
        call self%build_sparse_from_mask(mask)
        if( allocated(mask) ) deallocate(mask)
    end subroutine build_neigh_mask_from_prev_geom

    !===========================================================
    ! Build per-particle neighborhood mask from subspace peaks.
    ! Evaluates each subspace projection for every particle,
    ! picks the top npeak_use peaks, and unions their neighbor
    ! masks (plus the previous-orientation neighbor mask) to
    ! produce the logical mask(:,i). Evaluates all active states.
    !===========================================================
    subroutine build_neigh_mask_from_subspace_peaks(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%build_neigh_mask_from_subspace_peaks_impl()
    end subroutine build_neigh_mask_from_subspace_peaks

    !===========================================================
    ! Single-state variant of build_neigh_mask_from_subspace_peaks.
    ! Evaluates subspace projections for every particle using only
    ! the supplied state, picks the top npeak_use peaks, and unions
    ! their neighbor masks (plus the previous-orientation neighbor
    ! mask) to produce the logical mask(:,i).
    !===========================================================
    subroutine build_neigh_mask_from_subspace_peaks_one_state(self, state)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: state
        ! validate state index
        if( state < 1 .or. state > self%p_ptr%nstates )then
            THROW_HARD('build_neigh_mask_from_subspace_peaks_one_state: state index out of range')
        endif
        if( .not. self%state_exists(state) )then
            THROW_HARD('build_neigh_mask_from_subspace_peaks_one_state: requested state does not exist')
        endif
        call self%build_neigh_mask_from_subspace_peaks_impl(state)
    end subroutine build_neigh_mask_from_subspace_peaks_one_state

    !===========================================================
    ! Build per-particle neighborhood mask using the tree-based
    ! coarse search logic from simple_strategy3D_ptree. Subspace
    ! representatives are scored first, best score per tree is kept,
    ! then the top unique trees are selected and only the local
    ! extrema encountered during the stochastic tree descent are
    ! inserted into the sparse candidate set.
    !===========================================================
    subroutine build_neigh_mask_from_ptree_srch(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        call self%build_neigh_mask_from_ptree_srch_impl()
    end subroutine build_neigh_mask_from_ptree_srch

    subroutine build_neigh_mask_from_ptree_srch_one_state(self, state)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer,                   intent(in)    :: state
        if( state < 1 .or. state > self%p_ptr%nstates )then
            THROW_HARD('build_neigh_mask_from_ptree_srch_one_state: state index out of range')
        endif
        if( .not. self%state_exists(state) )then
            THROW_HARD('build_neigh_mask_from_ptree_srch_one_state: requested state does not exist')
        endif
        call self%build_neigh_mask_from_ptree_srch_impl(state)
    end subroutine build_neigh_mask_from_ptree_srch_one_state

    subroutine build_neigh_mask_from_ptree_srch_impl(self, state_opt)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, optional,         intent(in)    :: state_opt
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori) :: o_prev
        logical, allocatable :: mask(:,:)
        integer, allocatable :: peak_trees(:,:)
        real,    allocatable :: inpl_corrs(:,:), tree_best_corrs(:,:), peak_tree_corrs(:,:)
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_shift(2)
        logical :: do_shift_first
        integer :: nrots, nspace_sub, npeak_use, ntrees
        integer :: i, isub, ithr, iptcl, istate, iproj_full, irot, iref_start, iref_prev, itree, state_i
        integer :: npeak_trees, ipeak, prev_iproj
        real    :: corr_best
        call self%check_subspace_prereqs()
        nspace_sub     = self%p_ptr%nspace_sub
        ntrees         = self%b_ptr%block_tree%get_n_trees()
        if( ntrees <= 0 )then
            THROW_HARD('build_neigh_mask_from_ptree_srch_impl: block_tree has no trees')
        endif
        npeak_use      = max(1, min(self%p_ptr%npeaks, ntrees))
        nrots          = self%b_ptr%pftc%get_nrots()
        do_shift_first = self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift
        allocate(mask(self%p_ptr%nspace, self%nptcls), source=.false.)
        allocate(inpl_corrs(nrots, nthr_glob), tree_best_corrs(ntrees, nthr_glob), &
                 peak_trees(npeak_use, nthr_glob), peak_tree_corrs(npeak_use, nthr_glob))
        if( do_shift_first )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
        endif
        !$omp parallel do default(shared) private(i,iptcl,ithr,cxy_shift,o_prev,prev_iproj,istate,iref_start,iref_prev,irot,cxy,isub,iproj_full,itree,corr_best,state_i,npeak_trees,ipeak) &
        !$omp proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            cxy_shift = [0.,0.]
            if( do_shift_first )then
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                prev_iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)
                istate     = o_prev%get_state()
                if( present(state_opt) ) istate = state_opt
                if( istate >= 1 .and. istate <= self%p_ptr%nstates )then
                    if( self%state_exists(istate) .and. prev_iproj > 0 )then
                        if( self%proj_exists(prev_iproj, istate) )then
                            iref_start = (istate - 1) * self%p_ptr%nspace
                            iref_prev  = iref_start + prev_iproj
                            irot       = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                            call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                            cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                            cxy_shift = 0.
                            if( irot /= 0 ) cxy_shift = cxy(2:3)
                        endif
                    endif
                endif
                call o_prev%kill
            endif
            tree_best_corrs(:,ithr) = -huge(1.0)
            do isub = 1, nspace_sub
                iproj_full = self%b_ptr%subspace_inds(isub)
                if( iproj_full < 1 .or. iproj_full > self%p_ptr%nspace )then
                    THROW_HARD('build_neigh_mask_from_ptree_srch_impl: representative projection index out of range')
                endif
                itree = self%b_ptr%subspace_full2sub_map(iproj_full)
                if( itree < 1 .or. itree > ntrees )then
                    THROW_HARD('build_neigh_mask_from_ptree_srch_impl: tree index out of range')
                endif
                if( present(state_opt) )then
                    if( .not. self%proj_exists(iproj_full, state_opt) ) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((state_opt - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                    corr_best = maxval(inpl_corrs(:,ithr))
                    tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                else
                    do state_i = 1, self%nstates
                        istate = self%active_state_indices(state_i)
                        if( .not. self%proj_exists(iproj_full, istate) ) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_corrs(:,ithr))
                        corr_best = maxval(inpl_corrs(:,ithr))
                        tree_best_corrs(itree,ithr) = max(tree_best_corrs(itree,ithr), corr_best)
                    enddo
                endif
            enddo
            mask(:,i) = .false.
            peak_trees(:,ithr)      = 0
            peak_tree_corrs(:,ithr) = -huge(1.0)
            call select_peak_trees(tree_best_corrs(:,ithr), peak_trees(:,ithr), peak_tree_corrs(:,ithr), npeak_trees)
            do ipeak = 1, npeak_trees
                call trace_tree_prob(self%b_ptr%block_tree, peak_trees(ipeak,ithr), mark_ptree_ref_eval)
            enddo
            if( .not. any(mask(:,i) .and. (self%proj_active_state_count > 0)) )then
                mask(:,i) = (self%proj_active_state_count > 0)
            endif
        enddo
        !$omp end parallel do
        if( do_shift_first )then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call self%build_sparse_from_mask(mask)
        if( allocated(mask)            ) deallocate(mask)
        if( allocated(inpl_corrs)      ) deallocate(inpl_corrs)
        if( allocated(tree_best_corrs) ) deallocate(tree_best_corrs)
        if( allocated(peak_trees)      ) deallocate(peak_trees)
        if( allocated(peak_tree_corrs) ) deallocate(peak_tree_corrs)

    contains

        subroutine mark_ptree_ref_eval(ref_idx, best_corr)
            integer, intent(in)  :: ref_idx
            real,    intent(out) :: best_corr
            integer :: istate, state_i
            logical :: have_valid_ref
            real    :: inpl_corrs(nrots)
            best_corr      = -huge(1.0)
            have_valid_ref = .false.
            if( ref_idx < 1 .or. ref_idx > self%p_ptr%nspace ) return
            if( present(state_opt) )then
                if( .not. self%proj_exists(ref_idx, state_opt) ) return
                call self%b_ptr%pftc%gen_objfun_vals((state_opt - 1) * self%p_ptr%nspace + ref_idx, iptcl, cxy_shift, inpl_corrs)
                best_corr      = maxval(inpl_corrs)
                have_valid_ref = .true.
            else
                do state_i = 1, self%nstates
                    istate = self%active_state_indices(state_i)
                    if( .not. self%proj_exists(ref_idx, istate) ) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + ref_idx, iptcl, cxy_shift, inpl_corrs)
                    best_corr      = max(best_corr, maxval(inpl_corrs))
                    have_valid_ref = .true.
                enddo
            endif
            if( have_valid_ref ) mask(ref_idx,i) = .true.
        end subroutine mark_ptree_ref_eval

    end subroutine build_neigh_mask_from_ptree_srch_impl

    !===========================================================
    ! Implementation: unified routine for both multi-state and
    ! single-state cases. If state_opt is present, evaluates only
    ! that state; otherwise evaluates all active states.
    !===========================================================
    subroutine build_neigh_mask_from_subspace_peaks_impl(self, state_opt)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, optional,         intent(in)    :: state_opt
        type(eulspace_neigh_map) :: neigh_map
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori) :: o_prev
        logical, allocatable :: mask(:,:)
        integer, allocatable :: peak_sub_idxs(:,:), all_sub_idxs(:)
        real,    allocatable :: inpl_dists(:,:), coarse_best_dist(:,:)
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_shift(2), huge_dist
        logical :: do_shift_first, prev_mask(self%p_ptr%nspace), have_prev_mask
        integer :: nrots, nspace_sub, npeak_use, i, isub, ithr, iptcl, istate, iproj_full, irot, iref_start, iref_prev
        integer :: state_i, nvalid_sub, prev_sub
        call self%check_subspace_prereqs()
        nspace_sub     = self%p_ptr%nspace_sub
        npeak_use      = max(1, min(self%p_ptr%npeaks, nspace_sub))
        nrots          = self%b_ptr%pftc%get_nrots()
        huge_dist      = huge(1.)
        do_shift_first = self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, self%p_ptr%nspace_sub)
        allocate(mask(self%p_ptr%nspace, self%nptcls), source=.false.)
        allocate(inpl_dists(nrots, nthr_glob), coarse_best_dist(nspace_sub, nthr_glob), peak_sub_idxs(npeak_use, nthr_glob), all_sub_idxs(nspace_sub))
        all_sub_idxs = (/(isub, isub=1, nspace_sub)/)
        if( do_shift_first )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
        endif
        !$omp parallel do default(shared) private(i,ithr,iptcl,o_prev,istate,isub,iproj_full,irot,iref_start,iref_prev,cxy,cxy_shift,prev_mask,have_prev_mask,state_i,nvalid_sub,prev_sub) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            ithr  = omp_get_thread_num() + 1
            call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
            cxy_shift      = [0.,0.]
            prev_mask      = .false.
            have_prev_mask = .false.
            iproj_full     = self%b_ptr%eulspace%find_closest_proj(o_prev)
            if( iproj_full >= 1 .and. iproj_full <= self%p_ptr%nspace )then
                prev_sub = self%b_ptr%subspace_full2sub_map(iproj_full)
                if( prev_sub >= 1 .and. prev_sub <= nspace_sub )then
                    call neigh_map%get_neighbors_mask(prev_sub, prev_mask)
                    have_prev_mask = .true.
                endif
            endif
            if( present(state_opt) ) then
                ! single-state path: skip state validation (done at wrapper), use provided state
                istate = state_opt
                if( do_shift_first )then
                    if( iproj_full >= 1 .and. iproj_full <= self%p_ptr%nspace )then
                        if( self%proj_exists(iproj_full, istate) )then
                            iref_start = (istate - 1) * self%p_ptr%nspace
                            iref_prev  = iref_start + iproj_full
                            irot       = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                            call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                            cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                            cxy_shift = 0.
                            if( irot /= 0 ) cxy_shift = cxy(2:3) ! minimize may reset irot to 0 on failure
                        endif
                    endif
                endif
            else
                ! multi-state path: get state from particle, validate
                istate = o_prev%get_state()
                if( do_shift_first )then
                    if( istate >= 1 .and. istate <= self%p_ptr%nstates )then
                        if( self%state_exists(istate) )then
                            if( iproj_full >= 1 .and. iproj_full <= self%p_ptr%nspace )then
                                if( self%proj_exists(iproj_full,istate) )then
                                    iref_start = (istate - 1) * self%p_ptr%nspace
                                    iref_prev  = iref_start + iproj_full
                                    irot       = self%b_ptr%pftc%get_roind(360. - o_prev%e3get())
                                    call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                                    cxy_shift = 0.
                                    if( irot /= 0 ) cxy_shift = cxy(2:3) ! minimize may reset irot to 0 on failure
                                endif
                            endif
                        endif
                    endif
                endif
            endif
            coarse_best_dist(:,ithr) = huge_dist
            if( present(state_opt) ) then
                ! single-state: evaluate only provided state
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    if( iproj_full < 1 .or. iproj_full > self%p_ptr%nspace )then
                        THROW_HARD('build_neigh_mask_from_subspace_peaks_impl: representative projection index out of range')
                    endif
                    if( .not. self%proj_exists(iproj_full, state_opt) ) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((state_opt - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                    inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), self%p_ptr%cc_objfun)
                    irot = minloc(inpl_dists(:,ithr), dim=1)
                    coarse_best_dist(isub,ithr) = inpl_dists(irot,ithr)
                enddo
            else
                ! multi-state: evaluate all active states
                do isub = 1, nspace_sub
                    iproj_full = self%b_ptr%subspace_inds(isub)
                    if( iproj_full < 1 .or. iproj_full > self%p_ptr%nspace )then
                        THROW_HARD('build_neigh_mask_from_subspace_peaks_impl: representative projection index out of range')
                    endif
                    do state_i = 1, self%nstates
                        istate = self%active_state_indices(state_i)
                        if( .not. self%proj_exists(iproj_full,istate) ) cycle
                        call self%b_ptr%pftc%gen_objfun_vals((istate - 1) * self%p_ptr%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                        inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), self%p_ptr%cc_objfun)
                        irot = minloc(inpl_dists(:,ithr), dim=1)
                        coarse_best_dist(isub,ithr) = min(coarse_best_dist(isub,ithr), inpl_dists(irot,ithr))
                    enddo
                enddo
            endif
            mask(:,i) = .false.
            nvalid_sub = count(coarse_best_dist(:,ithr) < huge_dist / 2.)
            if( nvalid_sub > 0 )then
                if( nvalid_sub <= npeak_use )then
                    peak_sub_idxs(1:nvalid_sub,ithr) = pack(all_sub_idxs, mask=coarse_best_dist(:,ithr) < huge_dist / 2.)
                    call neigh_map%get_neighbors_mask_pooled(peak_sub_idxs(1:nvalid_sub,ithr), mask(:,i))
                else
                    peak_sub_idxs(:,ithr) = minnloc(coarse_best_dist(:,ithr), npeak_use)
                    call neigh_map%get_neighbors_mask_pooled(peak_sub_idxs(:,ithr), mask(:,i))
                endif
            endif
            if( have_prev_mask ) mask(:,i) = mask(:,i) .or. prev_mask
            if( .not. any(mask(:,i) .and. (self%proj_active_state_count > 0)) )then
                if( have_prev_mask .and. any(prev_mask .and. (self%proj_active_state_count > 0)) )then
                    mask(:,i) = prev_mask
                else
                    mask(:,i) = (self%proj_active_state_count > 0)
                endif
            endif
            call o_prev%kill
        enddo
        !$omp end parallel do
        if( do_shift_first )then
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%kill
            enddo
        endif
        call neigh_map%kill
        call self%build_sparse_from_mask(mask)
        if( allocated(mask)             ) deallocate(mask)
        if( allocated(inpl_dists)       ) deallocate(inpl_dists)
        if( allocated(coarse_best_dist) ) deallocate(coarse_best_dist)
        if( allocated(peak_sub_idxs)    ) deallocate(peak_sub_idxs)
        if( allocated(all_sub_idxs)     ) deallocate(all_sub_idxs)
    end subroutine build_neigh_mask_from_subspace_peaks_impl

    !===========================================================
    ! Convert neigh_mask to sparse edges.
    !
    ! Interpretation:
    !   for particle i (iptcl=pinds(i)):
    !     allowed iproj are those with neigh_mask(iproj,i)=.true.
    !     allowed refs are (state_i, iproj), filtered by proj_exists/state_exists.
    !
    ! A particle must always have at least one valid neighbor.
    ! If any particle has zero valid neighbors after filtering, fail hard.
    !===========================================================
    subroutine build_sparse_from_mask(self, neigh_mask)
        class(eul_prob_tab_neigh), intent(inout) :: self
        logical, intent(in) :: neigh_mask(:,:)         ! (nspace,nptcls)
        integer :: i, iproj, state_list_idx, istate, ref_idx, edge_write_pos, edge_idx, iptcl
        real :: x
        if( size(neigh_mask,1) /= self%p_ptr%nspace .or. size(neigh_mask,2) /= self%nptcls )then
            THROW_HARD('build_sparse_from_mask: neighborhood mask shape mismatch')
        endif
        if( allocated(self%ptcl_off)       ) deallocate(self%ptcl_off)
        if( allocated(self%edge_ref_index) ) deallocate(self%edge_ref_index)
        if( allocated(self%edge_ptcl)      ) deallocate(self%edge_ptcl)
        if( allocated(self%edge_val)       ) deallocate(self%edge_val)
        allocate(self%ptcl_off(self%nptcls+1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            self%ptcl_off(i+1) = self%ptcl_off(i) + sum(self%proj_active_state_count, mask=neigh_mask(:,i))
            if( self%ptcl_off(i+1) == self%ptcl_off(i) )then
                THROW_HARD('Particle has zero valid neighbors in neighborhood mask')
            endif
        enddo
        self%nedges = self%ptcl_off(self%nptcls+1) - 1
        self%maxdeg_ptcl = maxval(self%ptcl_off(2:self%nptcls+1) - self%ptcl_off(1:self%nptcls))
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            edge_write_pos = self%ptcl_off(i)
            do iproj = 1, self%p_ptr%nspace
                if (.not. neigh_mask(iproj,i)) cycle
                do state_list_idx = 1, self%nstates
                    istate = self%active_state_indices(state_list_idx)
                    if (.not. self%proj_exists(iproj,istate)) cycle
                    ref_idx = self%ref_index_map(iproj,istate)
                    if (ref_idx <= 0) cycle
                    edge_idx = edge_write_pos
                    self%edge_ref_index(edge_idx)  = ref_idx
                    self%edge_ptcl(edge_idx)       = i
                    self%edge_val(edge_idx)%pind   = iptcl
                    self%edge_val(edge_idx)%istate = istate
                    self%edge_val(edge_idx)%iproj  = iproj
                    self%edge_val(edge_idx)%inpl   = 0
                    self%edge_val(edge_idx)%dist   = huge(x)
                    self%edge_val(edge_idx)%x      = 0.
                    self%edge_val(edge_idx)%y      = 0.
                    self%edge_val(edge_idx)%has_sh = .false.
                    edge_write_pos = edge_write_pos + 1
                enddo
            enddo
        enddo
    end subroutine build_sparse_from_mask

    !===========================================================
    ! Build reference-major adjacency (CSR grouped by reference index)
    !===========================================================
    !
    ! build_ref_adjacency does 3 simple chores:
    ! (1) Count how many connections each reference has. It goes through the whole pile once and says:
    !         Reference 1 has 5 connections
    !         Reference 2 has 2
    !         Reference 3 has 7
    !         and so on.
    ! (2) Reserve shelf space for each reference. Using those counts, it computes where each reference’s block
    !     starts and ends inside one big array. That is what ref_edge_offsets stores:
    !         start of ref 1’s block
    !         start of ref 2’s block
    !         start of ref 3’s block
    !         and so on.
    !         (Like page numbers in an index.)
    ! (3) Fill the blocks with the actual connection IDs
    !         It walks the pile again and drops each edge index into the correct reference block.
    !         That filled list is ref_edge_indices.
    !         After this, we can quickly access all edges connected to reference 2 by looking at ref_edge_indices
    !         (self%ref_edge_offsets(2):self%ref_edge_offsets(3)-1)
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

    !===========================================================
    ! Fill sparse table: evaluate only candidate edges
    !===========================================================
    subroutine fill_tab_sparse(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori) :: o_prev
        integer :: nrots, i, e, ithr, iptcl, ri, istate, iproj, irot, iref_full, iref_start
        integer :: n_candidates, candidate_i, n_refine, state_i, refine_rank
        integer :: n_refs_to_refine, n_samples
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rotmat(2,2), rot_xy(2)
        real,    allocatable :: inpl_angle_thres(:) ! (nstates_total)
        real,    allocatable :: dists_inpl(:,:), dists_inpl_sorted(:,:), candidate_dist(:,:) 
        integer, allocatable :: inds_sorted(:,:), candidate_edge(:,:) 
        call seed_rnd
        nrots = self%b_ptr%pftc%get_nrots()
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
        ! Phase A: shift-first evaluation (if enabled), then optional neighborhood shift refinement.
        if (self%p_ptr%l_sh_first .and. self%p_ptr%l_doshift) then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                               maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            enddo
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,istate,iproj,irot,iref_start,cxy,rotmat,rot_xy,n_candidates,e,ri,candidate_i,n_refine,refine_rank,cxy_prob,iref_full) &
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! shift-first on previous orientation's closest ref (if exists)
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
            ! no shift-first branch: evaluate only candidate edges (optionally refine shifts if enabled)
            if (self%p_ptr%l_prob_sh .and. self%p_ptr%l_doshift) then
                lims(:,1)      = -self%p_ptr%trs
                lims(:,2)      =  self%p_ptr%trs
                lims_init(:,1) = -SHC_INPL_TRSHWDTH
                lims_init(:,2) =  SHC_INPL_TRSHWDTH
                do ithr = 1, nthr_glob
                    call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier, &
                                                   maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
                enddo
                !$omp parallel do default(shared) private(i,iptcl,ithr,n_candidates,e,ri,istate,iproj,irot,candidate_i,n_refine,refine_rank,cxy,iref_full,iref_start) &
                !$omp proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    n_candidates = self%ptcl_off(i+1) - self%ptcl_off(i)
                    candidate_i = 0
                    do e = self%ptcl_off(i), self%ptcl_off(i+1)-1
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
                !$omp parallel do default(shared) private(i,iptcl,ithr,e,ri,istate,iproj,irot,iref_full,iref_start) &
                !$omp proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    do e = self%ptcl_off(i), self%ptcl_off(i+1)-1
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

    !===========================================================
    ! Sparse normalization: per particle normalize over its candidate edges,
    ! then min/max normalize across all edges.
    !===========================================================
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

    !===========================================================
    ! Sort each reference’s candidate edge list by ascending dist (in place)
    !===========================================================
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

    ! In-place heapsort of edge index list by key edge_val(edge_idx)%dist (ascending)
    subroutine sort_eidx_by_dist(self, edge_idx_list)
        class(eul_prob_tab_neigh), intent(in) :: self
        integer, intent(inout) :: edge_idx_list(:)
        integer :: n_edges, root_idx, end_idx, tmp_edge_idx
        n_edges = size(edge_idx_list)
        if (n_edges <= 1) return
        do root_idx = n_edges/2, 1, -1
            call self%sift_down(edge_idx_list, root_idx, n_edges)
        enddo
        do end_idx = n_edges, 2, -1
            tmp_edge_idx          = edge_idx_list(1)
            edge_idx_list(1)      = edge_idx_list(end_idx)
            edge_idx_list(end_idx)= tmp_edge_idx
            call self%sift_down(edge_idx_list, 1, end_idx-1)
        enddo
    end subroutine sort_eidx_by_dist

    ! heapsort helper: sift down root of edge_idx_list to restore heap property, assuming subtrees are already heaps.
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
            tmp_edge_idx               = edge_idx_list(scan_idx)
            edge_idx_list(scan_idx)    = edge_idx_list(swap_idx)
            edge_idx_list(swap_idx)    = tmp_edge_idx
            scan_idx = swap_idx
        enddo
    end subroutine sift_down

    !===========================================================
    ! Globally coupled sparse assignment.
    ! Mirrors dense ref_assign: choose a reference from its current best available
    ! particle distance, assign that particle, then update only affected references.
    !===========================================================
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
            ref_frontier_edge(ref_idx)   = self%ref_edge_offsets(ref_idx)
            best_dist_by_ref(ref_idx) = huge(best_dist_by_ref(ref_idx))
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
            sampled_active_pos = angle_sampling(active_best_dist(1:n_active_refs), sorted_dist_scratch(1:n_active_refs),&
                                &sorted_idx_scratch(1:n_active_refs), projs_athres, self%p_ptr%prob_athres)
            ref_idx  = active_refs(sampled_active_pos)
            edge_idx = self%ref_edge_indices(ref_frontier_edge(ref_idx))
            ptcl_local_idx = self%edge_ptcl(edge_idx)
            if (particle_available(ptcl_local_idx)) then
                particle_available(ptcl_local_idx) = .false.
                self%assgn_map(ptcl_local_idx) = self%edge_val(edge_idx)
                n_assigned = n_assigned + 1
            endif
            ! update only incident refs of this particle
            do neighbor_edge = self%ptcl_off(ptcl_local_idx), self%ptcl_off(ptcl_local_idx+1)-1
                call update_ref_frontier(self%edge_ref_index(neighbor_edge))
            enddo
        enddo
        deallocate(particle_available, ref_frontier_edge, best_dist_by_ref, active_pos_by_ref, active_refs, active_best_dist, sorted_dist_scratch, sorted_idx_scratch)
    
    contains
        
        ! Refresh the current frontier edge for one reference based on still-available particles.
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
                    active_refs(active_pos)  = swapped_ref_idx
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

    !===========================================================
    ! Sparse table I/O.
    ! Each partition writes its particle-major CSR rows. The driver
    ! reads all partition files, rebuilds the global sparse graph,
    ! and only then performs the globally coupled assignment.
    !===========================================================
    subroutine write_tab(self, binfname)
        class(eul_prob_tab_neigh), intent(in) :: self
        class(string), intent(in) :: binfname
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

    subroutine read_tabs_to_glob(self, fbody, nparts, numlen)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string), intent(in) :: fbody
        integer, intent(in) :: nparts, numlen
        type(string) :: binfname
        type(ptcl_ref), allocatable :: edge_val_loc(:)
        integer, allocatable :: degree_by_ptcl(:), fill_ptr(:), pind_to_local(:), pinds_loc(:), ptcl_off_loc(:)
        integer :: funit, addr, io_stat, file_header(3), header_bytes
        integer :: ipart, nrefs_loc, nptcls_loc, nedges_loc, local_i, global_i, local_e, global_e, ref_idx
        if( nparts < 1 )then
            THROW_HARD('read_tabs_to_glob: nparts must be >= 1')
        endif
        if( allocated(self%ptcl_off)         ) deallocate(self%ptcl_off)
        if( allocated(self%edge_ref_index)   ) deallocate(self%edge_ref_index)
        if( allocated(self%edge_ptcl)        ) deallocate(self%edge_ptcl)
        if( allocated(self%edge_val)         ) deallocate(self%edge_val)
        if( allocated(self%ref_edge_offsets) ) deallocate(self%ref_edge_offsets)
        if( allocated(self%ref_edge_indices) ) deallocate(self%ref_edge_indices)
        allocate(degree_by_ptcl(self%nptcls), source=0)
        allocate(pind_to_local(self%p_ptr%nptcls), source=0)
        do global_i = 1, self%nptcls
            if( self%pinds(global_i) < 1 .or. self%pinds(global_i) > self%p_ptr%nptcls )then
                THROW_HARD('read_tabs_to_glob: global particle index out of range')
            endif
            if( pind_to_local(self%pinds(global_i)) /= 0 )then
                THROW_HARD('read_tabs_to_glob: duplicate particle index in global sampled set')
            endif
            pind_to_local(self%pinds(global_i)) = global_i
        enddo
        header_bytes = size(file_header) * sizeof(file_header(1))
        do ipart = 1, nparts
            binfname = fbody//int2str_pad(ipart, numlen)//'.dat'
            if( .not. file_exists(binfname) )then
                THROW_HARD('file '//binfname%to_char()//' does not exists!')
            else
                call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            endif
            call fileiochk('simple_eul_prob_tab_neigh; read_tabs_to_glob; file: '//binfname%to_char(), io_stat)
            read(unit=funit, pos=1) file_header
            nrefs_loc  = file_header(1)
            nptcls_loc = file_header(2)
            nedges_loc = file_header(3)
            if( nrefs_loc /= self%nrefs )then
                THROW_HARD('read_tabs_to_glob: nrefs mismatch between sparse partition file and global object')
            endif
            allocate(pinds_loc(nptcls_loc), ptcl_off_loc(nptcls_loc + 1))
            addr = header_bytes + 1
            read(unit=funit, pos=addr) pinds_loc
            addr = addr + nptcls_loc * sizeof(pinds_loc(1))
            read(unit=funit, pos=addr) ptcl_off_loc
            call fclose(funit)
            if( ptcl_off_loc(1) /= 1 )then
                THROW_HARD('read_tabs_to_glob: invalid sparse partition offsets')
            endif
            if( ptcl_off_loc(nptcls_loc + 1) - 1 /= nedges_loc )then
                THROW_HARD('read_tabs_to_glob: sparse partition header/offset mismatch')
            endif
            do local_i = 1, nptcls_loc
                if( pinds_loc(local_i) < 1 .or. pinds_loc(local_i) > self%p_ptr%nptcls )then
                    THROW_HARD('read_tabs_to_glob: partition file contains particle index out of range')
                endif
                global_i = pind_to_local(pinds_loc(local_i))
                if( global_i <= 0 )then
                    THROW_HARD('read_tabs_to_glob: partition file contains particle outside the sampled set')
                endif
                if( degree_by_ptcl(global_i) /= 0 )then
                    THROW_HARD('read_tabs_to_glob: duplicate particle across sparse partition files')
                endif
                degree_by_ptcl(global_i) = ptcl_off_loc(local_i + 1) - ptcl_off_loc(local_i)
                if( degree_by_ptcl(global_i) <= 0 )then
                    THROW_HARD('read_tabs_to_glob: sparse partition row has no edges')
                endif
            enddo
            deallocate(pinds_loc, ptcl_off_loc)
        enddo
        if( any(degree_by_ptcl <= 0) )then
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
            if( .not. file_exists(binfname) )then
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
                    if( edge_val_loc(local_e)%iproj < 1 .or. edge_val_loc(local_e)%iproj > self%p_ptr%nspace )then
                        THROW_HARD('read_tabs_to_glob: sparse partition edge has iproj out of range')
                    endif
                    if( edge_val_loc(local_e)%istate < 1 .or. edge_val_loc(local_e)%istate > self%p_ptr%nstates )then
                        THROW_HARD('read_tabs_to_glob: sparse partition edge has istate out of range')
                    endif
                    ref_idx = self%ref_index_map(edge_val_loc(local_e)%iproj, edge_val_loc(local_e)%istate)
                    if( ref_idx <= 0 )then
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
        if( any(fill_ptr /= self%ptcl_off(2:self%nptcls + 1)) )then
            THROW_HARD('read_tabs_to_glob: failed to reconstruct the sparse global table consistently')
        endif
        call self%build_ref_adjacency()
        deallocate(degree_by_ptcl, fill_ptr, pind_to_local)
    end subroutine read_tabs_to_glob

    !===========================================================
    ! Assignment I/O (same format as dense write_assignment/read_assignment)
    !===========================================================

    subroutine write_assignment(self, binfname)
        class(eul_prob_tab_neigh), intent(in) :: self
        class(string), intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1)          self%nptcls
        write(unit=funit, pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    subroutine read_assignment(self, binfname)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string), intent(in) :: binfname
        type(ptcl_ref), allocatable :: assgn_glob(:)
        integer, allocatable :: pind_to_local(:)
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

    !===========================================================
    ! Destructor
    !===========================================================
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

end module simple_eul_prob_tab_neigh