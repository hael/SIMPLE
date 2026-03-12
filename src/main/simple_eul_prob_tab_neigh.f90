!> Neighborhood sparse probabilistic ref table + globally coupled assignment
!> Separate module (does not modify simple_eul_prob_tab).
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,          only: builder
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_eul_prob_tab,     only: angle_sampling, calc_num2sample, calc_athres, eulprob_dist_switch
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
    procedure          :: new        => new_from_mask
    procedure          :: fill_tab   => fill_tab_sparse
    procedure          :: ref_assign => ref_assign_sparse
    procedure          :: write_assignment
    procedure          :: read_assignment
    procedure          :: kill
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
    ! Constructor: build sparse edges from per-particle neigh_mask(nspace,nptcls)
    ! neigh_mask(:,i) flags which projections (iproj) to evaluate for particle pinds(i).
    ! Those projections are interpreted in state = ptcl_state(i) (if provided),
    ! otherwise state is taken from previous orientation in spproj_field.
    !===========================================================
    subroutine new_from_mask(self, params, build, pinds, neigh_mask, ptcl_state, empty_okay)
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical,                   intent(in)    :: neigh_mask(:,:) ! (nspace, nptcls)
        integer, optional,         intent(in)    :: ptcl_state(:)   ! (nptcls)
        logical, optional,         intent(in)    :: empty_okay
        integer, parameter :: MIN_POP = 5
        logical :: l_empty
        integer :: nspace
        call self%kill
        self%p_ptr => params
        self%b_ptr => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        nspace = self%p_ptr%nspace
        if (size(neigh_mask,1) /= nspace .or. size(neigh_mask,2) /= self%nptcls) then
            THROW_HARD('neigh_mask must have shape (nspace, nptcls) matching params%nspace and size(pinds)')
        endif
        if( present(ptcl_state) )then
            if (size(ptcl_state) /= self%nptcls) then
                THROW_HARD('ptcl_state must have size nptcls')
            endif
        endif
        ! initialize assignment map
        call init_assgn_map(self)
        ! existence filtering (mirrors dense new_1)
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
        call self%build_ref_lists_and_map()
        ! Build sparse edge lists from neigh_mask + state per particle
        call self%build_sparse_from_mask(neigh_mask, ptcl_state)
        ! Build reference-major adjacency once
        call self%build_ref_adjacency()
    end subroutine new_from_mask

    subroutine init_assgn_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: ptcl_local_idx, iptcl
        real :: x
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
    end subroutine init_assgn_map

    subroutine build_ref_lists_and_map(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer :: state_list_idx, ref_idx, istate, iproj
        allocate(self%active_state_indices(self%nstates), self%ref_proj_indices(self%nrefs), self%ref_state_indices(self%nrefs))
        allocate(self%ref_index_map(self%p_ptr%nspace, self%p_ptr%nstates), source=0)
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
            enddo
        enddo
    end subroutine build_ref_lists_and_map

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
    subroutine build_sparse_from_mask(self, neigh_mask, ptcl_state)
        class(eul_prob_tab_neigh), intent(inout) :: self
        logical,           intent(in) :: neigh_mask(:,:)         ! (nspace,nptcls)
        integer, optional, intent(in) :: ptcl_state(:)           ! (nptcls)
        integer, allocatable :: ptcl_degree(:), state_for_ptcl(:)
        integer :: i, iproj, istate, ref_idx, edge_write_pos, edge_idx, iptcl
        type(ori) :: o_prev
        real :: x
        allocate(ptcl_degree(self%nptcls), state_for_ptcl(self%nptcls), source=0)
        ! Determine state per particle (either provided or from previous ori)
        if (present(ptcl_state)) then
            state_for_ptcl = ptcl_state
        else
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                state_for_ptcl(i) = o_prev%get_state()
            enddo
            call o_prev%kill
        endif
        ! Pass 1: count valid neighbors per particle and enforce at least 1 candidate
        do i = 1, self%nptcls
            istate = state_for_ptcl(i)
            ptcl_degree(i) = 0
            if (istate >= 1 .and. istate <= self%p_ptr%nstates .and. self%state_exists(istate)) then
                ! count only those projections that both mask and exist in cavgs
                ptcl_degree(i) = count(neigh_mask(:,i) .and. self%proj_exists(:,istate))
            endif
            if (ptcl_degree(i) == 0) then
                THROW_HARD('Particle has zero valid neighbors in neighborhood mask')
            endif
        enddo
        self%maxdeg_ptcl = maxval(ptcl_degree)
        allocate(self%ptcl_off(self%nptcls+1))
        self%ptcl_off(1) = 1
        do i = 1, self%nptcls
            self%ptcl_off(i+1) = self%ptcl_off(i) + ptcl_degree(i)
        enddo
        self%nedges = self%ptcl_off(self%nptcls+1) - 1
        allocate(self%edge_ref_index(self%nedges), self%edge_ptcl(self%nedges), self%edge_val(self%nedges))
        ! Pass 2: fill edges by scanning the mask
        do i = 1, self%nptcls
            iptcl  = self%pinds(i)
            istate = state_for_ptcl(i)
            edge_write_pos = self%ptcl_off(i)
            do iproj = 1, self%p_ptr%nspace
                if (.not. neigh_mask(iproj,i)) cycle
                if (.not. self%proj_exists(iproj,istate)) cycle
                ref_idx = self%ref_index_map(iproj,istate)
                if (ref_idx <= 0) cycle
                edge_idx = edge_write_pos
                self%edge_ref_index(edge_idx)         = ref_idx
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
        deallocate(ptcl_degree, state_for_ptcl)
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
                if (istate >= 1 .and. istate <= self%p_ptr%nstates) then
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
        call o_prev%kill
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
        integer :: funit, io_stat, nptcls_glob, headsz, local_idx, global_idx
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
        !$omp parallel do default(shared) private(local_idx,global_idx) proc_bind(close) schedule(static)
        do local_idx = 1, self%nptcls
            do global_idx = 1, nptcls_glob
                if (self%assgn_map(local_idx)%pind == assgn_glob(global_idx)%pind) then
                    self%assgn_map(local_idx) = assgn_glob(global_idx)
                    exit
                endif
            enddo
        enddo
        !$omp end parallel do
        deallocate(assgn_glob)
    end subroutine read_assignment

    !===========================================================
    ! Destructor
    !===========================================================
    subroutine kill(self)
        class(eul_prob_tab_neigh), intent(inout) :: self
        if (allocated(self%pinds))                deallocate(self%pinds)
        if (allocated(self%assgn_map))            deallocate(self%assgn_map)
        if (allocated(self%proj_exists))          deallocate(self%proj_exists)
        if (allocated(self%state_exists))         deallocate(self%state_exists)
        if (allocated(self%active_state_indices)) deallocate(self%active_state_indices)
        if (allocated(self%ref_proj_indices))     deallocate(self%ref_proj_indices)
        if (allocated(self%ref_state_indices))    deallocate(self%ref_state_indices)
        if (allocated(self%ref_index_map))        deallocate(self%ref_index_map)
        if (allocated(self%ptcl_off))             deallocate(self%ptcl_off)
        if (allocated(self%edge_ref_index))       deallocate(self%edge_ref_index)
        if (allocated(self%edge_ptcl))            deallocate(self%edge_ptcl)
        if (allocated(self%edge_val))             deallocate(self%edge_val)
        if (allocated(self%ref_edge_offsets))     deallocate(self%ref_edge_offsets)
        if (allocated(self%ref_edge_indices))     deallocate(self%ref_edge_indices)
        self%nptcls      = 0
        self%nstates     = 0
        self%nrefs       = 0
        self%nedges      = 0
        self%maxdeg_ptcl = 0
        self%b_ptr       => null()
        self%p_ptr       => null()
    end subroutine kill

end module simple_eul_prob_tab_neigh