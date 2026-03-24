!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is not just geometrically sparse, it is also state-sparse around the previous assignment only.
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab, angle_sampling, calc_num2sample, calc_athres, eulprob_dist_switch
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_ori,                only: ori
use simple_eulspace_neigh_map, only: eulspace_neigh_map
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    integer, allocatable :: eval_touched_refs(:,:)
    integer, allocatable :: eval_touched_counts(:)
    integer              :: eval_max_touched = 0
contains
    procedure :: new_neigh
    procedure :: fill_tab      => fill_tab_neigh
    procedure :: ref_normalize => ref_normalize_neigh
    procedure :: ref_assign    => ref_assign_neigh
    procedure :: kill          => kill_neigh
    procedure :: read_tabs_to_glob
end type eul_prob_tab_neigh

contains

    subroutine new_neigh( self, params, build, pinds, empty_okay )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer, allocatable :: sub_counts(:)
        integer :: nsubs, iproj, isub
        call self%kill
        call self%eul_prob_tab%new(params, build, pinds, empty_okay)
        nsubs = size(self%b_ptr%subspace_inds)
        if( nsubs < 1 )then
            THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; empty subspace indices')
        endif
        allocate(sub_counts(nsubs), source=0)
        do iproj = 1, size(self%b_ptr%subspace_full2sub_map)
            isub = self%b_ptr%subspace_full2sub_map(iproj)
            if( isub >= 1 .and. isub <= nsubs ) sub_counts(isub) = sub_counts(isub) + 1
        enddo
        self%eval_max_touched = max(1, maxval(sub_counts))
        deallocate(sub_counts)
        allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        allocate(self%eval_touched_counts(self%nptcls),                     source=0)
    end subroutine new_neigh

    subroutine fill_tab_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(eulspace_neigh_map) :: neigh_map
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori)              :: o_prev, osym
        type(ori), allocatable :: sub_oris(:)
        integer,   allocatable :: inds_sorted(:,:), locn(:,:), eval_refs(:,:)
        logical,   allocatable :: neigh_proj_mask(:,:)
        real,      allocatable :: inpl_athres(:), eval_dists(:,:), dists_inpl(:,:), dists_inpl_sorted(:,:)
        integer :: i, ri, iproj, isub, irot, istate, prev_state, ithr, projs_ns, j
        integer :: neval, seln, ri_eval, iref_start, iptcl, nsubs, best_sub,  n, si, k
        real    :: dtmp, inplrotdist, min_sub_dist
        real    :: rotmat(2,2), lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rot_xy(2)
        call seed_rnd
        nsubs = size(self%b_ptr%subspace_inds)
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, nsubs)
        if( self%eval_max_touched < 1 ) self%eval_max_touched = 1
        if( .not. allocated(self%eval_touched_refs) ) allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        if( .not. allocated(self%eval_touched_counts) ) allocate(self%eval_touched_counts(self%nptcls), source=0)
        !$omp parallel do default(shared) private(i,k,ri) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            do k = 1, self%eval_touched_counts(i)
                ri = self%eval_touched_refs(k,i)
                if( ri > 0 )then
                    self%loc_tab(ri,i)%dist   = huge(1.0)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                endif
            enddo
            self%eval_touched_counts(i) = 0
        enddo
        !$omp end parallel do
        allocate(sub_oris(nsubs))
        do isub = 1, nsubs
            iproj = self%b_ptr%subspace_inds(isub)
            call self%b_ptr%eulspace%get_ori(iproj, sub_oris(isub))
        enddo
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(self%b_ptr%pftc%get_nrots(),nthr_glob),&
        &dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob),&
        &inds_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob))
        ! determine max number of neighbors to refine per active state
        projs_ns = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, state=istate)
            projs_ns = max(projs_ns, n)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), eval_refs(self%nrefs,nthr_glob), source=0)
        allocate(neigh_proj_mask(self%p_ptr%nspace,nthr_glob),              source=.false.)
        allocate(eval_dists(self%nrefs,nthr_glob),                          source=huge(1.0))
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,osym,istate,irot,iproj,iref_start,cxy,ri,j,cxy_prob,rot_xy,rotmat,isub,neval,seln,ri_eval,dtmp,inplrotdist,best_sub,min_sub_dist)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! Geometry-based neighborhood around previous orientation (single closest subspace),
                ! with state fixed to the previous assignment.
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                prev_state   = o_prev%get_state()
                best_sub     = 1
                min_sub_dist = huge(1.0)
                do isub = 1, nsubs
                    call self%b_ptr%pgrpsyms%sym_dists(o_prev, sub_oris(isub), osym, dtmp, inplrotdist)
                    if( dtmp < min_sub_dist )then
                        min_sub_dist = dtmp
                        best_sub     = isub
                    endif
                    call osym%kill
                enddo
                call neigh_map%get_neighbors_mask(best_sub, neigh_proj_mask(:,ithr))
                istate = prev_state
                irot   = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                iproj  = self%b_ptr%eulspace%find_closest_proj(o_prev)
                if( self%state_exists(istate) .and. self%proj_exists(iproj,istate) )then
                    iref_start = (istate-1)*self%p_ptr%nspace
                    call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                    if( irot == 0 ) cxy(2:3) = 0.
                else
                    cxy(2:3) = 0.
                endif
                neval = 0
                ! Evaluate references in geometric neighborhood for the fixed state only.
                do ri = 1, self%nrefs
                    istate = self%sinds(ri)
                    if( istate /= prev_state ) cycle
                    iproj  = self%jinds(ri)
                    if( .not. neigh_proj_mask(iproj,ithr) ) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, cxy(2:3), dists_inpl(:,ithr))
                    dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                    irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres(istate), self%p_ptr%prob_athres)
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                    rot_xy                    = matmul(cxy(2:3), rotmat)
                    self%loc_tab(ri,i)%dist   = dists_inpl(irot,ithr)
                    self%loc_tab(ri,i)%inpl   = irot
                    self%loc_tab(ri,i)%x      = rot_xy(1)
                    self%loc_tab(ri,i)%y      = rot_xy(2)
                    self%loc_tab(ri,i)%has_sh = .true.
                    if( self%eval_touched_counts(i) >= self%eval_max_touched ) THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh; eval_touched overflow')
                    self%eval_touched_counts(i) = self%eval_touched_counts(i) + 1
                    self%eval_touched_refs(self%eval_touched_counts(i),i) = ri
                    neval                     = neval + 1
                    eval_refs(neval,ithr)     = ri
                    eval_dists(neval,ithr)    = dists_inpl(irot,ithr)
                enddo
                ! Dense-style refinement over best references in the evaluated neighborhood.
                locn(:,ithr) = 0
                if( neval > 0 )then
                    seln = min(projs_ns, neval)
                    locn(1:seln,ithr) = minnloc(eval_dists(1:neval,ithr), seln)
                else
                    seln = 0
                endif
                do j = 1, seln
                    ri_eval = locn(j,ithr)
                    if( ri_eval < 1 ) cycle
                    ri     = eval_refs(ri_eval,ithr)
                    if( ri < 1 ) cycle
                    istate = self%sinds(ri)
                    iproj  = self%jinds(ri)
                    call grad_shsrch_obj(ithr)%set_indices((istate-1)*self%p_ptr%nspace + iproj, iptcl)
                    irot     = self%loc_tab(ri,i)%inpl
                    cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=cxy(2:3))
                    if( irot > 0 )then
                        self%loc_tab(ri,i)%inpl   = irot
                        self%loc_tab(ri,i)%dist   = eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun)
                        self%loc_tab(ri,i)%x      = cxy_prob(2)
                        self%loc_tab(ri,i)%y      = cxy_prob(3)
                        self%loc_tab(ri,i)%has_sh = .true.
                    endif
                end do
                call o_prev%kill
            enddo
            !$omp end parallel do
        else
            ! no shift search - evaluate geometric neighborhoods with zero shift only
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,osym,ri,istate,iproj,irot,isub,dtmp,inplrotdist,best_sub,min_sub_dist)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                prev_state   = o_prev%get_state()
                best_sub     = 1
                min_sub_dist = huge(1.0)
                do isub = 1, nsubs
                    call self%b_ptr%pgrpsyms%sym_dists(o_prev, sub_oris(isub), osym, dtmp, inplrotdist)
                    if( dtmp < min_sub_dist )then
                        min_sub_dist = dtmp
                        best_sub     = isub
                    endif
                    call osym%kill
                enddo
                call neigh_map%get_neighbors_mask(best_sub, neigh_proj_mask(:,ithr))
                do ri = 1, self%nrefs
                    istate = self%sinds(ri)
                    if( istate /= prev_state ) cycle
                    iproj  = self%jinds(ri)
                    if( .not. neigh_proj_mask(iproj,ithr) ) cycle
                    call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, [0.,0.], dists_inpl(:,ithr))
                    dists_inpl(:,ithr)      = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                    irot                    = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres(istate), self%p_ptr%prob_athres)
                    self%loc_tab(ri,i)%dist = dists_inpl(irot,ithr)
                    self%loc_tab(ri,i)%inpl = irot
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                    if( self%eval_touched_counts(i) >= self%eval_max_touched ) THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh; eval_touched overflow')
                    self%eval_touched_counts(i) = self%eval_touched_counts(i) + 1
                    self%eval_touched_refs(self%eval_touched_counts(i),i) = ri
                enddo
                call o_prev%kill
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        do isub = 1, nsubs
            call sub_oris(isub)%kill
        enddo
        deallocate(sub_oris)
        call neigh_map%kill
        deallocate(eval_dists, eval_refs, neigh_proj_mask, locn, inds_sorted,&
        &dists_inpl_sorted, dists_inpl, inpl_athres)
    end subroutine fill_tab_neigh

    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(string) :: fname
        integer :: ipart
        do ipart = 1, nparts
            fname = fbody//int2str_pad(ipart,numlen)//'.dat'
            call self%read_tab_to_glob(fname)
        enddo
    end subroutine read_tabs_to_glob

    subroutine ref_normalize_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist, huge_val
        integer :: i, iref, neval
        huge_val = huge(1.0)
        ! normalize only over evaluated refs (dist < huge)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all,iref,neval)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%loc_tab(:,i)%dist, mask=(self%loc_tab(:,i)%dist < huge_val))
            if( sum_dist_all < TINY )then
                ! deterministic fallback: uniform weights over evaluated refs
                neval = count(self%loc_tab(:,i)%dist < huge_val)
                if( neval > 0 )then
                    do iref = 1, self%nrefs
                        if( self%loc_tab(iref,i)%dist < huge_val ) self%loc_tab(iref,i)%dist = 1. / real(neval)
                    enddo
                endif
            else
                ! divide only evaluated refs
                do iref = 1, self%nrefs
                    if( self%loc_tab(iref,i)%dist < huge_val )&
                    &self%loc_tab(iref,i)%dist = self%loc_tab(iref,i)%dist / sum_dist_all
                enddo
            endif
        enddo
        !$omp end parallel do
        if( .not. any(self%loc_tab(:,:)%dist < huge_val) ) return
        ! min/max normalization over evaluated refs only
        min_dist = minval(self%loc_tab(:,:)%dist, mask=(self%loc_tab(:,:)%dist < huge_val))
        max_dist = maxval(self%loc_tab(:,:)%dist, mask=(self%loc_tab(:,:)%dist < huge_val))
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab_neigh normalize')
        else
            where( self%loc_tab(:,:)%dist < huge_val )
                self%loc_tab(:,:)%dist = (self%loc_tab(:,:)%dist - min_dist) / (max_dist - min_dist)
            endwhere
        endif
    end subroutine ref_normalize_neigh

    ! sparse graph traversal on-the-fly for reference assignment
    ! with robust fallback to previous orientation if sparse graph leaves particles without valid candidates
    subroutine ref_assign_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
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
        type(ori) :: o_prev
        integer   :: i, iref, assigned_iref, assigned_ptcl, istate, iproj, irot, cand_ptcl, fallback_ref
        integer   :: k, idx, nactive, total, m, start, maxref, nleft, assigned_idx, nsel, pos, last_ref
        real      :: projs_athres
        real      :: huge_val
        huge_val = huge(1.0)
        ! normalization using neighborhood-specific logic (excludes unevaluated refs)
        call self%ref_normalize()
        allocate(graph%ref_counts(self%nrefs), graph%ref_offsets(self%nrefs+1), graph%ref_fill(self%nrefs), graph%ref_pos(self%nrefs),&
        &frontier%iref_dist(self%nrefs), frontier%dists_sorted(self%nrefs), frontier%inds_sorted(self%nrefs), frontier%ptcl_avail(self%nptcls),&
        &graph%ptcl_counts(self%nptcls), graph%ptcl_offsets(self%nptcls+1), graph%ptcl_fill(self%nptcls))
        graph%ref_counts = 0
        graph%ptcl_counts = 0
        do i = 1, self%nptcls
            do k = 1, self%eval_touched_counts(i)
                iref = self%eval_touched_refs(k,i)
                if( iref < 1 .or. iref > self%nrefs       ) cycle
                if( self%loc_tab(iref,i)%dist >= huge_val ) cycle
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
                if( self%loc_tab(iref,i)%dist >= huge_val ) cycle
                if( self%loc_tab(iref,i)%inpl <= 0        ) cycle
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
        projs_athres = 0.
        do istate = 1, self%nstates
            projs_athres = max(projs_athres, calc_athres(self%b_ptr%spproj_field, 'dist', self%p_ptr%prob_athres, state=istate))
        enddo
        ! current best per-reference candidate among available particles
        graph%ref_pos        = 1
        frontier%iref_dist   = huge_val
        frontier%ptcl_avail  = .true.
        nleft          = self%nptcls
        nsel           = 0
        do idx = 1, nactive
            iref  = graph%active_refs(idx)
            call advance_ref_head(iref)
            call sync_frontier_ref(iref)
        enddo
        do while( nleft > 0 )
            if( nsel == 0 ) exit
            ! sample only over currently valid active references
            assigned_idx = angle_sampling(frontier%sel_dists(1:nsel), frontier%sel_dists_sorted(1:nsel), frontier%inds_sorted(1:nsel), projs_athres, self%p_ptr%prob_athres)
            assigned_iref = frontier%sel_refs(assigned_idx)
            assigned_ptcl = graph%ref_list(graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1)
            frontier%ptcl_avail(assigned_ptcl) = .false.
            nleft                       = nleft - 1
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
            ! update only references that include the newly assigned particle
            do idx = graph%ptcl_offsets(assigned_ptcl), graph%ptcl_offsets(assigned_ptcl+1)-1
                iref  = graph%ptcl_refs(idx)
                m     = graph%ref_counts(iref)
                if( graph%ref_pos(iref) > m ) cycle
                start = graph%ref_offsets(iref)
                if( graph%ref_list(start + graph%ref_pos(iref) - 1) /= assigned_ptcl ) cycle
                call advance_ref_head(iref)
                call sync_frontier_ref(iref)
            enddo
        enddo
        ! robust fallback: if sparse graph left particles without valid candidates, seed from previous orientation
        do i = 1, self%nptcls
            if( .not. frontier%ptcl_avail(i) ) cycle
            call self%b_ptr%spproj_field%get_ori(self%pinds(i), o_prev)
            istate = o_prev%get_state()
            if( istate < 1 .or. istate > self%p_ptr%nstates ) istate = 1
            iproj  = self%b_ptr%eulspace%find_closest_proj(o_prev)
            iproj  = max(1, min(self%p_ptr%nspace, iproj))
            irot   = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
            if( irot < 1 .or. irot > self%b_ptr%pftc%get_nrots() ) irot = 1
            fallback_ref             = (istate-1) * self%p_ptr%nspace + iproj
            self%assgn_map(i)        = self%loc_tab(fallback_ref,i)
            self%assgn_map(i)%inpl   = irot
            self%assgn_map(i)%has_sh = .false.
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
        enddo
        call o_prev%kill

    contains

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
                    frontier%iref_dist(iref) = graph%ref_dists(sloc + graph%ref_pos(iref) - 1)
                    return
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

    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        if( allocated(self%eval_touched_refs)   ) deallocate(self%eval_touched_refs)
        if( allocated(self%eval_touched_counts) ) deallocate(self%eval_touched_counts)
        self%eval_max_touched = 0
        call self%eul_prob_tab%kill
    end subroutine kill_neigh

end module simple_eul_prob_tab_neigh
