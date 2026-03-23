!@descr: neighborhood extension of probabilistic 3D search table
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_eul_prob_tab,      only: eul_prob_tab, angle_sampling, calc_athres, eulprob_dist_switch
use simple_pftc_shsrch_grad,  only: pftc_shsrch_grad
use simple_ori,               only: ori
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    integer                 :: nrefs_sub = 0
    integer, allocatable    :: sub_ref_inds(:)
    integer, allocatable    :: sub_ref_offsets(:)
    integer, allocatable    :: sub_ref_list(:)
contains
    procedure :: new_neigh
    procedure :: fill_tab      => fill_tab_neigh
    procedure :: ref_normalize => ref_normalize_neigh
    procedure :: ref_assign    => ref_assign_neigh
    procedure :: kill          => kill_neigh
    procedure :: read_tabs_to_glob
end type eul_prob_tab_neigh

contains

    subroutine new_neigh( self, params, build, pinds, empty_okay, build_sparse_graph )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        logical, optional,         intent(in)    :: build_sparse_graph
        integer, allocatable :: sub_counts(:), sub_fill(:)
        integer :: isub, istate, iproj, ri, iref, nsubs, nfull2sub
        logical :: l_build
        l_build = .false.
        if( present(build_sparse_graph) ) l_build = build_sparse_graph
        if( trim(params%neigh_type) /= 'subspace_srch' )then
            THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; only neigh_type=subspace_srch is supported')
        endif
        call self%eul_prob_tab%new(params, build, pinds, empty_okay)
        nsubs = size(self%b_ptr%subspace_inds)
        nfull2sub = size(self%b_ptr%subspace_full2sub_map)
        self%nrefs_sub = nsubs * self%p_ptr%nstates
        if( allocated(self%sub_ref_inds) ) deallocate(self%sub_ref_inds)
        allocate(self%sub_ref_inds(self%nrefs_sub))
        ri = 0
        do istate = 1, self%p_ptr%nstates
            do isub = 1, nsubs
                ri    = ri + 1
                iproj = self%b_ptr%subspace_inds(isub)
                if( iproj < 1 .or. iproj > self%p_ptr%nspace )then
                    THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; subspace proj index out of bounds')
                endif
                self%sub_ref_inds(ri) = 0
                do iref = 1, self%nrefs
                    if( self%sinds(iref) /= istate ) cycle
                    if( self%jinds(iref) /= iproj  ) cycle
                    self%sub_ref_inds(ri) = iref
                    exit
                enddo
            enddo
        enddo
        ! Precompute all refs that belong to each subspace id for pooled refinement.
        if( allocated(self%sub_ref_offsets) ) deallocate(self%sub_ref_offsets)
        if( allocated(self%sub_ref_list) ) deallocate(self%sub_ref_list)
        allocate(sub_counts(nsubs), source=0)
        do iref = 1, self%nrefs
            iproj = self%jinds(iref)
            if( iproj < 1 .or. iproj > nfull2sub ) cycle
            isub = self%b_ptr%subspace_full2sub_map(iproj)
            if( isub < 1 .or. isub > nsubs ) cycle
            sub_counts(isub) = sub_counts(isub) + 1
        enddo
        allocate(self%sub_ref_offsets(nsubs+1), source=0)
        self%sub_ref_offsets(1) = 1
        do isub = 1, nsubs
            self%sub_ref_offsets(isub+1) = self%sub_ref_offsets(isub) + sub_counts(isub)
        enddo
        allocate(self%sub_ref_list(self%sub_ref_offsets(nsubs+1)-1))
        allocate(sub_fill(nsubs))
        sub_fill = self%sub_ref_offsets(1:nsubs)
        do iref = 1, self%nrefs
            iproj = self%jinds(iref)
            if( iproj < 1 .or. iproj > nfull2sub ) cycle
            isub = self%b_ptr%subspace_full2sub_map(iproj)
            if( isub < 1 .or. isub > nsubs ) cycle
            self%sub_ref_list(sub_fill(isub)) = iref
            sub_fill(isub) = sub_fill(isub) + 1
        enddo
        deallocate(sub_fill, sub_counts)
        if( l_build ) call self%fill_tab
    end subroutine new_neigh

    subroutine fill_tab_neigh( self )
        use simple_math, only: hpsort
        class(eul_prob_tab_neigh), intent(inout) :: self
        integer, allocatable :: chosen(:), inds_sorted(:), locn(:,:)
        integer, allocatable :: eval_refs(:,:)
        logical, allocatable :: sub_chosen(:,:)
        real,    allocatable :: dsub(:), inpl_athres(:), dists_inpl(:), dists_inpl_sorted(:), eval_dists(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(ori)              :: o_prev
        integer :: i, k, ri, iproj, isub, nchoose, irot, istate, ithr, projs_ns, j, iref_start, iptcl, nsubs, nfull2sub
        integer :: neval, seln, ri_eval
        real    :: rotmat(2,2), lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rot_xy(2)
        call seed_rnd
        nsubs = size(self%b_ptr%subspace_inds)
        nfull2sub = size(self%b_ptr%subspace_full2sub_map)
        nchoose = max(1, min(self%p_ptr%npeaks, self%nrefs_sub))
        allocate(chosen(nchoose), dsub(self%nrefs_sub))
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(self%b_ptr%pftc%get_nrots()),&
        &dists_inpl_sorted(self%b_ptr%pftc%get_nrots()),&
        &inds_sorted(self%b_ptr%pftc%get_nrots()))
        ! determine max number of neighbors to refine per state
        projs_ns = 0
        do istate = 1, self%p_ptr%nstates
            projs_ns = max(projs_ns, nchoose)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), source=0)
        allocate(eval_refs(self%nrefs,nthr_glob), source=0)
        allocate(eval_dists(self%nrefs,nthr_glob), source=huge(1.0))
        allocate(sub_chosen(nsubs,nthr_glob), source=.false.)
        if( self%p_ptr%l_doshift )then
            ! make shift search objects (shift-first + probabilistic refinement)
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,istate,irot,iproj,iref_start,cxy,ri,j,cxy_prob,rot_xy,rotmat,k,isub,neval,seln,ri_eval)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! reset all refs for this particle
                do ri = 1, self%nrefs
                    self%loc_tab(ri,i)%dist   = huge(1.0)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                enddo
                ! identify shifts using previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)
                istate     = o_prev%get_state()
                irot       = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                iproj      = self%b_ptr%eulspace%find_closest_proj(o_prev)
                if( self%state_exists(istate) .and. self%proj_exists(iproj,istate) )then
                    iref_start = (istate-1)*self%p_ptr%nspace
                    call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                    if( irot == 0 ) cxy(2:3) = 0.
                else
                    cxy(2:3) = 0.
                endif
                ! coarse pass on subspace refs
                do k = 1, self%nrefs_sub
                    ri = self%sub_ref_inds(k)
                    if( ri > 0 )then
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, cxy(2:3), dists_inpl)
                        dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                        irot       = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted, inpl_athres(istate), self%p_ptr%prob_athres)
                        dsub(k)    = dists_inpl(irot)
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                        rot_xy                    = matmul(cxy(2:3), rotmat)
                        self%loc_tab(ri,i)%dist   = dsub(k)
                        self%loc_tab(ri,i)%inpl   = irot
                        self%loc_tab(ri,i)%x      = rot_xy(1)
                        self%loc_tab(ri,i)%y      = rot_xy(2)
                        self%loc_tab(ri,i)%has_sh = .true.
                    else
                        dsub(k) = huge(1.0)
                    endif
                enddo
                chosen = minnloc(dsub, nchoose)
                sub_chosen(:,ithr) = .false.
                do k = 1, nchoose
                    ri    = self%sub_ref_inds(chosen(k))
                    if( ri < 1 ) cycle
                    iproj = self%jinds(ri)
                    if( iproj < 1 .or. iproj > nfull2sub ) cycle
                    isub  = self%b_ptr%subspace_full2sub_map(iproj)
                    if( isub < 1 .or. isub > nsubs ) cycle
                    sub_chosen(isub,ithr) = .true.
                enddo
                neval = 0
                ! evaluate pooled-neighborhood refs from selected subspace lists only
                do k = 1, nchoose
                    ri    = self%sub_ref_inds(chosen(k))
                    if( ri < 1 ) cycle
                    iproj = self%jinds(ri)
                    if( iproj < 1 .or. iproj > nfull2sub ) cycle
                    isub = self%b_ptr%subspace_full2sub_map(iproj)
                    if( isub < 1 .or. isub > nsubs ) cycle
                    if( .not. sub_chosen(isub,ithr) ) cycle
                    sub_chosen(isub,ithr) = .false.
                    do j = self%sub_ref_offsets(isub), self%sub_ref_offsets(isub+1)-1
                        ri     = self%sub_ref_list(j)
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, cxy(2:3), dists_inpl)
                        dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                        irot = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted, inpl_athres(istate), self%p_ptr%prob_athres)
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                        rot_xy                    = matmul(cxy(2:3), rotmat)
                        self%loc_tab(ri,i)%dist   = dists_inpl(irot)
                        self%loc_tab(ri,i)%inpl   = irot
                        self%loc_tab(ri,i)%x      = rot_xy(1)
                        self%loc_tab(ri,i)%y      = rot_xy(2)
                        self%loc_tab(ri,i)%has_sh = .true.
                        neval                     = neval + 1
                        eval_refs(neval,ithr)     = ri
                        eval_dists(neval,ithr)    = dists_inpl(irot)
                    enddo
                enddo
                ! refine shifts over evaluated pooled neighborhood refs only
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
            enddo
            !$omp end parallel do
        else
            ! no shift search - evaluate pooled neighborhoods with zero shift only
            !$omp parallel do default(shared) private(i,iptcl,ithr,ri,istate,iproj,irot,k,isub)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! reset all refs
                do ri = 1, self%nrefs
                    self%loc_tab(ri,i)%dist   = huge(1.0)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                enddo
                ! coarse pass on subspace refs
                do k = 1, self%nrefs_sub
                    ri = self%sub_ref_inds(k)
                    if( ri > 0 )then
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, [0.,0.], dists_inpl)
                        dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                        irot       = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted, inpl_athres(istate), self%p_ptr%prob_athres)
                        dsub(k)    = dists_inpl(irot)
                        self%loc_tab(ri,i)%dist   = dsub(k)
                        self%loc_tab(ri,i)%inpl   = irot
                    else
                        dsub(k) = huge(1.0)
                    endif
                enddo
                chosen = minnloc(dsub, nchoose)
                sub_chosen(:,ithr) = .false.
                do k = 1, nchoose
                    ri    = self%sub_ref_inds(chosen(k))
                    if( ri < 1 ) cycle
                    iproj = self%jinds(ri)
                    if( iproj < 1 .or. iproj > nfull2sub ) cycle
                    isub  = self%b_ptr%subspace_full2sub_map(iproj)
                    if( isub < 1 .or. isub > nsubs ) cycle
                    sub_chosen(isub,ithr) = .true.
                enddo
                ! evaluate pooled-neighborhood refs from selected subspace lists only
                do k = 1, nchoose
                    ri    = self%sub_ref_inds(chosen(k))
                    if( ri < 1 ) cycle
                    iproj = self%jinds(ri)
                    if( iproj < 1 .or. iproj > nfull2sub ) cycle
                    isub = self%b_ptr%subspace_full2sub_map(iproj)
                    if( isub < 1 .or. isub > nsubs ) cycle
                    if( .not. sub_chosen(isub,ithr) ) cycle
                    sub_chosen(isub,ithr) = .false.
                    do j = self%sub_ref_offsets(isub), self%sub_ref_offsets(isub+1)-1
                        ri     = self%sub_ref_list(j)
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, [0.,0.], dists_inpl)
                        dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                        irot       = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted, inpl_athres(istate), self%p_ptr%prob_athres)
                        self%loc_tab(ri,i)%dist   = dists_inpl(irot)
                        self%loc_tab(ri,i)%inpl   = irot
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
        deallocate(sub_chosen, eval_dists, eval_refs, locn, inds_sorted, dists_inpl_sorted,&
        &dists_inpl, inpl_athres, dsub, chosen)
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
        integer :: i, iref
        huge_val = huge(1.0)
        ! normalize only over evaluated refs (those with dist < huge_val)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%loc_tab(:,i)%dist, mask=(self%loc_tab(:,i)%dist < huge_val))
            if( sum_dist_all < TINY )then
                ! if no refs evaluated, randomize for stochasticity
                !$omp critical
                do iref = 1, self%nrefs
                    if( self%loc_tab(iref,i)%dist < huge_val )&
                    &self%loc_tab(iref,i)%dist = ran3()
                enddo
                !$omp end critical
            else
                ! divide only evaluated refs
                do iref = 1, self%nrefs
                    if( self%loc_tab(iref,i)%dist < huge_val )&
                    &self%loc_tab(iref,i)%dist = self%loc_tab(iref,i)%dist / sum_dist_all
                enddo
            endif
        enddo
        !$omp end parallel do
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

    subroutine ref_assign_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        type(ori) :: o_prev
        integer :: i, iref, assigned_iref, assigned_ptcl, istate,&
                    &iproj, irot, cand_ptcl, cand_ind, fallback_ref,&
                    &stab_inds(self%nptcls, self%nrefs), inds_sorted(self%nrefs), iref_dist_inds(self%nrefs)
        real    :: sorted_tab(self%nptcls, self%nrefs), projs_athres,iref_dist(self%nrefs), dists_sorted(self%nrefs)
        logical :: ptcl_avail(self%nptcls)
        real    :: huge_val
        huge_val = huge(1.0)
        ! normalization using neighborhood-specific logic (excludes unevaluated refs)
        call self%ref_normalize()
        ! sorting each column considering only evaluated refs (unevaluated stay at huge_val)
        sorted_tab = transpose(self%loc_tab(:,:)%dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i)
        do iref = 1, self%nrefs
            stab_inds(:,iref) = (/(i,i=1,self%nptcls)/)
            call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        projs_athres = 0.
        do istate = 1, self%nstates
            projs_athres = max(projs_athres, calc_athres(self%b_ptr%spproj_field, 'dist', self%p_ptr%prob_athres, state=istate))
        enddo
        ! first row is the current best reference distribution (now sparse)
        iref_dist_inds = 1
        iref_dist      = huge_val
        ptcl_avail     = .true.
        do iref = 1, self%nrefs
            do
                cand_ind = iref_dist_inds(iref)
                if( cand_ind > self%nptcls ) exit
                cand_ptcl = stab_inds(cand_ind, iref)
                if( ptcl_avail(cand_ptcl) .and. sorted_tab(cand_ind,iref) < huge_val .and. self%loc_tab(iref,cand_ptcl)%inpl > 0 )then
                    iref_dist(iref) = sorted_tab(cand_ind,iref)
                    exit
                endif
                if( cand_ind == self%nptcls )then
                    iref_dist(iref) = huge_val
                    exit
                endif
                iref_dist_inds(iref) = iref_dist_inds(iref) + 1
            enddo
        enddo
        do while( any(ptcl_avail) )
            if( minval(iref_dist) >= huge_val ) exit
            ! sampling the ref distribution to choose next iref to assign
            assigned_iref = angle_sampling(iref_dist, dists_sorted, inds_sorted, projs_athres, self%p_ptr%prob_athres)
            if( iref_dist(assigned_iref) >= huge_val ) exit
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
            ! update the iref_dist and iref_dist_inds, skipping unevaluated refs
            do iref = 1, self%nrefs
                do
                    cand_ind = iref_dist_inds(iref)
                    if( cand_ind > self%nptcls )then
                        iref_dist(iref) = huge_val
                        exit
                    endif
                    cand_ptcl = stab_inds(cand_ind, iref)
                    if( ptcl_avail(cand_ptcl) .and. sorted_tab(cand_ind,iref) < huge_val .and. self%loc_tab(iref,cand_ptcl)%inpl > 0 )then
                        iref_dist(iref) = sorted_tab(cand_ind,iref)
                        exit
                    endif
                    if( cand_ind == self%nptcls )then
                        iref_dist(iref) = huge_val
                        exit
                    endif
                    iref_dist_inds(iref) = iref_dist_inds(iref) + 1
                enddo
            enddo
        enddo
        ! robust fallback: if sparse graph left particles without valid candidates, seed from previous orientation
        do i = 1, self%nptcls
            if( .not. ptcl_avail(i) ) cycle
            call self%b_ptr%spproj_field%get_ori(self%pinds(i), o_prev)
            istate = o_prev%get_state()
            if( istate < 1 .or. istate > self%p_ptr%nstates ) istate = 1
            iproj  = self%b_ptr%eulspace%find_closest_proj(o_prev)
            iproj  = max(1, min(self%p_ptr%nspace, iproj))
            irot   = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
            if( irot < 1 .or. irot > self%b_ptr%pftc%get_nrots() ) irot = 1
            fallback_ref = (istate-1) * self%p_ptr%nspace + iproj
            self%assgn_map(i)        = self%loc_tab(fallback_ref,i)
            self%assgn_map(i)%inpl   = irot
            self%assgn_map(i)%has_sh = .false.
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
        enddo
        call o_prev%kill
    end subroutine ref_assign_neigh

    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        if( allocated(self%sub_ref_inds) ) deallocate(self%sub_ref_inds)
        if( allocated(self%sub_ref_offsets) ) deallocate(self%sub_ref_offsets)
        if( allocated(self%sub_ref_list) ) deallocate(self%sub_ref_list)
        call self%eul_prob_tab%kill
        self%nrefs_sub = 0
    end subroutine kill_neigh

end module simple_eul_prob_tab_neigh
