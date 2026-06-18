!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_pftc_srch_api
use simple_builder,          only: builder
use simple_eul_prob_tab_utils, only: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample,&
    &eulprob_dist_switch, materialize_seed_shift, read_seed_shift_table, write_seed_shift_table
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_type_defs,        only: OBJFUN_EUCLID
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"

integer, parameter :: PROB_TAB_IO_CHUNK = 1024

type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(ptcl_ref), allocatable :: loc_tab(:,:)    !< 2D search table (nspace*nstates, nptcls)
    type(ptcl_ref), allocatable :: state_tab(:,:)  !< 2D search table (nstates,        nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
    real,           allocatable :: seed_shifts(:,:) !< per-particle seeded shift       (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)   !< per-particle seeded shift flag  (nptcls)
    integer                     :: seed_nrots = 0   !< rotation grid used for deferred seeded shifts
    integer,        allocatable :: pinds(:)        !< particle indices for processing
    integer,        allocatable :: ssinds(:)       !< non-empty state indices
    integer,        allocatable :: sinds(:)        !< non-empty state indices of each ref
    integer,        allocatable :: jinds(:)        !< non-empty proj  indices of each ref
    logical,        allocatable :: proj_exists(:,:)
    logical,        allocatable :: state_exists(:)
    integer                     :: nptcls          !< size of pinds array
    integer                     :: nstates         !< states number
    integer                     :: nrefs           !< reference number
    contains
    ! CONSTRUCTOR
    procedure :: new
    procedure :: new_assignment
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: fill_tab
    procedure :: fill_tab_state_only
    procedure :: write_tab
    procedure :: write_state_tab
    procedure :: read_assignment
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: read_state_tab
    procedure :: read_tab_to_glob
    procedure :: ref_assign
    procedure :: write_assignment
    procedure :: state_assign
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: ref_normalize
    procedure, private :: state_normalize
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new( self, params, build, pinds, empty_okay, state_only )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        logical, optional,         intent(in)    :: state_only
        integer, parameter :: MIN_POP = 5   ! ignoring cavgs with less than 5 particles
        integer :: i, iproj, iptcl, istate, si, ri
        real    :: x
        logical :: l_empty, l_state_only
        call self%kill
        l_empty = (trim(params%empty3Dcavgs) .eq. 'yes')
        if( present(empty_okay) ) l_empty = empty_okay
        l_state_only = .false.
        if( present(state_only) ) l_state_only = state_only
        self%p_ptr => params
        self%b_ptr  => build
        self%nptcls       = size(pinds)
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        if( l_empty )then
            allocate(self%proj_exists(self%p_ptr%nspace,self%p_ptr%nstates), source=.false.)
            do istate = 1,self%p_ptr%nstates
                if( self%state_exists(istate) ) self%proj_exists(:,istate) = .true.
            enddo
        else
            self%proj_exists = self%b_ptr%spproj_field%projs_exist(self%p_ptr%nstates,self%p_ptr%nspace, thres=MIN_POP)
            do istate = 1,self%p_ptr%nstates
                if( .not. self%state_exists(istate) ) self%proj_exists(:,istate) = .false.
            enddo
        endif
        self%nrefs = count(self%proj_exists .eqv. .true.)
        allocate(self%ssinds(self%nstates),self%jinds(self%nrefs),self%sinds(self%nrefs))
        si = 0
        ri = 0
        do istate = 1, self%p_ptr%nstates
            if( .not. self%state_exists(istate) )cycle
            si              = si + 1
            self%ssinds(si) = istate
            do iproj = 1,self%p_ptr%nspace
                if( .not. self%proj_exists(iproj,istate) )cycle
                ri             = ri + 1
                self%jinds(ri) = iproj
                self%sinds(ri) = istate
            enddo
        enddo
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls), self%state_tab(self%nstates,self%nptcls))
        if( .not. l_state_only ) allocate(self%loc_tab(self%nrefs,self%nptcls))
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
        !$omp parallel do default(shared) private(i,iptcl,si,ri) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            do si = 1,self%nstates
                self%state_tab(si,i)%pind   = iptcl
                self%state_tab(si,i)%istate = self%ssinds(si)
                self%state_tab(si,i)%iproj  = 0
                self%state_tab(si,i)%inpl   = 0
                self%state_tab(si,i)%dist   = huge(x)
                self%state_tab(si,i)%x      = 0.
                self%state_tab(si,i)%y      = 0.
                self%state_tab(si,i)%has_sh = .false.
            enddo
            if( allocated(self%loc_tab) )then
                do ri = 1, self%nrefs
                    self%loc_tab(ri,i)%pind   = iptcl
                    self%loc_tab(ri,i)%istate = self%sinds(ri)
                    self%loc_tab(ri,i)%iproj  = self%jinds(ri)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%dist   = huge(x)
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine new

    subroutine new_assignment( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: i, iptcl
        real    :: x
        call self%kill
        self%p_ptr  => params
        self%b_ptr  => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        do i = 1,self%nptcls
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
        enddo
    end subroutine new_assignment

    ! partition-wise table filling, used only in shared-memory commander 'exec_prob_tab'
    subroutine fill_tab( self )
        class(eul_prob_tab), intent(inout) :: self
        integer, allocatable   :: locn(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
        type(ori)              :: o_prev
        integer :: i, si, ri, j, iproj, iptcl, n, projs_ns, ithr, irot, inds_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob),&
                  &istate, iref_start
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), inpl_athres(self%p_ptr%nstates)
        real    :: dists_inpl(self%b_ptr%pftc%get_nrots(),nthr_glob), dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), dists_refs(self%nrefs,nthr_glob)
        logical :: l_prob_objfun
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        l_prob_objfun   = (self%p_ptr%cc_objfun == OBJFUN_EUCLID)
        call seed_rnd
        projs_ns = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, state=istate)
            projs_ns            = max(projs_ns, n)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), source=0)
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
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,istate,irot,iproj,iref_start,cxy,ri,j,cxy_prob)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! (1) identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)        ! previous ori
                istate     = o_prev%get_state()
                irot       = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj      = self%b_ptr%eulspace%find_closest_proj(o_prev) ! previous projection direction
                cxy = 0.
                if( istate >= 1 .and. istate <= self%p_ptr%nstates .and. iproj >= 1 .and. iproj <= self%p_ptr%nspace )then
                    if( self%state_exists(istate) .and. self%proj_exists(iproj,istate) )then
                        iref_start = (istate-1)*self%p_ptr%nspace
                        ! BFGS over shifts
                        call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        if( irot == 0 ) cxy(2:3) = 0.
                    endif
                endif
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
                ! (2) search projection directions using those shifts for all references
                do ri = 1, self%nrefs
                    istate = self%sinds(ri)
                    iproj  = self%jinds(ri)
                    if( l_prob_objfun )then
                        call self%b_ptr%pftc%gen_prob_objfun_val((istate-1)*self%p_ptr%nspace + iproj, iptcl, cxy(2:3),&
                            &inpl_athres(istate), self%p_ptr%prob_athres, self%loc_tab(ri,i)%dist, irot,&
                            &dists_inpl_sorted(:,ithr), inds_sorted(:,ithr))
                    else
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, cxy(2:3), dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                        irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr),&
                            &inpl_athres(istate), self%p_ptr%prob_athres)
                        self%loc_tab(ri,i)%dist = dists_inpl(irot,ithr)
                    endif
                    dists_refs(ri,ithr)     = self%loc_tab(ri,i)%dist
                    self%loc_tab(ri,i)%inpl = irot
                enddo
                ! (3) see if we can refine the shifts by re-searching them for individual references in the 
                !     identified probabilistic neighborhood
                locn(:,ithr) = minnloc(dists_refs(:,ithr), projs_ns)
                do j = 1,projs_ns
                    ri     = locn(j,ithr)
                    istate = self%sinds(ri)
                    iproj  = self%jinds(ri)
                    ! BFGS over shifts
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
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,ri,istate,iproj,irot) proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                do ri = 1, self%nrefs
                    istate = self%sinds(ri)
                    iproj  = self%jinds(ri)
                    if( l_prob_objfun )then
                        call self%b_ptr%pftc%gen_prob_objfun_val((istate-1)*self%p_ptr%nspace + iproj, iptcl, [0.,0.],&
                            &inpl_athres(istate), self%p_ptr%prob_athres, self%loc_tab(ri,i)%dist, irot,&
                            &dists_inpl_sorted(:,ithr), inds_sorted(:,ithr))
                    else
                        call self%b_ptr%pftc%gen_objfun_vals((istate-1)*self%p_ptr%nspace + iproj, iptcl, [0.,0.], dists_inpl(:,ithr))
                        dists_inpl(:,ithr)      = eulprob_dist_switch(dists_inpl(:,ithr), self%p_ptr%cc_objfun)
                        irot                    = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr),&
                            &inpl_athres(istate), self%p_ptr%prob_athres)
                        self%loc_tab(ri,i)%dist = dists_inpl(irot,ithr)
                    endif
                    self%loc_tab(ri,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab

    subroutine fill_tab_state_only( self )
        class(eul_prob_tab), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)  !< origin shift search object, L-BFGS with gradient
        type(ori)               :: o_prev
        integer :: i, iproj, iptcl, ithr, irot, istate, iref_start, is
        real    :: lims(2,2), lims_init(2,2), cxy(3)
        call seed_rnd
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,iproj,is,istate,irot,iref_start,cxy)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                iproj = self%b_ptr%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate     = self%ssinds(is)
                    irot       = self%b_ptr%pftc%get_roind(360.-o_prev%e3get()) ! in-plane angle index
                    iref_start = (istate-1)*self%p_ptr%nspace
                    ! BFGS over shifts
                    call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        irot     = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                        cxy(1)   = real(self%b_ptr%pftc%gen_corr_for_rot_8(iref_start+iproj, iptcl, irot))
                        cxy(2:3) = 0.
                    endif
                    self%state_tab(is,i)%dist   = eulprob_dist_switch(cxy(1), self%p_ptr%cc_objfun)
                    self%state_tab(is,i)%iproj  = iproj
                    self%state_tab(is,i)%inpl   = irot
                    self%state_tab(is,i)%x      = cxy(2)
                    self%state_tab(is,i)%y      = cxy(3)
                    self%state_tab(is,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        else
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,o_prev,irot,iproj,is,istate,iref_start)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ! identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                irot  = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj = self%b_ptr%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate      = self%ssinds(is)
                    iref_start  = (istate-1)*self%p_ptr%nspace
                    self%state_tab(is,i)%dist   = eulprob_dist_switch(&
                        &real(self%b_ptr%pftc%gen_corr_for_rot_8(iref_start+iproj, iptcl, irot)), self%p_ptr%cc_objfun)
                    self%state_tab(is,i)%iproj  = iproj
                    self%state_tab(is,i)%inpl   = irot
                    self%state_tab(is,i)%x      = 0.
                    self%state_tab(is,i)%y      = 0.
                    self%state_tab(is,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab_state_only

    ! Legacy in-place normalization retained for subclasses that override this
    ! binding. Plain prob assignment uses compact score vectors instead,
    ! preserving raw distances for ASSIGNMENT.dat and downstream score reporting.
    subroutine ref_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, iref
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%loc_tab(:,i)%dist)
            if( sum_dist_all < TINY )then
                self%loc_tab(:,i)%dist = 0.
            else
                self%loc_tab(:,i)%dist = self%loc_tab(:,i)%dist / sum_dist_all
            endif
        enddo
        !$omp end parallel do
        !$omp parallel workshare proc_bind(close)
        min_dist = minval(self%loc_tab(:,:)%dist)
        max_dist = maxval(self%loc_tab(:,:)%dist)
        !$omp end parallel workshare
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab')
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    self%loc_tab(iref,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            self%loc_tab(:,:)%dist = (self%loc_tab(:,:)%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine ref_normalize

    ! ptcl -> (proj, state) assignment using the global normalized dist value table
    subroutine ref_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer, allocatable :: stab_inds(:,:), inds_sorted(:), iref_dist_inds(:)
        integer, allocatable :: greedy_state(:), active_refs(:), active_inds(:), score_mode(:)
        real,    allocatable :: score_work(:,:), score_min(:), score_spread(:)
        real,    allocatable :: iref_dist(:), dists_sorted(:)
        real,    allocatable :: state_projs_athres(:), active_dists(:), active_dists_sorted(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: i, iref, assigned_iref, assigned_ptcl, istate, si, active_idx, nactive, ithr
        integer :: alloc_stat
        real    :: projs_athres
        character(len=256) :: alloc_msg
        allocate(stab_inds(self%nptcls, self%nrefs), score_work(self%nptcls,nthr_glob),&
            &score_min(self%nptcls), score_spread(self%nptcls), score_mode(self%nptcls),&
            &inds_sorted(self%nrefs), iref_dist_inds(self%nrefs), greedy_state(self%nptcls),&
            &active_refs(self%nrefs), active_inds(self%nrefs), iref_dist(self%nrefs),&
            &dists_sorted(self%nrefs), state_projs_athres(self%p_ptr%nstates), active_dists(self%nrefs),&
            &active_dists_sorted(self%nrefs), ptcl_avail(self%nptcls), stat=alloc_stat, errmsg=alloc_msg)
        if( alloc_stat /= 0 )then
            write(logfhandle,*) 'eul_prob_tab%ref_assign allocation failed for nptcls/nrefs: ', self%nptcls, self%nrefs
            write(logfhandle,*) trim(alloc_msg)
            THROW_HARD('failed allocating probability assignment work arrays')
        endif
        call prepare_ref_score_vectors
        ! Sort each reference using one scratch score column per worker thread.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i,ithr)
        do iref = 1, self%nrefs
            ithr = omp_get_thread_num() + 1
            do i = 1,self%nptcls
                stab_inds(i,iref) = i
                score_work(i,ithr) = ref_score(iref, i)
            enddo
            call hpsort(score_work(:,ithr), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        projs_athres = 0.
        state_projs_athres = 0.
        do si = 1, self%nstates
            istate = self%ssinds(si)
            state_projs_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist',&
                &self%p_ptr%prob_athres, state=istate)
            projs_athres = max(projs_athres, state_projs_athres(istate))
        enddo
        greedy_state = 0
        if( self%nstates > 1 )then
            call assign_greedy_state_labels()
            do si = 1,self%nstates
                call assign_refs_for_state(self%ssinds(si))
            enddo
        else
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                assigned_iref = angle_sampling(iref_dist, dists_sorted, inds_sorted, projs_athres, self%p_ptr%prob_athres)
                call assign_current_ref()
                do iref = 1,self%nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
        endif
        if( allocated(stab_inds) ) deallocate(stab_inds, inds_sorted, iref_dist_inds, greedy_state, active_refs, active_inds,&
            &score_work, score_min, score_spread, score_mode, iref_dist, dists_sorted, state_projs_athres, active_dists,&
            &active_dists_sorted, ptcl_avail)

    contains

        subroutine prepare_ref_score_vectors()
            real    :: min_dist, max_dist, dist_val, spread, invalid_dist
            integer :: i_loc, iref_loc, nvalid, nflat, nbad
            invalid_dist = 0.1 * huge(invalid_dist)
            nflat = 0
            nbad  = 0
            score_min    = 0.
            score_spread = 0.
            score_mode   = 0
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i_loc,iref_loc,min_dist,max_dist,dist_val,spread,nvalid) reduction(+:nflat,nbad)
            do i_loc = 1, self%nptcls
                min_dist = huge(min_dist)
                max_dist = -huge(max_dist)
                nvalid   = 0
                do iref_loc = 1,self%nrefs
                    dist_val = self%loc_tab(iref_loc,i_loc)%dist
                    if( ieee_is_finite(dist_val) .and. dist_val < invalid_dist )then
                        min_dist = min(min_dist, dist_val)
                        max_dist = max(max_dist, dist_val)
                        nvalid   = nvalid + 1
                    endif
                enddo
                if( nvalid == 0 )then
                    nbad = nbad + 1
                    cycle
                endif
                spread = max_dist - min_dist
                score_min(i_loc)    = min_dist
                score_spread(i_loc) = spread
                if( spread <= 0. )then
                    score_mode(i_loc) = 1
                    nflat = nflat + 1
                else
                    score_mode(i_loc) = 2
                endif
            enddo
            !$omp end parallel do
            if( nbad == self%nptcls ) THROW_HARD('all probability-table reference distances are invalid')
            if( nflat == self%nptcls ) THROW_HARD('all probability-table reference distances are flat')
            if( nbad > 0 ) write(logfhandle,*) 'WARNING: particles with invalid probability-table distances: ', nbad
            if( nflat > 0 ) write(logfhandle,*) 'WARNING: particles with flat probability-table reference distances: ', nflat
        end subroutine prepare_ref_score_vectors

        real function ref_score( iref_loc, iptcl_loc ) result(score)
            integer, intent(in) :: iref_loc, iptcl_loc
            real :: dist_val, invalid_dist
            invalid_dist = 0.1 * huge(invalid_dist)
            dist_val = self%loc_tab(iref_loc,iptcl_loc)%dist
            score = 1.
            if( .not.(ieee_is_finite(dist_val) .and. dist_val < invalid_dist) ) return
            select case(score_mode(iptcl_loc))
                case(1)
                    score = 0.5
                case(2)
                    score = (dist_val - score_min(iptcl_loc)) / score_spread(iptcl_loc)
            end select
        end function ref_score

        subroutine reset_ref_frontier()
            do iref = 1,self%nrefs
                iref_dist(iref) = ref_score(iref, stab_inds(1,iref))
            enddo
            iref_dist_inds = 1
            ptcl_avail     = .true.
        end subroutine reset_ref_frontier

        subroutine advance_ref_head( iref_current, state_filter )
            integer, intent(in) :: iref_current, state_filter
            integer :: iptcl_current
            do while( iref_dist_inds(iref_current) <= self%nptcls )
                iptcl_current = stab_inds(iref_dist_inds(iref_current), iref_current)
                if( ptcl_avail(iptcl_current) .and. &
                    &(state_filter == 0 .or. greedy_state(iptcl_current) == state_filter) )then
                    iref_dist(iref_current) = ref_score(iref_current, iptcl_current)
                    return
                endif
                iref_dist_inds(iref_current) = iref_dist_inds(iref_current) + 1
            enddo
            iref_dist(iref_current) = huge(iref_dist(iref_current))
        end subroutine advance_ref_head

        subroutine assign_current_ref()
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
        end subroutine assign_current_ref

        subroutine assign_greedy_state_labels()
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                assigned_iref  = minloc(iref_dist, dim=1)
                assigned_ptcl  = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
                greedy_state(assigned_ptcl) = self%sinds(assigned_iref)
                ptcl_avail(assigned_ptcl)  = .false.
                do iref = 1,self%nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
        end subroutine assign_greedy_state_labels

        subroutine gather_active_refs( state_filter )
            integer, intent(in) :: state_filter
            nactive = 0
            do iref = 1,self%nrefs
                if( self%sinds(iref) /= state_filter ) cycle
                if( iref_dist_inds(iref) > self%nptcls ) cycle
                nactive = nactive + 1
                active_refs(nactive)  = iref
                active_dists(nactive) = iref_dist(iref)
            enddo
        end subroutine gather_active_refs

        subroutine assign_refs_for_state( state_filter )
            integer, intent(in) :: state_filter
            call reset_ref_frontier()
            ptcl_avail = greedy_state == state_filter
            do iref = 1,self%nrefs
                if( self%sinds(iref) == state_filter ) call advance_ref_head(iref, state_filter)
            enddo
            do while( any(ptcl_avail) )
                call gather_active_refs(state_filter)
                if( nactive == 0 ) THROW_HARD('no active refs left in multi-state probability assignment')
                if( nactive == 1 )then
                    active_idx = 1
                else
                    active_idx = angle_sampling(active_dists(1:nactive), active_dists_sorted(1:nactive),&
                        &active_inds(1:nactive), state_projs_athres(state_filter), self%p_ptr%prob_athres)
                endif
                assigned_iref = active_refs(active_idx)
                call assign_current_ref()
                do iref = 1,self%nrefs
                    if( self%sinds(iref) == state_filter ) call advance_ref_head(iref, state_filter)
                enddo
            enddo
        end subroutine assign_refs_for_state
    end subroutine ref_assign

    ! state normalization (same energy) of the state_tab
    ! [0,1] normalization of the whole table
    subroutine state_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, istate
        ! normalize so prob of each ptcl is between [0,1] for all states
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%state_tab(:,i)%dist)
            if( sum_dist_all < TINY )then
                self%state_tab(:,i)%dist = 0.
            else
                self%state_tab(:,i)%dist = self%state_tab(:,i)%dist / sum_dist_all
            endif
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        max_dist = 0.
        min_dist = huge(min_dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)&
        !$omp reduction(min:min_dist) reduction(max:max_dist)
        do i = 1, self%nptcls
            max_dist = max(max_dist, maxval(self%state_tab(:,i)%dist, dim=1))
            min_dist = min(min_dist, minval(self%state_tab(:,i)%dist, dim=1))
        enddo
        !$omp end parallel do
        ! special case of numerical unstability of dist values
        if( (max_dist - min_dist) < TINY )then
            ! randomize dist so the assignment is stochastic
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(istate,i)
            do istate = 1, self%nstates
                do i = 1, self%nptcls
                    self%state_tab(istate,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            self%state_tab%dist = (self%state_tab%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine state_normalize

    ! ptcl -> state (using assigned iproj or previous iproj) assignment
    subroutine state_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, istate, assigned_istate, assigned_ptcl, state_dist_inds(self%nstates),&
                    &stab_inds(self%nptcls, self%nstates)
        real    :: sorted_tab(self%nptcls, self%nstates), state_dist(self%nstates)
        logical :: ptcl_avail(self%nptcls)
        if( self%nstates == 1 )then
            self%assgn_map = self%state_tab(1,:)
            return
        endif
        call self%state_normalize
        ! sorting each columns
        sorted_tab = transpose(self%state_tab%dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(istate,i)
        do istate = 1, self%nstates
            stab_inds(:,istate) = (/(i,i=1,self%nptcls)/)
            call hpsort(sorted_tab(:,istate), stab_inds(:,istate))
        enddo
        !$omp end parallel do
        ! first row is the current best state distribution
        state_dist_inds = 1
        state_dist      = sorted_tab(1,:)
        ptcl_avail      = .true.
        do while( any(ptcl_avail) )
            ! choose next istate to assign !!! SHOULD DO PROBABILISTIC SAMPLING HERE
            assigned_istate = minloc(state_dist, dim=1)
            assigned_ptcl   = stab_inds(state_dist_inds(assigned_istate), assigned_istate)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%state_tab(assigned_istate,assigned_ptcl)
            ! update the state_dist and state_dist_inds
            do istate = 1, self%nstates
                do while( state_dist_inds(istate) < self%nptcls .and. .not.(ptcl_avail(stab_inds(state_dist_inds(istate), istate))))
                    state_dist_inds(istate) = state_dist_inds(istate) + 1
                    state_dist(istate)      = sorted_tab(state_dist_inds(istate), istate)
                enddo
            enddo
        enddo
    end subroutine state_assign

    ! FILE IO

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_tab( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        type(ptcl_ref) :: ptcl_ref_sample
        integer         :: funit, addr, io_stat, file_header(2)
        integer         :: first_ptcl, last_ptcl
        integer(kind=8) :: addr8, expected_bytes, file_bytes, elem_bytes
        file_header(1) = self%nrefs
        file_header(2) = self%nptcls
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        call fileiochk('simple_eul_prob_tab; write_tab; file: '//binfname%to_char(), io_stat)
        write(unit=funit,pos=1,iostat=io_stat) file_header
        call fileiochk('simple_eul_prob_tab; write_tab header; file: '//binfname%to_char(), io_stat)
        addr = sizeof(file_header) + 1
        call write_seed_shift_table(funit, addr, self%seed_nrots, self%seed_shifts, self%seed_has_sh)
        addr8 = int(addr,kind=8)
        elem_bytes = int(sizeof(ptcl_ref_sample), kind=8)
        do first_ptcl = 1, self%nptcls, PROB_TAB_IO_CHUNK
            last_ptcl = min(first_ptcl + PROB_TAB_IO_CHUNK - 1, self%nptcls)
            write(funit,pos=addr8,iostat=io_stat) self%loc_tab(:,first_ptcl:last_ptcl)
            call fileiochk('simple_eul_prob_tab; write_tab loc_tab; file: '//binfname%to_char(), io_stat)
            addr8 = addr8 + int(self%nrefs,kind=8) * int(last_ptcl-first_ptcl+1,kind=8) * elem_bytes
        enddo
        expected_bytes = addr8 - 1_8
        inquire(unit=funit, size=file_bytes, iostat=io_stat)
        if( io_stat == 0 .and. file_bytes < expected_bytes )then
            write(logfhandle,*) 'prob_tab write size check failed: ', binfname%to_char()
            write(logfhandle,*) 'actual/expected bytes: ', file_bytes, expected_bytes
            THROW_HARD('probability table write truncated: '//binfname%to_char())
        endif
        call fclose(funit)
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global reg object's dist value table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: mat_chunk(:,:)
        type(ptcl_ref)                     :: ptcl_ref_sample
        real,                allocatable   :: seed_shifts_loc(:,:)
        logical,             allocatable   :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, nrefs_loc, i_loc, i_glob, pind, max_pind, seed_nrots_loc
        integer :: first_loc, last_loc, nchunk, i
        integer(kind=8) :: addr8, file_bytes, expected_bytes, tab_bytes, elem_bytes
        if( file_exists(binfname) )then
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        ! reading header and the nprojs/nptcls in this partition file
        read(unit=funit,pos=1) file_header
        nrefs_loc  = file_header(1)
        nptcls_loc = file_header(2)
        if( nrefs_loc .ne. self%nrefs ) THROW_HARD( 'nrefs should be the same as nrefs in this partition file!' )
        allocate(seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        ! read partition information
        addr = sizeof(file_header) + 1
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        addr8 = int(addr,kind=8)
        elem_bytes     = int(sizeof(ptcl_ref_sample),kind=8)
        tab_bytes      = int(nrefs_loc,kind=8) * int(nptcls_loc,kind=8) * elem_bytes
        expected_bytes = addr8 + tab_bytes - 1_8
        file_bytes = -1_8
        inquire(file=binfname%to_char(), size=file_bytes, iostat=io_stat)
        if( io_stat /= 0 .or. file_bytes < expected_bytes )then
            write(logfhandle,*) 'prob_tab read size check failed: ', binfname%to_char()
            write(logfhandle,*) 'nrefs/nptcls/actual/expected bytes: ', nrefs_loc, nptcls_loc, file_bytes, expected_bytes
            THROW_HARD('probability table is missing or truncated: '//binfname%to_char())
        endif
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab%read_tab_to_glob')
        max_pind = maxval(self%pinds)
        if( max_pind < 1 )then
            call fclose(funit)
            deallocate(seed_shifts_loc, seed_has_sh_loc)
            return
        endif
        allocate(pind2glob(max_pind), source=0)
        do i = 1, size(self%pinds)
            pind = self%pinds(i)
            if( pind > 0 .and. pind <= max_pind ) pind2glob(pind) = i
        enddo
        allocate(mat_chunk(nrefs_loc, min(PROB_TAB_IO_CHUNK, nptcls_loc)))
        do first_loc = 1, nptcls_loc, PROB_TAB_IO_CHUNK
            last_loc = min(first_loc + PROB_TAB_IO_CHUNK - 1, nptcls_loc)
            nchunk   = last_loc - first_loc + 1
            read(unit=funit,pos=addr8,iostat=io_stat) mat_chunk(:,1:nchunk)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob loc_tab; file: '//binfname%to_char(), io_stat)
            addr8 = addr8 + int(nrefs_loc,kind=8) * int(nchunk,kind=8) * elem_bytes
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
            do i_loc = 1, nchunk
                pind = mat_chunk(1,i_loc)%pind
                if( pind < 1 .or. pind > max_pind ) cycle
                i_glob = pind2glob(pind)
                if( i_glob > 0 )then
                    self%loc_tab(:,i_glob)    = mat_chunk(:,i_loc)
                    self%seed_shifts(:,i_glob)= seed_shifts_loc(:,first_loc+i_loc-1)
                    self%seed_has_sh(i_glob)  = seed_has_sh_loc(first_loc+i_loc-1)
                endif
            end do
            !$omp end parallel do
        enddo
        call fclose(funit)
        deallocate(mat_chunk, seed_shifts_loc, seed_has_sh_loc, pind2glob)
    end subroutine read_tab_to_glob

    subroutine write_state_tab( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%state_tab
        call fclose(funit)
    end subroutine write_state_tab

    subroutine read_state_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: state_tab_glob(:,:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_state_tab; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(state_tab_glob(self%nstates,nptcls_glob))
        read(unit=funit,pos=headsz + 1) state_tab_glob
        call fclose(funit)
        pinds_glob = state_tab_glob(1,:)%pind
        call build_pind_lookup(pinds_glob, self%pinds, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(state_tab_glob, pind2glob, pinds_glob)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, self%nptcls
            pind = self%pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 ) self%state_tab(:,i_loc) = state_tab_glob(:,i_glob)
        end do
        !$omp end parallel do
        deallocate(state_tab_glob, pind2glob, pinds_glob)
    end subroutine read_state_tab

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: assgn_glob(:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit,pos=headsz + 1) assgn_glob
        call fclose(funit)
        pinds_glob = assgn_glob(:)%pind
        call build_pind_lookup(pinds_glob, self%pinds, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(assgn_glob, pind2glob, pinds_glob)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, self%nptcls
            pind = self%pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 ) self%assgn_map(i_loc) = assgn_glob(i_glob)
        end do
        !$omp end parallel do
        deallocate(assgn_glob, pind2glob, pinds_glob)
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        if( allocated(self%loc_tab)      ) deallocate(self%loc_tab)
        if( allocated(self%state_tab)    ) deallocate(self%state_tab)
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%ssinds)       ) deallocate(self%ssinds)
        if( allocated(self%sinds)        ) deallocate(self%sinds)
        if( allocated(self%jinds)        ) deallocate(self%jinds)
        if( allocated(self%state_exists) ) deallocate(self%state_exists)
        if( allocated(self%proj_exists)  ) deallocate(self%proj_exists)
        self%seed_nrots = 0
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

end module simple_eul_prob_tab
