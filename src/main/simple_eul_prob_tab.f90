! orientation eul_prob_tab, used in refine3D
module simple_eul_prob_tab
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
implicit none

public :: eul_prob_tab
public :: calc_num2sample, calc_athres, eulprob_dist_switch, eulprob_corr_switch, angle_sampling
private
#include "simple_local_flags.inc"

interface angle_sampling
    module procedure angle_sampling_1
    module procedure angle_sampling_2
end interface

type :: eul_prob_tab
    type(ptcl_ref), allocatable :: loc_tab(:,:)    !< 2D search table (nspace*nstates, nptcls)
    type(ptcl_ref), allocatable :: state_tab(:,:)  !< 2D search table (nstates,        nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
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
    procedure, private :: new_1, new_2
    generic            :: new => new_1, new_2
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

    subroutine new_1( self, pinds, empty_okay )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: pinds(:)
        logical, optional,   intent(in)    :: empty_okay
        integer, parameter :: MIN_POP = 5   ! ignoring cavgs with less than 5 particles
        integer :: i, iproj, iptcl, istate, si, ri
        real    :: x
        logical :: l_empty
        l_empty = (trim(params_glob%empty3Dcavgs) .eq. 'yes')
        if( present(empty_okay) ) l_empty = empty_okay
        call self%kill
        self%nptcls       = size(pinds)
        self%state_exists = build_glob%spproj_field%states_exist(params_glob%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        if( l_empty )then
            allocate(self%proj_exists(params_glob%nspace,params_glob%nstates), source=.true.)
        else
            self%proj_exists = build_glob%spproj_field%projs_exist(params_glob%nstates,params_glob%nspace, thres=MIN_POP)
        endif
        self%nrefs = count(self%proj_exists .eqv. .true.)
        allocate(self%ssinds(self%nstates),self%jinds(self%nrefs),self%sinds(self%nrefs))
        si = 0
        ri = 0
        do istate = 1, params_glob%nstates
            if( .not. self%state_exists(istate) )cycle
            si              = si + 1
            self%ssinds(si) = istate
            do iproj = 1,params_glob%nspace
                if( .not. self%proj_exists(iproj,istate) )cycle
                ri             = ri + 1
                self%jinds(ri) = iproj
                self%sinds(ri) = istate
            enddo
        enddo
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%loc_tab(self%nrefs,self%nptcls), self%assgn_map(self%nptcls),self%state_tab(self%nstates,self%nptcls))
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
        end do
        !$omp end parallel do
    end subroutine new_1

    subroutine new_2( self, os )
        class(eul_prob_tab), intent(inout) :: self
        class(oris),         intent(in)    :: os
        integer, allocatable :: states(:)
        integer :: i, iproj, iptcl, istate, n, iref
        real    :: x
        call self%kill
        self%nptcls  = os%get_noris(consider_state=.true.)
        self%nstates = params_glob%nstates
        n = os%get_noris()
        allocate(states(n), source=nint(os%get_all('state')))
        self%pinds = pack((/(i,i=1,n)/), mask=states > 0)
        deallocate(states)
        allocate(self%loc_tab(self%nrefs,self%nptcls), self%assgn_map(self%nptcls),self%state_tab(self%nstates,self%nptcls))
        !$omp parallel do default(shared) private(i,iptcl,istate,iproj,iref) proc_bind(close) schedule(static)
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
            do istate = 1,self%nstates
                self%state_tab(istate,i)%pind   = iptcl
                self%state_tab(istate,i)%istate = istate
                self%state_tab(istate,i)%iproj  = 0
                self%state_tab(istate,i)%inpl   = 0
                self%state_tab(istate,i)%dist   = huge(x)
                self%state_tab(istate,i)%x      = 0.
                self%state_tab(istate,i)%y      = 0.
                self%state_tab(istate,i)%has_sh = .false.
                do iproj = 1,params_glob%nspace
                    iref = (istate-1)*params_glob%nspace + iproj
                    self%loc_tab(iref,i)%pind   = iptcl
                    self%loc_tab(iref,i)%istate = istate
                    self%loc_tab(iref,i)%iproj  = iproj
                    self%loc_tab(iref,i)%inpl   = 0
                    self%loc_tab(iref,i)%dist   = huge(x)
                    self%loc_tab(iref,i)%x      = 0.
                    self%loc_tab(iref,i)%y      = 0.
                    self%loc_tab(iref,i)%has_sh = .false.
                end do
            end do
        end do
        !$omp end parallel do
        if( file_exists(DIST_FBODY//'.dat') )then
            call self%read_tab_to_glob(string(DIST_FBODY//'.dat'))
        endif
    end subroutine new_2

    ! partition-wise table filling, used only in shared-memory commander 'exec_prob_tab'
    subroutine fill_tab( self, pftc )
        use simple_polarft_calc,  only: polarft_calc
        use simple_pftc_shsrch_grad, only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
        class(eul_prob_tab),     intent(inout) :: self
        class(polarft_calc), intent(inout) :: pftc
        integer,                 allocatable   :: locn(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
        type(ori)               :: o_prev
        integer :: i, si, ri, j, iproj, iptcl, n, projs_ns, ithr, irot, inds_sorted(pftc%get_nrots(),nthr_glob),&
                  &istate, iref_start
        logical :: l_doshift
        real    :: rotmat(2,2), lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rot_xy(2), inpl_athres(params_glob%nstates)
        real    :: dists_inpl(pftc%get_nrots(),nthr_glob), dists_inpl_sorted(pftc%get_nrots(),nthr_glob), dists_refs(self%nrefs,nthr_glob)
        call seed_rnd
        projs_ns = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(params_glob%nspace, 'dist', n, state=istate)
            projs_ns            = max(projs_ns, n)
            inpl_athres(istate) = calc_athres('dist_inpl', state=istate)
        enddo
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), source=0)
        if( params_glob%l_sh_first .and. params_glob%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                    &maxits=params_glob%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,istate,irot,iproj,iref_start,cxy,ri,j,cxy_prob,rot_xy,rotmat)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! (1) identify shifts using the previously assigned best reference
                call build_glob%spproj_field%get_ori(iptcl, o_prev)        ! previous ori
                istate     = o_prev%get_state()
                irot       = pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj      = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
                if( self%state_exists(istate) .and. self%proj_exists(iproj,istate) )then
                    iref_start = (istate-1)*params_glob%nspace
                    ! BFGS over shifts
                    call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                    if( irot == 0 ) cxy(2:3) = 0.
                else
                    cxy(2:3) = 0.
                endif
                ! (2) search projection directions using those shifts for all references
                do ri = 1, self%nrefs
                    istate = self%sinds(ri)
                    iproj  = self%jinds(ri)
                    call pftc%gen_objfun_vals((istate-1)*params_glob%nspace + iproj, iptcl, cxy(2:3), dists_inpl(:,ithr))
                    dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                    irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres(istate))
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(pftc%get_rot(irot), rotmat)
                    rot_xy                    = matmul(cxy(2:3), rotmat)
                    self%loc_tab(ri,i)%dist   = dists_inpl(irot,ithr)
                    dists_refs(  ri,ithr)     = dists_inpl(irot,ithr)
                    self%loc_tab(ri,i)%inpl   = irot
                    self%loc_tab(ri,i)%x      = rot_xy(1)
                    self%loc_tab(ri,i)%y      = rot_xy(2)
                    self%loc_tab(ri,i)%has_sh = .true.
                enddo
                ! (3) see if we can refine the shifts by re-searching them for individual references in the 
                !     identified probabilistic neighborhood
                if( params_glob%l_prob_sh )then
                    locn(:,ithr) = minnloc(dists_refs(:,ithr), projs_ns)
                    do j = 1,projs_ns
                        ri     = locn(j,ithr)
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        ! BFGS over shifts
                        call grad_shsrch_obj(ithr)%set_indices((istate-1)*params_glob%nspace + iproj, iptcl)
                        irot     = self%loc_tab(ri,i)%inpl
                        cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=cxy(2:3))
                        if( irot > 0 )then
                            self%loc_tab(ri,i)%inpl   = irot
                            self%loc_tab(ri,i)%dist   = eulprob_dist_switch(cxy_prob(1))
                            self%loc_tab(ri,i)%x      = cxy_prob(2)
                            self%loc_tab(ri,i)%y      = cxy_prob(3)
                            self%loc_tab(ri,i)%has_sh = .true.
                        endif
                    end do
                endif
            enddo
            !$omp end parallel do
        else
            l_doshift = params_glob%l_prob_sh .and. params_glob%l_doshift
            if( l_doshift )then
                ! make shift search objects
                lims(:,1)      = -params_glob%trs
                lims(:,2)      =  params_glob%trs
                lims_init(:,1) = -SHC_INPL_TRSHWDTH
                lims_init(:,2) =  SHC_INPL_TRSHWDTH
                do ithr = 1,nthr_glob
                    call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                        &maxits=params_glob%maxits_sh, opt_angle=.true.)
                end do
                ! fill the table
                !$omp parallel do default(shared) private(i,iptcl,ithr,ri,istate,iproj,irot,j,cxy) proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    do ri = 1, self%nrefs
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        call pftc%gen_objfun_vals((istate-1)*params_glob%nspace + iproj, iptcl, [0.,0.], dists_inpl(:,ithr))
                        dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                        irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres(istate))
                        self%loc_tab(ri,i)%dist = dists_inpl(irot,ithr)
                        self%loc_tab(ri,i)%inpl = irot
                        dists_refs(  ri,ithr)   = dists_inpl(irot,ithr)
                    enddo
                    locn(:,ithr) = minnloc(dists_refs(:,ithr), projs_ns)
                    do j = 1,projs_ns
                        ri     = locn(j,ithr)
                        istate = self%sinds(ri)
                        iproj  = self%jinds(ri)
                        ! BFGS over shifts
                        call grad_shsrch_obj(ithr)%set_indices((istate-1)*params_glob%nspace + iproj, iptcl)
                        irot = self%loc_tab(ri,i)%inpl
                        cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                        if( irot > 0 )then
                            self%loc_tab(ri,i)%inpl = irot
                            self%loc_tab(ri,i)%dist = eulprob_dist_switch(cxy(1))
                            self%loc_tab(ri,i)%x    = cxy(2)
                            self%loc_tab(ri,i)%y    = cxy(3)
                        endif
                        self%loc_tab(ri,i)%has_sh = .true.
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
                        call pftc%gen_objfun_vals((istate-1)*params_glob%nspace + iproj, iptcl, [0.,0.], dists_inpl(:,ithr))
                        dists_inpl(:,ithr)      = eulprob_dist_switch(dists_inpl(:,ithr))
                        irot                    = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres(istate))
                        self%loc_tab(ri,i)%dist = dists_inpl(irot,ithr)
                        self%loc_tab(ri,i)%inpl = irot
                    enddo
                enddo
                !$omp end parallel do
            endif
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab

    subroutine fill_tab_state_only( self, pftc )
        use simple_polarft_calc,  only: polarft_calc
        use simple_pftc_shsrch_grad, only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
        class(eul_prob_tab),     intent(inout) :: self
        class(polarft_calc), intent(inout) :: pftc
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)  !< origin shift search object, L-BFGS with gradient
        type(ori)               :: o_prev
        integer :: i, iproj, iptcl, ithr, irot, istate, iref_start, is
        real    :: lims(2,2), lims_init(2,2), cxy(3)
        call seed_rnd
        if( params_glob%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                    &maxits=params_glob%maxits_sh, opt_angle=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,iproj,is,istate,irot,iref_start,cxy)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! identify shifts using the previously assigned best reference
                call build_glob%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                iproj = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate     = self%ssinds(is)
                    irot       = pftc%get_roind(360.-o_prev%e3get()) ! in-plane angle index
                    iref_start = (istate-1)*params_glob%nspace
                    ! BFGS over shifts
                    call grad_shsrch_obj(ithr)%set_indices(iref_start + iproj, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        irot     = pftc%get_roind(360.-o_prev%e3get())
                        cxy(1)   = real(pftc%gen_corr_for_rot_8(iref_start+iproj, iptcl, irot))
                        cxy(2:3) = 0.
                    endif
                    self%state_tab(istate,i)%dist   = eulprob_dist_switch(cxy(1))
                    self%state_tab(istate,i)%iproj  = iproj
                    self%state_tab(istate,i)%inpl   = irot
                    self%state_tab(istate,i)%x      = cxy(2)
                    self%state_tab(istate,i)%y      = cxy(3)
                    self%state_tab(istate,i)%has_sh = .true.
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
                call build_glob%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                irot  = pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate      = self%ssinds(is)
                    iref_start  = (istate-1)*params_glob%nspace
                    self%state_tab(istate,i)%dist   = eulprob_dist_switch(real(pftc%gen_corr_for_rot_8(iref_start+iproj, iptcl, irot)))
                    self%state_tab(istate,i)%iproj  = iproj
                    self%state_tab(istate,i)%inpl   = irot
                    self%state_tab(istate,i)%x      = 0.
                    self%state_tab(istate,i)%y      = 0.
                    self%state_tab(istate,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab_state_only

    ! reference normalization (same energy) of the loc_tab
    ! [0,1] normalization
    subroutine ref_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, iref
        ! normalize so prob of each ptcl is between [0,1] for all projs
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
        ! min/max normalization to obtain values between 0 and 1
        !$omp parallel workshare proc_bind(close)
        min_dist = minval(self%loc_tab(:,:)%dist)
        max_dist = maxval(self%loc_tab(:,:)%dist)
        !$omp end parallel workshare
        ! special case of numerical unstability of dist values
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab')
            ! randomize dist so the assignment is stochastic
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
        integer :: i, iref, assigned_iref, assigned_ptcl, istate,&
                    &stab_inds(self%nptcls, self%nrefs), inds_sorted(self%nrefs), iref_dist_inds(self%nrefs)
        real    :: sorted_tab(self%nptcls, self%nrefs), projs_athres,iref_dist(self%nrefs), dists_sorted(self%nrefs)
        logical :: ptcl_avail(self%nptcls)
        ! normalization
        call self%ref_normalize
        ! sorting each columns
        sorted_tab = transpose(self%loc_tab(:,:)%dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i)
        do iref = 1, self%nrefs
            stab_inds(:,iref) = (/(i,i=1,self%nptcls)/)
            call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        projs_athres = 0.
        do istate = 1, self%nstates
            projs_athres = max(projs_athres, calc_athres('dist', state=istate))
        enddo
        ! first row is the current best reference distribution
        iref_dist_inds = 1
        iref_dist      = sorted_tab(1,:)
        ptcl_avail     = .true.
        do while( any(ptcl_avail) )
            ! sampling the ref distribution to choose next iref to assign
            assigned_iref = angle_sampling(iref_dist, dists_sorted, inds_sorted, projs_athres)
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
            ! update the iref_dist and iref_dist_inds
            do iref = 1, self%nrefs
                do while( iref_dist_inds(iref) < self%nptcls .and. .not.(ptcl_avail(stab_inds(iref_dist_inds(iref),iref))))
                    iref_dist_inds(iref) = iref_dist_inds(iref) + 1
                    iref_dist(     iref) = sorted_tab(iref_dist_inds(iref), iref)
                enddo
            enddo
        enddo
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
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = self%nrefs
        file_header(2) = self%nptcls
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit,pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global reg object's dist value table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: mat_loc(:,:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, nrefs_loc, i_loc, i_glob
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
        allocate(mat_loc(nrefs_loc, nptcls_loc))
        ! read partition information
        addr = sizeof(file_header) + 1
        read(unit=funit,pos=addr) mat_loc
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_glob = 1, self%nptcls
            do i_loc = 1, nptcls_loc
                if( mat_loc(1,i_loc)%pind == self%loc_tab(1,i_glob)%pind )then
                    self%loc_tab(:,i_glob) = mat_loc(:,i_loc)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        deallocate(mat_loc)
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
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
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
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( self%state_tab(1,i_loc)%pind == state_tab_glob(1,i_glob)%pind )then
                    self%state_tab(:,i_loc) = state_tab_glob(:,i_glob)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        deallocate(state_tab_glob)
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
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
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
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( self%assgn_map(i_loc)%pind == assgn_glob(i_glob)%pind )then
                    self%assgn_map(i_loc) = assgn_glob(i_glob)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        if( allocated(self%loc_tab)      ) deallocate(self%loc_tab)
        if( allocated(self%state_tab)    ) deallocate(self%state_tab)
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%ssinds)       ) deallocate(self%ssinds)
        if( allocated(self%sinds)        ) deallocate(self%sinds)
        if( allocated(self%jinds)        ) deallocate(self%jinds)
        if( allocated(self%state_exists) ) deallocate(self%state_exists)
        if( allocated(self%proj_exists)  ) deallocate(self%proj_exists)
    end subroutine kill

    ! PUBLIC UTILITITES

    subroutine calc_num2sample( num_all, field_str, num_smpl, state)
        integer,           intent(in)  :: num_all
        character(len=*),  intent(in)  :: field_str
        integer,           intent(out) :: num_smpl
        integer, optional, intent(in)  :: state
        real :: athres
        athres   = calc_athres(field_str, state)
        num_smpl = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_num2sample

    function calc_athres( field_str, state ) result( athres )
        character(len=*),  intent(in) :: field_str
        integer, optional, intent(in) :: state
        real, allocatable :: vals(:)
        real :: athres, dist_thres
        vals       = build_glob%spproj_field%get_all_sampled(trim(field_str), state=state)
        dist_thres = sum(vals) / real(size(vals))
        athres     = params_glob%prob_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
    end function calc_athres

    ! switch corr in [0,1] to [0, infinity) to do greedy_sampling
    elemental function eulprob_dist_switch( corr ) result(dist)
        real, intent(in) :: corr
        real :: dist
        dist = corr
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                if( corr < 0. )then
                    dist = 0.
                else
                    dist = corr
                endif
                dist = 1. - dist
            case(OBJFUN_EUCLID)
                if( corr < TINY )then
                    dist = huge(dist)
                else
                    dist = - log(corr)
                endif
        end select
    end function eulprob_dist_switch

    ! switch corr in [0,1] to [0, infinity) to do greedy_sampling
    elemental function eulprob_corr_switch( dist ) result(corr)
        real, intent(in) :: dist
        real :: corr
        corr = dist
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                corr = 1 - dist
            case(OBJFUN_EUCLID)
                corr = exp(-dist)
        end select
    end function eulprob_corr_switch

    function angle_sampling_1( pvec, athres_ub_in ) result( which )
        real,    intent(in)  :: pvec(:)        !< probabilities
        real,    intent(in)  :: athres_ub_in
        real,    allocatable :: pvec_sorted(:)
        integer, allocatable :: sorted_inds(:)
        integer :: which, n
        n = size(pvec)
        allocate(pvec_sorted(n),sorted_inds(n))
        which = angle_sampling_2(pvec, pvec_sorted, sorted_inds, athres_ub_in)
    end function angle_sampling_1

    function angle_sampling_2( pvec, pvec_sorted, sorted_inds, athres_ub_in ) result( which )
        real,    intent(in)    :: pvec(:)        !< probabilities
        real,    intent(inout) :: pvec_sorted(:) !< sorted probabilities
        integer, intent(inout) :: sorted_inds(:)
        real,    intent(in)    :: athres_ub_in
        integer :: which, num_lb, num_ub, n
        real    :: athres_ub, athres_lb
        n         = size(pvec)
        athres_ub = min(params_glob%prob_athres, athres_ub_in)
        athres_lb = min(athres_ub / 10., 1.) ! athres lower bound is 1/10 of athres upper bound, max at 1 degree
        num_ub    = min(n,max(1,int(athres_ub * real(n) / 180.)))
        num_lb    = 1 + floor(athres_lb / athres_ub * num_ub)
        which     = greedy_sampling(pvec, pvec_sorted, sorted_inds, num_ub, num_lb)
    end function angle_sampling_2

end module simple_eul_prob_tab
