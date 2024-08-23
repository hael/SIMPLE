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

integer, parameter :: SHIFT_NUM = 50
integer, parameter :: NSHIFTS   = SHIFT_NUM**2    !< number of discretized shift space

type :: eul_prob_tab
    type(ptcl_ref), allocatable :: loc_tab(:,:,:) !< 3D search table (nspace,  nptcls, nstates)
    type(ptcl_ref), allocatable :: state_tab(:,:) !< 2D search table (nstates, nptcls)
    type(ptcl_ref), allocatable :: shift_tab(:,:) !< 2D search table (nshifts, nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)   !< assignment map           (nptcls)
    integer,        allocatable :: pinds(:)       !< particle indices for processing
    integer                     :: nptcls         !< size of pinds array
    integer                     :: nstates        !< states number
    contains
    ! CONSTRUCTOR
    procedure, private :: new_1, new_2
    generic            :: new => new_1, new_2
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: fill_tab
    procedure :: write_tab, write_shift
    procedure :: read_assignment, read_shift
    procedure :: trim_tab
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: read_tab_to_glob
    procedure :: prob_assign
    procedure :: fill_shift_tab
    procedure :: shift_assign
    procedure :: write_assignment
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: proj_normalize
    procedure, private :: proj_assign
    procedure, private :: state_normalize
    procedure, private :: state_assign
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new_1( self, pinds )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: pinds(:)
        integer :: i, iproj, iptcl, istate, ishift
        real    :: x
        call self%kill
        self%nptcls  = size(pinds)
        self%nstates = params_glob%nstates
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%loc_tab(params_glob%nspace,self%nptcls,self%nstates), self%assgn_map(self%nptcls),&
                    &self%state_tab(self%nstates,self%nptcls), self%shift_tab(NSHIFTS,self%nptcls))
        !$omp parallel do default(shared) private(i,iptcl,istate,iproj,ishift) proc_bind(close) schedule(static)
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
                    self%loc_tab(iproj,i,istate)%pind   = iptcl
                    self%loc_tab(iproj,i,istate)%istate = istate
                    self%loc_tab(iproj,i,istate)%iproj  = iproj
                    self%loc_tab(iproj,i,istate)%inpl   = 0
                    self%loc_tab(iproj,i,istate)%dist   = huge(x)
                    self%loc_tab(iproj,i,istate)%x      = 0.
                    self%loc_tab(iproj,i,istate)%y      = 0.
                    self%loc_tab(iproj,i,istate)%has_sh = .false.
                end do
            end do
            do ishift = 1,NSHIFTS
                self%shift_tab(ishift,i)%pind   = iptcl
                self%shift_tab(ishift,i)%istate = 0
                self%shift_tab(ishift,i)%iproj  = 0
                self%shift_tab(ishift,i)%inpl   = 0
                self%shift_tab(ishift,i)%dist   = huge(x)
                self%shift_tab(ishift,i)%x      = 0.
                self%shift_tab(ishift,i)%y      = 0.
                self%shift_tab(ishift,i)%has_sh = .false.
            end do
        end do
        !$omp end parallel do 
    end subroutine new_1

    subroutine new_2( self, os )
        class(eul_prob_tab), intent(inout) :: self
        class(oris),         intent(in)    :: os
        integer, allocatable :: states(:)
        integer :: i, iproj, iptcl, istate, n
        real    :: x
        call self%kill
        self%nptcls  = os%get_noris(consider_state=.true.)
        self%nstates = params_glob%nstates
        n = os%get_noris()
        allocate(states(n), source=nint(os%get_all('state')))
        self%pinds = pack((/(i,i=1,n)/), mask=states > 0)
        deallocate(states)
        allocate(self%loc_tab(params_glob%nspace,self%nptcls,self%nstates), self%assgn_map(self%nptcls),&
                    &self%state_tab(self%nstates,self%nptcls))
        !$omp parallel do default(shared) private(i,iptcl,istate,iproj) proc_bind(close) schedule(static)
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
                    self%loc_tab(iproj,i,istate)%pind   = iptcl
                    self%loc_tab(iproj,i,istate)%istate = istate
                    self%loc_tab(iproj,i,istate)%iproj  = iproj
                    self%loc_tab(iproj,i,istate)%inpl   = 0
                    self%loc_tab(iproj,i,istate)%dist   = huge(x)
                    self%loc_tab(iproj,i,istate)%x      = 0.
                    self%loc_tab(iproj,i,istate)%y      = 0.
                    self%loc_tab(iproj,i,istate)%has_sh = .false.
                end do
            end do
        end do
        !$omp end parallel do
        if( file_exists(trim(DIST_FBODY)//'.dat') )then
            call self%read_tab_to_glob(trim(DIST_FBODY)//'.dat')
        endif
    end subroutine new_2

    ! partition-wise table filling, used only in shared-memory commander 'exec_prob_tab'
    subroutine fill_tab( self, pftcc )
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        class(eul_prob_tab),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 allocatable   :: locn(:)
        type(pftcc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
        type(ori)               :: o_prev
        integer :: i, j, iproj, iptcl, projs_ns, ithr, irot, inds_sorted(pftcc%nrots,nthr_glob), istate, iref
        logical :: l_doshift
        real    :: dists_inpl(pftcc%nrots,nthr_glob), dists_inpl_sorted(pftcc%nrots,nthr_glob), rotmat(2,2)
        real    :: dists_projs(params_glob%nspace,nthr_glob), lims(2,2), lims_init(2,2), cxy(3), rot_xy(2), inpl_athres
        call seed_rnd
        if( trim(params_glob%sh_first).eq.'yes' )then
            ! make shift search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                    &maxits=params_glob%maxits_sh, opt_angle=(trim(params_glob%sh_opt_angle).eq.'yes'), coarse_init=.true.)
            end do
            ! fill the table
            do istate = 1, self%nstates
                iref        = (istate-1)*params_glob%nspace
                inpl_athres = calc_athres('dist_inpl', state=istate)
                call calc_num2sample(params_glob%nspace, 'dist', projs_ns, state=istate)
                if( allocated(locn) ) deallocate(locn)
                allocate(locn(projs_ns), source=0)
                !$omp parallel do default(shared) private(i,j,iptcl,ithr,o_prev,iproj,irot,cxy,rot_xy,rotmat,locn) proc_bind(close) schedule(static)
                do i = 1, self%nptcls
                    iptcl = self%pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    if( trim(params_glob%sh_glob) .eq. 'yes' )then  ! retrieve ptcl shift from assign_map
                        cxy(2:3) = [self%assgn_map(i)%x, self%assgn_map(i)%y]
                    else                                            ! using previous ori for shift search
                        call build_glob%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                        irot  = pftcc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                        iproj = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
                        ! BFGS over shifts
                        call grad_shsrch_obj(ithr)%set_indices(iref + iproj, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        if( irot < TINY ) cxy(2:3) = 0.
                    endif
                    if( params_glob%l_prob_sh )then
                        do iproj = 1, params_glob%nspace
                            call pftcc%gencorrs(iref + iproj, iptcl, cxy(2:3), dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres)
                            self%loc_tab(iproj,i,istate)%dist = dists_inpl(irot,ithr)
                            self%loc_tab(iproj,i,istate)%inpl = irot
                            dists_projs(iproj,ithr) = dists_inpl(irot,ithr)
                        enddo
                        locn = minnloc(dists_projs(:,ithr), projs_ns)
                        do j = 1,projs_ns
                            iproj = locn(j)
                            ! BFGS over shifts
                            call grad_shsrch_obj(ithr)%set_indices(iref + iproj, iptcl)
                            irot = self%loc_tab(iproj,i,istate)%inpl
                            cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                            if( irot > 0 )then
                                self%loc_tab(iproj,i,istate)%inpl = irot
                                self%loc_tab(iproj,i,istate)%dist = eulprob_dist_switch(cxy(1))
                                self%loc_tab(iproj,i,istate)%x    = cxy(2)
                                self%loc_tab(iproj,i,istate)%y    = cxy(3)
                            endif
                            self%loc_tab(iproj,i,istate)%has_sh = .true.
                        end do
                    else
                        do iproj = 1, params_glob%nspace
                            call pftcc%gencorrs(iref + iproj, iptcl, cxy(2:3), dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres)
                            ! rotate the shift vector to the frame of reference
                            call rotmat2d(pftcc%get_rot(irot), rotmat)
                            rot_xy = matmul(cxy(2:3), rotmat)
                            self%loc_tab(iproj,i,istate)%dist   = dists_inpl(irot,ithr)
                            self%loc_tab(iproj,i,istate)%inpl   = irot
                            self%loc_tab(iproj,i,istate)%x      = rot_xy(1)
                            self%loc_tab(iproj,i,istate)%y      = rot_xy(2)
                            self%loc_tab(iproj,i,istate)%has_sh = .true.
                        enddo
                    endif
                enddo
                !$omp end parallel do
            enddo
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
                        &maxits=params_glob%maxits_sh, opt_angle=(trim(params_glob%sh_opt_angle).eq.'yes'))
                end do
                ! fill the table
                do istate = 1, self%nstates
                    iref        = (istate-1)*params_glob%nspace
                    inpl_athres = calc_athres('dist_inpl', state=istate)
                    call calc_num2sample(params_glob%nspace, 'dist', projs_ns, state=istate)
                    if( allocated(locn) ) deallocate(locn)
                    allocate(locn(projs_ns), source=0)
                    !$omp parallel do default(shared) private(i,j,iptcl,ithr,iproj,irot,cxy,locn) proc_bind(close) schedule(static)
                    do i = 1, self%nptcls
                        iptcl = self%pinds(i)
                        ithr  = omp_get_thread_num() + 1
                        do iproj = 1, params_glob%nspace
                            call pftcc%gencorrs(iref + iproj, iptcl, dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres)
                            self%loc_tab(iproj,i,istate)%dist = dists_inpl(irot,ithr)
                            self%loc_tab(iproj,i,istate)%inpl = irot
                            dists_projs(iproj,ithr) = dists_inpl(irot,ithr)
                        enddo
                        locn = minnloc(dists_projs(:,ithr), projs_ns)
                        do j = 1,projs_ns
                            iproj = locn(j)
                            ! BFGS over shifts
                            call grad_shsrch_obj(ithr)%set_indices(iref + iproj, iptcl)
                            irot = self%loc_tab(iproj,i,istate)%inpl
                            cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                            if( irot > 0 )then
                                self%loc_tab(iproj,i,istate)%inpl = irot
                                self%loc_tab(iproj,i,istate)%dist = eulprob_dist_switch(cxy(1))
                                self%loc_tab(iproj,i,istate)%x    = cxy(2)
                                self%loc_tab(iproj,i,istate)%y    = cxy(3)
                            endif
                            self%loc_tab(iproj,i,istate)%has_sh = .true.
                        end do
                    enddo
                    !$omp end parallel do
                enddo
            else
                ! fill the table
                do istate = 1, self%nstates
                    iref        = (istate-1)*params_glob%nspace
                    inpl_athres = calc_athres('dist_inpl', state=istate)
                    !$omp parallel do default(shared) private(i,iptcl,ithr,iproj,irot) proc_bind(close) schedule(static)
                    do i = 1, self%nptcls
                        iptcl = self%pinds(i)
                        ithr  = omp_get_thread_num() + 1
                        do iproj = 1, params_glob%nspace
                            call pftcc%gencorrs(iref + iproj, iptcl, dists_inpl(:,ithr))
                            dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                            irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres)
                            self%loc_tab(iproj,i,istate)%dist = dists_inpl(irot,ithr)
                            self%loc_tab(iproj,i,istate)%inpl = irot
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
            endif
        endif
    end subroutine fill_tab

    ! ptcl -> (proj, state) assignment used in the global prob_align commander, in 'exec_prob_align'
    subroutine prob_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        call self%proj_normalize
        call self%proj_assign
        if( self%nstates > 1 )then
            call self%state_normalize
            call self%state_assign
        else
            self%assgn_map = self%state_tab(1,:)
        endif
    end subroutine prob_assign

    subroutine fill_shift_tab( self, pftcc )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(eul_prob_tab),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        type(ori) :: o_prev
        integer   :: i, iptcl, iproj, ithr, ix, iy, ishift, irot, inds_sorted(pftcc%nrots,nthr_glob)
        real      :: lims(2,2), stepx, stepy, x, y, xy(2), inpl_athres,&
                    &dists_inpl(pftcc%nrots,nthr_glob), dists_inpl_sorted(pftcc%nrots,nthr_glob),&
                    &shift_points(2,NSHIFTS)
        lims(:,1)   = -params_glob%trs
        lims(:,2)   =  params_glob%trs
        stepx       = real(lims(1,2) - lims(1,1), dp) / real(SHIFT_NUM, dp)
        stepy       = real(lims(2,2) - lims(2,1), dp) / real(SHIFT_NUM, dp)
        inpl_athres = calc_athres('dist_inpl', state=1)
        ! generating discretized shifts
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(ix,iy,x,y,ishift)
        do ix = 1,SHIFT_NUM
            x = lims(1,1) + stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,SHIFT_NUM
                ishift = (ix-1)*SHIFT_NUM + iy
                y      = lims(2,1) + stepy/2. + real(iy-1,dp)*stepy
                shift_points(:,ishift) = [x,y]
            end do
        end do
        !$omp end parallel do
        ! filling the tab
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,iptcl,ithr,iproj,ishift,irot,xy,o_prev)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            call build_glob%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
            iproj = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
            ithr  = omp_get_thread_num() + 1
            do ishift = 1, NSHIFTS
                xy = shift_points(:,ishift)
                call pftcc%gencorrs(iproj, iptcl, xy, dists_inpl(:,ithr))
                dists_inpl(:,ithr) = eulprob_dist_switch(dists_inpl(:,ithr))
                irot = angle_sampling(dists_inpl(:,ithr), dists_inpl_sorted(:,ithr), inds_sorted(:,ithr), inpl_athres)
                self%shift_tab(ishift,i)%x      = xy(1)
                self%shift_tab(ishift,i)%y      = xy(2)
                self%shift_tab(ishift,i)%inpl   = irot
                self%shift_tab(ishift,i)%has_sh = .true.
                self%shift_tab(ishift,i)%dist   = dists_inpl(irot,ithr)
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_shift_tab

    subroutine shift_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, ishift, dist_inds(NSHIFTS), stab_inds(self%nptcls, NSHIFTS), assigned_ishift, assigned_ptcl
        real    :: sum_dist_all, min_dist, max_dist, sorted_tab(self%nptcls, NSHIFTS), shift_dist(NSHIFTS)
        logical :: ptcl_avail(self%nptcls)
        ! normalization
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%shift_tab(:,i)%dist)
            if( sum_dist_all < TINY )then
                self%shift_tab(:,i)%dist = 0.
            else
                self%shift_tab(:,i)%dist = self%shift_tab(:,i)%dist / sum_dist_all
            endif
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        max_dist = 0.
        min_dist = huge(min_dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)&
        !$omp reduction(min:min_dist) reduction(max:max_dist)
        do i = 1, self%nptcls
            max_dist = max(max_dist, maxval(self%shift_tab(:,i)%dist, dim=1))
            min_dist = min(min_dist, minval(self%shift_tab(:,i)%dist, dim=1))
        enddo
        !$omp end parallel do
        if( (max_dist - min_dist) < TINY )then
            self%shift_tab%dist = 0.
        else
            self%shift_tab%dist = (self%shift_tab%dist - min_dist) / (max_dist - min_dist)
        endif
        ! assigning shift points to iptcl:
        ! sorting each columns
        sorted_tab = transpose(self%shift_tab%dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ishift,i)
        do ishift = 1, NSHIFTS
            stab_inds(:,ishift) = (/(i,i=1,self%nptcls)/)
            call hpsort(sorted_tab(:,ishift), stab_inds(:,ishift))
        enddo
        !$omp end parallel do
        ! first row is the current best state distribution
        dist_inds  = 1
        shift_dist = sorted_tab(1,:)
        ptcl_avail = .true.
        do while( any(ptcl_avail) )
            ! choose next ishift to assign !!! SHOULD DO PROBABILISTIC SAMPLING HERE
            assigned_ishift = minloc(shift_dist, dim=1)
            assigned_ptcl   = stab_inds(dist_inds(assigned_ishift), assigned_ishift)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%shift_tab(assigned_ishift,assigned_ptcl)
            ! update the shift_dist and dist_inds
            do ishift = 1, NSHIFTS
                do while( dist_inds(ishift) < self%nptcls .and. .not.(ptcl_avail(stab_inds(dist_inds(ishift), ishift))))
                    dist_inds(ishift)  = dist_inds(ishift) + 1
                    shift_dist(ishift) = sorted_tab(dist_inds(ishift), ishift)
                enddo
            enddo
        enddo
    end subroutine shift_assign

    ! projection normalization (same energy) of the 3D loc_tab (for each state)
    ! [0,1] normalization for each state
    subroutine proj_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, istate
        do istate = 1,self%nstates
            ! normalize so prob of each ptcl is between [0,1] for all projs
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
            do i = 1, self%nptcls
                sum_dist_all = sum(self%loc_tab(:,i,istate)%dist)
                if( sum_dist_all < TINY )then
                    self%loc_tab(:,i,istate)%dist = 0.
                else
                    self%loc_tab(:,i,istate)%dist = self%loc_tab(:,i,istate)%dist / sum_dist_all
                endif
            enddo
            !$omp end parallel do
            ! min/max normalization to obtain values between 0 and 1
            !$omp parallel workshare proc_bind(close)
            min_dist = minval(self%loc_tab(:,:,istate)%dist)
            max_dist = maxval(self%loc_tab(:,:,istate)%dist)
            !$omp end parallel workshare
            if( (max_dist - min_dist) < TINY )then
                self%loc_tab(:,:,istate)%dist = 0.
            else
                self%loc_tab(:,:,istate)%dist = (self%loc_tab(:,:,istate)%dist - min_dist) / (max_dist - min_dist)
            endif
        enddo
    end subroutine proj_normalize

    ! (for each state) ptcl -> proj assignment using the global normalized dist value table
    subroutine proj_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, iproj, istate, assigned_iproj, assigned_ptcl, proj_dist_inds(params_glob%nspace, self%nstates),&
                    &stab_inds(self%nptcls, params_glob%nspace, self%nstates), inds_sorted(params_glob%nspace, self%nstates)
        real    :: sorted_tab(self%nptcls, params_glob%nspace, self%nstates), projs_athres,&
                    &proj_dist(params_glob%nspace, self%nstates), dists_sorted(params_glob%nspace, self%nstates)
        logical :: ptcl_avail(self%nptcls, self%nstates)
        ! sorting each columns
        do istate = 1, self%nstates
            sorted_tab(:,:,istate) = transpose(self%loc_tab(:,:,istate)%dist)
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iproj,i)
            do iproj = 1, params_glob%nspace
                stab_inds(:,iproj,istate) = (/(i,i=1,self%nptcls)/)
                call hpsort(sorted_tab(:,iproj,istate), stab_inds(:,iproj,istate))
            enddo
            !$omp end parallel do
        enddo
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(istate,projs_athres,assigned_iproj,assigned_ptcl,iproj)
        do istate = 1, self%nstates
            projs_athres             = calc_athres('dist', state=istate)
            ! first row is the current best proj distribution
            proj_dist_inds(:,istate) = 1
            proj_dist(     :,istate) = sorted_tab(1,:,istate)
            ptcl_avail(    :,istate) = .true.
            do while( any(ptcl_avail(:,istate)) )
                ! sampling the proj distribution to choose next iproj to assign
                assigned_iproj = angle_sampling(proj_dist(:,istate), dists_sorted(:,istate), inds_sorted(:,istate), projs_athres)
                assigned_ptcl  = stab_inds(proj_dist_inds(assigned_iproj,istate), assigned_iproj, istate)
                ptcl_avail(assigned_ptcl,istate)     = .false.
                self%state_tab(istate,assigned_ptcl) = self%loc_tab(assigned_iproj,assigned_ptcl,istate)
                ! update the proj_dist and proj_dist_inds
                do iproj = 1, params_glob%nspace
                    do while( proj_dist_inds(iproj,istate) < self%nptcls .and. &
                                .not.(ptcl_avail(stab_inds(proj_dist_inds(iproj,istate),iproj,istate),istate)))
                        proj_dist_inds(iproj,istate) = proj_dist_inds(iproj,istate) + 1
                        proj_dist(     iproj,istate) = sorted_tab(proj_dist_inds(iproj,istate), iproj, istate)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine proj_assign

    ! state normalization (same energy) of the state_tab
    ! [0,1] normalization of the whole table
    subroutine state_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i
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
        if( (max_dist - min_dist) < TINY )then
            self%state_tab%dist = 0.
        else
            self%state_tab%dist = (self%state_tab%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine state_normalize

    ! ptcl -> state assignment
    subroutine state_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, istate, assigned_istate, assigned_ptcl, state_dist_inds(self%nstates),&
                    &stab_inds(self%nptcls, self%nstates)
        real    :: sorted_tab(self%nptcls, self%nstates), state_dist(self%nstates)
        logical :: ptcl_avail(self%nptcls)
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

    subroutine trim_tab( self, os )
        class(eul_prob_tab), intent(inout) :: self
        class(oris),         intent(in)    :: os
        integer, allocatable :: states(:), sampled(:), inds(:)
        logical, allocatable :: mask(:)
        integer :: i, n, ntot
        ntot = os%get_noris()
        allocate(sampled(ntot), source=nint(os%get_all('sampled')))
        allocate(states(ntot), source=nint(os%get_all('state')))
        mask = sampled > 0 .and. states > 0
        n    = count(mask)
        deallocate(states,sampled)
        if( n == self%nptcls )return
        inds           = pack((/(i,i=1,size(self%pinds))/), mask=mask(self%pinds(:)))
        self%loc_tab   = self%loc_tab(:,inds(:),:)
        self%assgn_map = self%assgn_map(inds(:))
        self%state_tab = self%state_tab(:,inds(:))
        self%pinds     = self%pinds(inds(:))
        self%nptcls    = n
        deallocate(mask,inds)
    end subroutine trim_tab

    ! FILE IO

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_shift( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        character(len=*),    intent(in) :: binfname
        integer :: funit, addr, io_stat
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)                     self%nptcls
        write(unit=funit,pos=sizeof(self%nptcls)+1) self%shift_tab
        call fclose(funit)
    end subroutine write_shift

    subroutine read_shift( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: shift_tab_glob(:,:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('read_tab_to_glob; read_shift; file: '//trim(binfname), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(shift_tab_glob(NSHIFTS,nptcls_glob))
        read(unit=funit,pos=headsz + 1) shift_tab_glob
        call fclose(funit)
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( shift_tab_glob(1,i_glob)%pind == self%shift_tab(1,i_loc)%pind ) self%shift_tab(:,i_loc) = shift_tab_glob(:,i_glob)
            end do
        end do
        !$omp end parallel do
    end subroutine read_shift

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_tab( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        character(len=*),    intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = params_glob%nspace
        file_header(2) = self%nptcls
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit,pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global reg object's dist value table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: mat_loc(:,:,:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, nprojs_loc, i_loc, i_glob
        if( file_exists(trim(binfname)) )then
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob; file: '//trim(binfname), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        ! reading header and the nprojs/nptcls in this partition file
        read(unit=funit,pos=1) file_header
        nprojs_loc = file_header(1)
        nptcls_loc = file_header(2)
        allocate(mat_loc(nprojs_loc, nptcls_loc, self%nstates))
        if( nprojs_loc .ne. params_glob%nspace ) THROW_HARD( 'npsace should be the same as nprojs in this partition file!' )
        ! read partition information
        addr = sizeof(file_header) + 1
        read(unit=funit,pos=addr) mat_loc
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_glob = 1, self%nptcls
            do i_loc = 1, nptcls_loc
                if( mat_loc(1,i_loc,1)%pind == self%loc_tab(1,i_glob,1)%pind )then
                    self%loc_tab(:,i_glob,:) = mat_loc(:,i_loc,:)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine read_tab_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        character(len=*),    intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        character(len=*),    intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: assgn_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('read_tab_to_glob; read_assignment; file: '//trim(binfname), io_stat)
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
        if( allocated(self%loc_tab)   ) deallocate(self%loc_tab)
        if( allocated(self%state_tab) ) deallocate(self%state_tab)
        if( allocated(self%assgn_map) ) deallocate(self%assgn_map)
        if( allocated(self%shift_tab) ) deallocate(self%shift_tab)
        if( allocated(self%pinds)     ) deallocate(self%pinds)
    end subroutine kill

    ! PUBLIC UTILITITES

    subroutine calc_num2sample( num_all, field_str, num_smpl, state)
        integer,           intent(in)  :: num_all
        character(len=*),  intent(in)  :: field_str
        integer,           intent(out) :: num_smpl
        integer, optional, intent(in)  :: state
        real :: athres
        athres   = calc_athres( field_str, state )
        num_smpl = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_num2sample

    function calc_athres( field_str, state ) result( athres )
        character(len=*),  intent(in) :: field_str
        integer, optional, intent(in) :: state
        real, allocatable :: vals(:)
        real :: athres, dist_thres
        if( params_glob%l_batchfrac )then
            vals = build_glob%spproj_field%get_all_sampled(trim(field_str), state=state, lowerbound=0.5)
        else
            vals = build_glob%spproj_field%get_all_sampled(trim(field_str), state=state)
        endif
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
