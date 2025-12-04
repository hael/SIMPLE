! class probability table eul_prob_tab2D, used in abinitio2D
module simple_eul_prob_tab2D
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_polarft_calc,  only: pftc_glob
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_eul_prob_tab,      only: eulprob_dist_switch, eulprob_corr_switch
use simple_decay_funs,        only: extremal_decay2D
implicit none

public :: eul_prob_tab2D, squared_sampling, power_sampling, neighfrac2nsmpl
private
#include "simple_local_flags.inc"

type :: eul_prob_tab2D
    type(ptcl_rec), allocatable :: loc_tab(:,:)   !< 2D search table (ncls x nptcls)
    type(ptcl_rec), allocatable :: assgn_map(:)   !< assignment map  (nptcls)
    integer,        allocatable :: pinds(:)       !< particle indices
    integer,        allocatable :: clsinds(:)     !< non-empty class indices
    integer,        allocatable :: indcls(:)      !< reverse lookup non-empty class indices
    logical,        allocatable :: populated(:)   !< nonempty classes mask
    integer                     :: nptcls         !< size of pinds array
    integer                     :: neffcls        !< # of non-empty classes
    integer                     :: ncls           !< # of classes
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! TABLE FILLING
    procedure          :: fill_table_greedy
    procedure          :: fill_table_smpl
    procedure          :: fill_table_smpl_stream
    ! ASSIGNMENT FROM TABLE
    procedure          :: normalize_table
    procedure, private :: select_top_ptcls
    procedure, private :: prevcls_withdrawal
    procedure          :: assign_greedy
    procedure          :: assign_shc
    procedure          :: assign_smpl
    procedure          :: assign_smpl_stream
    procedure          :: assign_prob
    ! I/O
    procedure          :: write_table
    procedure          :: read_table_parts_to_glob
    procedure          :: write_assignment
    procedure          :: read_assignment
    ! DESTRUCTOR
    procedure          :: kill
end type eul_prob_tab2D

! particle class pair record
type ptcl_rec
    integer :: ptcl = 0, cls = 0, inpl = 0
    real    :: dist = 0., x = 0., y = 0.
    logical :: has_sh = .false., incl=.true.
  contains
    procedure :: set
end type ptcl_rec

contains

    ! CONSTRUCTORS

    elemental subroutine set( self, ptcl, cls )
        class(ptcl_rec), intent(inout) :: self
        integer, intent(in) :: ptcl, cls
        self%ptcl   = ptcl
        self%cls    = cls
        self%inpl   = 0
        self%dist   = huge(self%dist)
        self%x      = 0.
        self%y      = 0.
        self%has_sh = .false.
        self%incl   = .true.
    end subroutine set

    subroutine new( self, pinds )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: pinds(:)
        integer, parameter   :: MIN_POP = 2   ! ignoring classes with one particle
        integer, allocatable :: pops(:)
        integer :: i, iptcl, icls
        call self%kill
        call seed_rnd
        self%nptcls = size(pinds)
        self%ncls   = params_glob%ncls
        allocate(self%loc_tab(self%ncls,self%nptcls), self%assgn_map(self%nptcls),self%pinds(self%nptcls))
        ! Particles
        !$omp parallel do default(shared) private(i,iptcl,icls) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            iptcl = pinds(i)
            self%pinds(i) = iptcl
            call self%assgn_map(i)%set(iptcl, 0)
            call self%loc_tab(:,i)%set(iptcl,(/(icls,icls=1,self%ncls)/))
        end do
        !$omp end parallel do
        ! Classes (similar to strategy2D_alloc)
        if( build_glob%spproj%os_cls2D%get_noris() == 0 )then
            if( build_glob%spproj_field%isthere('class') )then
                call build_glob%spproj%os_ptcl2D%get_pops(pops, 'class', maxn=self%ncls)
            else
                allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
            endif
        else
            if( build_glob%spproj_field%isthere('class') )then
                if( build_glob%spproj%os_cls2D%get_noris() /= self%ncls )then
                    ! to be able to restart after having run cleanup with fewer classes
                    allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
                else
                    pops = nint(build_glob%spproj%os_cls2D%get_all('pop'))
                    where( pops < MIN_POP ) pops = 0 ! ignoring classes with one particle
                endif
            else
                allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
            endif
        endif
        if( all(pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        ! non-empty classses
        self%populated = pops > 0
        self%neffcls   = count(self%populated)
        allocate(self%clsinds(self%neffcls),self%indcls(self%ncls),source=0)
        i = 0
        do icls = 1,self%ncls
            if( self%populated(icls) )then
                i = i + 1
                self%clsinds(i)   = icls
                self%indcls(icls) = i
            endif
        enddo
    end subroutine new

    ! TABLE

    ! Fill the probability table taking the best in-plane angle
    subroutine fill_table_greedy( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real    :: scores(pftc_glob%get_nrots())
        real    :: lims(2,2), lims_init(2,2), cxy(3), best_score
        integer :: i, j, iptcl, ithr, irot, icls, best_rot
        if( params_glob%l_doshift )then
            ! search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                    &maxits=params_glob%maxits_sh, opt_angle=.true.)
            end do
            ! search
            !$omp parallel do default(shared) private(i,iptcl,ithr,icls,irot,best_rot,best_score,scores,cxy)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                do icls = 1, self%ncls
                    if( .not.self%populated(icls) ) cycle
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    irot     = maxloc(scores, dim=1)
                    best_rot = irot
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        best_score = scores(best_rot)
                        cxy(2:3)   = 0.
                    else
                        best_rot   = irot
                        best_score = cxy(1)
                    endif
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(best_score)
                    self%loc_tab(icls,i)%inpl   = best_rot
                    self%loc_tab(icls,i)%x      = cxy(2)
                    self%loc_tab(icls,i)%y      = cxy(3)
                    self%loc_tab(icls,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,j,iptcl,icls,irot,scores)&
            !$omp proc_bind(close) schedule(static) collapse(2)
            do i = 1, self%nptcls
                do j = 1, self%neffcls
                    iptcl = self%pinds(i)
                    icls  = self%clsinds(j)
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    irot = maxloc(scores, dim=1)
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(scores(irot))
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_table_greedy

    ! Fill the probability table choosing an in-plane angle stochastically
    ! after shift search against a per-particle subset of top ranking classes
    subroutine fill_table_smpl( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real,       allocatable :: sorted_scores(:)
        integer,    allocatable :: sorted_inds(:)
        real    :: scores(pftc_glob%get_nrots()),lims(2,2),lims_init(2,2),cxy(3),P,score,neigh_frac
        integer :: vec(pftc_glob%get_nrots())
        integer :: nrots, i, j, iptcl, ithr, irot, icls, jrot, rank, ninpl_smpl, ncls_smpl
        nrots = pftc_glob%get_nrots()
        ! power of sampling distribution
        P = EXTR_POWER
        if( params_glob%extr_iter > params_glob%extr_lim ) P = POST_EXTR_POWER
        ! size of stochastic neighborhood (# of in-plane angles to draw from)
        neigh_frac = extremal_decay2D(params_glob%extr_iter, params_glob%extr_lim)
        ninpl_smpl = neighfrac2nsmpl(neigh_frac, nrots)
        ncls_smpl  = neighfrac2nsmpl(neigh_frac, self%neffcls)
        ! Fork
        if( params_glob%l_doshift )then
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init,&
                    &shbarrier=params_glob%shbarrier, maxits=params_glob%maxits_sh,&
                    &opt_angle=.true., coarse_init=.false.)
            end do
            allocate(sorted_scores(self%neffcls),sorted_inds(self%neffcls))
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i,iptcl,ithr,icls,j,irot,jrot,score,scores,sorted_scores,sorted_inds,vec,cxy,rank)
            do i = 1,self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! exhaustive evaluation without shifts
                do j = 1,self%neffcls
                    icls  = self%clsinds(j)
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    call power_sampling(P, nrots, scores, vec, ninpl_smpl, irot, rank, score)
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl   = irot
                    self%loc_tab(icls,i)%x      = 0.
                    self%loc_tab(icls,i)%y      = 0.
                    self%loc_tab(icls,i)%has_sh = .true.
                enddo
                ! subset
                sorted_scores = self%loc_tab(self%clsinds(:),i)%dist
                sorted_inds   = (/(j,j=1,self%neffcls)/)
                call hpsort(sorted_scores, sorted_inds)
                ! shift search of top-ranking subset
                do j = 1,ncls_smpl
                    icls = self%clsinds(sorted_inds(j))
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    irot = self%loc_tab(icls,i)%inpl
                    jrot = irot
                    cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                    if( irot == 0 )then
                        irot = jrot
                        cxy  = [real(pftc_glob%gen_corr_for_rot_8(icls, iptcl, irot)), 0.,0.]
                    endif
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(cxy(1))
                    self%loc_tab(icls,i)%inpl = irot
                    self%loc_tab(icls,i)%x    = cxy(2)
                    self%loc_tab(icls,i)%y    = cxy(3)
                enddo
            enddo
            !$omp end parallel do
            deallocate(sorted_scores,sorted_inds)
        else
            !$omp parallel do default(shared) private(i,j,iptcl,icls,irot,scores,vec,rank,score)&
            !$omp proc_bind(close) schedule(static) collapse(2)
            do i = 1, self%nptcls
                do j = 1, self%neffcls
                    iptcl = self%pinds(i)
                    icls  = self%clsinds(j)
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    call power_sampling(P, nrots, scores, vec, ninpl_smpl, irot, rank, score)
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_table_smpl

    ! Fill the probability table choosing an in-plane angle stochastically
    ! after shift search against a per-particle subset of top ranking classes
    ! New particles are searched in a greedy fashion
    subroutine fill_table_smpl_stream( self, os )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real,       allocatable :: sorted_scores(:)
        integer,    allocatable :: sorted_inds(:)
        real    :: scores(pftc_glob%get_nrots()),lims(2,2),lims_init(2,2),cxy(3),P,score
        integer :: vec(pftc_glob%get_nrots())
        integer :: nrots, i, j, iptcl, ithr, irot, icls, jrot, rank, ninpl_smpl, ncls_smpl
        logical :: greedy
        nrots = pftc_glob%get_nrots()
        ! power of sampling distribution
        P = EXTR_POWER
        ! size of stochastic neighborhood (# of in-plane angles to draw from)
        ninpl_smpl = neighfrac2nsmpl(0., nrots)
        ncls_smpl  = neighfrac2nsmpl(0., self%neffcls)
        ! Fork
        if( params_glob%l_doshift )then
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init,&
                    &shbarrier=params_glob%shbarrier, maxits=params_glob%maxits_sh,&
                    &opt_angle=.true., coarse_init=.false.)
            end do
            allocate(sorted_scores(self%neffcls),sorted_inds(self%neffcls))
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i,iptcl,ithr,icls,j,irot,jrot,score,scores,sorted_scores,sorted_inds,vec,cxy,rank,greedy)
            do i = 1,self%nptcls
                iptcl  = self%pinds(i)
                ithr   = omp_get_thread_num() + 1
                greedy = (os%get_updatecnt(iptcl)==1) .or. (.not.os%has_been_searched(iptcl))
                ! exhaustive evaluation without shifts
                do j = 1,self%neffcls
                    icls  = self%clsinds(j)
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    if( greedy )then
                        ! greedy in-plane
                        irot  = maxloc(scores,dim=1)
                        score = scores(irot)
                    else
                        ! stochastic in-plane
                        call power_sampling(P, nrots, scores, vec, ninpl_smpl, irot, rank, score)
                    endif
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl   = irot
                    self%loc_tab(icls,i)%x      = 0.
                    self%loc_tab(icls,i)%y      = 0.
                    self%loc_tab(icls,i)%has_sh = .true.
                enddo
                ! subset
                sorted_scores = self%loc_tab(self%clsinds(:),i)%dist
                sorted_inds   = (/(j,j=1,self%neffcls)/)
                call hpsort(sorted_scores, sorted_inds)
                ! shift search of top-ranking subset
                do j = 1,ncls_smpl
                    icls = self%clsinds(sorted_inds(j))
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    irot = self%loc_tab(icls,i)%inpl
                    jrot = irot
                    cxy  = grad_shsrch_obj(ithr)%minimize(irot=irot)
                    if( irot == 0 )then
                        irot = jrot
                        cxy  = [real(pftc_glob%gen_corr_for_rot_8(icls, iptcl, irot)), 0.,0.]
                    endif
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(cxy(1))
                    self%loc_tab(icls,i)%inpl = irot
                    self%loc_tab(icls,i)%x    = cxy(2)
                    self%loc_tab(icls,i)%y    = cxy(3)
                enddo
            enddo
            !$omp end parallel do
            deallocate(sorted_scores,sorted_inds)
        else
            !$omp parallel do default(shared) private(i,j,iptcl,icls,irot,scores,vec,rank,score,greedy)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl  = self%pinds(i)
                greedy = (os%get_updatecnt(iptcl)==1).or.(.not.os%has_been_searched(iptcl))
                do j = 1, self%neffcls
                    icls  = self%clsinds(j)
                    call pftc_glob%gen_objfun_vals(icls, iptcl, [0.,0.], scores)
                    if( greedy )then
                        ! new particle
                        irot  = maxloc(scores,dim=1)
                        score = scores(irot)
                    else
                        call power_sampling(P, nrots, scores, vec, ninpl_smpl, irot, rank, score)
                    endif
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_table_smpl_stream

    ! Normalization of loc_tab of each class such that prob of each ptcl is in [0,1]
    subroutine normalize_table( self )
        class(eul_prob_tab2D), intent(inout) :: self
        real(dp) :: sumdist
        real     :: mindist, maxdist
        integer  :: i
        ! normalize each ptcl for all non-empty classes
        mindist = huge(mindist)
        maxdist = -1.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sumdist)&
        !$omp reduction(min:mindist) reduction(max:maxdist)
        do i = 1, self%nptcls
            sumdist = sum(real(self%loc_tab(:,i)%dist,dp), mask=self%populated)
            if( sumdist < DTINY )then
                self%loc_tab(:,i)%dist = 0.0
            else
                where( self%populated )
                    self%loc_tab(:,i)%dist = self%loc_tab(:,i)%dist / real(sumdist)
                else where
                    self%loc_tab(:,i)%dist = huge(mindist)
                end where
            endif
            mindist = min(mindist, minval(self%loc_tab(:,i)%dist, mask=self%populated))
            maxdist = max(maxdist, maxval(self%loc_tab(:,i)%dist, mask=self%populated))
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        if( (maxdist - mindist) < TINY )then
            !$omp parallel workshare proc_bind(close)
            self%loc_tab(:,:)%dist = 0.
            !$omp end parallel workshare
        else
            !$omp parallel workshare proc_bind(close)
            self%loc_tab(:,:)%dist = (self%loc_tab(:,:)%dist - mindist) / (maxdist - mindist)
            !$omp end parallel workshare
        endif
    end subroutine normalize_table

    ! set state/weights to zero for worst ranking particles per class
    subroutine select_top_ptcls( self, pops )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,     optional, intent(in)    :: pops(self%ncls)
        real,    allocatable :: dists(:)
        integer, allocatable :: inds(:)
        real    :: t
        integer :: counts(self%ncls), i, j, icls
        if( present(pops) )then
            counts = pops
        else
            counts = 0
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(icls)
            do icls = 1,self%ncls
                if( self%populated(icls) ) counts(icls) = count(self%assgn_map(:)%cls == icls)
            enddo
            !$omp end parallel do
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,j,icls,inds,dists,t)
        do i = 1,self%neffcls
            icls = self%clsinds(i)
            if( counts(icls) > params_glob%maxpop )then
                inds  = pack((/(j,j=1,self%nptcls)/), mask=self%assgn_map(:)%cls==icls)
                dists = self%assgn_map(inds(:))%dist
                call hpsort(dists)
                t = dists(params_glob%maxpop)               ! threshold within the class
                where( self%assgn_map(inds(:))%dist > t )
                    self%assgn_map(inds(:))%incl = .false.  ! this will set the restoration weight to ZERO
                else where
                    self%assgn_map(inds(:))%incl = .true.
                end where
                counts(icls) = count(self%assgn_map(inds(:))%incl)
            endif
        enddo
        !$omp end parallel do
    end subroutine select_top_ptcls

    ! Self-withdrawal, sets distance to previous class to maximum
    subroutine prevcls_withdrawal( self, os )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        real    :: a
        integer :: i,icls,iptcl
        if( params_glob%extr_iter <= params_glob%extr_lim )then
            a = huge(a)
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,iptcl,icls)
            do i = 1,self%nptcls
                iptcl = self%pinds(i)
                if( os%get_updatecnt(iptcl) == 1 ) cycle ! no withdrawal when first searched
                icls = os%get_class(iptcl)
                if( icls > 0 ) self%loc_tab(icls,i)%dist = a
            enddo
            !$omp end parallel do
        endif
    end subroutine prevcls_withdrawal

    ! Assigns best class to all particles
    subroutine assign_greedy( self, rank_ptcls )
        class(eul_prob_tab2D), intent(inout) :: self
        logical,               intent(in)    :: rank_ptcls
        integer :: pops(self%ncls), i, icls
        pops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,icls) reduction(+:pops)
        do i = 1,self%nptcls
            icls = minloc(self%loc_tab(:,i)%dist, dim=1, mask=self%populated)
            self%assgn_map(i) = self%loc_tab(icls,i)
            pops(icls) = pops(icls) + 1
        enddo
        !$omp end parallel do
        if( rank_ptcls ) call self%select_top_ptcls(pops)
    end subroutine assign_greedy

    ! Assign class to all particles stochastically with self-withdrawal
    subroutine assign_smpl( self, os, rank_ptcls )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        logical,               intent(in)    :: rank_ptcls
        real    :: pdists(self%neffcls), P, score, neigh_frac
        integer :: pops(self%ncls), vec(self%neffcls), i, icls, ind, rank, ncls_smpl
        call self%prevcls_withdrawal(os)
        ! size of stochastic neighborhood (# of classes to draw from)
        neigh_frac = extremal_decay2D(params_glob%extr_iter, params_glob%extr_lim)
        ncls_smpl  = neighfrac2nsmpl(neigh_frac, self%neffcls)
        P = EXTR_POWER
        if( params_glob%extr_iter > params_glob%extr_lim ) P = POST_EXTR_POWER
        ! select class stochastically
        pops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,pdists,vec,ind,rank,score,icls) reduction(+:pops)
        do i = 1,self%nptcls
            pdists = eulprob_corr_switch(self%loc_tab(self%clsinds(:),i)%dist)  ! distances to score
            call power_sampling(P, self%neffcls, pdists, vec, ncls_smpl, ind, rank, score) ! stochastic sampling
            icls = self%clsinds(ind)                                            ! class ID
            self%assgn_map(i) = self%loc_tab(icls,i)                            ! updates assignement
            pops(icls) = pops(icls) + 1                                         ! updates population
        enddo
        !$omp end parallel do
        ! toss worst particles
        if( rank_ptcls ) call self%select_top_ptcls(pops)
    end subroutine assign_smpl

    ! Assign class to old particles stochastically without self-withdrawal
    ! New particles are assigned greedily
    subroutine assign_smpl_stream( self, os, rank_ptcls )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        logical,               intent(in)    :: rank_ptcls
        real    :: pdists(self%neffcls), P, score
        integer :: pops(self%ncls), vec(self%neffcls)
        integer :: iptcl, i, icls, ind, rank, ncls_smpl
        logical :: greedy
        ! power of sampling distribution
        P = EXTR_POWER
        ! size of stochastic neighborhood (# of classes to draw from)
        ncls_smpl = neighfrac2nsmpl(0., self%neffcls)
        pops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,pdists,vec,ind,rank,score,icls,iptcl,greedy) reduction(+:pops)
        do i = 1,self%nptcls
            iptcl  = self%pinds(i)
            greedy = (os%get_updatecnt(iptcl)==1).or.(.not.os%has_been_searched(iptcl))
            if( greedy )then
                pdists = self%loc_tab(self%clsinds(:),i)%dist
                ind    = minloc(pdists,dim=1)                                   ! greedy sampling
            else
                pdists = eulprob_corr_switch(self%loc_tab(self%clsinds(:),i)%dist)             ! distances to score
                call power_sampling(P, self%neffcls, pdists, vec, ncls_smpl, ind, rank, score) ! stochastic sampling
            endif
            icls              = self%clsinds(ind)                               ! class ID
            self%assgn_map(i) = self%loc_tab(icls,i)                            ! updates assignement
            pops(icls)        = pops(icls) + 1                                  ! updates population
        enddo
        !$omp end parallel do
        ! ignore worst particles
        if( rank_ptcls ) call self%select_top_ptcls(pops)
    end subroutine assign_smpl_stream

    ! Assignment based on SHC
    subroutine assign_shc( self, os, rank_ptcls )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        logical,               intent(in)    :: rank_ptcls
        real    :: scores(self%neffcls)
        integer :: pops(self%ncls), i, icls, prev_cls
        pops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(i,prev_cls,icls,scores) reduction(+:pops)
        do i = 1,self%nptcls
            prev_cls = os%get_class(self%pinds(i))
            if( prev_cls <= 0 )then
                icls = minloc(self%loc_tab(self%clsinds(:),i)%dist,dim=1)
            else
                scores = eulprob_corr_switch(self%loc_tab(self%clsinds(:),i)%dist)
                icls   = shcloc(self%neffcls, scores, scores(prev_cls))
                if( icls == 0 ) icls = maxloc(scores,dim=1)
            endif
            icls              = self%clsinds(icls)
            self%assgn_map(i) = self%loc_tab(icls,i)
            pops(icls)        = pops(icls) + 1
        enddo
        !$omp end parallel do
        ! toss worst particles
        if( rank_ptcls ) call self%select_top_ptcls(pops)
    end subroutine assign_shc

    ! Same assignment as 3D (cf eul_prob_tab)
    subroutine assign_prob( self, os, rank_ptcls )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        logical,               intent(in)    :: rank_ptcls
        real,    allocatable :: sorted_tab(:,:)
        integer, allocatable :: stab_inds(:,:)
        integer :: inds_sorted(self%ncls), cls_dist_inds(self%ncls)
        integer :: i,j,icls, cls, ptcl, ncls_smpl, r, iptcl
        real    :: cls_dist(self%ncls), tmp(self%ncls), P, neigh_frac, s
        logical :: ptcl_avail(self%nptcls)
        call self%prevcls_withdrawal(os)
        ! column sorting
        sorted_tab = transpose(self%loc_tab(:,:)%dist)
        allocate(stab_inds(self%nptcls, self%ncls),source=0)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(icls,i,j)
        do i = 1,self%neffcls
            icls = self%clsinds(i)
            stab_inds(:,icls) = (/(j,j=1,self%nptcls)/)
            call hpsort(sorted_tab(:,icls), stab_inds(:,icls))
        enddo
        !$omp end parallel do
        neigh_frac = extremal_decay2D(params_glob%extr_iter, params_glob%extr_lim)
        ncls_smpl = neighfrac2nsmpl(neigh_frac, self%neffcls)
        P = EXTR_POWER
        if( params_glob%extr_iter > params_glob%extr_lim ) P = POST_EXTR_POWER
        ! first row is the current best reference distribution
        cls_dist_inds = 1
        cls_dist      = sorted_tab(1,:)
        ptcl_avail    = .true.
        ! assignment
        do j = 1,self%nptcls
            tmp = eulprob_corr_switch(cls_dist)
            call power_sampling(P, self%ncls, tmp, inds_sorted, ncls_smpl, cls, r, s)
            ptcl                 = stab_inds(cls_dist_inds(cls), cls)
            ptcl_avail(ptcl)     = .false.
            self%assgn_map(ptcl) = self%loc_tab(cls,ptcl)
            do i = 1, self%neffcls
                icls  = self%clsinds(i)
                iptcl = cls_dist_inds(icls)
                do while( iptcl < self%nptcls .and. .not.(ptcl_avail(stab_inds(iptcl,icls))))
                    iptcl               = iptcl + 1
                    cls_dist(icls)      = sorted_tab(iptcl, icls)
                    cls_dist_inds(icls) = iptcl
                enddo
            enddo
        enddo
        deallocate(sorted_tab,stab_inds)
        ! toss worst particles
        if( rank_ptcls ) call self%select_top_ptcls
    end subroutine assign_prob

    ! FILE I/O

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_table( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = self%ncls
        file_header(2) = self%nptcls
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit,pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_table

    ! read the partition records binary files to the global table
    ! particles indices are assumed in increasing order
    subroutine read_table_parts_to_glob( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(string) :: fname
        integer :: funit, addr, io_stat, file_header(2), nptcls, ncls, ipart, istart, iend
        istart = 1
        do ipart = 1,params_glob%nparts
            fname = DIST_FBODY//int2str_pad(ipart,params_glob%numlen)//'.dat'
            if( file_exists(fname) )then
                call fopen(funit,fname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
                call fileiochk('simple_eul_prob_tab2D; read_table_parts_to_glob; file: '//fname%to_char(), io_stat)
            else
                THROW_HARD('Missing file: '//fname%to_char())
            endif
            ! reading header
            read(unit=funit,pos=1) file_header
            ncls   = file_header(1)
            nptcls = file_header(2)
            if( ncls .ne. params_glob%ncls ) THROW_HARD( 'NCLS should be the same in this partition file!' )
            ! read partition information
            iend = istart + nptcls - 1
            if( iend > self%nptcls ) THROW_HARD('More particle records than required!')
            addr = sizeof(file_header) + 1
            read(unit=funit,pos=addr) self%loc_tab(:,istart:iend)
            call fclose(funit)
            istart = istart + nptcls
        enddo
    end subroutine read_table_parts_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(ptcl_rec), allocatable :: assgn_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab2D; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit,pos=headsz + 1) assgn_glob
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( self%assgn_map(i_loc)%ptcl == assgn_glob(i_glob)%ptcl )then
                    self%assgn_map(i_loc) = assgn_glob(i_glob)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        deallocate(assgn_glob)
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        self%neffcls = 0
        self%ncls    = 0
        self%nptcls  = 0
        if( allocated(self%loc_tab)   )then
            deallocate(self%loc_tab, self%assgn_map, self%pinds,&
            &self%populated, self%clsinds, self%indcls)
        endif
    end subroutine kill

    ! PUBLIC FUNCTIONS

    subroutine power_sampling( P, n, corrs, order, nb, ind, rank, cc )
        real,    intent(in)    :: P
        integer, intent(in)    :: n, nb
        real,    intent(inout) :: corrs(n), cc
        integer, intent(inout) :: order(n), ind, rank
        real    :: cdf(nb), r
        integer :: i
        if( nb == 1 )then
            rank = n
            ind  = maxloc(corrs,dim=1)
            cc   = corrs(ind)
            return
        endif
        order = (/(i,i=1,n)/)
        call hpsort(corrs, order)
        cdf = corrs(n-nb+1:n)
        if( all(cdf<TINY) )then
            rank = n
            ind  = order(rank)
            cc   = corrs(rank)
            return
        endif
        where( cdf < TINY ) cdf = 0.
        do i = 2,nb
            cdf(i) = cdf(i) + cdf(i-1)
        enddo
        cdf = cdf / sum(cdf)
        r   = ran3()
        r   = min(1.,max(0.,1.-r**P))
        rank = 0
        do i = 1,nb
            if( cdf(i) > r )then
                rank = i
                exit
            endif
        enddo
        if( rank == 0 ) rank = nb
        rank = n - nb + rank ! rank of selected value
        ind  = order(rank)   ! index
        cc   = corrs(rank)   ! value
    end subroutine power_sampling

    subroutine squared_sampling( n, corrs, order, nb, ind, rank, cc )
        integer, intent(in)    :: n, nb
        real,    intent(inout) :: corrs(n), cc
        integer, intent(inout) :: order(n), ind, rank
        call power_sampling( 2.0, n, corrs, order, nb, ind, rank, cc)
    end subroutine squared_sampling

    pure integer function neighfrac2nsmpl( neighfrac, n )
        real,    intent(in) :: neighfrac
        integer, intent(in) :: n
        neighfrac2nsmpl = nint(real(n)*(neighfrac/3.))
        neighfrac2nsmpl = max(2,min(neighfrac2nsmpl, n))
    end function neighfrac2nsmpl

end module simple_eul_prob_tab2D
