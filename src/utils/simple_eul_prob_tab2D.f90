! class probability table eul_prob_tab2D, used in abinitio2D
module simple_eul_prob_tab2D
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_polarft_corrcalc,  only: pftcc_glob
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad
use simple_eul_prob_tab,      only: eulprob_dist_switch , eulprob_corr_switch
use simple_decay_funs,        only: extremal_decay2D
implicit none

public :: eul_prob_tab2D, squared_sampling
private
#include "simple_local_flags.inc"

type :: eul_prob_tab2D
    type(ptcl_ref), allocatable :: loc_tab(:,:)   !< 2D search table (ncls x nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)   !< assignment map  (nptcls)
    integer,        allocatable :: pinds(:)       !< particle indices
    integer,        allocatable :: clsinds(:)     !< non-empty class indices
    logical,        allocatable :: populated(:)   !< nonempty classes mask
    integer                     :: nptcls         !< size of pinds array
    integer                     :: neffcls        !< # of non-empty classes
    integer                     :: ncls           !< # of classes
    contains
    ! CONSTRUCTOR
    procedure :: new
    ! TABLE FILLING
    procedure :: fill_table_greedy_inpl
    procedure :: fill_table_stoch_inpl
    ! ASSIGNMENT FROM TABLE
    procedure :: normalize_table
    procedure :: assign_cls_stoch
    procedure :: assign_cls_greedy
    ! I/O
    procedure :: write_table
    procedure :: read_table_to_glob
    procedure :: write_assignment
    procedure :: read_assignment
    ! DESTRUCTOR
    procedure :: kill
end type eul_prob_tab2D

contains

    ! CONSTRUCTORS

    subroutine new( self, pinds )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: pinds(:)
        integer, allocatable :: pops(:)
        integer :: i, iptcl, icls
        real    :: x
        call self%kill
        call seed_rnd
        self%nptcls = size(pinds)
        self%ncls   = params_glob%ncls
        allocate(self%loc_tab(self%ncls,self%nptcls), self%assgn_map(self%nptcls),self%pinds(self%nptcls))
        ! Particles
        !$omp parallel do default(shared) private(i,iptcl,icls) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            self%pinds(i) = pinds(i)
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 1
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            do icls = 1,self%ncls
                self%loc_tab(icls,i)%pind   = iptcl
                self%loc_tab(icls,i)%istate = 1
                self%loc_tab(icls,i)%iproj  = icls
                self%loc_tab(icls,i)%inpl   = 0
                self%loc_tab(icls,i)%dist   = huge(x)
                self%loc_tab(icls,i)%x      = 0.
                self%loc_tab(icls,i)%y      = 0.
                self%loc_tab(icls,i)%has_sh = .false.
            end do
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
                    where( pops < 2 ) pops = 0 ! ignoring classes with one particle
                endif
            else
                allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
            endif
        endif
        if( all(pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        self%populated = pops > 0
        self%neffcls   = count(self%populated)
        self%clsinds   = pack((/(i,i=1,self%ncls)/),mask=self%populated)
    end subroutine new

    ! TABLE

    ! Fill the probability table taking the best in-plane angle
    subroutine fill_table_greedy_inpl( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftcc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real    :: scores(pftcc_glob%get_nrots())
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
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
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
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
                    irot = maxloc(scores, dim=1)
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(scores(irot))
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_table_greedy_inpl

    ! Fill the probability table choosing an in-plane angle stochastically
    subroutine fill_table_stoch_inpl( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftcc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real    :: scores(pftcc_glob%get_nrots()),lims(2,2),lims_init(2,2),cxy(3),score,neigh_frac
        integer :: vec(pftcc_glob%get_nrots())
        integer :: nrots, i, j, iptcl, ithr, irot, icls, jrot, rank, ninpl_stoch
        nrots = pftcc_glob%get_nrots()
        ! size of stochastic neighborhood (# of in-plane angles to draw from)
        neigh_frac  = extremal_decay2D(params_glob%extr_iter, params_glob%extr_lim)
        ninpl_stoch = nint(real(pftcc_glob%get_nrots())*(neigh_frac/3.))
        ninpl_stoch = max(2,min(ninpl_stoch,pftcc_glob%get_nrots()))
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
            !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2)&
            !$omp private(i,iptcl,ithr,icls,j,irot,jrot,score,scores,vec,cxy,rank)
            do i = 1,self%nptcls
                do j = 1,self%neffcls
                    iptcl = self%pinds(i)
                    icls  = self%clsinds(j)
                    ithr  = omp_get_thread_num() + 1
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
                    call squared_sampling(nrots, scores, vec, ninpl_stoch, irot, rank, score)
                    jrot = irot
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        cxy(2:3) = 0.
                    else
                        jrot  = irot
                        score = cxy(1)
                    endif
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl   = jrot
                    self%loc_tab(icls,i)%x      = cxy(2)
                    self%loc_tab(icls,i)%y      = cxy(3)
                    self%loc_tab(icls,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,j,iptcl,icls,irot,scores,vec,rank,score)&
            !$omp proc_bind(close) schedule(static) collapse(2)
            do i = 1, self%nptcls
                do j = 1, self%neffcls
                    iptcl = self%pinds(i)
                    icls  = self%clsinds(j)
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
                    call squared_sampling(nrots, scores, vec, ninpl_stoch, irot, rank, score)
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(score)
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_table_stoch_inpl

    ! Normalization of loc_tab of each class such that prob of each ptcl is in [0,1]
    subroutine normalize_table( self )
        class(eul_prob_tab2D), intent(inout) :: self
        real(dp) :: sumdist
        real     :: mindist, maxdist
        integer :: i
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

    ! Assigns best class to all particles
    subroutine assign_cls_greedy( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: i, icls
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,icls)
        do i = 1,self%nptcls
            icls = minloc(self%loc_tab(:,i)%dist, dim=1, mask=self%populated)
            self%assgn_map(i) = self%loc_tab(icls,i)
        enddo
        !$omp end parallel do
    end subroutine assign_cls_greedy

    ! Assign class to all particles stochastically with self-subtraction
    subroutine assign_cls_stoch( self, os )
        class(eul_prob_tab2D), intent(inout) :: self
        class(oris),           intent(in)    :: os
        real,    allocatable :: dists(:)
        integer, allocatable :: inds(:), prev_cls(:)
        real    :: tmp(self%ncls), pdists(self%neffcls), score, t, neigh_frac
        integer :: pops(self%ncls), vec(self%neffcls), i, icls, ind, rank, ncls_stoch, prev_ref
        logical :: l_self_subtr
        ! previous class
        prev_cls = os%get_all_asint('class')
        ! size of stochastic neighborhood (# of classes to draw from)
        l_self_subtr = params_glob%extr_iter <= params_glob%extr_lim ! self-subtraction
        neigh_frac   = extremal_decay2D(params_glob%extr_iter, params_glob%extr_lim)
        ncls_stoch   = nint(real(self%ncls)*(neigh_frac/3.))
        ncls_stoch   = max(2,min(ncls_stoch,self%ncls))
        ! select class stochastically
        pops = 0
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(i,pdists,vec,ind,rank,score,icls,tmp,prev_ref,inds,dists,t)
        !$omp do schedule(static) reduction(+:pops)
        do i = 1,self%nptcls
            tmp = self%loc_tab(:,i)%dist
            if( l_self_subtr )then
                prev_ref = prev_cls(self%pinds(i))
                if( prev_ref > 0 ) tmp(prev_ref) = 1.           ! previous class set to maximal distance
            endif
            pdists = pack(tmp,mask=self%populated)              ! distances for non-empty classes
            pdists = eulprob_corr_switch(pdists)                ! distances back to score
            call squared_sampling(self%neffcls, pdists, vec, ncls_stoch, ind, rank, score) ! stochastic sampling
            icls = self%clsinds(ind)                            ! class ID
            self%assgn_map(i) = self%loc_tab(icls,i)            ! updates assignement
            pops(icls) = pops(icls) + 1                         ! updates population
        enddo
        !$omp end do
        ! taking the best maxpop particles
        !$omp do schedule(static)
        do i = 1,self%neffcls
            icls = self%clsinds(i)
            if( pops(icls) > params_glob%maxpop )then
                inds  = pack((/(ind,ind=1,self%nptcls)/), mask=self%assgn_map(:)%iproj==icls)
                dists = self%assgn_map(inds(:))%dist
                call hpsort(dists)
                t = dists(params_glob%maxpop)           ! threshold within the class
                where( self%assgn_map(inds(:))%dist > t )
                    self%assgn_map(inds(:))%istate = 0  ! this will set the restoration weight to ZERO
                else where
                    self%assgn_map(inds(:))%istate = 1
                end where
            endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine assign_cls_stoch

    ! FILE I/O

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_table( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        character(len=*),      intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = self%ncls
        file_header(2) = self%nptcls
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit,pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_table

    ! read the partition-wise dist value binary file to global object's table
    subroutine read_table_to_glob( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        character(len=*),      intent(in)    :: binfname
        type(ptcl_ref), allocatable :: mat_loc(:,:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, ncls_loc, i_loc, i_glob
        if( file_exists(trim(binfname)) )then
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab2D; read_table_to_glob; file: '//trim(binfname), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        ! reading header and the ncls/nptcls in this partition file
        read(unit=funit,pos=1) file_header
        ncls_loc = file_header(1)
        nptcls_loc = file_header(2)
        if( ncls_loc .ne. params_glob%ncls ) THROW_HARD( 'NCLS should be the same in this partition file!' )
        allocate(mat_loc(ncls_loc, nptcls_loc))
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
    end subroutine read_table_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        character(len=*),      intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        character(len=*),      intent(in)    :: binfname
        type(ptcl_ref), allocatable :: assgn_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab2D; read_assignment; file: '//trim(binfname), io_stat)
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
            &self%populated, self%clsinds)
        endif
    end subroutine kill

    ! PUBLIC FUNCTION

    subroutine squared_sampling( n, corrs, order, nb, ind, rank, cc )
        integer, intent(in)    :: n, nb
        real,    intent(inout) :: corrs(n), cc
        integer, intent(inout) :: order(n), ind, rank
        integer, parameter :: P=2
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
        r   = 1.-r**P
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
    end subroutine squared_sampling

end module simple_eul_prob_tab2D
