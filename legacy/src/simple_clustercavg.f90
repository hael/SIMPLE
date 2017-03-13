!==Class simple_clustercavg
!
module simple_clustercavg
use simple_defs
use simple_build,        only: build
use simple_oris,         only: oris
use simple_ori,          only: ori
use simple_jiffys,       only: alloc_err, progress
use simple_filehandling, only: get_fileunit
use simple_image,        only: image
use simple_params,       only: params
use simple_shc_cluster,  only: shc_cluster
use simple_cavgppca,     only: cavgppca

implicit none

public :: clustercavg
private

type clustercavg
    private
    class(params), pointer, private   :: pp=>null()
    class(build),  pointer, private   :: pb=>null()
    type(oris)               :: oris             ! private copy of oris of current class
    type(image), allocatable :: ptcls(:)     ! Particles stack
    type(image), allocatable :: ptcls_ref(:)     ! Particles stack
    integer, allocatable     :: ptcls_inds(:)    ! Particles indices to originating oris
    real, allocatable        :: latent_vecs(:,:) ! Particles projection in feature space
    real, allocatable        :: simmat(:,:)      ! Similarity matrix
    integer                  :: n      = 0       ! number of ptcls
    integer                  :: nfeats = 30      ! number of latent vectors for PPCA
    integer                  :: ncls   = 0       ! 
    integer                  :: pop_thresh = 10  ! 
    character(len=STDLEN)    :: stk
    character(len=STDLEN)    :: metric
    logical                  :: usemsk = .false. ! Whether to mask ptcls
    logical                  :: exists = .false.
  contains
    ! CONSTRUCTORS
    procedure          :: new
    procedure          :: exec
    procedure, private :: init_simmat
    procedure, private :: build_simmat_realcc
    procedure, private :: build_simmat_feat
    procedure, private :: build_simmat_pftcc
    procedure, private :: write_simmat
    procedure, private :: make_cavgs
    procedure, private :: mask_ptcls
    procedure, private :: print_pop
    procedure, private :: exec_shc
    procedure, private :: exec_2dshc
    procedure, private :: exec_shc_refine
    procedure, private :: exec_bishc
    procedure, private :: write_stk
    procedure, private :: get_minpop
    procedure, private :: get_maxpopind
    procedure, private :: prepptcls
    procedure, private :: refine_cls
    procedure, private :: denoise_cls
    procedure, private :: cls_mixing
    procedure, private :: get_feature_dist
    procedure, private :: simmat2weights
    procedure          :: get_ncls
    ! DESTRUCTOR
    procedure          :: kill
end type

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new( self, b, p, class, stk, nfeats, msk )
        class(clustercavg), intent(inout)     :: self
        class(build), intent(inout), target  :: b
        class(params), intent(in),target     :: p
        integer, intent(in)                  :: class, nfeats 
        logical, intent(in), optional        :: msk    
        integer                              :: ind, alloc_stat, i
        character(len=STDLEN)                :: stk
        call self%kill
        if( present(msk) )self%usemsk=msk
        self%pp         => p
        self%pb         => b
        self%nfeats     = nfeats
        self%n          = b%a%get_clspop( class )
        self%stk        = stk
        self%ptcls_inds = b%a%get_cls_pinds( class )
        self%ncls       = 1
        !self%pop_thresh = 2*p%minp
        self%pop_thresh = 2*p%minp
        call self%oris%new( self%n )
        ! Stack
        allocate(self%ptcls(self%n),self%ptcls_ref(self%n), stat=alloc_stat)
        call alloc_err('In: new; simple_clustercavg 1', alloc_stat)
        do i=1,self%n
            call self%ptcls(i)%new([self%pp%box,self%pp%box,1], self%pp%smpd)
            call self%ptcls(i)%read( self%stk,self%ptcls_inds(i) )
            call self%ptcls_ref(i)%new([self%pp%box,self%pp%box,1], self%pp%smpd)
            call self%ptcls_ref(i)%copy( self%ptcls(i) )
        enddo
        ! oris
        do i=1,self%n
            ind = self%ptcls_inds(i)
            call self%oris%set_ori(i, b%a%get_ori(ind))
        enddo
        ! init class label
        do i=1,self%n
            call self%oris%set(i, 'class', 1.)
        enddo
        self%exists = .true.
    end subroutine new

    !>   \brief  is the master routine
    subroutine exec( self, metric, clust_mode )
        class(clustercavg), intent(inout)      :: self
        type(cavgppca)                         :: cavg_ppca        ! Probabilistic PCA object
        character(len=*), intent(in), optional :: metric, clust_mode
        integer                                :: i,ind
        character(len=STDLEN) :: fname ! debug
        integer               :: icls  ! debug
        ! Parameters
        if( metric.ne.'realcc' .and. metric.ne.'feat' )then
            stop 'Unknown similarity mode; simple_clustercavg%exec'
        else
            self%metric = metric
        endif
        if( clust_mode.ne.'shc'  .and. clust_mode.ne.'bishc' .and. clust_mode.ne.'2dshc'  ) &
            stop 'Unknown clustering mode; simple_clustercavg%exec'
        ! translation, rotation
        call self%prepptcls
        if( clust_mode.eq.'2dshc' )then
        ! PRIME2D-like CLUSTERING
            call self%exec_2dshc()
        else
            ! PPCA & SHC based clustering
            ! PPCA
            ! debug
            icls = nint(self%oris%get(1,'class'))
            ! write(fname,'(I4.4)')icls
            ! fname = 'stk_'//trim(fname)//'.mrc'
            ! call self%write_stk(fname)
            ! end debug
            if( self%usemsk )then
                call cavg_ppca%new( self%pp,self%pb,self%oris,self%ptcls,self%nfeats,self%metric,.true. )
            else
                call cavg_ppca%new( self%pp,self%pb,self%oris,self%ptcls,self%nfeats,self%metric,.false. )
            endif
            call cavg_ppca%exec
            if( self%metric.eq.'feat' )call cavg_ppca%get_latentvecs( self%latent_vecs )
            call cavg_ppca%kill
            ! debug
            write(fname,'(I4.4)')icls
            fname = 'stk_denoise_'//trim(fname)//'.mrc'
            call self%write_stk(fname)
            !end debug
            ! SIMILARITY MATRIX
            if( self%metric.eq.'realcc' )then
                ! Real space cross-correlation
                call self%build_simmat_realcc
                !call self%write_simmat
            elseif( self%metric.eq.'feat' )then
                ! Features similarity
                call self%build_simmat_feat
            endif
            ! CLUSTERING
            if( clust_mode.eq.'shc'.and.self%pp%refine.eq.'yes' )then
                call self%exec_shc_refine  ! Iteratice refinment+feature extraction+SHC
            elseif( clust_mode.eq.'shc' )then
                call self%exec_shc         ! SHC
            elseif( clust_mode.eq.'bishc' )then
                call self%exec_bishc       ! bisective SHC
            endif
            deallocate( self%simmat )
        endif
        ! copy new classes of builder oris
        do i=1,self%n
            ind = self%ptcls_inds(i)
            call self%pb%a%set( ind,'class',self%oris%get(i,'class') )
        enddo
    end subroutine exec

    subroutine exec_bishc( self )
        use simple_cluster_valid 
        class(clustercavg), intent(inout)      :: self
        type(cluster_valid)                    :: cvalid
        type(shc_cluster)                      :: shcc
        type(ori)                              :: o
        type(oris)                             :: os, os_prev, ios
        integer, allocatable       :: inds(:)
        real, allocatable          :: iS(:,:)
        real           :: sim, sim_prev, elm_ri,overlap, bisim
        integer        :: i,ncls, minpop,ipop,maxpop_thresh,jcls,jptcl,kptcl,jind
        integer        :: icls, incls, nrepeats, iind, iptcl
        logical        :: converged
        call os%new( self%n )
        os        = self%oris
        sim       = calc_sim( os )
        overlap   = 0.
        converged = .false.
        i         = 0
        nrepeats  = 0
        do while( .not.converged )
            sim_prev = sim
            i        = i+1
            ncls     = os%get_ncls()
            if( i>1 )then
                icls = get_worstsim_ind( os )
            else
                icls = 1
            endif
            if( icls == 0 )then
                converged = .true.
                exit
            endif
            ipop = os%get_clspop( icls )
            print *,'selected',icls,ipop
            if( ipop>self%pop_thresh )then
                ! loop oris
                inds = os%get_cls_pinds( icls )
                call ios%new( ipop )
                do iptcl=1,ipop
                    o = os%get_ori( inds(iptcl) )
                    call o%set( 'class',1. )
                    call ios%set_ori( iptcl,o )
                enddo
                ! loop simmat
                allocate( iS(ipop,ipop) )
                do jptcl=1,ipop
                    jind = inds(jptcl)
                    do kptcl=jptcl+1,ipop
                        iS(jptcl,kptcl) = self%simmat( jind,inds(kptcl) )
                        iS(kptcl,jptcl) = iS(jptcl,kptcl)
                    enddo
                enddo
                ! Actual bisective SHC clustering
                call shcc%new( ipop, 2, iS, ios )
                call shcc%shc( .false.,'class',bisim)
                call shcc%kill
                incls = ios%get_ncls()
                if( incls==1 )then
                    ! unsuccessful split
                    nrepeats = nrepeats+1
                    if( nrepeats>2 )converged = .true. ! largest cluster did not get split after 3 attempts
                else
                    ! successful split
                    nrepeats = 0
                    ! updates clusters numbering
                    do jptcl=1,self%n
                        jcls = nint( os%get(jptcl,'class') )
                        if( jcls>icls )call os%set( jptcl,'class',real(jcls+1) )
                    enddo
                    do iptcl=1,ipop
                        iind = inds( iptcl )
                        if ( ios%get(iptcl,'class')==2. ) &
                            call os%set( iind,'class',real(icls+1) )
                    enddo
                    ! global convergergence criterion
                    sim    = calc_sim( os )
                    ncls   = os%get_ncls()
                    call cvalid%new( os,ncls,'class',S=self%simmat )
                    elm_ri = cvalid%ratio_index()
                    call cvalid%kill
                    if( sim_prev/sim>.9 )converged = .true.
                    if( i>50 )             converged = .true.
                    print *,ncls,sim_prev/sim,elm_ri, sim_prev,sim,bisim
                endif
                call ios%kill
                deallocate( inds,iS )
            else
                converged = .true.
            endif
        enddo
        self%oris = os
        call self%print_pop( self%oris )
        self%ncls = self%oris%get_ncls()
        write(*,'(A,I3)')'>>> SHC CLUSTERS FOUND=', self%ncls
        call os%kill
        call os_prev%kill

        contains 
            function calc_sim( localoris )result( s )
                type( oris ), intent(inout) :: localoris
                integer :: i,j,pop,n
                real    :: s,si,icls
                n = localoris%get_noris()
                s = 0.
                do i=1,n
                    si   = 0.
                    pop  = 0
                    icls = localoris%get(i,'class')
                    do j=1,n
                        if( i==j ) cycle
                        if( icls == localoris%get(j,'class') )then
                            si = si+self%simmat(i,j)
                            pop = pop+1
                        endif
                    enddo
                    si = si/real(pop)
                    s  = s+si
                enddo
                s = s/real(n)
            end function

            function get_worstsim_ind( localoris )result( ind )
                type( oris ), intent(inout) :: localoris
                real,allocatable  :: sims(:)
                integer,allocatable :: pops(:)
                integer :: i,j,pop,n,ncls,ind,icls
                real    :: min_sim
                ncls = localoris%get_ncls()
                n    = localoris%get_noris()
                allocate( sims(ncls),pops(ncls) )
                sims = 0.
                pops = 0
                do i=1,n
                    icls       = nint(localoris%get(i,'class'))
                    pops(icls) = pops(icls)+1
                    do j=1,n
                        if( i==j ) cycle
                        if( icls == nint(localoris%get(j,'class')) )then
                            sims(icls) = sims(icls)+self%simmat(i,j)
                        endif
                    enddo
                enddo
                min_sim = huge(min_sim)
                ind     = 0
                do icls=1,ncls
                    pop        = pops(icls)*(pops(icls)-1)
                    sims(icls) = sims(icls)/real(pop)
                    !print *,'i,sim,pop',icls,sims(icls),pops(icls)
                    if( sims(icls)<min_sim .and. pops(icls)>self%pop_thresh )then
                        ind     = icls
                        min_sim = sims(icls)
                    endif
                enddo
                deallocate( sims,pops )
            end function
    end subroutine exec_bishc

    subroutine exec_shc_refine( self )
        !use simple_cluster_valid 
        class(clustercavg), intent(inout)  :: self
        type(shc_cluster)                  :: shcc
        real           :: sim, sim_prev
        integer        :: i,ncls, minpop, icls
        logical        :: converged
        ncls         = 1
        minpop    = self%n
        sim       = 0.
        do i=1,self%n
            sim = sim+sum(self%simmat(i,1:i-1))/real(self%n-1)
            sim = sim+sum(self%simmat(i,i+1:self%n))/real(self%n-1)
        enddo
        sim       = sim/real(self%n)
        converged = .false.
        i         = 0
        print *,'i,cls,sim:',i,ncls,sim
        do while( .not.converged )
            i         = i+1
            sim_prev  = sim
            ncls      = ncls+1
            ! shc clustering
            call shcc%new( self%n, ncls, self%simmat, self%oris )
            call shcc%shc(.false., 'class', sim)
            call shcc%kill
            minpop = self%get_minpop( self%oris )
            ! cluster metrics
            print *,'i,cls,sim,sim_prev/simpop,minpop:',i,ncls,sim,sim_prev/sim,minpop ! debug
            ! Convergence
            ! Minimum population criterion
            if( minpop<=self%pop_thresh )converged=.true.
            ! Similarity improvement criterion
            if( sim_prev/sim>.9 )   converged=.true.
            ! Maximum iterations
            if( ncls>50 )             converged=.true.
            ! refinment
            do icls=1,ncls
                call self%refine_cls( icls )
            enddo
            ! denoising/feature extraction
            call self%prepptcls
            do icls=1,ncls
                call self%denoise_cls( icls )
            enddo
            ! updates similarity matrix
            call self%build_simmat_realcc
            sim = calc_sim( self%oris )
            print *,'i,cls,sim new',i,ncls,sim ! debug
        enddo
        call self%print_pop( self%oris )
        self%ncls = self%oris%get_ncls()
        write(*,'(A,I3)')'>>> SHC CLUSTERS FOUND=', self%ncls

        contains
            function calc_sim( localoris )result( s )
                type( oris ), intent(inout) :: localoris
                integer :: i,j,pop,n
                real    :: s,si,icls
                n = localoris%get_noris()
                s = 0.
                do i=1,n
                    si   = 0.
                    pop  = 0
                    icls = localoris%get(i,'class')
                    do j=1,n
                        if( i==j ) cycle
                        if( icls == nint(localoris%get(j,'class')) )then
                            si  = si+self%simmat(i,j)
                            pop = pop+1
                        endif
                    enddo
                    s  = s+si/real(pop)
                enddo
                s = s/real(n)
            end function
    end subroutine exec_shc_refine

    subroutine exec_2dshc( self )
        use simple_cluster_valid
        ! use simple_shccluster_cavg
        class(clustercavg), intent(inout)  :: self
        ! type(oris) :: os_copy
        ! type(shccluster_cavg)              :: shcc_cavg
        ! real           :: sim, sim_prev
        ! integer        :: step,ncls, minpop
        ! logical        :: converged
        ! ncls   = 1
        ! call self%mask_ptcls
        ! call self%build_simmat_pftcc
        ! sim = calc_sim( self%oris )
        ! converged = .false.
        ! step      = 0
        ! os_copy   = self%oris
        ! print *,'step,cls,sim:',step,ncls,sim
        ! do while( .not.converged )
        !     step      = step+1
        !     sim_prev  = sim
        !     ncls      = ncls+1
        !     self%oris = os_copy
        !     ! shc clustering
        !     call shcc_cavg%new( self%pb,self%pp,ncls,self%oris,self%ptcls_ref,doprint=.true.,rotsh_srch=.true. )
        !     call shcc_cavg%exec
        !     call shcc_cavg%kill
        !     minpop = self%get_minpop( self%oris )
        !     ! simmilarity uodate
        !     !call self%prepptcls
        !     !call self%mask_ptcls
        !     !call self%build_simmat_pftcc
        !     sim = calc_sim( self%oris )
        !     ! Convergence
        !     print *,'step,cls,sim,sim_prev/simpop,minpop:',step,ncls,sim,sim_prev/sim,minpop ! debug
        !     if( minpop<=self%pop_thresh )converged=.true.
        !     if( sim_prev/sim>.98 )       converged=.true.
        !     if( step>10 )                converged=.true.
        ! enddo
        ! self%ncls = self%oris%get_ncls()
        ! write(*,'(A,I3)')'>>> SHC CLUSTERS FOUND=', self%ncls

        ! contains 
        !     function calc_sim( localoris )result( s )
        !         type( oris ), intent(inout) :: localoris
        !         integer :: i,j,pop,n
        !         real    :: s,si,icls
        !         n = localoris%get_noris()
        !         s = 0.
        !         do i=1,n
        !             si   = 0.
        !             pop  = 0
        !             icls = localoris%get(i,'class')
        !             do j=1,n
        !                 if( i==j ) cycle
        !                 if( icls == localoris%get(j,'class') )then
        !                     si = si+self%simmat(i,j)
        !                     pop = pop+1
        !                 endif
        !             enddo
        !             si = si/real(pop)
        !             s  = s+si
        !         enddo
        !         s = s/real(n)
        !     end function
    end subroutine exec_2dshc

    subroutine mask_ptcls( self )
        class(clustercavg), intent(inout)  :: self
        integer :: i
        do i=1,self%n
            call self%ptcls(i)%mask(self%pp%msk,'soft')
        enddo
    end subroutine

    subroutine exec_shc( self )
        !use simple_cluster_valid 
        class(clustercavg), intent(inout)  :: self
        type(shc_cluster)                  :: shcc
        real           :: sim, sim_prev
        integer        :: i,step,ncls, minpop
        logical        :: converged
        ! Init
        ncls      = 1
        minpop    = self%n
        converged = .false.
        step      = 0
        ! Initial similarity
        sim    = 0.
        do i=1,self%n
            sim = sim+sum(self%simmat(i,1:i-1))/real(self%n-1)
            sim = sim+sum(self%simmat(i,i+1:self%n))/real(self%n-1)
        enddo
        sim = sim/real(self%n)
        print *,'i,cls,sim:',i,ncls,sim
        do while( .not.converged )
            step     = step+1
            sim_prev = sim
            ncls     = ncls+1
            ! shc clustering
            call shcc%new( self%n, ncls, self%simmat, self%oris )
            call shcc%shc( .false., 'class', sim )
            call shcc%kill
            minpop = self%get_minpop( self%oris )
            ! cluster metrics
            print *,'i,cls,sim,sim_prev/simpop,minpop:',step,ncls,sim,sim_prev/sim,minpop ! debug
            ! Convergence
            ! Minimum population criterion
            if( minpop<=self%pop_thresh )converged=.true.
            if( sim_prev/sim>.95 )       converged=.true.
            if( step>50 )                converged=.true.
        enddo
        call self%print_pop( self%oris ) ! debug
        self%ncls = self%oris%get_ncls()
        write(*,'(A,I3)')'>>> SHC CLUSTERS FOUND=', self%ncls
    end subroutine exec_shc

    function get_feature_dist( self, os )result( dist )
        class(clustercavg), intent(inout) :: self
        type(oris), intent(inout) :: os
        real, allocatable         :: center(:),ini_center(:),vec(:)
        integer, allocatable      :: inds(:)
        real    :: dist
        integer :: i,ind,j,nfeat,alloc_stat,ncls,pop
        ncls  = os%get_ncls()
        nfeat = size( self%latent_vecs,2 )
        allocate( vec(nfeat),stat=alloc_stat )
        allocate( ini_center(nfeat),stat=alloc_stat )
        allocate( center(nfeat),stat=alloc_stat )
        ini_center = 0.
        dist       = 0.
        do i=1,nfeat
            ini_center(i) = sum(self%latent_vecs(:,i))/real(self%n)
        enddo
        do i=1,ncls
            pop = os%get_clspop( i )
            if( pop>1 )then
                center = 0.
                inds   = os%get_cls_pinds( i )
                do j=1,pop
                    ind    = inds(j)
                    center = center+self%latent_vecs(ind,:)/real(pop)
                enddo
                vec  = center-ini_center
                dist = dist+real(pop-1)*sum( vec**2 )
                deallocate( vec,inds )
            endif
        enddo
        deallocate( center,ini_center )
    end function get_feature_dist

    function cls_mixing( self, os, os_prev )result( mix )
        class(clustercavg), intent(inout)      :: self
        type(oris), intent(inout)              :: os, os_prev
        integer                                :: i, iptcl, ipop, pop_prev, ncls, ncls_prev
        integer, allocatable                   :: cls_arr(:,:)
        real                                   :: mix
        integer                                :: alloc_stat, cls_ind,loc(1), membership
        integer  :: icls, icls_prev
        ncls      = os%get_ncls()
        ncls_prev = os_prev%get_ncls()
        allocate( cls_arr(ncls,ncls_prev), stat=alloc_stat)
        call alloc_err('In: new; simple_clustercavg; cls_mixing', alloc_stat)
        cls_arr = 0
        do iptcl=1,os%get_noris()
            icls      = nint( os%get(iptcl,'class') )
            icls_prev = nint( os_prev%get(iptcl,'class') )
            cls_arr( icls,icls_prev ) = cls_arr( icls,icls_prev )+1
        enddo
        mix = 0.
        do icls=1,ncls
            ipop = sum( cls_arr(icls,:) )
            mix  = mix+real( maxval(cls_arr(icls,:)) )/real(ipop)
        enddo
        mix = mix/real(ncls)
    end function cls_mixing

    function get_minpop( self, os )result( minpop )
        class(clustercavg), intent(inout)      :: self
        type(oris), intent(inout) :: os
        integer                   :: i,ncls,minpop,pop
        ncls    = os%get_ncls()
        minpop  = os%get_clspop(1)
        do i=2,ncls
            pop = os%get_clspop(i)
            if( pop<minpop )minpop=pop
        enddo
    end function get_minpop

    function get_maxpopind( self, os )result( ind )
        class(clustercavg), intent(inout)      :: self
        type(oris), intent(inout) :: os
        integer                   :: i,ncls,maxpop,pop,ind
        ncls   = os%get_ncls()
        maxpop = 0
        ind    = 0
        do i=1,ncls
            pop = os%get_clspop(i)
            if( pop>maxpop )then
                maxpop = pop
                ind    = i
            endif
        enddo
    end function get_maxpopind

    function get_ncls( self )result( ncls )
        class(clustercavg), intent(inout)   :: self
        integer                             :: ncls
        ncls = self%oris%get_ncls()
    end function get_ncls

    subroutine print_pop( self, os )
        class(clustercavg), intent(inout)  :: self
        type(oris), intent(inout)          :: os
        integer                            :: i,n
        n = os%get_ncls()
        do i=1,n
            write(*,*)i,os%get_clspop(i)
        enddo
    end subroutine print_pop

    subroutine write_stk( self, fname )
        class(clustercavg), intent(inout) :: self
        character(len=STDLEN),intent(in) :: fname
        type(ori)   :: o
        real    :: x,y
        integer :: i
        do i=1,self%n
            call self%ptcls(i)%write( fname,i )
        enddo
    end subroutine write_stk

    !>  \brief  simmilatrity matrix based on feature euclidian distance
    subroutine build_simmat_feat( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(clustercavg), intent(inout) :: self
        integer                           :: i,j
        call self%init_simmat( 1.)
        !$omp parallel do default(shared) private(i,j)
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(i,j) = sqrt(sum( (self%latent_vecs(i,:)-self%latent_vecs(j,:))**2 ))
            enddo
        enddo
        !$omp end parallel do
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(j,i) = self%simmat(i,j)
            enddo
        enddo
        self%simmat = maxval(self%simmat)-self%simmat
    end subroutine

    subroutine build_simmat_realcc( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(clustercavg), intent(inout) :: self
        type(image)                       :: mask
        integer                           :: i,j
        if( .not.self%usemsk )stop 'Real space correlation must be used with a mask'
        call self%init_simmat( 1.)
        call mask%new( self%ptcls(1)%get_ldim(), self%pp%smpd )
        mask = 1.
        call mask%mask( self%pp%msk, 'hard' )
        !$omp parallel do default(shared) private(i,j)
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(i,j) = self%ptcls(i)%real_corr(self%ptcls(j))
            enddo
        enddo
        !$omp end parallel do
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(j,i) = self%simmat(i,j)
            enddo
        enddo
    end subroutine

    subroutine simmat2weights( self )
        use simple_stat, only: corrs2weights
        class(clustercavg), intent(inout) :: self
        real, allocatable                 :: weights(:), corrs(:)
        integer                           :: i,j,cnt,n,alloc_stat
        n = self%n*(self%n-1)/2
        allocate( corrs(n),stat=alloc_stat )
        call alloc_err('In: simple_clustercavg%simmat2weights', alloc_stat)
        cnt = 0
        do i=1,self%n-1
            do j=i+1,self%n
                cnt          = cnt+1
                corrs( cnt ) = self%simmat( i,j )
            enddo
        enddo
        weights = corrs2weights( corrs )
        cnt = 0
        do i=1,self%n-1
            do j=i+1,self%n
                cnt              = cnt+1
                self%simmat(i,j) = weights( cnt )
                self%simmat(j,i) = weights( cnt )
            enddo
        enddo
        deallocate( weights, corrs )
    end subroutine simmat2weights

    subroutine build_simmat_pftcc( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(clustercavg), intent(inout) :: self
        type(polarft_corrcalc)            :: pftcc
        integer                           :: i,j,roind
        if( .not.self%usemsk )stop 'Polar FT correlation must be used with a mask'
        ! PFTCC
        call pftcc%new(self%n, [1,self%n],self%ptcls(1)%get_ldim(),self%pp%kfromto, self%pp%ring2, 'no')
        do i=1,self%n
            call self%pb%proj%img2polarft( i,self%ptcls(i),pftcc,isptcl=.true.)
            call self%pb%proj%img2polarft( i,self%ptcls(i),pftcc,isptcl=.false.)
        enddo
        ! Similarity matrix
        roind = pftcc%get_roind(0.)
        call self%init_simmat(1.)
        !$omp parallel do default(shared) private(i,j)
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(i,j) = pftcc%corr( i,j,roind )
            enddo
        enddo
        !$omp end parallel do
        do i=1,self%n-1
            do j=i+1,self%n
                self%simmat(j,i) = self%simmat(i,j)
            enddo
        enddo
        call pftcc%kill
    end subroutine build_simmat_pftcc

    !>  \brief  builds class averages and masks
    subroutine make_cavgs( self )
        class(clustercavg),intent(inout) :: self
        type(image),allocatable          :: cavgs(:)
        integer                          :: i, icls, pop, alloc_stat
        allocate( cavgs(self%ncls),stat=alloc_stat )
        call alloc_err('In: simple_clustercavg%make_cavgs;', alloc_stat)
        ! Init CAvgs
        do i=1,self%ncls
            call cavgs(i)%new([self%pp%box,self%pp%box,1], self%pp%smpd)
            cavgs(i) = 0.
        enddo
        ! Summation
        do i=1,self%n
            icls = nint(self%oris%get(i,'class'))
            call cavgs( icls )%add( self%ptcls(i) )
        enddo
        ! Averaging
        do i=1,self%ncls
            pop = self%oris%get_clspop(i)
            if( pop>1 )then
                call cavgs(i)%div( real(pop) )    
                ! noise norm & Mask 
                call cavgs(i)%noise_norm( self%pp%msk )
                call cavgs(i)%mask( self%pp%msk,'soft' )
            endif
            !call cavgs(i)%write( 'final_cavgs.mrc', i)
        enddo
        deallocate( cavgs )
    end subroutine make_cavgs

    subroutine init_simmat( self, diag )
        class(clustercavg), intent(inout)  :: self
        real, intent(in)                   :: diag
        integer                            :: i, alloc_stat
        if( .not.allocated(self%simmat) )then
            allocate( self%simmat(self%n,self%n), stat=alloc_stat)
            call alloc_err('In: simple_clustercavg, init_simmat', alloc_stat)
        endif
        self%simmat = 0.
        do i=1,self%n
            self%simmat(i,i) = diag
        enddo
    end subroutine init_simmat

    subroutine write_simmat( self )
        class(clustercavg), intent(inout)  :: self
        integer :: i,j,filnum
        character(len=STDLEN) :: str
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file='simmat.txt')
        do i=1,self%n
            str = ''
            do j=1,self%n-1
                write(filnum,'(F8.3)',advance='no') self%simmat(i,j)
            enddo
            write(filnum,'(F8.3)') self%simmat(i,self%n)
        enddo
        close(filnum)
    end subroutine write_simmat

    !>  \brief  performs rotation, translation, low pass
    subroutine prepptcls( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(clustercavg), intent(inout)   :: self
        real                  :: e3s(self%n)
        type(ori)             :: o
        integer               :: i
        do i=1,self%n
            call self%ptcls(i)%copy( self%ptcls_ref(i) )
            ! translation
            o      = self%oris%get_ori( i )
            e3s(i) = o%e3get()
            call self%ptcls(i)%fwd_ft
            call self%ptcls(i)%shift( -o%get('x'),-o%get('y') )
            call self%ptcls(i)%bwd_ft
        enddo
        ! openmp does crash here
        !!!$omp parallel do default(shared) private(i)
        do i=1,self%n
            call self%ptcls(i)%rtsq(-e3s(i), 0., 0.)
        enddo
        !!!$omp end parallel do
    end subroutine prepptcls

    !>  \brief  refines the selected cluster
    subroutine refine_cls( self, icls )
        use simple_refinecluster
        class(clustercavg), intent(inout)   :: self
        integer, intent(in)                 :: icls
        type( image ),allocatable           :: tmp_imgs(:)
        integer, allocatable                :: inds(:)
        type(refinecluster)                 :: refinec
        type( oris )                        :: tmp_os
        type( ori )                         :: o
        integer                             :: i,ind,pop,alloc_stat
        pop = self%oris%get_clspop( icls )
        if( pop>1 )then
            allocate( tmp_imgs(pop), stat=alloc_stat)
            call alloc_err('In: new; simple_biclust_cavg%refine_cls', alloc_stat)
            call tmp_os%new( pop )
            inds = self%oris%get_cls_pinds( icls )
            do i=1,pop
                ind = inds(i)
                o = self%oris%get_ori( ind )
                call tmp_os%set_ori( i,o )
                call tmp_imgs( i )%copy( self%ptcls_ref(ind) )
            enddo
            call refinec%new( self%pb, self%pp, os=tmp_os, imgs=tmp_imgs, doprint=.false. )
            call refinec%exec
            call refinec%kill
            do i=1,pop
                o = tmp_os%get_ori( i )
                call self%oris%set_ori( inds(i),o )
            enddo
            deallocate( inds )
        endif
    end subroutine refine_cls

    !>  \brief  ppca-based feature extraction the selected cluster
    subroutine denoise_cls( self, icls )
        use simple_cavgppca
        class(clustercavg), intent(inout)   :: self
        integer, intent(in)                 :: icls
        type( image ),allocatable :: tmp_imgs(:)
        integer,      allocatable :: inds(:)
        type( cavgppca )          :: cavg_ppca
        type( oris )              :: tmp_os
        type( ori )               :: o
        integer                   :: i,ind,pop,alloc_stat
        pop = self%oris%get_clspop( icls )
        if( pop>self%pop_thresh )then
            allocate( tmp_imgs(pop), stat=alloc_stat)
            call alloc_err('In: new; simple_biclust_cavg%refine_cls', alloc_stat)
            call tmp_os%new( pop )
            inds = self%oris%get_cls_pinds( icls )
            do i=1,pop
                ind = inds(i)
                o = self%oris%get_ori( ind )
                call tmp_os%set_ori( i,o )
                call tmp_imgs( i )%copy( self%ptcls(ind) )
            enddo
            if( self%usemsk )then
                call cavg_ppca%new( self%pp,self%pb,tmp_os,tmp_imgs,self%nfeats,'realcc',.true. )
            else
                call cavg_ppca%new( self%pp,self%pb,tmp_os,tmp_imgs,self%nfeats,'realcc',.false. )
            endif
            call cavg_ppca%exec
            call cavg_ppca%kill
            do i=1,pop
                ind = inds(i)
                call self%ptcls( ind )%copy( tmp_imgs(i) )
            enddo
            deallocate( inds )
        endif
    end subroutine denoise_cls

    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(clustercavg), intent(inout)   :: self
        self%pp    => null()
        self%pb    => null()
        if( allocated(self%ptcls_inds)  )deallocate(self%ptcls_inds)
        if( allocated(self%ptcls)   )deallocate(self%ptcls)
        if( allocated(self%ptcls_ref)   )deallocate(self%ptcls_ref)
        if( allocated(self%latent_vecs) )deallocate(self%latent_vecs)
        if( allocated(self%simmat) )deallocate(self%simmat)
        call self%oris%kill
        self%exists = .false.
    end subroutine kill

end module simple_clustercavg
