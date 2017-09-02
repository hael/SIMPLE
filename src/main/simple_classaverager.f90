module simple_classaverager
use simple_build,   only: build
use simple_params,  only: params
use simple_syslib,  only: alloc_errchk
use simple_image,   only: image
use simple_strings, only: int2str_pad
use simple_defs     ! use all in there
implicit none

public :: classaverager
private

type ptcl_record
    integer              :: pind   = 0
    integer              :: eo     = -1  ! even is 0, odd is 1, default is -1
    real                 :: pw     = 0.0
    real                 :: kv     = 0.0
    real                 :: cs     = 0.0
    real                 :: fraca  = 0.0
    real                 :: dfx    = 0.0
    real                 :: dfy    = 0.0 
    real                 :: angast = 0.0
    integer, allocatable :: classes(:)
    integer, allocatable :: states(:)
    real,    allocatable :: ows(:)
    real,    allocatable :: e3s(:)
    real,    allocatable :: shifts(:,:)
end type ptcl_record

type classaverager
    private
    class(build),      pointer     :: bp => null()          !< pointer to build
    class(params),     pointer     :: pp => null()          !< pointer to params
    integer                        :: istart, iend          !< particle index range
    integer                        :: partsz                !< size of partition
    integer                        :: ncls                  !< number of classes
    integer                        :: nstates               !< number of states
    type(ptcl_record), allocatable :: precs(:)              !< particle records     
    type(image),       allocatable :: cavgs_even(:,:)       !< class averages
    type(image),       allocatable :: cavgs_odd(:,:)        !< -"-
    type(image),       allocatable :: cavgs_merged(:,:)     !< -"-
    type(image),       allocatable :: ctfsqsums_even(:,:)   !< CTF**2 sums for Wiener normalisation
    type(image),       allocatable :: ctfsqsums_odd(:,:)    !< -"-
    type(image),       allocatable :: ctfsqsums_merged(:,:) !< -"-
    logical                        :: l_hard_assign = .true.
    logical                        :: l_grid        = .true.
    logical                        :: exists        = .false.
  contains
    ! constructors
    procedure          :: new
    procedure          :: new_minimal
    ! setters/getters
    procedure, private :: init_cavgs_sums
    procedure, private :: get_indices
    procedure, private :: class_pop
    procedure          :: get_cavg
    ! calculators
    procedure          :: assemble_sums
    procedure, private :: apply_ctf_and_shift
    procedure          :: merge_eos_and_norm
    ! I/O
    procedure          :: write
    procedure          :: write_partial_sums
    procedure          :: read

    ! destructor
    procedure          :: kill
end type classaverager

integer, parameter :: BATCHTHRSZ = 20

contains

    !>  \brief  is a constructor
    !!          data is now managed so that all exclusions are taken care of here
    !!          which means properly balanced batches can be produced for both soft
    !!          and hard clustering solutions
    subroutine new( self, b, p, grid, prime3Dsrchobj )
        use simple_ori,          only: ori
        use simple_oris,         only: oris
        use simple_prime3D_srch, only: prime3D_srch
        class(classaverager),          intent(inout) :: self
        class(build),        target,   intent(inout) :: b
        class(params),       target,   intent(inout) :: p
        logical,             optional, intent(in)    :: grid
        class(prime3D_srch), optional, intent(inout) :: prime3Dsrchobj(p%fromp:p%top)
        real, allocatable :: ori_weights(:)
        type(ori)         :: orientation
        type(oris)        :: prime3D_oris, a_here
        integer           :: alloc_stat, cnt, n_incl, iori
        integer           :: cnt_ori, istate, icls, iptcl
        logical           :: l_reduce_projs
        ! destruct possibly pre-existing instance
        call self%kill
        ! set pointers
        self%bp => b
        self%pp => p
        ! set interpolation flag
        self%l_grid = .true.
        if( present(grid) ) self%l_grid = grid
        ! set nstates
        self%nstates = p%nstates
        ! work out range and partsz
        if( p%l_distr_exec )then
            self%istart = p%fromp
            self%iend   = p%top
            self%partsz = self%iend - self%istart + 1
        else
            self%istart = 1
            self%iend   = p%nptcls
            self%partsz = p%nptcls
        endif
        ! create the particle records
        allocate(self%precs(self%partsz), stat=alloc_stat)
        call alloc_errchk('new; simple_classaverager, self%precs', alloc_stat)
        ! create a copy of a that can be modified
        a_here = b%a
        cnt    = 0
        ! fetch data from a_here
        do iptcl=self%istart,self%iend
            cnt = cnt + 1
            ! exclusion condition
            if( nint(a_here%get(iptcl,'state')) == 0 .or.&
               &nint(a_here%get(iptcl,'state_balance')) == 0 .or.&
               &a_here%get(iptcl,'w') < TINY )then
                self%precs(cnt)%pind  = 0
                cycle
            endif
            ! parameter transfer
            self%precs(cnt)%pind  = iptcl
            self%precs(cnt)%eo    = nint(a_here%get(iptcl,'eo'))
            self%precs(cnt)%pw    = a_here%get(iptcl,'w')
            self%precs(cnt)%kv    = a_here%get(iptcl,'kv')
            self%precs(cnt)%cs    = a_here%get(iptcl,'cs')
            self%precs(cnt)%fraca = a_here%get(iptcl,'fraca')
            select case(p%tfplan%mode)
                case('astig') ! astigmatic CTF
                    self%precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    self%precs(cnt)%dfy    = a_here%get(iptcl,'dfy')
                    self%precs(cnt)%angast = a_here%get(iptcl,'angast')
                case('noastig') ! non-astigmatic CTF
                    self%precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    self%precs(cnt)%dfy    = self%precs(cnt)%dfx
                    self%precs(cnt)%angast = 0.
            end select
        end do
        self%l_hard_assign = .true.
        if( present(prime3Dsrchobj) )then
            l_reduce_projs = .false.
            if( p%nspace > NSPACE_BALANCE )then
                ! reduce the number of projection directions used for the class average representation
                call a_here%reduce_projs(NSPACE_BALANCE, p%nsym, p%eullims)
                l_reduce_projs = .true.
                self%ncls = NSPACE_BALANCE
            else
                self%ncls = p%nspace
            endif
            cnt = 0
            do iptcl=self%istart,self%iend
                cnt = cnt + 1
                ! inclusion condition
                if( self%precs(cnt)%pind > 0 )then
                    ! ori template
                    orientation = a_here%get_ori(iptcl)
                    if( p%npeaks > 1 )then
                        self%l_hard_assign = .false.
                        ! get orientation distribution
                        call prime3Dsrchobj(iptcl)%get_oris(prime3D_oris, orientation)
                        if( l_reduce_projs ) call prime3D_oris%reduce_projs(NSPACE_BALANCE, p%nsym, p%eullims)
                        ori_weights = prime3D_oris%get_all('ow')
                        n_incl = count(ori_weights > TINY)
                        if( n_incl >= 1 )then
                            ! allocate & set info in record
                            allocate( self%precs(cnt)%classes(n_incl),  self%precs(cnt)%states(n_incl),&
                                      self%precs(cnt)%ows(n_incl),      self%precs(cnt)%e3s(n_incl),&
                                      self%precs(cnt)%shifts(n_incl,2), stat=alloc_stat )
                            call alloc_errchk('new; simple_classaverager, record arrays', alloc_stat)
                            cnt_ori = 0
                            do iori=1,prime3D_oris%get_noris()
                                if( ori_weights(iori) > TINY )then
                                    cnt_ori = cnt_ori + 1
                                    self%precs(cnt)%classes(cnt_ori)  = nint(prime3D_oris%get(iori, 'proj'))
                                    self%precs(cnt)%states(cnt_ori)   = nint(prime3D_oris%get(iori, 'state'))
                                    self%precs(cnt)%ows(cnt_ori)      = prime3D_oris%get(iori, 'w')
                                    self%precs(cnt)%e3s(cnt_ori)      = prime3D_oris%e3get(iori)
                                    self%precs(cnt)%shifts(cnt_ori,1) = prime3D_oris%get(iori, 'x')
                                    self%precs(cnt)%shifts(cnt_ori,2) = prime3D_oris%get(iori, 'y')
                                endif
                            end do
                        else
                            self%precs(cnt)%pind = 0
                        endif
                        deallocate(ori_weights)
                        call prime3D_oris%kill
                    else
                        cycle
                    endif
                endif
            end do
        else
            self%ncls = p%ncls
            cnt       = 0
            do iptcl=self%istart,self%iend
                cnt = cnt + 1
                ! inclusion condition
                if( self%precs(cnt)%pind > 0 )then
                    ! allocate & set info in record
                    allocate( self%precs(cnt)%classes(1),  self%precs(cnt)%states(1),&
                              self%precs(cnt)%ows(1),      self%precs(cnt)%e3s(1),&
                              self%precs(cnt)%shifts(1,2), stat=alloc_stat )
                    call alloc_errchk('new; simple_classaverager, record arrays', alloc_stat)
                    self%precs(cnt)%classes(1)  = nint(a_here%get(iptcl, 'class'))
                    self%precs(cnt)%states(1)   = nint(a_here%get(iptcl, 'state'))
                    self%precs(cnt)%ows(1)      = a_here%get(iptcl, 'w')
                    self%precs(cnt)%e3s(1)      = a_here%e3get(iptcl)
                    self%precs(cnt)%shifts(1,1) = a_here%get(iptcl, 'x')
                    self%precs(cnt)%shifts(1,2) = a_here%get(iptcl, 'y')
                endif
            end do
        endif
        ! build cavg/CTF**2 arrays
        allocate(self%cavgs_even(p%nstates,self%ncls), self%cavgs_odd(p%nstates,self%ncls),&
        &self%cavgs_merged(p%nstates,self%ncls), self%ctfsqsums_even(p%nstates,self%ncls),&
        &self%ctfsqsums_odd(p%nstates,self%ncls), self%ctfsqsums_merged(p%nstates,self%ncls), stat=alloc_stat)
        call alloc_errchk('new; simple_classaverager, cavg/ctfsqsum arrays', alloc_stat)
        do istate=1,p%nstates
            do icls=1,self%ncls
                call self%cavgs_even(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%cavgs_odd(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%cavgs_merged(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_even(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_odd(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_merged(istate,icls)%new([p%box,p%box,1],p%smpd)
            end do
        end do
        call a_here%kill
        self%exists = .true.
    end subroutine new

    !>  \brief  is a minimal constructor
    subroutine new_minimal( self, b, p )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(classaverager),  intent(inout) :: self
        class(build),  target, intent(inout) :: b
        class(params), target, intent(inout) :: p
        integer :: alloc_stat, istate, icls
        ! set pointers
        self%bp => b
        self%pp => p
        ! set nstates
        self%nstates = p%nstates
        ! build cavg/CTF**2 arrays
        allocate(self%cavgs_even(p%nstates,self%ncls), self%cavgs_odd(p%nstates,self%ncls),&
        &self%cavgs_merged(p%nstates,self%ncls), self%ctfsqsums_even(p%nstates,self%ncls),&
        &self%ctfsqsums_odd(p%nstates,self%ncls), self%ctfsqsums_merged(p%nstates,self%ncls), stat=alloc_stat)
        call alloc_errchk('new; simple_classaverager, cavg/ctfsqsum arrays', alloc_stat)
        do istate=1,p%nstates
            do icls=1,self%ncls
                call self%cavgs_even(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%cavgs_odd(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%cavgs_merged(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_even(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_odd(istate,icls)%new([p%box,p%box,1],p%smpd)
                call self%ctfsqsums_merged(istate,icls)%new([p%box,p%box,1],p%smpd)
            end do
        end do
        self%exists = .true.
    end subroutine new_minimal

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                self%cavgs_even(istate,icls)       = 0.
                self%cavgs_odd(istate,icls)        = 0.
                self%cavgs_merged(istate,icls)     = 0.
                self%ctfsqsums_even(istate,icls)   = cmplx(0.,0.)
                self%ctfsqsums_odd(istate,icls)    = cmplx(0.,0.)
                self%ctfsqsums_merged(istate,icls) = cmplx(0.,0.)
            end do
        end do
    end subroutine init_cavgs_sums

    !>  \brief  is for getting allocatable arrays with particle/record/ori indices
    subroutine get_indices( self, state, class, pinds, iprecs, ioris )
        class(classaverager), intent(in)  :: self
        integer,              intent(in)  :: state, class
        integer, allocatable, intent(out) :: pinds(:)
        integer, allocatable, intent(out) :: iprecs(:)
        integer, allocatable, intent(out) :: ioris(:)
        integer :: pop, alloc_stat, i, sz, iprec, cnt
        logical, allocatable :: l_state_class(:)
        pop = self%class_pop(state, class)
        if( allocated(pinds) )  deallocate(pinds)
        if( allocated(iprecs) ) deallocate(iprecs)
        if( allocated(ioris)  ) deallocate(ioris)
        allocate(pinds(pop), iprecs(pop), ioris(pop), stat=alloc_stat)
        call alloc_errchk('get_iprecs_ioris; simple_classaverager', alloc_stat)
        cnt = 0
        do iprec=1,self%partsz
            sz = size(self%precs(iprec)%classes)
            allocate(l_state_class(sz))
            where( self%precs(iprec)%states .eq. state .and. self%precs(iprec)%classes .eq. class )
                l_state_class = .true.
            else where
                l_state_class = .false.
            endwhere
            if( any(l_state_class) )then
                do i=1,sz
                    if( l_state_class(i) )then
                        cnt = cnt + 1
                        pinds(cnt)  = self%precs(iprec)%pind
                        iprecs(cnt) = iprec
                        ioris(cnt)  = i
                    endif
                enddo
            endif
            deallocate(l_state_class)
        end do
    end subroutine get_indices

    !>  \brief  is for calculating class population
    function class_pop( self, state, class ) result( pop )
        class(classaverager), intent(in) :: self
        integer,              intent(in) :: state, class
        integer :: pop, iprec, sz
        logical, allocatable :: l_state_class(:)
        pop = 0
        do iprec=1,self%partsz
            sz = size(self%precs(iprec)%classes)
            allocate(l_state_class(sz))
            where( self%precs(iprec)%states .eq. state .and. self%precs(iprec)%classes .eq. class )
                l_state_class = .true.
            else where
                l_state_class = .false.
            endwhere
            pop = pop + count(l_state_class)
            deallocate(l_state_class)
        end do
    end function class_pop

    !>  \brief  is for getting a class average
    subroutine get_cavg( self, state, class, which, img )
        class(classaverager), intent(inout) :: self
        integer,              intent(in)    :: state, class
        character(len=*),     intent(in)    :: which
        class(image),         intent(out)   :: img
        select case(which)
            case('even')
                img = self%cavgs_even(state,class)
            case('odd')
                img = self%cavgs_odd(state,class)
            case('merged')
                img = self%cavgs_merged(state,class)
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: get_cavg'
        end select
    end subroutine get_cavg

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !!          using gridding interpolation in Fourier space or quadratic interpolation
    !!          in real-space
    subroutine assemble_sums( self )
        use simple_math,       only: cyci_1d, sqwin_2d, rotmat2d
        use simple_kbinterpol, only: kbinterpol
        use simple_gridding,   only: prep4cgrid
        use simple_projector,  only: projector
        use simple_map_reduce, only: split_nobjs_even
        class(classaverager), intent(inout)  :: self
        type(kbinterpol)             :: kbwin
        type(image)                  :: batch_imgsum_even, batch_imgsum_odd
        type(image),     allocatable :: batch_imgs(:)
        type(projector), allocatable :: padded_imgs(:)
        complex,         allocatable :: cmat_even(:,:), cmat_odd(:,:), comps(:,:)
        real,            allocatable :: w(:,:)
        integer,         allocatable :: ptcls_inds(:), batches(:,:), iprecs(:)
        integer,         allocatable :: ioris(:), cyc1(:), cyc2(:)
        complex   :: comp, zero
        real      :: loc(2), mat(2,2), winsz, pw
        integer   :: icls, iptcl, icls_pop, iprec, iori, cnt_progress
        integer   :: i, nbatches, batch, batchsz, cnt, lims(3,2), istate
        integer   :: ldim(3), ldim_pd(3), logi(3), phys(3), win(2,2)
        integer   :: cyc_lims(3,2), alloc_stat, wdim, incr, h, k, l, m
        if( .not. self%pp%l_distr_exec )then
            write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        endif
        ! init
        call self%init_cavgs_sums
        kbwin      = kbinterpol(KBWINSZ, KBALPHA)
        zero       = cmplx(0.,0.)
        winsz      = KBWINSZ
        ldim       = [self%pp%box,self%pp%box,1]
        ldim_pd    = nint(KBALPHA)*ldim
        ldim_pd(3) = 1
        wdim       = ceiling(KBALPHA*winsz) + 1
        ! state loop 
        cnt_progress = 0
        do istate=1,self%nstates
            ! class loop
            do icls=1,self%ncls
                cnt_progress = cnt_progress + 1
                call progress(cnt_progress, self%nstates * self%ncls)
                icls_pop = self%class_pop(istate, icls)
                if( icls_pop == 0 ) cycle
                call self%get_indices(istate, icls, ptcls_inds, iprecs, ioris)
                ! batch planning
                nbatches = ceiling(real(icls_pop)/real(self%pp%nthr*BATCHTHRSZ))
                batches  = split_nobjs_even(icls_pop, nbatches)
                ! batch loop, prep
                do batch=1,nbatches
                    ! prep batch
                    batchsz = batches(batch,2) - batches(batch,1) + 1
                    allocate(batch_imgs(batchsz), padded_imgs(batchsz))
                    ! batch particles loop
                    do i=1,batchsz
                        iptcl = self%istart - 1 + ptcls_inds(batches(batch,1) + i - 1)
                        iprec = iprecs(batches(batch,1) + i - 1)
                        iori  = ioris(batches(batch,1)  + i - 1)
                        ! stash images (this goes here or suffer bugs)
                        call read_img_from_stk( self%bp, self%pp, iptcl )
                        batch_imgs(i) = self%bp%img
                        ! CTF square sum & shift
                        call self%apply_ctf_and_shift(iprec, iori, batch_imgs(i))
                        ! create padded imgs
                        call padded_imgs(i)%new(ldim_pd, self%pp%smpd)
                        call prep4cgrid(batch_imgs(i), padded_imgs(i), self%pp%msk, kbwin)
                    enddo
                    if( self%l_grid )then ! Fourier rotation
                        lims     = padded_imgs(1)%loop_lims(2)
                        cyc_lims = padded_imgs(1)%loop_lims(3)
                        allocate(cmat_even(lims(1,1):lims(1,2), lims(2,1):lims(2,2)), cyc1(wdim), cyc2(wdim),&
                        &cmat_odd(lims(1,1):lims(1,2), lims(2,1):lims(2,2)), w(wdim, wdim), comps(wdim, wdim), stat=alloc_stat)
                        call alloc_errchk('assemble_sums; simple_classaverager', alloc_stat)
                        cmat_even = zero
                        cmat_odd = zero
                        !$omp parallel do default(shared) private(i,iprec,iori,h,k,l,m,loc,mat,logi,phys,cyc1,cyc2,w,comps,win,incr,pw)&
                        !$omp schedule(static) reduction(+:cmat_even,cmat_odd) proc_bind(close)
                        ! batch loop, convolution interpolation
                        do i=1,batchsz
                            iprec = iprecs(batches(batch,1) + i - 1)
                            iori  = ioris(batches(batch,1)  + i - 1)
                            ! prep weight
                            if( self%l_hard_assign )then
                                pw = self%precs(iprec)%pw
                            else
                                pw = self%precs(iprec)%pw * self%precs(iprec)%ows(iori)
                            endif
                            mat = rotmat2d( -self%precs(iprec)%e3s(iori) )
                            ! Fourier components loop
                            do h=lims(1,1),lims(1,2)
                                do k=lims(2,1),lims(2,2)
                                    loc   = matmul(real([h,k]),mat)
                                    win   = sqwin_2d(loc(1),loc(2), winsz)
                                    comps = zero
                                    w     = 1.
                                    do l=1,wdim
                                        incr = l - 1
                                        ! circular addresses
                                        cyc1(l) = cyci_1d(cyc_lims(1,:), win(1,1) + incr)
                                        cyc2(l) = cyci_1d(cyc_lims(2,:), win(2,1) + incr)
                                        ! interpolation kernel matrix
                                        w(l,:) = w(l,:) * kbwin%apod( real(win(1,1) + incr) - loc(1) )
                                        w(:,l) = w(:,l) * kbwin%apod( real(win(2,1) + incr) - loc(2) )
                                    enddo
                                    ! fetch fourier components
                                    do l=1,wdim
                                        do m=1,wdim
                                            if( w(l,m) == 0.) cycle
                                            logi       = [cyc1(l), cyc2(m), 0]
                                            phys       = padded_imgs(i)%comp_addr_phys(logi)
                                            comps(l,m) = padded_imgs(i)%get_fcomp(logi, phys)
                                        end do
                                    end do
                                    ! SUM( kernel x components )
                                    select case(self%precs(iprec)%eo)
                                        case(0)
                                            cmat_even(h,k) = cmat_even(h,k) + pw * sum(w * comps)
                                        case(1)
                                            cmat_odd(h,k)  = cmat_odd(h,k)  + pw * sum(w * comps)
                                    end select
                                    ! above is an optimized version of:
                                    ! cmat(h,k) = cmat(h,k) + padded_imgs(i)%extr_gridfcomp( [loc(1),loc(2),0.] )
                                end do
                            end do
                            ! cleanup
                            call padded_imgs(i)%kill
                        enddo
                        !$omp end parallel do
                        ! transfer to images
                        call batch_imgsum_even%new(ldim_pd, self%pp%smpd)
                        call batch_imgsum_odd%new(ldim_pd, self%pp%smpd)
                        call batch_imgsum_even%set_ft(.true.)
                        call batch_imgsum_odd%set_ft(.true.)
                        batch_imgsum_even = zero
                        batch_imgsum_odd  = zero
                        !$omp parallel do collapse(2) default(shared) private(h,k,logi,phys,comp)&
                        !$omp schedule(static) proc_bind(close)
                        do h=lims(1,1),lims(1,2)
                            do k=lims(2,1),lims(2,2)
                                logi = [h, k, 0]
                                phys = batch_imgsum_even%comp_addr_phys(logi)
                                comp = cmat_even(h, k)
                                call batch_imgsum_even%set_fcomp(logi, phys, comp)
                                comp = cmat_odd(h, k)
                                call batch_imgsum_odd%set_fcomp(logi, phys, comp)
                            end do 
                        end do
                        !$omp end parallel do
                        ! real space & clipping
                        call batch_imgsum_even%bwd_ft
                        call batch_imgsum_odd%bwd_ft
                        call batch_imgsum_even%clip_inplace(ldim)
                        call batch_imgsum_odd%clip_inplace(ldim)
                        ! cleanup
                        deallocate(padded_imgs,comps,w,cyc1,cyc2,cmat_even,cmat_odd)
                    else
                        ! real space rotation
                        call batch_imgsum_even%new([self%pp%box, self%pp%box, 1], self%pp%smpd)
                        call batch_imgsum_odd%new([self%pp%box, self%pp%box, 1], self%pp%smpd)
                        ! batch loop, quadratic interpolation
                        do i=1,batchsz
                            ! set indices
                            iptcl = self%istart - 1 + ptcls_inds(batches(batch,1) + i - 1)
                            iprec = iprecs(batches(batch,1) + i - 1)
                            iori  = ioris(batches(batch,1)  + i - 1)
                            ! rotate image
                            call batch_imgs(i)%rtsq( -self%precs(iprec)%e3s(iori), 0., 0. )
                            ! prep weight
                            if( self%l_hard_assign )then
                                pw = self%precs(iprec)%pw
                            else
                                pw = self%precs(iprec)%pw * self%precs(iprec)%ows(iori)
                            endif
                            ! add to sums
                            select case(self%precs(iprec)%eo)
                                case(0)
                                    call batch_imgsum_even%add(batch_imgs(i), pw)
                                case(1)
                                    call batch_imgsum_odd%add(batch_imgs(i), pw)
                            end select
                        enddo
                    endif
                    ! batch summation
                    call self%cavgs_even(istate,icls)%add( batch_imgsum_even )
                    call self%cavgs_odd(istate,icls)%add( batch_imgsum_odd )
                    ! batch cleanup
                    do i=1,batchsz
                        call batch_imgs(i)%kill
                    enddo
                    deallocate(batch_imgs)
                enddo
                ! class cleanup
                deallocate(ptcls_inds, iprecs, ioris)
            enddo
        enddo
        if( .not.self%pp%l_distr_exec ) call self%merge_eos_and_norm
    end subroutine assemble_sums

    !>  \brief  is for CTF application and shifting 
    !!          the class CTF**2 sum is updated and the image is FTed on exit
    subroutine apply_ctf_and_shift( self, iprec, iori, img )
        use simple_ctf,   only: ctf
        class(classaverager), intent(inout) :: self
        integer,              intent(in)    :: iprec, iori
        class(image),         intent(inout) :: img
        type(image) :: ctfsq
        type(ctf)   :: tfun
        real        :: pw
        integer     :: state, class
        call ctfsq%new(img%get_ldim(), self%pp%smpd)
        call ctfsq%set_ft(.true.)
        tfun = ctf(img%get_smpd(), self%precs(iprec)%kv, self%precs(iprec)%cs, self%precs(iprec)%fraca)
        call img%fwd_ft
        ! take care of the nominator
        select case(self%pp%tfplan%flag)
            case('yes')  ! multiply with CTF
                call tfun%apply_and_shift(img, ctfsq, self%precs(iprec)%shifts(iori,1),&
                &self%precs(iprec)%shifts(iori,2), self%precs(iprec)%dfx, 'ctf',&
                &self%precs(iprec)%dfy, self%precs(iprec)%angast)
            case('flip') ! multiply with abs(CTF)
                call tfun%apply_and_shift(img, ctfsq, self%precs(iprec)%shifts(iori,1),&
                &self%precs(iprec)%shifts(iori,2), self%precs(iprec)%dfx, 'abs',&
                &self%precs(iprec)%dfy, self%precs(iprec)%angast)
            case('mul','no')
                call tfun%apply_and_shift(img, ctfsq, self%precs(iprec)%shifts(iori,1),&
                &self%precs(iprec)%shifts(iori,2), self%precs(iprec)%dfx, '',&
                &self%precs(iprec)%dfy, self%precs(iprec)%angast)
        end select
        ! prep weight
        if( self%l_hard_assign )then
            pw = self%precs(iprec)%pw
        else
            pw = self%precs(iprec)%pw * self%precs(iprec)%ows(iori)
        endif
        ! prep indices
        state = self%precs(iprec)%states(iori)
        class = self%precs(iprec)%classes(iori)
        ! add to sums
        select case(self%precs(iprec)%eo)
            case(0)
                call self%ctfsqsums_even(state, class)%add(ctfsq, pw)
            case(1)
                call self%ctfsqsums_odd(state, class)%add(ctfsq, pw)
        end select
    end subroutine apply_ctf_and_shift

    !>  \brief  merges the even/odd pairs and normalises the sums
    subroutine merge_eos_and_norm( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                ! merging
                self%cavgs_merged(istate,icls) = 0.
                call self%cavgs_merged(istate,icls)%add(self%cavgs_even(istate,icls))
                call self%cavgs_merged(istate,icls)%add(self%cavgs_odd(istate,icls))
                self%ctfsqsums_merged(istate,icls) = cmplx(0.,0.)
                call self%ctfsqsums_merged(istate,icls)%add(self%ctfsqsums_even(istate,icls))
                call self%ctfsqsums_merged(istate,icls)%add(self%ctfsqsums_odd(istate,icls))
                ! (w*CTF)**2 density correction
                call self%cavgs_even(istate,icls)%fwd_ft
                call self%cavgs_even(istate,icls)%ctf_dens_correct(self%ctfsqsums_even(istate,icls))
                call self%cavgs_even(istate,icls)%bwd_ft
                call self%cavgs_odd(istate,icls)%fwd_ft
                call self%cavgs_odd(istate,icls)%ctf_dens_correct(self%ctfsqsums_odd(istate,icls))
                call self%cavgs_odd(istate,icls)%bwd_ft
                call self%cavgs_merged(istate,icls)%fwd_ft
                call self%cavgs_merged(istate,icls)%ctf_dens_correct(self%ctfsqsums_merged(istate,icls))
                call self%cavgs_merged(istate,icls)%bwd_ft
            end do
         end do
    end subroutine merge_eos_and_norm

    !>  \brief  writes class averages to disk
    subroutine write( self, which_iter, fname )
        use simple_fileio, only: add2fbody
        class(classaverager),       intent(inout) :: self
        integer,          optional, intent(in)    :: which_iter
        character(len=*), optional, intent(in)    :: fname
        integer :: istate, icls
        character(len=:), allocatable :: refname, refname_even, refname_odd
        if( present(which_iter) )then
            if( present(fname) ) stop &
            &'fname cannot be present together with which_iter; simple_classaverager :: write'
            self%pp%refs = 'cavgs_iter'//int2str_pad(which_iter,3)//self%pp%ext
        else
            if( present(fname) )then
                self%pp%refs = fname
            else
                self%pp%refs = 'startcavgs'//self%pp%ext
            endif
        endif
        if( self%pp%chunktag .ne. '' ) self%pp%refs = trim(self%pp%chunktag)//trim(self%pp%refs)
        do istate=1,self%nstates
            if( self%nstates > 1 )then
                refname = add2fbody(self%pp%refs, self%pp%ext, '_state'//int2str_pad(istate,2))
            else
                allocate(refname, source=trim(self%pp%refs))
            endif
            refname_even = add2fbody(refname, self%pp%ext, '_even')
            refname_odd  = add2fbody(refname, self%pp%ext, '_odd')
            do icls=1,self%ncls
                call self%cavgs_merged(istate, icls)%write(refname)
                call self%cavgs_even(istate, icls)%write(refname_even)
                call self%cavgs_odd(istate, icls)%write(refname_odd)
            end do
            deallocate(refname, refname_even, refname_odd)
        enddo
    end subroutine write

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine write_partial_sums( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        character(len=:), allocatable :: cae, cao, cte, cto
        do istate=1,self%nstates
            if( self%nstates > 1 )then
                allocate(cae, source='cavgs_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
            else
                allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
            endif
            do icls=1,self%ncls
                call self%cavgs_even(istate, icls)%write(cae)
                call self%cavgs_odd(istate, icls)%write(cao)
                call self%ctfsqsums_even(istate, icls)%write(cte)
                call self%ctfsqsums_odd(istate, icls)%write(cto)
            end do
            deallocate(cae, cao, cte, cto)
        end do
    end subroutine write_partial_sums

    !>  \brief  reads class averages from disk
    subroutine read( self )
        use simple_fileio, only: add2fbody
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        character(len=:), allocatable :: refname, refname_even, refname_odd
        do istate=1,self%nstates
            if( self%nstates > 1 )then
                refname = add2fbody(self%pp%refs, self%pp%ext, '_state'//int2str_pad(istate,2))
            else
                allocate(refname, source=trim(self%pp%refs))
            endif
            refname_even = add2fbody(refname, self%pp%ext, '_even')
            refname_odd  = add2fbody(refname, self%pp%ext, '_odd')
            do icls=1,self%ncls
                call self%cavgs_merged(istate, icls)%read(refname)
                call self%cavgs_even(istate, icls)%read(refname_even)
                call self%cavgs_odd(istate, icls)%read(refname_odd)
            end do
            deallocate(refname, refname_even, refname_odd)
        enddo
    end subroutine read

    !>  \brief  re-generates the object after distributed execution
    subroutine assemble_sums_from_parts( self )
        use simple_fileio, only: file_exists
        class(classaverager), intent(inout) :: self
        character(len=:), allocatable :: cae, cao, cte, cto
        integer :: ipart, istate, icls
        call self%init_cavgs_sums
        do istate=1,self%nstates
            if( self%nstates > 1 )then
                allocate(cae, source='cavgs_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
            else
                allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
            endif
            do ipart=1,self%pp%nparts
                if( file_exists(cae) )then
                    do icls=1,self%ncls
                        call self%bp%img%read(cae, icls)
                        call self%cavgs_even(istate,icls)%add(self%bp%img)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cae)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                if( file_exists(cao) )then
                    do icls=1,self%ncls
                        call self%bp%img%read(cao, icls)
                        call self%cavgs_odd(istate,icls)%add(self%bp%img)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cao)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                if( file_exists(cte) )then
                    do icls=1,self%ncls
                        call self%bp%img%read(cte, icls)
                        call self%ctfsqsums_even(istate,icls)%add(self%bp%img)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cte)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif

                if( file_exists(cto) )then
                    do icls=1,self%ncls
                        call self%bp%img%read(cto, icls)
                        call self%ctfsqsums_odd(istate,icls)%add(self%bp%img)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cto)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
            end do
        end do
        call self%merge_eos_and_norm()
    end subroutine assemble_sums_from_parts

    !>  \brief  is a destructor
    subroutine kill( self )
        class(classaverager),  intent(inout) :: self
        integer :: istate, icls, iprec
        if( self%exists )then
            self%bp => null()
            self%pp => null()
            do istate=1,self%nstates
                do icls=1,self%ncls
                    call self%cavgs_even(istate,icls)%kill
                    call self%cavgs_odd(istate,icls)%kill
                    call self%cavgs_merged(istate,icls)%kill
                    call self%ctfsqsums_even(istate,icls)%kill
                    call self%ctfsqsums_odd(istate,icls)%kill
                    call self%ctfsqsums_merged(istate,icls)%kill
                end do
            end do
            do iprec=1,self%partsz
                deallocate( self%precs(iprec)%classes, self%precs(iprec)%states,&
                &self%precs(iprec)%ows, self%precs(iprec)%e3s, self%precs(iprec)%shifts)
            end do
            deallocate(self%precs)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_classaverager
