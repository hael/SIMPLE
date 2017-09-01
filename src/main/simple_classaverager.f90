module simple_classaverager
use simple_build,  only: build
use simple_params, only: params
use simple_syslib, only: alloc_errchk
implicit none

public :: classaverager
private

type ptcl_record
    integer              :: pind   = 0
    integer              :: eo     = -1 ! even is 0, odd is 1, default is -1
    real                 :: pw     = 0.0
    real                 :: kv     = 0.0
    real                 :: cs     = 0.0
    real                 :: fraca  = 0.0
    real                 :: dfx    = 0.0
    real                 :: dfy    = 0.0 
    real                 :: angast = 0.0
    integer, allocatable :: classes(:)
    real,    allocatable :: ows(:)
    real,    allocatable :: e3s(:)
    real,    allocatable :: shifts(:,:)
end type ptcl_record

type classaverager
    private
    class(build),      pointer     :: bp => null ()         !< pointer to build
    class(params),     pointer     :: pp => null()          !< pointer to params
    integer                        :: istart, iend          !< particle index range
    integer                        :: partsz                !< size of partition
    integer                        :: ncls                  !< number of classes
    integer                        :: nstates               !< number of states
    type(ptcl_record), allocatable :: precs(:)              !< particle records     
    type(image),       allocatable :: cavgs_even(:,:)       !< class averages
    type(image),       allocatable :: cavgs_odd(:,:)        ! -"-
    type(image),       allocatable :: cavgs_merged(:,:)     ! -"-
    type(image),       allocatable :: ctfsqsums_even(:,:)   !< CTF**2 sums for Wiener normalisation
    type(image),       allocatable :: ctfsqsums_odd(:,:)    !< -"-
    type(image),       allocatable :: ctfsqsums_merged(:,:) !< -"-
    logical                        :: l_hard_assign = .false.
    logical                        :: l_grid        = .false.
    logical                        :: exists        = .false.
contains
    ! constructor
    procedure :: new
    ! setters/getters
    procedure :: init_cavgs_sums
    procedure :: get_indices
    procedure :: class_pop


end type classaverager

integer, parameter :: BATCHTHRSZ = 20

contains

    !>  \brief  is the constructor
    !!          data is now managed so that all exclusions are taken care of here
    !!          which means properly balanced batches can be produced for both soft
    !!          and hard clustering solutions
    subroutine new( self, b, p, grid, prime3Dsrchobj )
        use simple_ori,    only: ori
        use simple_oris,   only: oris
        class(classaverager),  intent(inout) :: self
        class(build),  target, intent(inout) :: b
        class(params), target, intent(inout) :: b
        real, allocatable :: ori_weights(:)
        type(ori)         :: orientation
        type(oris)        :: prime3D_oris, a_here
        integer           :: alloc_stat, cnt, n_incl, iori, cnt_ori
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
                        call primesrch3D(iptcl)%get_oris(prime3D_oris, orientation)
                        if( l_reduce_projs ) call prime3D_oris%reduce_projs(NSPACE_BALANCE, p%nsym, p%eullims)
                        ori_weights = prime3D_oris%get_all('ow')
                        n_incl = count(ori_weights > TINY)
                        if( n_incl >= 1 )then
                            ! allocate & set info in record
                            allocate( self%precs(cnt)%classes(n_incl), self%precs(cnt)%states(n_incl),&
                                      self%precs(cnt)%ows(n_incl), self%precs(cnt)%e3s(n_incl),&
                                      self%precs(cnt)%shifts(n_incl,2), stat=alloc_stat
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
            cnt = 0
            do iptcl=self%istart,self%iend
                cnt = cnt + 1
                ! inclusion condition
                if( self%precs(cnt)%pind > 0 )then
                    ! allocate & set info in record
                    allocate( self%precs(cnt)%classes(1), self%precs(cnt)%states(1),&
                              self%precs(cnt)%ows(1), self%precs(cnt)%e3s(1),&
                              self%precs(cnt)%shifts(1,2), stat=alloc_stat
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
        &self%ctfsqsums_odd(p%nstates,self%ncls), self%ctfsqsums_merged(p%nstates,self%ncls, stat=alloc_stat)
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

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums( self )
        class(classaverager), intent(inout) :: self
        do istate=1,self%nstates
            do icls=1,self%ncls
                self%cavgs_even(icls)       = 0.
                self%cavgs_odd(icls)        = 0.
                self%cavgs_merged(icls)     = 0.
                self%ctfsqsums_even(icls)   = cmplx(0.,0.)
                self%ctfsqsums_odd(icls)    = cmplx(0.,0.)
                self%ctfsqsums_merged(icls) = cmplx(0.,0.)
            end do
        end do
    end subroutine init_cavgs_sums

    !>  \brief  is for getting allocatable arrays with particle/record/ori indices
    function get_indices( self, state, class, pinds, iprecs, ioris )
        class(classaverager), intent(in)  :: self
        integer,              intent(in)  :: state, class
        integer, allocatable, intent(out) :: pinds(:)
        integer, allocatable, intent(out) :: iprecs(:)
        integer, allocatable, intent(out) :: ioris(:)
        integer :: pop, alloc_stat, i, sz, iprec
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
            endif
            if( any(l_state_class) )then
                do i=1,sz
                    if( l_state_class(i) )then
                        cnt = cnt + 1
                        pinds(cnt)  = self%precs(iprec)%pind
                        iprecs(cnt) = iprec
                        ioris(cnt)  = i
                    endif
                endif
            endif
            deallocate(l_state_class)
        end do
    end function get_indices

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
            endif
            pop = pop + count(l_state_class)
            deallocate(l_state_class)
        end do
    end function class_pop

    !>  \brief  is for CTF application and shifting 
    !!          the class CTF**2 sum is updated and the image is FTed on exit
    subroutine apply_ctf_and_shift( self, iprec, iori, img )
        use simple_image, only: image
        use simple_ctf,   only: ctf
        class(classaverager), intent(inout) :: self
        integer,              intent(in)    :: iprec, iori
        class(image),         intent(inout) :: img
        type(image) :: ctfsq
        type(ctf)   :: tfun
        real        :: w
        integer     :: state, class
        call ctfsq%new(img%get_ldim(), p%smpd)
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
            w = self%precs(iprec)%pw
        else
            w = self%precs(iprec)%pw * self%precs(iprec)%ows(iori)
        endif
        ! prep indices
        state = self%precs(iprec)%states(iori)
        class = self%precs(iprec)%classes(iori)
        ! add to sums
        select case(self%precs(iprec)%eo)
            case(0)
                call self%ctfsqsums_even(state, class)%add(ctfsq, w)
            case(1)
                call self%ctfsqsums_odd(state, class)%add(ctfsq, w)
        end select
    end subroutine apply_ctf_and_shift

    subroutine assemble_sums( b, p, grid )
        use simple_projector_hlev, only: rot_imgbatch
        use simple_map_reduce,     only: split_nobjs_even
        use simple_oris,           only: oris
        use simple_ctf,            only: ctf
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        logical, optional, intent(in)    :: grid
        type(oris)               :: a_here, batch_oris
        type(ori)                :: orientation
        type(image)              :: batch_imgsum, cls_imgsum
        type(image), allocatable :: batch_imgs(:) 
        integer,     allocatable :: ptcls_inds(:), batches(:,:), iprecs(:), ioris(:)
        real      :: w
        integer   :: icls, iptcl, icls_pop, iprec, iori, cnt_progress
        integer   :: i, nbatches, batch, batchsz, cnt

        
        


        if( .not. self%pp%l_distr_exec )then
            write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        endif
        call self%init_cavgs_sums
        




        ! state loop 
        cnt_progress = 0
        do istate=1,self%nstates
            ! class loop
            do icls=1,self%ncls
                cnt_progress = cnt_progress + 1
                call progress(cnt_progress, self%nstates * self%ncls)
                icls_pop = self%class_pop(istate, icls)
                if( icls_pop == 0 ) cycle
                call cls_imgsum%new([self%pp%box, self%pp%box, 1], self%pp%smpd)
                call self%get_indices(istate, icls, ptcls_inds, iprecs, ioris)
                ! batch planning
                nbatches = ceiling(real(icls_pop)/real(self%pp%nthr*BATCHTHRSZ))
                batches  = split_nobjs_even(icls_pop, nbatches)
                ! batch loop
                do batch=1,nbatches
                    ! prep batch
                    batchsz = batches(batch,2) - batches(batch,1) + 1
                    allocate(batch_imgs(batchsz))
                    ! batch particles loop
                    do i=1,batchsz
                        iptcl = self%istart - 1 + ptcls_inds(batches(batch,1) + i - 1)
                        iprec = iprecs(batches(batch,1) + i - 1)
                        iori  = ioris(batches(batch,1)  + i - 1)
                        ! stash images (this goes here or suffer bugs)
                        call read_img_from_stk( self%pb, self%pp, iptcl )
                        batch_imgs(i) = self%bp%img
                        ! CTF square sum & shift
                        call apply_ctf_and_shift(iprec, iori, batch_imgs(i))
                    enddo







                    if( self%l_grid )then
                        ! rotate batch by gridding
                        call rot_imgbatch(batch_imgs, batch_oris, batch_imgsum, p%msk)
                    else
                        ! real space rotation
                        call batch_imgsum%new([p%box, p%box, 1], p%smpd)
                        do i = 1,batchsz
                            iptcl = istart - 1 + ptcls_inds(batches(batch,1)+i-1)
                            orientation = b%a%get_ori(iptcl)
                            w = orientation%get('w')
                            call batch_imgs(i)%rtsq( -orientation%e3get(), 0., 0. )
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            !call batch_imgsum%add(batch_imgs(i), w)
                        enddo
                    endif
                    ! batch summation
                    call cls_imgsum%add( batch_imgsum )
                    ! batch cleanup
                    do i = 1, batchsz
                        call batch_imgs(i)%kill
                    enddo
                    deallocate(batch_imgs)
                enddo
                ! set class
                b%cavgs(icls) = cls_imgsum
                ! class cleanup
                deallocate(ptcls_inds)
            enddo
        enddo
        if( .not.p%l_distr_exec ) call prime2D_norm_sums( b, p )
    end subroutine assemble_sums





end module simple_classaverager
