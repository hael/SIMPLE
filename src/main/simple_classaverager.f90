module simple_classaverager
#include "simple_lib.f08"
use simple_ctf,     only: ctf
use simple_build,   only: build
use simple_params,  only: params
use simple_image,   only: image
use simple_timer    ! use all in there
implicit none

public :: classaverager
private

type ptcl_record
    type(ctf)            :: tfun          !< transfer function
    integer              :: pind    = 0   !< particle index in stack
    integer              :: eo      = -1  !< even is 0, odd is 1, default is -1
    real                 :: pw      = 0.0 !< particle weight
    real                 :: dfx     = 0.0 !< defocus in x (microns)
    real                 :: dfy     = 0.0 !< defocus in y (microns)
    real                 :: angast  = 0.0 !< angle of astigmatism (in degrees)
    real                 :: phshift = 0.0 !< additional phase shift from the Volta
    integer, allocatable :: classes(:)    !< class assignments (can be many per particle in 3D case, hence the array)
    integer, allocatable :: states(:)     !< state assignments
    integer, allocatable :: inpl_inds(:)  !< in-plane rotation indices
    real,    allocatable :: ows(:)        !< orientation weights
    real,    allocatable :: e3s(:)        !< in-plane rotations
    real,    allocatable :: shifts(:,:)   !< rotational origin shifts
end type ptcl_record

type classaverager
    private
    class(build),      pointer     :: bp => null()             !< pointer to build
    class(params),     pointer     :: pp => null()             !< pointer to params
    type(CTFFLAGTYPE)              :: ctf                      !< ctf flag <yes|no|mul|flip>
    integer                        :: istart     = 0, iend = 0 !< particle index range
    integer                        :: partsz     = 0           !< size of partition
    integer                        :: ncls       = 0           !< # classes
    integer                        :: nstates    = 0           !< # states
    integer                        :: nrots      = 0           !< # in-plane rotations
    integer                        :: ldim(3)    = [0,0,0]     !< logical dimension of image
    integer                        :: ldim_pd(3) = [0,0,0]     !< logical dimension of image, padded
    real                           :: smpd       = 0.          !< sampling distance
    type(ptcl_record), allocatable :: precs(:)                 !< particle records     
    type(image),       allocatable :: cavgs_even(:,:)          !< class averages
    type(image),       allocatable :: cavgs_odd(:,:)           !< -"-
    type(image),       allocatable :: cavgs_merged(:,:)        !< -"-
    type(image),       allocatable :: ctfsqsums_even(:,:)      !< CTF**2 sums for Wiener normalisation
    type(image),       allocatable :: ctfsqsums_odd(:,:)       !< -"-
    type(image),       allocatable :: ctfsqsums_merged(:,:)    !< -"-
    real,              allocatable :: inpl_rots(:)             !< in-plane rotations (sign shifted)
    logical                        :: phaseplate    = .false.  !< Volta phaseplate images or not
    logical                        :: use_kbmem     = .false.  !< use memoization of K-B kernel or not
    logical                        :: l_is_class    = .true.   !< for prime2D or not
    logical                        :: l_hard_assign = .true.   !< npeaks == 1 or not
    logical                        :: exists        = .false.  !< to flag instance existence 
  contains
    ! constructors
    procedure          :: new
    ! setters/getters
    procedure          :: transf_oridat
    procedure, private :: init_cavgs_sums
    procedure, private :: get_indices
    procedure, private :: class_pop
    procedure          :: get_cavg
    procedure          :: set_cavg
    ! calculators
    procedure          :: assemble_sums
    procedure, private :: calc_tfun_vals
    procedure          :: merge_eos_and_norm
    procedure          :: calc_and_write_frcs
    procedure          :: eoavg
    ! I/O
    procedure          :: write
    procedure          :: read
    procedure          :: write_partial_sums
    procedure          :: assemble_sums_from_parts
    ! destructor
    procedure          :: kill
end type classaverager

integer, parameter      :: BATCHTHRSZ = 50
logical, parameter      :: L_BENCH    = .true.
integer(timer_int_kind) :: t_batch_loop, t_gridding, t_tot
real(timer_int_kind)    :: rt_batch_loop, rt_gridding, rt_tot
character(len=STDLEN)   :: benchfname

contains

    !>  \brief  is a constructor
    !!          data is now managed so that all exclusions are taken care of here
    !!          which means properly balanced batches can be produced for both soft
    !!          and hard clustering solutions
    subroutine new( self, b, p, which, inpl_rots )
        class(classaverager),  intent(inout) :: self         !< instance
        class(build),  target, intent(inout) :: b            !< builder
        class(params), target, intent(inout) :: p            !< params
        character(len=*),      intent(in)    :: which        !< class/proj
        real, optional,        intent(in)    :: inpl_rots(:) !< in-plane rotations
        integer :: alloc_stat, istate, icls
        ! destruct possibly pre-existing instance
        call self%kill
        ! set pointers
        self%bp => b
        self%pp => p
        ! set nstates
        self%nstates = p%nstates
        ! take care of in-plane rots
        if( present(inpl_rots) )then
            self%use_kbmem = .true.
            self%nrots = size(inpl_rots)
            allocate(self%inpl_rots(self%nrots), source=inpl_rots)
        endif
        ! class or proj
        select case(which)
            case('class')
                self%l_is_class = .true.
                self%ncls       = p%ncls
            case('proj')
                self%l_is_class = .false.
                ! possible reduction of # projection directions used 
                ! for the class average representation
                self%ncls = min(NSPACE_BALANCE,p%nspace)
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: new'
        end select
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
        ! CTF logics
        select case(p%ctf)
            case('no')
                self%ctf%flag = CTFFLAG_NO
            case('yes')
                self%ctf%flag = CTFFLAG_YES
            case('mul')
                stop 'ERROR ctf=mul deprecated; simple_classaverager_dev :: new'
            case('flip')
                self%ctf%flag = CTFFLAG_FLIP
        end select
        self%phaseplate = p%tfplan%l_phaseplate
        ! smpd
        self%smpd       = p%smpd
        ! set ldims
        self%ldim       = [self%pp%box,self%pp%box,1]
        self%ldim_pd    = nint(KBALPHA) * self%ldim
        self%ldim_pd(3) = 1
        ! build arrays
        allocate(self%precs(self%partsz), self%cavgs_even(p%nstates,self%ncls), self%cavgs_odd(p%nstates,self%ncls),&
        &self%cavgs_merged(p%nstates,self%ncls), self%ctfsqsums_even(p%nstates,self%ncls),&
        &self%ctfsqsums_odd(p%nstates,self%ncls), self%ctfsqsums_merged(p%nstates,self%ncls), stat=alloc_stat)
        call alloc_errchk('new; simple_classaverager_dev', alloc_stat)
        do istate=1,p%nstates
            do icls=1,self%ncls
                call self%cavgs_even(istate,icls)%new(self%ldim_pd,p%smpd)
                call self%cavgs_odd(istate,icls)%new(self%ldim_pd,p%smpd)
                call self%cavgs_merged(istate,icls)%new(self%ldim_pd,p%smpd)
                call self%ctfsqsums_even(istate,icls)%new(self%ldim_pd,p%smpd)
                call self%ctfsqsums_odd(istate,icls)%new(self%ldim_pd,p%smpd)
                call self%ctfsqsums_merged(istate,icls)%new(self%ldim_pd,p%smpd)
            end do
        end do
        ! flag existence
        self%exists = .true.
    end subroutine new

    ! setters/getters

    !>  \brief  transfers orientation data to the instance
    subroutine transf_oridat( self, a, prime3Dsrchobj )
        use simple_ori,          only: ori
        use simple_oris,         only: oris
        use simple_prime3D_srch, only: prime3D_srch
        class(classaverager),          intent(inout) :: self  !< instance
        class(oris),                   intent(in)    :: a
        class(prime3D_srch), optional, intent(inout) :: prime3Dsrchobj(self%istart:self%iend)
        real, allocatable :: ori_weights(:)
        type(ori)         :: orientation
        type(oris)        :: prime3D_oris, a_here
        integer           :: alloc_stat, cnt, n_incl, iori
        integer           :: cnt_ori, istate, icls, iptcl
        logical           :: l_reduce_projs
        ! create a copy of a that can be modified
        a_here = a
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
            self%precs(cnt)%tfun  = ctf(self%pp%smpd, a_here%get(iptcl,'kv'), a_here%get(iptcl,'cs'), a_here%get(iptcl,'fraca'))
            select case(self%pp%tfplan%mode)
                case('astig') ! astigmatic CTF
                    self%precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    self%precs(cnt)%dfy    = a_here%get(iptcl,'dfy')
                    self%precs(cnt)%angast = a_here%get(iptcl,'angast')
                case('noastig') ! non-astigmatic CTF
                    self%precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    self%precs(cnt)%dfy    = self%precs(cnt)%dfx
                    self%precs(cnt)%angast = 0.
            end select
            self%precs(cnt)%phshift = 0.
            if( self%phaseplate ) self%precs(cnt)%phshift = a_here%get(iptcl,'phshift')
        end do
        self%l_hard_assign = .true.
        if( present(prime3Dsrchobj) )then
            l_reduce_projs = .false.
            if( self%pp%nspace > NSPACE_BALANCE )then
                ! reduce # projection directions 
                ! used for the class average representation
                call a_here%reduce_projs(NSPACE_BALANCE, self%pp%nsym, self%pp%eullims)
                l_reduce_projs = .true.
            endif
            cnt = 0
            do iptcl=self%istart,self%iend
                cnt = cnt + 1
                ! inclusion condition
                if( self%precs(cnt)%pind > 0 )then
                    ! ori template
                    orientation = a_here%get_ori(iptcl)
                    if( self%pp%npeaks > 1 )then
                        self%l_hard_assign = .false.
                        ! get orientation distribution
                        call prime3Dsrchobj(iptcl)%get_oris(prime3D_oris, orientation)
                        if( l_reduce_projs ) call prime3D_oris%reduce_projs(NSPACE_BALANCE, self%pp%nsym, self%pp%eullims)
                        ori_weights = prime3D_oris%get_all('ow')
                        n_incl = count(ori_weights > TINY)
                        if( n_incl >= 1 )then
                            ! allocate & set info in record
                            if( allocated(self%precs(cnt)%classes)  ) deallocate(self%precs(cnt)%classes)
                            if( allocated(self%precs(cnt)%inpl_inds)) deallocate(self%precs(cnt)%inpl_inds)
                            if( allocated(self%precs(cnt)%states)   ) deallocate(self%precs(cnt)%states)
                            if( allocated(self%precs(cnt)%ows)      ) deallocate(self%precs(cnt)%ows)
                            if( allocated(self%precs(cnt)%e3s)      ) deallocate(self%precs(cnt)%e3s)
                            if( allocated(self%precs(cnt)%shifts)   ) deallocate(self%precs(cnt)%shifts)
                            allocate( self%precs(cnt)%classes(n_incl),  self%precs(cnt)%states(n_incl),&
                                      self%precs(cnt)%ows(n_incl),      self%precs(cnt)%e3s(n_incl),&
                                      self%precs(cnt)%shifts(n_incl,2), self%precs(cnt)%inpl_inds(n_incl), stat=alloc_stat )
                            call alloc_errchk('new; simple_classaverager, record arrays', alloc_stat)
                            cnt_ori = 0
                            do iori=1,prime3D_oris%get_noris()
                                if( ori_weights(iori) > TINY )then
                                    cnt_ori = cnt_ori + 1
                                    self%precs(cnt)%classes(cnt_ori)   = nint(prime3D_oris%get(iori, 'proj'))
                                    self%precs(cnt)%inpl_inds(cnt_ori) = nint(prime3D_oris%get(iori, 'proj'))
                                    self%precs(cnt)%states(cnt_ori)    = nint(prime3D_oris%get(iori, 'inpl'))
                                    self%precs(cnt)%ows(cnt_ori)       = prime3D_oris%get(iori, 'w')
                                    self%precs(cnt)%e3s(cnt_ori)       = prime3D_oris%e3get(iori)
                                    self%precs(cnt)%shifts(cnt_ori,1)  = prime3D_oris%get(iori, 'x')
                                    self%precs(cnt)%shifts(cnt_ori,2)  = prime3D_oris%get(iori, 'y')
                                endif
                            end do
                        else
                            self%precs(cnt)%pind = 0
                        endif
                        deallocate(ori_weights)
                        call prime3D_oris%kill
                    endif
                endif
            end do
        else
            cnt = 0
            do iptcl=self%istart,self%iend
                cnt = cnt + 1
                ! inclusion condition
                if( self%precs(cnt)%pind > 0 )then
                    ! allocate & set info in record
                    if( allocated(self%precs(cnt)%classes) )  deallocate(self%precs(cnt)%classes)
                    if( allocated(self%precs(cnt)%inpl_inds)) deallocate(self%precs(cnt)%inpl_inds)
                    if( allocated(self%precs(cnt)%states)  )  deallocate(self%precs(cnt)%states)
                    if( allocated(self%precs(cnt)%ows)     )  deallocate(self%precs(cnt)%ows)
                    if( allocated(self%precs(cnt)%e3s)     )  deallocate(self%precs(cnt)%e3s)
                    if( allocated(self%precs(cnt)%shifts)  )  deallocate(self%precs(cnt)%shifts)
                    allocate( self%precs(cnt)%classes(1),  self%precs(cnt)%states(1),&
                              self%precs(cnt)%ows(1),      self%precs(cnt)%e3s(1),&
                              self%precs(cnt)%shifts(1,2), self%precs(cnt)%inpl_inds(1), stat=alloc_stat )
                    call alloc_errchk('new; simple_classaverager, record arrays', alloc_stat)
                    self%precs(cnt)%classes(1)   = nint(a_here%get(iptcl, 'class'))
                    self%precs(cnt)%inpl_inds(1) = nint(a_here%get(iptcl, 'inpl'))
                    self%precs(cnt)%states(1)    = nint(a_here%get(iptcl, 'state'))
                    self%precs(cnt)%ows(1)       = a_here%get(iptcl, 'w')
                    self%precs(cnt)%e3s(1)       = a_here%e3get(iptcl)
                    self%precs(cnt)%shifts(1,1)  = a_here%get(iptcl, 'x')
                    self%precs(cnt)%shifts(1,2)  = a_here%get(iptcl, 'y')
                endif
            end do
        endif
        call a_here%kill
    end subroutine transf_oridat

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                call self%cavgs_even(istate,icls)%new(self%ldim_pd,self%smpd)
                call self%cavgs_odd(istate,icls)%new(self%ldim_pd,self%smpd)
                call self%cavgs_merged(istate,icls)%new(self%ldim_pd,self%smpd)
                self%cavgs_even(istate,icls)       = cmplx(0.,0.)
                self%cavgs_odd(istate,icls)        = cmplx(0.,0.)
                self%cavgs_merged(istate,icls)     = cmplx(0.,0.)
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
            if( allocated(self%precs(iprec)%classes) )then
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
            endif
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
            if( allocated(self%precs(iprec)%classes) )then
                sz = size(self%precs(iprec)%classes)
                allocate(l_state_class(sz))
                where( self%precs(iprec)%states .eq. state .and. self%precs(iprec)%classes .eq. class )
                    l_state_class = .true.
                else where
                    l_state_class = .false.
                endwhere
                pop = pop + count(l_state_class)
                deallocate(l_state_class)
            endif
        end do
    end function class_pop

    !>  \brief  is for getting a class average
    subroutine get_cavg( self, class, which, img, state )
        class(classaverager), intent(inout) :: self
        integer,              intent(in)    :: class
        character(len=*),     intent(in)    :: which
        class(image),         intent(inout) :: img
        integer, optional,    intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                img = self%cavgs_even(sstate,class)
            case('odd')
                img = self%cavgs_odd(sstate,class)
            case('merged')
                img = self%cavgs_merged(sstate,class)
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: get_cavg'
        end select
    end subroutine get_cavg

    !>  \brief  is for setting a class average
    subroutine set_cavg( self, class, which, img, state )
        class(classaverager), intent(inout) :: self
        integer,              intent(in)    :: class
        character(len=*),     intent(in)    :: which
        class(image),         intent(in)    :: img
        integer, optional,    intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                self%cavgs_even(sstate,class)   = img
            case('odd')
                self%cavgs_odd(sstate,class)    = img
            case('merged')
                self%cavgs_merged(sstate,class) = img
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: set_cavg'
        end select
    end subroutine set_cavg

    ! calculators

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !!          using gridding interpolation in Fourier space or quadratic 
    !!          interpolation in real-space      
    subroutine assemble_sums( self )
        use simple_kbinterpol,      only: kbinterpol
        use simple_prep4cgrid,      only: prep4cgrid
        use simple_map_reduce,      only: split_nobjs_even
        use simple_hadamard_common, only: read_img_from_stk
        use simple_kbinterpol_mem,  only: kbinterpol_mem
        class(classaverager), intent(inout) :: self
        type(kbinterpol)         :: kbwin
        type(kbinterpol_mem)     :: kbmem
        type(prep4cgrid)         :: gridprep
        type(image)              :: batch_imgsum_even, batch_imgsum_odd
        type(image)              :: batch_rhosum_even, batch_rhosum_odd
        type(image), allocatable :: padded_imgs(:), batch_imgs(:)
        complex,     allocatable :: cmat_even(:,:), cmat_odd(:,:)
        real,        allocatable :: w(:,:), rho_even(:,:), rho_odd(:,:)
        integer,     allocatable :: ptcls_inds(:), batches(:,:), iprecs(:)
        integer,     allocatable :: ioris(:)
        complex   :: comp, zero, oshift
        real      :: loc(2), mat(2,2), winsz, pw, tval, tvalsq, vec(2)
        integer   :: cnt_progress, nbatches, batch, icls_pop, iprec, iori, i, batchsz, fnr, icls_popmax
        integer   :: lims(3,2), istate, sh, nyq, logi(3), phys(3), win(2,2), win4w(2,2), win_tst(2,2)
        integer   :: cyc_lims(3,2), alloc_stat, wdim, h, k, hh, kk, icls, iptcl, batchsz_max
        if( .not. self%pp%l_distr_exec ) write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! init
        call self%init_cavgs_sums
        cyc_lims  = self%ctfsqsums_even(1,1)%loop_lims(3)
        lims      = cyc_lims
        lims(1,1) = 0
        nyq       = self%ctfsqsums_even(1,1)%get_lfny(1)
        kbwin     = kbinterpol(KBWINSZ, KBALPHA)
        zero      = cmplx(0.,0.)
        winsz     = KBWINSZ
        wdim      = ceiling(KBALPHA*winsz) + 1
        allocate(cmat_even(cyc_lims(1,1):cyc_lims(1,2),cyc_lims(2,1):cyc_lims(2,2)),&
                &cmat_odd(cyc_lims(1,1):cyc_lims(1,2),cyc_lims(2,1):cyc_lims(2,2)),&
                &rho_even(cyc_lims(1,1):cyc_lims(1,2),cyc_lims(2,1):cyc_lims(2,2)),&
                &rho_odd(cyc_lims(1,1):cyc_lims(1,2),cyc_lims(2,1):cyc_lims(2,2)),&
                &w(wdim, wdim), stat=alloc_stat)
        call alloc_errchk('assemble_sums; simple_classaverager', alloc_stat)
        call batch_imgsum_even%new(self%ldim_pd, self%pp%smpd)
        call batch_imgsum_odd%new(self%ldim_pd, self%pp%smpd)
        call batch_rhosum_even%new(self%ldim_pd, self%pp%smpd)
        call batch_rhosum_odd%new(self%ldim_pd, self%pp%smpd)
        call gridprep%new(self%bp%img, kbwin)
        if( self%use_kbmem )then
            ! memoize K-B kernel
            call kbmem%new(self%inpl_rots, cyc_lims(1:2,1:2))
            call kbmem%memoize_kb(kbwin)
        endif
        if( L_BENCH )then
            rt_batch_loop = 0.
            rt_gridding   = 0.
            rt_tot        = 0.
            t_tot         = tic()
        endif
        ! find maximum batch size
        batchsz_max = 0
        do istate=1,self%nstates
            ! class loop
            do icls=1,self%ncls
                ! batch planning
                nbatches = ceiling(real(icls_pop)/real(self%pp%nthr*BATCHTHRSZ))
                batches  = split_nobjs_even(icls_pop, nbatches)
                ! batch loop, prep
                do batch=1,nbatches
                    ! prep batch
                    batchsz = batches(batch,2) - batches(batch,1) + 1
                    if( batchsz > batchsz_max ) batchsz_max = batchsz
                end do
            end do
        end do
        ! pre-create images based on maximum batchsz
        allocate(padded_imgs(batchsz_max), batch_imgs(batchsz_max))
        do i=1,batchsz_max
            call batch_imgs(i)%new(self%ldim, self%pp%smpd)
            call padded_imgs(i)%new(self%ldim_pd, self%pp%smpd, wthreads=.false.)
        end do
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
                    ! batch particles loop
                    if( L_BENCH ) t_batch_loop = tic()
                    do i=1,batchsz
                        iptcl = ptcls_inds(batches(batch,1) + i - 1)
                        call read_img_from_stk( self%bp, self%pp, iptcl )
                        batch_imgs(i)  = self%bp%img
                        padded_imgs(i) = 0.
                    enddo
                    if( L_BENCH ) rt_batch_loop = rt_batch_loop + toc(t_batch_loop)
                    if( L_BENCH ) t_gridding = tic()
                    cmat_even = zero
                    cmat_odd  = zero
                    rho_odd   = 0.
                    rho_even  = 0.
                    !$omp parallel do default(shared) schedule(static) reduction(+:cmat_even,cmat_odd,rho_even,rho_odd) proc_bind(close)&
                    !$omp private(sh,i,iprec,iori,h,k,loc,mat,logi,phys,w,win4w,win,pw,tval,tvalsq,vec,hh,kk,oshift,comp,win_tst)
                    ! batch loop, convolution interpolation
                    do i=1,batchsz
                        call gridprep%prep(batch_imgs(i), padded_imgs(i))
                        iprec = iprecs(batches(batch,1) + i - 1)
                        iori  = ioris(batches(batch,1)  + i - 1)
                        ! prep weight
                        if( self%l_hard_assign )then
                            pw = self%precs(iprec)%pw
                        else
                            pw = self%precs(iprec)%pw * self%precs(iprec)%ows(iori)
                        endif
                        mat = rotmat2d( self%precs(iprec)%e3s(iori) )
                        ! Fourier components loop
                        do h=cyc_lims(1,1),cyc_lims(1,2)
                            do k=cyc_lims(2,1),cyc_lims(2,2)
                                sh = nint(hyp(real(h),real(k)))
                                if( sh > nyq + 1 )cycle
                                ! fetch component & calc shift
                                logi     = [h,k,0]
                                phys     = padded_imgs(i)%comp_addr_phys(logi)
                                comp     = padded_imgs(i)%get_fcomp(logi, phys)
                                oshift   = padded_imgs(i)%oshift(logi, [-self%precs(iprec)%shifts(iori,1),-self%precs(iprec)%shifts(iori,2),0.])
                                ! rotation
                                vec      = real([h,k])
                                loc      = matmul(vec,mat)
                                ! kernel limits
                                win = sqwin_2d(loc(1),loc(2), winsz, cyc_lims(1:2,1:2))
                                ! evaluate the transfer function
                                call self%calc_tfun_vals(iprec, vec, tval, tvalsq)
                                ! kernel
                                win4w(1,:) = win(1,:) - win(1,1) + 1
                                win4w(2,:) = win(2,:) - win(2,1) + 1
                                if( self%use_kbmem )then
                                    call kbmem%fetch(self%precs(iprec)%inpl_inds(iori), h, k, w)
                                    w = w * pw
                                else
                                    w = pw
                                    do hh=win4w(1,1),win4w(1,2)
                                        w(hh,:) = w(hh,:) * kbwin%apod( real(win(1,1)-1 + hh)-loc(1) )
                                    end do
                                    do kk=win4w(2,1),win4w(2,2)
                                        w(:,kk) = w(:,kk) * kbwin%apod( real(win(2,1)-1 + kk)-loc(2) )
                                    end do
                                endif
                                ! summation
                                select case(self%precs(iprec)%eo)
                                case(0,-1)
                                    cmat_even(win(1,1):win(1,2),win(2,1):win(2,2)) = cmat_even(win(1,1):win(1,2),win(2,1):win(2,2)) + &
                                        &(w(win4w(1,1):win4w(1,2),win4w(2,1):win4w(2,2))*tval*comp)*oshift
                                    rho_even(win(1,1):win(1,2),win(2,1):win(2,2)) = rho_even(win(1,1):win(1,2),win(2,1):win(2,2)) + &
                                        &w(win4w(1,1):win4w(1,2),win4w(2,1):win4w(2,2))*tvalsq
                                case(1)
                                    cmat_odd(win(1,1):win(1,2),win(2,1):win(2,2)) = cmat_odd(win(1,1):win(1,2),win(2,1):win(2,2)) + &
                                        &(w(win4w(1,1):win4w(1,2),win4w(2,1):win4w(2,2))*tval*comp)*oshift
                                    rho_odd(win(1,1):win(1,2),win(2,1):win(2,2)) = rho_odd(win(1,1):win(1,2),win(2,1):win(2,2)) + &
                                        &w(win4w(1,1):win4w(1,2),win4w(2,1):win4w(2,2))*tvalsq
                                end select
                            end do
                        end do
                    enddo
                    !$omp end parallel do
                    if( L_BENCH ) rt_gridding = rt_gridding + toc(t_gridding)
                    ! compression
                    batch_imgsum_even = zero
                    batch_imgsum_odd  = zero
                    batch_rhosum_even = zero
                    batch_rhosum_odd  = zero
                    !$omp parallel do collapse(2) default(shared) private(h,k,logi,phys)&
                    !$omp schedule(static) proc_bind(close)                    
                    do h=lims(1,1),lims(1,2)
                        do k=lims(2,1),lims(2,2)
                            logi = [h,k,0]
                            phys = batch_imgsum_even%comp_addr_phys([h,k,0])
                            call batch_imgsum_even%set_fcomp(logi, phys, cmat_even(h,k))
                            call batch_imgsum_odd%set_fcomp( logi, phys, cmat_odd(h,k))
                            call batch_rhosum_even%set_fcomp(logi, phys, cmplx(rho_even(h,k),0.))
                            call batch_rhosum_odd%set_fcomp( logi, phys, cmplx(rho_odd(h,k), 0.))
                        enddo
                    enddo
                    !$omp end parallel do

                    ! THIS CANNOT BE THREADED THIS WAY BECAUSE IT SCREWS UP THE FRIEDEL SYMMETRY
                    ! $omp parallel do collapse(2) default(shared) private(h,k,phys)&
                    ! $omp schedule(static) proc_bind(close)                    
                    ! do h=cyc_lims(1,1),cyc_lims(1,2)
                    !     do k=cyc_lims(2,1),cyc_lims(2,2)
                    !         phys = batch_imgsum_even%comp_addr_phys([h,k,0])
                    !         call batch_imgsum_even%set_cmat_at(phys, cmat_even(h,k))
                    !         call batch_imgsum_odd%set_cmat_at(phys,  cmat_odd(h,k))
                    !         call batch_rhosum_even%set_cmat_at(phys, cmplx(rho_even(h,k),0.))
                    !         call batch_rhosum_odd%set_cmat_at(phys,  cmplx(rho_odd(h,k),0.))
                    !     enddo
                    ! enddo
                    ! $omp end parallel do

                    ! batch summation
                    call self%cavgs_even(istate,icls)%add( batch_imgsum_even)
                    call self%cavgs_odd(istate,icls)%add( batch_imgsum_odd)
                    call self%ctfsqsums_even(istate,icls)%add( batch_rhosum_even)
                    call self%ctfsqsums_odd(istate,icls)%add( batch_rhosum_odd)
                    
                enddo
                ! class cleanup
                deallocate(ptcls_inds, batches, iprecs, ioris)
            enddo
        enddo
        ! batch cleanup
        do i=1,batchsz_max
            call padded_imgs(i)%kill
            call batch_imgs(i)%kill
        enddo
        deallocate(padded_imgs, batch_imgs)
        if( self%use_kbmem ) call kbmem%kill
        call batch_imgsum_even%kill
        call batch_imgsum_odd%kill
        call batch_rhosum_even%kill
        call batch_rhosum_odd%kill
        deallocate(cmat_even, cmat_odd, rho_even, rho_odd, w)
        if( .not. self%pp%l_distr_exec ) call self%merge_eos_and_norm
        if( L_BENCH )then
            rt_tot = rt_tot + toc(t_tot)
            benchfname = 'CLASSAVERAGER_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'batch loop : ', rt_batch_loop
            write(fnr,'(a,1x,f9.2)') 'gridding   : ', rt_gridding
            write(fnr,'(a,1x,f9.2)') 'total time : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** REATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'batch loop : ', (rt_batch_loop/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding   : ', (rt_gridding/rt_tot)   * 100.
            call fclose(fnr)
        endif
    end subroutine assemble_sums

    !>  \brief  merges the even/odd pairs and normalises the sums
    subroutine merge_eos_and_norm( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                ! merging
                self%cavgs_merged(istate,icls) = cmplx(0.,0.)
                call self%cavgs_merged(istate,icls)%add(self%cavgs_even(istate,icls))
                call self%cavgs_merged(istate,icls)%add(self%cavgs_odd(istate,icls))
                self%ctfsqsums_merged(istate,icls) = cmplx(0.,0.)
                call self%ctfsqsums_merged(istate,icls)%add(self%ctfsqsums_even(istate,icls))
                call self%ctfsqsums_merged(istate,icls)%add(self%ctfsqsums_odd(istate,icls))
                ! (w*CTF)**2 density correction
                call self%cavgs_even(istate,icls)%ctf_dens_correct(self%ctfsqsums_even(istate,icls))
                call self%cavgs_even(istate,icls)%bwd_ft
                call self%cavgs_even(istate,icls)%clip_inplace(self%ldim)
                call self%cavgs_odd(istate,icls)%ctf_dens_correct(self%ctfsqsums_odd(istate,icls))
                call self%cavgs_odd(istate,icls)%bwd_ft
                call self%cavgs_odd(istate,icls)%clip_inplace(self%ldim)
                call self%cavgs_merged(istate,icls)%ctf_dens_correct(self%ctfsqsums_merged(istate,icls))
                call self%cavgs_merged(istate,icls)%bwd_ft
                call self%cavgs_merged(istate,icls)%clip_inplace(self%ldim)
            end do
         end do
    end subroutine merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    subroutine calc_and_write_frcs( self, fname )
        class(classaverager), intent(inout) :: self
        character(len=*),     intent(in)    :: fname
        type(image)       :: even_img, odd_img
        real, allocatable :: res(:), frc(:)
        integer           :: istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                even_img = self%cavgs_even(istate,icls)
                odd_img  = self%cavgs_odd(istate,icls)
                call even_img%norm
                call odd_img%norm
                if( self%pp%l_innermsk )then
                    call even_img%mask(self%pp%msk, 'soft', inner=self%pp%inner, width=self%pp%width)
                    call odd_img%mask(self%pp%msk, 'soft', inner=self%pp%inner, width=self%pp%width)
                else
                    call even_img%mask(self%pp%msk, 'soft')
                    call odd_img%mask(self%pp%msk, 'soft')
                endif
                call even_img%fwd_ft
                call odd_img%fwd_ft
                call even_img%fsc(odd_img, res, frc)
                call self%bp%projfrcs%set_frc(icls, frc, istate)
                call even_img%kill
                call odd_img%kill
            end do
        end do
        call self%bp%projfrcs%write(fname)
    end subroutine calc_and_write_frcs

    subroutine calc_tfun_vals( self, iprec, vec, tval, tvalsq )
        class(classaverager), intent(inout) :: self         !< instance  
        integer,              intent(in)    :: iprec        !< particle record
        real,                 intent(in)    :: vec(2)       !< nonuniform sampling location
        real,                 intent(out)   :: tval, tvalsq !< CTF and CTF**2.
        real :: sqSpatFreq, ang, inv1, inv2
        if( self%ctf%flag /= CTFFLAG_NO )then
            inv1       = vec(1)*(1./real(self%ldim_pd(1)))
            inv2       = vec(2)*(1./real(self%ldim_pd(2)))
            sqSpatFreq = inv1*inv1+inv2*inv2
            ang        = atan2(vec(2), vec(1))
            ! calculate CTF and CTF**2 values
            if( self%phaseplate )then
                tval = self%precs(iprec)%tfun%eval(sqSpatFreq, self%precs(iprec)%dfx,&
                    &self%precs(iprec)%dfy, self%precs(iprec)%angast, ang, self%precs(iprec)%phshift)
            else
                tval = self%precs(iprec)%tfun%eval(sqSpatFreq, self%precs(iprec)%dfx,&
                    &self%precs(iprec)%dfy, self%precs(iprec)%angast, ang)
            endif
            tvalsq = tval * tval
            if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
        else
            tval   = 1.
            tvalsq = tval
        endif
    end subroutine calc_tfun_vals

    !> \brief average low-resolution info between eo pairs
    subroutine eoavg( self )
        class(classaverager), intent(inout) :: self
        integer :: find, istate, icls
        do istate=1,self%nstates
            do icls=1,self%ncls
                find = self%bp%projfrcs%estimate_find_for_eoavg(icls, istate)
                call self%cavgs_merged(istate,icls)%fwd_ft
                call self%cavgs_even(istate,icls)%fwd_ft
                call self%cavgs_odd(istate,icls)%fwd_ft
                call self%cavgs_even(istate,icls)%insert_lowres(self%cavgs_merged(istate,icls), find)
                call self%cavgs_odd(istate,icls)%insert_lowres(self%cavgs_merged(istate,icls), find)
                call self%cavgs_merged(istate,icls)%bwd_ft
                call self%cavgs_even(istate,icls)%bwd_ft
                call self%cavgs_odd(istate,icls)%bwd_ft
            end do
        end do
    end subroutine eoavg

    ! I/O

    !>  \brief  writes class averages to disk
    subroutine write( self, fname, which, state )
        class(classaverager), intent(inout) :: self
        character(len=*),     intent(in)    :: fname, which
        integer, optional,    intent(in)    :: state
        integer :: sstate, icls
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                do icls=1,self%ncls
                    call self%cavgs_even(sstate, icls)%write(fname, icls)
                end do
            case('odd')
                do icls=1,self%ncls
                    call self%cavgs_odd(sstate, icls)%write(fname, icls)
                end do
            case('merged')
                 do icls=1,self%ncls
                    call self%cavgs_merged(sstate, icls)%write(fname, icls)
                end do
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: get_cavg'
        end select
    end subroutine write

    !>  \brief  reads class averages from disk
    subroutine read( self, fname, which, state )
        class(classaverager), intent(inout) :: self
        character(len=*),     intent(in)    :: fname, which
        integer, optional,    intent(in)    :: state
        integer :: sstate, icls
        if( .not. file_exists(fname) )then
            write(*,*) 'file does not exist in cwd: ', trim(fname)
            stop 'simple_classaverager :: read'
        endif
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                do icls=1,self%ncls
                    call self%cavgs_even(sstate,icls)%new(self%ldim,self%smpd)
                    call self%cavgs_even(sstate, icls)%read(fname, icls)
                end do
            case('odd')
                do icls=1,self%ncls
                call self%cavgs_odd(sstate,icls)%new(self%ldim,self%smpd)
                    call self%cavgs_odd(sstate, icls)%read(fname, icls)
                end do
            case('merged')
                 do icls=1,self%ncls
                    call self%cavgs_merged(sstate,icls)%new(self%ldim,self%smpd)
                    call self%cavgs_merged(sstate, icls)%read(fname, icls)
                end do
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: get_cavg'
        end select
    end subroutine read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine write_partial_sums( self )
        class(classaverager), intent(inout) :: self
        integer :: istate, icls
        character(len=:), allocatable :: cae, cao, cte, cto
        do istate=1,self%nstates
            if( self%nstates > 1 )then
                allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
            else
                allocate(cae, source='cavgs_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cao, source='cavgs_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cte, source='ctfsqsums_even_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
                allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(self%pp%part,self%pp%numlen)//self%pp%ext)
            endif
            do icls=1,self%ncls
                call self%cavgs_even(istate, icls)%write(cae, icls)
                call self%cavgs_odd(istate, icls)%write(cao, icls)
                call self%ctfsqsums_even(istate, icls)%write(cte, icls)
                call self%ctfsqsums_odd(istate, icls)%write(cto, icls)
            end do
            deallocate(cae, cao, cte, cto)
        end do
    end subroutine write_partial_sums

    !>  \brief  re-generates the object after distributed execution
    subroutine assemble_sums_from_parts( self )
        class(classaverager), intent(inout) :: self
        type(image)                   :: img4read
        character(len=:), allocatable :: cae, cao, cte, cto
        integer :: ipart, istate, icls
        call self%init_cavgs_sums
        call img4read%new(self%ldim_pd, self%smpd)
        call img4read%set_ft(.true.)
        do istate=1,self%nstates
            do ipart=1,self%pp%nparts
                if( self%nstates > 1 )then
                    allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                else
                    allocate(cae, source='cavgs_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cao, source='cavgs_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cte, source='ctfsqsums_even_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                    allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(ipart,self%pp%numlen)//self%pp%ext)
                endif
                if( file_exists(cae) )then
                    do icls=1,self%ncls
                        call img4read%read(cae, icls)
                        call self%cavgs_even(istate,icls)%add(img4read)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cae)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                if( file_exists(cao) )then
                    do icls=1,self%ncls
                        call img4read%read(cao, icls)
                        call self%cavgs_odd(istate,icls)%add(img4read)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cao)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                if( file_exists(cte) )then
                    do icls=1,self%ncls
                        call img4read%read(cte, icls)
                        call self%ctfsqsums_even(istate,icls)%add(img4read)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cte)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                if( file_exists(cto) )then
                    do icls=1,self%ncls
                        call img4read%read(cto, icls)
                        call self%ctfsqsums_odd(istate,icls)%add(img4read)
                    end do
                else
                    write(*,*) 'File does not exists: ', trim(cto)
                    stop 'In: simple_classaverager :: assemble_sums_from_parts'
                endif
                deallocate(cae, cao, cte, cto )
            end do
        end do
        call img4read%kill
        call self%merge_eos_and_norm()
    end subroutine assemble_sums_from_parts

    ! destructor

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
            deallocate( self%cavgs_even, self%cavgs_odd, self%cavgs_merged,&
            &self%ctfsqsums_even, self%ctfsqsums_odd, self%ctfsqsums_merged)
            do iprec=1,self%partsz
                if( allocated(self%precs(iprec)%classes) ) deallocate(self%precs(iprec)%classes)
                if( allocated(self%precs(iprec)%states)  ) deallocate(self%precs(iprec)%states)
                if( allocated(self%precs(iprec)%ows)     ) deallocate(self%precs(iprec)%ows)
                if( allocated(self%precs(iprec)%e3s)     ) deallocate(self%precs(iprec)%e3s)
                if( allocated(self%precs(iprec)%shifts)  ) deallocate(self%precs(iprec)%shifts)
            end do
            deallocate(self%precs)
            if( allocated(self%inpl_rots) ) deallocate(self%inpl_rots)
            self%istart        = 0
            self%iend          = 0
            self%partsz        = 0
            self%ncls          = 0
            self%nstates       = 0
            self%l_is_class    = .true.
            self%l_hard_assign = .true.
            self%exists        = .false.
        endif
    end subroutine kill

end module simple_classaverager
