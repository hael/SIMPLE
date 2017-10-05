! 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
!$ use omp_lib
!$ use omp_lib_kinds
use, intrinsic :: iso_c_binding
#include "simple_lib.f08"
!! import functions
use simple_fftw3,      only: fftwf_alloc_real, fftwf_free
use simple_imghead,    only: find_ldim_nptcls
use simple_timer,      only: tic, toc, timer_int_kind
use simple_gridding,   only: prep4cgrid
!! import classes
use simple_ctf,        only: ctf
use simple_ori,        only: ori
use simple_oris,       only: oris
use simple_params,     only: params
use simple_sym,        only: sym
use simple_kbinterpol, only: kbinterpol
use simple_kbfast,     only: kbfast
use simple_image,      only: image
implicit none

public :: reconstructor
private

!! CTF flag type
enum, bind(c) 
    enumerator :: CTFFLAG_NO = 0, CTFFLAG_YES = 1, CTFFLAG_MUL = 2,  CTFFLAG_FLIP = 3
end enum

type :: CTFFLAGTYPE
    private
    integer(kind(CTFFLAG_NO)) :: flag=CTFFLAG_NO
end type CTFFLAGTYPE

type, extends(image) :: reconstructor
    private
    type(kbinterpol)            :: kbwin                        !< window function object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    type(ctf)                   :: tfun                         !< CTF object
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable        :: cmat_exp(:,:,:)              !< Fourier components of expanded reconstructor
    real,    allocatable        :: rho_exp(:,:,:)               !< sampling+CTF**2 density of expanded reconstructor
    real                        :: winsz         = 1.           !< window half-width
    real                        :: alpha         = 2.           !< oversampling ratio
    integer                     :: wdim          = 0            !< dim of interpolation matrix
    integer                     :: lfny          = 0            !< Nyqvist Fourier index
    integer                     :: ldim_img(3)   = 0            !< logical dimension of the original image
    type(CTFFLAGTYPE)           :: ctf                          !< ctf flag <yes|no|mul|flip>
    logical                     :: tfneg              = .false. !< invert contrast or not
    logical                     :: tfastig            = .false. !< astigmatic CTF or not
    logical                     :: rho_allocated      = .false. !< existence of rho matrix
    logical                     :: rho_exp_allocated  = .false. !< existence of rho expanded matrix
    logical                     :: cmat_exp_allocated = .false. !< existence of Fourier components expanded matrix
  contains
    ! CONSTRUCTORS 
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    ! GETTER
    procedure          :: get_kbwin
    ! I/O
    procedure          :: write_rho
    procedure          :: read_rho
    ! INTERPOLATION
    procedure, private :: inout_fcomp
    procedure, private :: calc_tfun_vals
    procedure          :: inout_fplane
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    ! SUMMATION
    procedure          :: sum
    ! RECONSTRUCTION
    procedure          :: rec
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

real            :: dfx=0., dfy=0., angast=0.
real, parameter :: SHTHRESH=0.0001
#include "simple_local_flags.inc"

contains

    ! CONSTRUCTOR

    subroutine alloc_rho( self, p, expand )
        class(reconstructor), intent(inout) :: self   !< this instance
        class(params),        intent(in)    :: p      !< parameters object
        logical, optional,    intent(in)    :: expand !< expand flag
        integer :: rho_shape(3), lims(3,2), rho_exp_lims(3,2), ldim_exp(3,2)
        integer :: dim, maxlims
        logical :: l_expand
        l_expand = .true.
        if(.not. self%exists()   ) stop 'construct image before allocating rho; alloc_rho; simple_reconstructor'
        if(      self%is_2d()    ) stop 'only for volumes; alloc_rho; simple_reconstructor'
        if( present(expand) )l_expand = expand
        self%ldim_img = self%get_ldim()
        if( self%ldim_img(3) < 2 ) stop 'reconstructor need to be 3D 4 now; alloc_rho; simple_reconstructor'
        call self%dealloc_rho
        self%winsz = p%winsz
        self%wdim  = 2*ceiling(self%winsz) + 1
        self%alpha = p%alpha
        select case(p%ctf)
            case('no')
                self%ctf%flag = CTFFLAG_NO
            case('yes')
                self%ctf%flag = CTFFLAG_YES
            case('mul')
                self%ctf%flag = CTFFLAG_MUL
            case('flip')
                self%ctf%flag = CTFFLAG_FLIP
        end select
        self%tfneg = .false.
        if( p%neg .eq. 'yes' ) self%tfneg = .true.
        self%tfastig = .false.
        if( p%tfplan%mode .eq. 'astig' ) self%tfastig = .true.
        self%kbwin = kbinterpol(self%winsz,self%alpha)
        ! Work out dimensions of the rho array
        rho_shape(1)   = fdim(self%ldim_img(1))
        rho_shape(2:3) = self%ldim_img(2:3)
        ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
        self%kp = fftwf_alloc_real(int(product(rho_shape),c_size_t))
        ! Set up the rho array which will point at the allocated memory
        call c_f_pointer(self%kp,self%rho,rho_shape)
        self%rho_allocated = .true.
        if( l_expand )then
            ! setup expanded matrices
            lims    = self%loop_lims(2)
            maxlims = maxval(abs(lims))
            dim     = ceiling(sqrt(2.)*maxlims + self%winsz)
            ldim_exp(1,:) = [-dim, dim]
            ldim_exp(2,:) = [-dim, dim]
            ldim_exp(3,:) = [-dim, dim]
            rho_exp_lims  = ldim_exp
            allocate(self%cmat_exp( ldim_exp(1,1):ldim_exp(1,2),&
                &ldim_exp(2,1):ldim_exp(2,2),&
                &ldim_exp(3,1):ldim_exp(3,2)), stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor cmat_exp")
            allocate(self%rho_exp( rho_exp_lims(1,1):rho_exp_lims(1,2),&
                &rho_exp_lims(2,1):rho_exp_lims(2,2),&
                &rho_exp_lims(3,1):rho_exp_lims(3,2)), stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor rho_exp")
            self%cmat_exp           = cmplx(0.,0.)
            self%rho_exp            = 0.
            self%cmat_exp_allocated = .true.
            self%rho_exp_allocated  = .true.
        end if       
        call self%reset
    end subroutine alloc_rho

    ! SETTERS

    ! resets the reconstructor object before reconstruction
    subroutine reset( self )
        class(reconstructor), intent(inout) :: self !< this instance
        self     = cmplx(0.,0.)
        self%rho = 0.
    end subroutine reset

    ! GETTERS
    !> get the kbintpol window
    function get_kbwin( self ) result( wf )
        class(reconstructor), intent(inout) :: self !< this instance
        type(kbinterpol) :: wf   !< return kbintpol window
        wf = kbinterpol(self%winsz,self%alpha)
    end function get_kbwin

    ! I/O
    !>Write reconstructed image
    subroutine write_rho( self, kernam )
        class(reconstructor), intent(in) :: self   !< this instance
        character(len=*),     intent(in) :: kernam !< kernel name
        integer :: filnum, ierr
        call del_file(trim(kernam))
        call fopen(filnum, trim(kernam), status='NEW', action='WRITE', access='STREAM', iostat=ierr)
        call fileio_errmsg( 'simple_reconstructor ; write rho '//trim(kernam), ierr)
        write(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileio_errmsg('read_rho; simple_reconstructor writing '//trim(kernam), ierr)
        call fclose(filnum,errmsg='simple_reconstructor ; write rho  fclose ')
    end subroutine write_rho

    !> Read sampling density matrix
    subroutine read_rho( self, kernam )
        class(reconstructor), intent(inout) :: self !< this instance
        character(len=*),     intent(in)    :: kernam !< kernel name
        integer :: filnum, ierr
        call fopen(filnum, file=trim(kernam), status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileio_errmsg('read_rho; simple_reconstructor opening '//trim(kernam), ierr)
        read(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileio_errmsg('simple_reconstructor::read_rho; simple_reconstructor reading '&
            &// trim(kernam), ierr)
        call fclose(filnum,errmsg='read_rho; simple_reconstructor closing '//trim(kernam))
    end subroutine read_rho

    ! INTERPOLATION

    subroutine inout_fcomp( self, h, k, e, inoutmode, comp, oshift, pwght)
        class(reconstructor), intent(inout) :: self      !< this object
        integer,              intent(in)    :: h, k      !< Fourier indices
        class(ori),           intent(inout) :: e         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        complex,              intent(in)    :: comp      !< input component, if not given only sampling density calculation
        complex,              intent(in)    :: oshift    !< origin shift
        real,                 intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        real    :: w(self%wdim,self%wdim,self%wdim)
        real    :: vec(3), loc(3), tval, tvalsq
        integer :: i, win(3,2)
        ! calculate non-uniform sampling location
        vec = [real(h), real(k), 0.]
        loc = matmul(vec, e%get_mat())
        ! evaluate the transfer function
        call self%calc_tfun_vals(vec, tval, tvalsq)
        ! initiate kernel matrix
        win = sqwin_3d(loc(1), loc(2), loc(3), self%winsz)
        ! (weighted) kernel values
        w = pwght
        do i=1,self%wdim
            w(i,:,:) = w(i,:,:) * self%kbwin%apod(real(win(1,1)+i-1)-loc(1))
            w(:,i,:) = w(:,i,:) * self%kbwin%apod(real(win(2,1)+i-1)-loc(2))
            w(:,:,i) = w(:,:,i) * self%kbwin%apod(real(win(3,1)+i-1)-loc(3))
        enddo
        ! expanded matrices update
        if( inoutmode )then
            ! addition
            ! CTF and w modulates the component before origin shift
            self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + (comp*tval*w)*oshift
            self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + tvalsq*w
        else
            ! substraction
            ! CTF and w modulates the component before origin shift
            self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - (comp*tval*w)*oshift
            self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - tvalsq*w
        endif
    end subroutine inout_fcomp

    subroutine calc_tfun_vals( self, vec, tval, tvalsq )
        class(reconstructor), intent(inout) :: self         !< instance
        real,                 intent(in)    :: vec(3)       !< nonuniform sampling location
        real,                 intent(out)   :: tval, tvalsq !< CTF and CTF**2.
        real    :: sqSpatFreq,ang,inv1,inv2
        integer :: ldim(3)
        if( self%ctf%flag /= CTFFLAG_NO )then
            ldim       = self%get_ldim()
            inv1       = vec(1)*(1./real(ldim(1)))
            inv2       = vec(2)*(1./real(ldim(2)))
            sqSpatFreq = inv1*inv1+inv2*inv2
            ang        = atan2(vec(2), vec(1))
            ! calculate CTF and CTF**2 values
            ! dfx, dfy, angast are class vars that are set by inout fplane
            tval   = self%tfun%eval(sqSpatFreq, dfx, dfy, angast, ang) ! no bfactor 4 now
            tvalsq = tval * tval
            if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
            if( self%ctf%flag == CTFFLAG_MUL  ) tval = 1.
            if( self%tfneg ) tval = -tval
        else
            tval   = 1.
            tvalsq = tval
            if( self%tfneg ) tval = -tval
        endif
    end subroutine calc_tfun_vals

    !> insert or uninsert Fourier plane
    subroutine inout_fplane( self, o, inoutmode, fpl, pwght )
        class(reconstructor), intent(inout) :: self      !< instance
        class(ori),           intent(inout) :: o         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        class(image),         intent(inout) :: fpl       !< Fourier plane
        real,                 intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        integer :: h, k, lims(3,2), logi(3), phys(3)
        complex :: oshift
        real    :: x, y, xtmp
        if( self%ctf%flag /= CTFFLAG_NO )then ! make CTF object & get CTF info
            self%tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            dfx  = o%get('dfx')
            if( self%tfastig )then ! astigmatic CTF model
                dfy = o%get('dfy')
                angast = o%get('angast')
            else                   ! non-astigmatic CTF model
                dfy = dfx
                angast = 0.
            endif
        endif
        lims = self%loop_lims(2)
        x    = o%get('x')
        y    = o%get('y')
        !$omp parallel do collapse(2) default(shared) schedule(static)&
        !$omp private(h,k,oshift,logi,phys) proc_bind(close)
        do k=lims(1,1),lims(1,2)
            do h=lims(1,1),lims(1,2)
                logi   = [h,k,0]
                oshift = fpl%oshift(logi, [-x,-y,0.])
                phys   = fpl%comp_addr_phys(logi)
                call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp(logi,phys),oshift,pwght)
            end do
        end do
        !$omp end parallel do
    end subroutine inout_fplane

    !>  is for uneven distribution of orientations correction 
    !>  from Pipe & Menon 1999
    subroutine sampl_dens_correct( self, maxits )
        use simple_gridding,   only: mul_w_instr
        class(reconstructor), intent(inout) :: self
        integer,    optional, intent(in)    :: maxits
        type(kbinterpol)     :: kbwin 
        type(image)          :: W_img, Wprev_img
        real                 :: Wnorm, val_prev, val, Wnorm_prev, invrho
        integer              :: h, k, m, lims(3,2),  phys(3), iter
        integer              :: maxits_here
        complex, parameter   :: zero = cmplx(0.,0.), one = cmplx(1.,0.)
        real,    parameter   :: winsz  = 2.
        maxits_here = GRIDCORR_MAXITS
        if( present(maxits) )maxits_here = maxits
        lims = self%loop_lims(2)
        call W_img%new(self%ldim_img, self%get_smpd())
        call Wprev_img%new(self%ldim_img, self%get_smpd())
        call W_img%set_ft(.true.)
        call Wprev_img%set_ft(.true.)
        ! kernel
        kbwin = kbinterpol(winsz, self%alpha)
        ! weights init to 1.
        W_img = one
        Wnorm = SMALL
        do iter = 1, maxits_here
            Wnorm_prev = Wnorm
            Wprev_img  = W_img 
            ! W <- W * rho
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys) proc_bind(close)
            do h = lims(1,1),lims(1,2)
                do k = lims(2,1),lims(2,2)
                    do m = lims(3,1),lims(3,2)
                        phys  = W_img%comp_addr_phys([h,k,m])
                        call W_img%mul_cmat_at(self%rho(phys(1),phys(2),phys(3)), phys)
                    end do
                end do
            end do
            !$omp end parallel do
            ! W <- (W / rho) x kernel
            call W_img%bwd_ft
            call mul_w_instr(W_img, kbwin)
            call W_img%fwd_ft
            ! W <- Wprev / ((W/ rho) x kernel)
            Wnorm = 0.
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys,val,val_prev) proc_bind(close)&
            !$omp reduction(+:Wnorm)
            do h = lims(1,1),lims(1,2)
                do k = lims(2,1),lims(2,2)
                    do m = lims(3,1),lims(3,2)
                        phys     = W_img%comp_addr_phys([h, k, m])
                        val      = mycabs(W_img%get_cmat_at(phys))   !! ||C|| == ||C*||
                        val_prev = real(Wprev_img%get_cmat_at(phys)) !! Real(C) == Real(C*)
                        val      = min(val_prev/val, 1.e20)
                        call W_img%set_cmat_at( phys, cmplx(val, 0.)) 
                        Wnorm    = Wnorm + val ! values are positive, no need to square (numerical stability)
                    end do
                end do
            end do
            !$omp end parallel do
            Wnorm = log10(1. + Wnorm / sqrt(real(product(lims(:,2) - lims(:,1)))))
            if( Wnorm/Wnorm_prev < 1.05 ) exit
        enddo
        call Wprev_img%kill
        ! Fourier comps / rho
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,invrho) proc_bind(close)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do m = lims(3,1),lims(3,2)
                    phys   = W_img%comp_addr_phys([h, k, m])
                    invrho = real(W_img%get_cmat_at(phys)) !! Real(C) == Real(C*)
                    call self%mul_cmat_at(invrho,phys)
                end do
            end do
        end do
        !$omp end parallel do
        ! cleanup
        call W_img%kill
    end subroutine sampl_dens_correct

    subroutine compress_exp( self )
        class(reconstructor), intent(inout) :: self
        complex :: comp
        integer :: lims(3,2), phys(3), h, k, m
        if(.not. self%cmat_exp_allocated .or. .not.self%rho_allocated)then
            stop 'expanded complex or rho matrices do not exist; simple_reconstructor::compress_exp'
        endif
        call self%reset
        lims = self%loop_lims(3)
        ! Fourier components & rho matrices compression
        !$omp parallel do collapse(3) private(h,k,m,phys,comp) schedule(static) default(shared) proc_bind(close)
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                do m = lims(3,1),lims(3,2)
                    comp = self%cmat_exp(h,k,m)
                    if(abs(comp) < TINY) cycle
                    !  addition because FC and its Friedel mate must be summed
                    !  as the expansion updates one or the other
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + self%ldim_img(2) * MERGE(1,0,k < 0)
                        phys(3) = m + 1 + self%ldim_img(3) * MERGE(1,0,m < 0)
                        call self%add2_cmat_at(phys, comp)
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + self%ldim_img(2) * MERGE(1,0,-k < 0)
                        phys(3) = -m + 1 + self%ldim_img(3) * MERGE(1,0,-m < 0)
                        call self%add2_cmat_at(phys, conjg(comp))
                    endif
                    self%rho(phys(1),phys(2),phys(3)) = &
                        &self%rho(phys(1),phys(2),phys(3)) + self%rho_exp(h,k,m)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compress_exp

    ! SUMMATION

    !> for summing reconstructors generated by parallel execution
    subroutine sum( self, self_in )
         class(reconstructor), intent(inout) :: self !< this instance
         class(reconstructor), intent(in)    :: self_in !< other instance
         call self%add(self_in)
         !$omp parallel workshare proc_bind(close)
         self%rho = self%rho+self_in%rho
         !$omp end parallel workshare
    end subroutine sum

    ! RECONSTRUCTION
    !> reconstruction routine
    subroutine rec( self, fname, p, o, se, state, mul, part )
        class(reconstructor), intent(inout) :: self      !< this object
        character(len=*),     intent(inout) :: fname     !< spider/MRC stack filename
        class(params),        intent(in)    :: p         !< parameters
        class(oris),          intent(inout) :: o         !< orientations
        class(sym),           intent(inout) :: se        !< symmetry element
        integer,              intent(in)    :: state     !< state to reconstruct
        real,    optional,    intent(in)    :: mul       !< shift multiplication factor
        integer, optional,    intent(in)    :: part      !< partition (4 parallel rec)
        type(image) :: img, img_pd
        real        :: skewness
        integer     :: statecnt(p%nstates), i, cnt, n, ldim(3)
        integer     :: state_here, state_glob
        call find_ldim_nptcls(fname, ldim, n)
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; rec; simple_reconstructor'
        ! stash global state index
        state_glob = state
        ! make the images
        call img_pd%new([p%boxpd,p%boxpd,1],self%get_smpd())
        call img%new([p%box,p%box,1],self%get_smpd())
        ! population balancing logics
        if( p%balance > 0 )then
            call o%balance( p%balance, NSPACE_BALANCE, p%nsym, p%eullims, skewness )
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call o%set_all2single('state_balance', 1.0)
        endif
        ! zero the Fourier volume and rho
        call self%reset
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt+1
                state_here = nint(o%get(i,'state'))
                if( state_here > 0 .and. (state_here == state) )then
                    statecnt(state) = statecnt(state)+1
                    call rec_dens
                endif
            endif
        end do
        if( present(part) )then
            return
        else
            write(*,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
            call self%compress_exp
            call self%sampl_dens_correct
        endif
        call self%bwd_ft
        call img%kill
        call img_pd%kill
        ! report how many particles were used to reconstruct each state
        if( p%nstates > 1 )then
            write(*,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif
        contains

            !> \brief  the density reconstruction functionality
            subroutine rec_dens
                character(len=:), allocatable :: stkname
                type(ori) :: orientation, o_sym
                integer   :: j, state_balance, ind
                real      :: pw
                state_balance = nint(o%get(i, 'state_balance'))
                if( state_balance == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    if( p%l_stktab_input )then
                        call p%stkhandle%get_stkname_and_ind(i, stkname, ind)
                        call img%read(stkname, ind)
                    else
                        if( p%l_distr_exec )then
                            call img%read(p%stk_part, cnt)
                        else
                            call img%read(fname, i)
                        endif
                    endif
                    call prep4cgrid(img, img_pd, p%msk, self%kbwin)
                    if( p%pgrp == 'c1' )then
                        call self%inout_fplane(orientation, .true., img_pd, pwght=pw)
                    else
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            call self%inout_fplane(o_sym, .true., img_pd, pwght=pw)
                        end do
                    endif
                endif
            end subroutine rec_dens

    end subroutine rec

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine dealloc_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        if( allocated(self%rho_exp) ) deallocate(self%rho_exp)
        if( allocated(self%cmat_exp) ) deallocate(self%cmat_exp)
        self%rho_exp_allocated  = .false.
        self%cmat_exp_allocated = .false.
    end subroutine dealloc_exp

    !>  \brief  is a destructor
    subroutine dealloc_rho( self )
        class(reconstructor), intent(inout) :: self !< this instance
        call self%dealloc_exp
        if( self%rho_allocated )then
            call fftwf_free(self%kp)
            self%rho => null()
            self%rho_allocated = .false.
        endif
    end subroutine dealloc_rho

end module simple_reconstructor
