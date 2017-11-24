! 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
!$ use omp_lib
!$ use omp_lib_kinds
use, intrinsic :: iso_c_binding
#include "simple_lib.f08"
!! import functions
use simple_fftw3,      only: fftwf_alloc_real, fftwf_free
use simple_timer,      only: tic, toc, timer_int_kind
!! import classes
use simple_ctf,        only: ctf
use simple_ori,        only: ori
use simple_oris,       only: oris
use simple_params,     only: params
use simple_sym,        only: sym
use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
implicit none

public :: reconstructor
private

type, extends(image) :: reconstructor
    private
    type(kbinterpol)            :: kbwin                        !< window function object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    type(ctf)                   :: tfun                         !< CTF object
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable        :: cmat_exp(:,:,:)              !< Fourier components of expanded reconstructor
    real,    allocatable        :: rho_exp(:,:,:)               !< sampling+CTF**2 density of expanded reconstructor
    real                        :: winsz          = RECWINSZ    !< window half-width
    real                        :: alpha          = KBALPHA     !< oversampling ratio
    real                        :: dfx=0., dfy=0., angast=0.    !< CTF params
    real                        :: phshift        = 0.          !< additional phase shift from the Volta 
    real                        :: shconst_rec(3) = [0.,0.,0.]  !< memoized constants for origin shifting 
    integer                     :: wdim           = 0           !< dim of interpolation matrix
    integer                     :: nyq            = 0           !< Nyqvist Fourier index
    integer                     :: ldim_img(3)    = 0           !< logical dimension of the original image
    integer                     :: lims(3,2)      = 0           !< Friedel limits
    integer                     :: cyc_lims(3,2)  = 0           !< redundant limits
    type(CTFFLAGTYPE)           :: ctf                          !< ctf flag <yes|no|mul|flip>
    logical                     :: tfastig        = .false.     !< astigmatic CTF or not
    logical                     :: phaseplate     = .false.     !< Volta phaseplate images or not
    logical                     :: rho_allocated  = .false.     !< existence of rho matrix
  contains
    ! CONSTRUCTORS 
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    procedure          :: reset_exp
    procedure          :: apply_weight
    ! GETTER
    procedure          :: get_kbwin
    ! I/O
    ! procedure          :: write_rho
    ! procedure          :: read_rho
    procedure          :: write_rho_exp
    procedure          :: read_rho_exp
    procedure          :: write_cmat_exp
    procedure          :: read_cmat_exp
    ! INTERPOLATION
    procedure, private :: inout_fcomp
    procedure, private :: inout_fplane_1
    procedure, private :: inout_fplane_2
    generic            :: inout_fplane => inout_fplane_1, inout_fplane_2
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    ! SUMMATION
    procedure          :: sum
    procedure          :: sum_exp
    ! RECONSTRUCTION
    procedure          :: rec
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

real, parameter :: SHTHRESH=0.0001
#include "simple_local_flags.inc"

contains

    ! CONSTRUCTOR

    subroutine alloc_rho( self, p, expand )
        class(reconstructor), intent(inout) :: self   !< this instance
        class(params),        intent(in)    :: p      !< parameters object
        logical, optional,    intent(in)    :: expand !< expand flag
        integer :: rho_shape(3), ldim_exp(3,2), dim
        logical :: l_expand
        if(.not. self%exists()) stop 'construct image before allocating rho; alloc_rho; simple_reconstructor'
        if(      self%is_2d() ) stop 'only for volumes; alloc_rho; simple_reconstructor'
        call self%dealloc_rho
        l_expand = .true.
        if( present(expand) ) l_expand = expand
        self%ldim_img = self%get_ldim()
        self%nyq      = self%get_lfny(1)
        self%winsz    = p%winsz
        self%alpha    = p%alpha
        select case(p%ctf)
            case('no')
                self%ctf%flag = CTFFLAG_NO
            case('yes')
                self%ctf%flag = CTFFLAG_YES
            case('mul')
                stop 'ERROR! ctf=mul deprecated; simple_reconstructor :: alloc_rho'
            case('flip')
                self%ctf%flag = CTFFLAG_FLIP
        end select
        self%tfastig = .false.
        if( p%tfplan%mode .eq. 'astig' ) self%tfastig = .true.
        self%phaseplate  = p%tfplan%l_phaseplate
        self%kbwin       = kbinterpol(self%winsz,self%alpha)
        self%wdim        = self%kbwin%get_wdim()
        self%lims        = self%loop_lims(2)
        self%cyc_lims    = self%loop_lims(3)
        self%shconst_rec = self%get_shconst() 
        ! Work out dimensions of the rho array
        rho_shape(1)    = fdim(self%ldim_img(1))
        rho_shape(2:3)  = self%ldim_img(2:3)
        ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
        self%kp = fftwf_alloc_real(int(product(rho_shape),c_size_t))
        ! Set up the rho array which will point at the allocated memory
        call c_f_pointer(self%kp,self%rho,rho_shape)
        self%rho_allocated = .true.
        if( l_expand )then
            ! setup expanded matrices
            dim  = maxval(abs(self%lims)) + ceiling(KBWINSZ)
            ldim_exp(1,:) = [self%lims(1,1)-self%wdim, dim]
            ldim_exp(2,:) = [-dim, dim]
            ldim_exp(3,:) = [-dim, dim]
            allocate(self%cmat_exp( ldim_exp(1,1):ldim_exp(1,2),ldim_exp(2,1):ldim_exp(2,2),&
                &ldim_exp(3,1):ldim_exp(3,2)), source=cmplx(0.,0.), stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor cmat_exp")
            allocate(self%rho_exp( ldim_exp(1,1):ldim_exp(1,2),ldim_exp(2,1):ldim_exp(2,2),&
                &ldim_exp(3,1):ldim_exp(3,2)), source=0., stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor rho_exp")
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

    ! resets the reconstructor expanded matrices before reconstruction
    subroutine reset_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        if(allocated(self%cmat_exp))self%cmat_exp = cmplx(0.,0.)
        if(allocated(self%rho_exp)) self%rho_exp  = 0.
    end subroutine reset_exp

    subroutine apply_weight( self, w )
        class(reconstructor), intent(inout) :: self !< this instance
        real,                 intent(in)    :: w
        if(allocated(self%cmat_exp)) self%cmat_exp = self%cmat_exp * w
        if(allocated(self%rho_exp))  self%rho_exp  = self%rho_exp  * w
    end subroutine apply_weight

    ! GETTERS
    !> get the kbintpol window
    function get_kbwin( self ) result( wf )
        class(reconstructor), intent(inout) :: self !< this instance
        type(kbinterpol) :: wf                      !< return kbintpol window
        wf = kbinterpol(self%winsz,self%alpha)
    end function get_kbwin

    ! I/O
    !>Write reconstructed image
    ! subroutine write_rho( self, kernam )
    !     class(reconstructor), intent(in) :: self   !< this instance
    !     character(len=*),     intent(in) :: kernam !< kernel name
    !     integer :: filnum
    !     call fopen(filnum, trim(kernam), status='REPLACE', action='WRITE', access='STREAM')
    !     write(filnum, pos=1) self%rho
    !     call fclose(filnum)
    ! end subroutine write_rho

    !> Read sampling density matrix
    ! subroutine read_rho( self, kernam )
    !     class(reconstructor), intent(inout) :: self   !< this instance
    !     character(len=*),     intent(in)    :: kernam !< kernel name
    !     integer :: filnum
    !     call fopen(filnum, file=trim(kernam), status='OLD', action='READ', access='STREAM')
    !     read(filnum, pos=1) self%rho
    !     call fclose(filnum)
    ! end subroutine read_rho

    !>Write reconstructed image
    subroutine write_rho_exp( self, kernam )
        class(reconstructor), intent(in) :: self   !< this instance
        character(len=*),     intent(in) :: kernam !< kernel name
        integer :: filnum
        call fopen(filnum, trim(kernam), status='REPLACE', action='WRITE', access='STREAM')
        write(filnum, pos=1) self%rho_exp
        call fclose(filnum)
    end subroutine write_rho_exp

    !> Read sampling density matrix
    subroutine read_rho_exp( self, kernam )
        class(reconstructor), intent(inout) :: self   !< this instance
        character(len=*),     intent(in)    :: kernam !< kernel name
        integer :: filnum
        call fopen(filnum, file=trim(kernam), status='OLD', action='READ', access='STREAM')
        read(filnum, pos=1) self%rho_exp
        call fclose(filnum)
    end subroutine read_rho_exp

    !>Write reconstructed image
    subroutine write_cmat_exp( self, kernam )
        class(reconstructor), intent(in) :: self   !< this instance
        character(len=*),     intent(in) :: kernam !< kernel name
        integer :: filnum
        call fopen(filnum, trim(kernam), status='REPLACE', action='WRITE', access='STREAM')
        write(filnum, pos=1) self%cmat_exp
        call fclose(filnum)
    end subroutine write_cmat_exp

    !> Read sampling density matrix
    subroutine read_cmat_exp( self, kernam )
        class(reconstructor), intent(inout) :: self   !< this instance
        character(len=*),     intent(in)    :: kernam !< kernel name
        integer :: filnum
        call fopen(filnum, file=trim(kernam), status='OLD', action='READ', access='STREAM')
        read(filnum, pos=1) self%cmat_exp
        call fclose(filnum)
    end subroutine read_cmat_exp

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
        real    :: sqSpatFreq, ang, inv1, inv2
        integer :: i, win(3,2)
        ! calculate non-uniform sampling location
        vec = [real(h), real(k), 0.]
        loc = matmul(vec, e%get_mat())
        ! initiate kernel matrix
        call sqwin_3d(loc(1), loc(2), loc(3), self%winsz, win)
        ! no need to update outside the non-redundant Friedel limits
        ! consistently with compress_exp
        if( win(1,2) < self%lims(1,1) )return
        ! evaluate the transfer function
        if( self%ctf%flag /= CTFFLAG_NO )then
            inv1       = vec(1)*(1./real(self%ldim_img(1)))
            inv2       = vec(2)*(1./real(self%ldim_img(2)))
            sqSpatFreq = inv1 * inv1 + inv2 * inv2
            ang        = atan2(vec(2), vec(1))
            ! calculate CTF and CTF**2 values
            if( self%phaseplate )then
                tval = self%tfun%eval(sqSpatFreq, self%dfx, self%dfy, self%angast, ang, self%phshift)
            else
                tval = self%tfun%eval(sqSpatFreq, self%dfx, self%dfy, self%angast, ang)
            endif
            tvalsq = tval * tval
            if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
        else
            tval   = 1.
            tvalsq = tval
        endif
        ! (weighted) kernel values
        w = pwght
        do i=1,self%wdim
            w(i,:,:) = w(i,:,:) * self%kbwin%apod(real(win(1,1) + i - 1) - loc(1))
            w(:,i,:) = w(:,i,:) * self%kbwin%apod(real(win(2,1) + i - 1) - loc(2))
            w(:,:,i) = w(:,:,i) * self%kbwin%apod(real(win(3,1) + i - 1) - loc(3))
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
            ! subtraction
            ! CTF and w modulates the component before origin shift
            self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - (comp*tval*w)*oshift
            self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - tvalsq*w
        endif
    end subroutine inout_fcomp

    !> insert or uninsert Fourier plane
    subroutine inout_fplane_1( self, o, inoutmode, fpl, pwght )
        class(reconstructor), intent(inout) :: self      !< instance
        class(ori),           intent(inout) :: o         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        class(image),         intent(inout) :: fpl       !< Fourier plane
        real,                 intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        integer :: h, k, lims(3,2), logi(3), phys(3), sh
        complex :: oshift
        real    :: x, y
        if( self%ctf%flag /= CTFFLAG_NO )then ! make CTF object & get CTF info
            self%tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            self%dfx  = o%get('dfx')
            if( self%tfastig )then            ! astigmatic CTF model
                self%dfy = o%get('dfy')
                self%angast = o%get('angast')
            else                              ! non-astigmatic CTF model
                self%dfy = self%dfx
                self%angast = 0.
            endif
            ! additional phase shift from the Volta
            self%phshift = 0.
            if( self%phaseplate ) self%phshift = o%get('phshift')
        endif
        lims = fpl%loop_lims(3)
        x    = o%get('x')
        y    = o%get('y')
        !$omp parallel do collapse(2) default(shared) schedule(static)&
        !$omp private(h,k,sh,oshift,logi,phys) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sh = nint(hyp(real(h),real(k)))
                if( sh > self%nyq + 1 )cycle
                logi   = [h,k,0]
                oshift = fpl%oshift(logi, [-x,-y,0.])
                phys   = fpl%comp_addr_phys(logi)
                call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp(logi,phys),oshift,pwght)
            end do
        end do
        !$omp end parallel do
    end subroutine inout_fplane_1

    !> insert or uninsert Fourier plane
    subroutine inout_fplane_2( self, o, inoutmode, cmat, pwght )
        class(reconstructor), intent(inout) :: self      !< instance
        class(ori),           intent(inout) :: o         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        complex(sp),          intent(in)    :: cmat(self%cyc_lims(1,1):self%cyc_lims(1,2),self%cyc_lims(2,1):self%cyc_lims(2,2))
        real,                 intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        integer :: h, k, sh
        complex :: oshift
        real    :: x, y, arg
        if( self%ctf%flag /= CTFFLAG_NO )then ! make CTF object & get CTF info
            self%tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            self%dfx  = o%get('dfx')
            if( self%tfastig )then            ! astigmatic CTF model
                self%dfy = o%get('dfy')
                self%angast = o%get('angast')
            else                              ! non-astigmatic CTF model
                self%dfy = self%dfx
                self%angast = 0.
            endif
            ! additional phase shift from the Volta
            self%phshift = 0.
            if( self%phaseplate ) self%phshift = o%get('phshift')
        endif
        x = o%get('x')
        y = o%get('y')
        !$omp parallel do collapse(2) default(shared) schedule(static)&
        !$omp private(h,k,sh,arg,oshift) proc_bind(close)
        do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
            do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
                sh = nint(hyp(real(h),real(k)))
                if( sh > self%nyq + 1 )cycle
                arg    = real(h) * (-x) * self%shconst_rec(1) + real(k) * (-y) * self%shconst_rec(2)
                oshift = cmplx(cos(arg),sin(arg))
                call self%inout_fcomp(h,k,o,inoutmode,cmat(h,k),oshift,pwght)
            end do
        end do
        !$omp end parallel do
    end subroutine inout_fplane_2

    !>  is for uneven distribution of orientations correction 
    !>  from Pipe & Menon 1999
    subroutine sampl_dens_correct( self, maxits )
        use simple_gridding, only: mul_w_instr
        class(reconstructor), intent(inout) :: self
        integer,    optional, intent(in)    :: maxits
        type(kbinterpol)     :: kbwin 
        type(image)          :: W_img, Wprev_img
        real                 :: val_prev, val, invrho
        integer              :: h, k, m, phys(3), iter
        integer              :: maxits_here
        complex, parameter   :: one = cmplx(1.,0.)
        real,    parameter   :: winsz  = 2.
        maxits_here = GRIDCORR_MAXITS
        if( present(maxits) )maxits_here = maxits
        call W_img%new(self%ldim_img, self%get_smpd())
        call Wprev_img%new(self%ldim_img, self%get_smpd())
        call W_img%set_ft(.true.)
        call Wprev_img%set_ft(.true.)
        ! kernel
        kbwin = kbinterpol(winsz, self%alpha)
        ! weights init to 1.
        W_img = one
        do iter = 1, maxits_here
            Wprev_img  = W_img 
            ! W <- W * rho
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys) proc_bind(close)
            do h = self%lims(1,1),self%lims(1,2)
                do k = self%lims(2,1),self%lims(2,2)
                    do m = self%lims(3,1),self%lims(3,2)
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
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys,val,val_prev) proc_bind(close)
            do h = self%lims(1,1),self%lims(1,2)
                do k = self%lims(2,1),self%lims(2,2)
                    do m = self%lims(3,1),self%lims(3,2)
                        phys     = W_img%comp_addr_phys([h, k, m])
                        val      = mycabs(W_img%get_cmat_at(phys))   !! ||C|| == ||C*||
                        val_prev = real(Wprev_img%get_cmat_at(phys)) !! Real(C) == Real(C*)
                        val      = min(val_prev/val, 1.e20)
                        call W_img%set_cmat_at( phys, cmplx(val, 0.)) 
                    end do
                end do
            end do
            !$omp end parallel do
        enddo
        call Wprev_img%kill
        ! Fourier comps / rho
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,invrho) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
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
        integer :: phys(3), h, k, m, logi(3)
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            stop 'expanded complex or rho matrices do not exist; simple_reconstructor::compress_exp'
        endif
        call self%reset
        ! Fourier components & rho matrices compression
        !$omp parallel do collapse(3) private(h,k,m,phys,logi) schedule(static) default(shared) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    if(abs(self%cmat_exp(h,k,m)) < TINY) cycle
                    logi = [h,k,m]
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(self%ldim_img(2),0,k < 0)
                        phys(3) = m + 1 + MERGE(self%ldim_img(3),0,m < 0)
                        call self%set_cmat_at(phys, self%cmat_exp(h,k,m))
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(self%ldim_img(2),0,-k < 0)
                        phys(3) = -m + 1 + MERGE(self%ldim_img(3),0,-m < 0)
                        call self%set_cmat_at(phys, conjg(self%cmat_exp(h,k,m)))
                    endif
                    self%rho(phys(1),phys(2),phys(3)) = self%rho_exp(h,k,m)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compress_exp

    ! SUMMATION

    !> for summing reconstructors generated by parallel execution
    subroutine sum( self, self_to_add )
        class(reconstructor), intent(inout) :: self        !< this instance
        class(reconstructor), intent(in)    :: self_to_add !< other instance
        call self%add_workshare(self_to_add, self%rho, self_to_add%rho)
    end subroutine sum

    !> for summing reconstructors generated by parallel execution
    subroutine sum_exp( self, self_to_add )
        class(reconstructor), intent(inout) :: self        !< this instance
        class(reconstructor), intent(in)    :: self_to_add !< other instance
        !$omp parallel workshare proc_bind(close)
        self%cmat_exp = self%cmat_exp + self_to_add%cmat_exp
        self%rho_exp  = self%rho_exp  + self_to_add%rho_exp
        !$omp end parallel workshare
    end subroutine sum_exp

    ! RECONSTRUCTION
    !> reconstruction routine
    subroutine rec( self, p, o, se, state, mul, part )
        use simple_prep4cgrid, only: prep4cgrid
        class(reconstructor), intent(inout) :: self  !< this object
        class(params),        intent(in)    :: p     !< parameters
        class(oris),          intent(inout) :: o     !< orientations
        class(sym),           intent(inout) :: se    !< symmetry element
        integer,              intent(in)    :: state !< state to reconstruct
        real,    optional,    intent(in)    :: mul   !< shift multiplication factor
        integer, optional,    intent(in)    :: part  !< partition (4 parallel rec)
        type(image)      :: img, img_pad
        type(prep4cgrid) :: gridprep
        real             :: skewness
        integer          :: statecnt(p%nstates), i, cnt, state_here, state_glob
        ! stash global state index
        state_glob = state
        ! make the images
        call img%new([p%box,p%box,1], p%smpd)
        call img_pad%new([p%boxpd,p%boxpd,1], p%smpd)
        ! make the gridding prepper
        call gridprep%new(img, self%kbwin, [p%boxpd,p%boxpd,1])
        ! population balancing logics
        if( p%balance > 0 )then
            call o%balance( p%balance, NSPACE_BALANCE, p%nsym, p%eullims, skewness )
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call o%set_all2single('state_balance', 1.0)
        endif
        ! zero the Fourier volume and rho
        call self%reset
        call self%reset_exp
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt + 1
                state_here = nint(o%get(i,'state'))
                if( state_here > 0 .and. (state_here == state) )then
                    statecnt(state) = statecnt(state) + 1
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
        call img_pad%kill
        ! report how many particles were used to reconstruct each state
        if( p%nstates > 1 )then
            write(*,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif
        contains

            !> \brief  the density reconstruction functionality
            subroutine rec_dens
                character(len=:), allocatable :: stkname
                type(ori) :: orientation, o_sym
                integer   :: j, state, state_balance, ind
                real      :: pw
                state         = nint(o%get(i, 'state'))
                state_balance = nint(o%get(i, 'state_balance'))
                if( state == 0 .or. state_balance == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    if( p%l_stktab_input )then
                        call p%stkhandle%get_stkname_and_ind(i, stkname, ind)
                    else
                        if( p%l_distr_exec )then
                            ind = cnt
                            allocate(stkname, source=trim(p%stk_part))
                        else
                            ind = i
                            allocate(stkname, source=trim(p%stk))
                        endif
                    endif
                    call img%read(stkname, ind)
                    call gridprep%prep(img, img_pad)
                    if( p%pgrp == 'c1' )then
                        call self%inout_fplane(orientation, .true., img_pad, pwght=pw)
                    else
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            call self%inout_fplane(o_sym, .true., img_pad, pwght=pw)
                        end do
                    endif
                    deallocate(stkname)
                endif
            end subroutine rec_dens

    end subroutine rec

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine dealloc_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        if( allocated(self%rho_exp) ) deallocate(self%rho_exp)
        if( allocated(self%cmat_exp) ) deallocate(self%cmat_exp)
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
