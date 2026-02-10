!@descr: 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
use simple_core_module_api
use simple_fftw3
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none

public :: reconstructor
private
#include "simple_local_flags.inc"

type, extends(image) :: reconstructor
    private
    type(kbinterpol)            :: kbwin                        !< window function object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable        :: cmat_exp(:,:,:)              !< Fourier components of expanded reconstructor
    real,    allocatable        :: rho_exp(:,:,:)               !< sampling+CTF**2 density of expanded reconstructor
    real                        :: shconst_rec(3) = 0.          !< memoized constants for origin shifting
    integer                     :: wdim           = 0           !< dim of interpolation matrix
    integer                     :: nyq            = 0           !< Nyqvist Fourier index
    integer                     :: sh_lim         = 0           !< limit for shell-limited reconstruction
    integer                     :: ldim_img(3)    = 0           !< logical dimension of the original image
    integer                     :: ldim_exp(3,2)  = 0           !< logical dimension of the expanded complex matrix
    integer                     :: lims(3,2)      = 0           !< Friedel limits
    integer                     :: rho_shape(3)   = 0           !< shape of sampling density matrix
    integer                     :: cyc_lims(3,2)  = 0           !< redundant limits of the 2D image
    integer(kind(ENUM_CTFFLAG)) :: ctfflag                      !< ctf flag <yes=1|no=0|flip=2>
    logical                     :: phaseplate     = .false.     !< Volta phaseplate images or not
    logical                     :: rho_allocated  = .false.     !< existence of rho matrix
  contains
    ! CONSTRUCTORS
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    procedure          :: reset_exp
    procedure          :: apply_weight
    procedure          :: set_sh_lim
    procedure          :: pad_with_zeros
    ! GETTERS
    procedure          :: get_kbwin
    ! I/O
    procedure          :: write_rho, write_rho_as_mrc, write_absfc_as_mrc
    procedure          :: read_rho, read_raw_rho
    ! CONVOLUTION INTERPOLATION
    procedure          :: insert_plane_strided
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    procedure          :: expand_exp
    ! SUMMATION
    procedure          :: sum_reduce
    procedure          :: add_invtausq2rho
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

contains

    ! CONSTRUCTORS

    subroutine alloc_rho( self, spproj, expand )
        use simple_sp_project, only: sp_project
        class(reconstructor), intent(inout) :: self           !< this instance
        class(sp_project),    intent(inout) :: spproj         !< project description
        logical,              intent(in)    :: expand         !< expand flag
        integer :: dim
        logical :: l_expand
        if(.not. self%exists() ) THROW_HARD('construct image before allocating rho; alloc_rho')
        if(      self%is_2d()  ) THROW_HARD('only for volumes; alloc_rho')
        call self%dealloc_rho
        l_expand            = expand
        self%ldim_img       = self%get_ldim()
        self%nyq            = self%get_lfny(1)
        self%sh_lim         = self%nyq
        self%ctfflag        = spproj%get_ctfflag_type(params_glob%oritype)
        self%phaseplate     = spproj%has_phaseplate(params_glob%oritype)
        self%kbwin          = kbinterpol(KBWINSZ,KBALPHA)
        self%wdim           = self%kbwin%get_wdim()
        self%lims           = self%loop_lims(2)
        self%cyc_lims       = self%loop_lims(3)
        self%shconst_rec    = self%get_shconst()
        ! Work out dimensions of the rho array
        self%rho_shape(1)   = fdim(self%ldim_img(1))
        self%rho_shape(2:3) = self%ldim_img(2:3)
        ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
        self%kp = fftwf_alloc_real(int(product(self%rho_shape),c_size_t))
        ! Set up the rho array which will point at the allocated memory
        call c_f_pointer(self%kp,self%rho,self%rho_shape)
        self%rho_allocated = .true.
        if( l_expand )then
            ! setup expanded matrices
            dim = maxval(abs(self%lims)) + ceiling(KBWINSZ)
            self%ldim_exp(1,:) = [self%lims(1,1)-self%wdim, dim]
            self%ldim_exp(2,:) = [-dim, dim]
            self%ldim_exp(3,:) = [-dim, dim]
            allocate(self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),self%ldim_exp(2,1):self%ldim_exp(2,2),&
                &self%ldim_exp(3,1):self%ldim_exp(3,2)), source=cmplx(0.,0.))
            allocate(self%rho_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),self%ldim_exp(2,1):self%ldim_exp(2,2),&
                &self%ldim_exp(3,1):self%ldim_exp(3,2)), source=0.)
        end if
        call self%reset
    end subroutine alloc_rho

    ! SETTERS

    ! Resets the reconstructor object before reconstruction.
    ! the workshare pragma is faster than a parallel do
    subroutine reset( self )
        class(reconstructor), intent(inout) :: self !< this instance
        call self%set_ft(.true.)
        call self%reset_mats(self%rho)
    end subroutine reset

    ! resets the reconstructor expanded matrices before reconstruction
    subroutine reset_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        if(allocated(self%cmat_exp) .and. allocated(self%rho_exp) )then
            !$omp parallel workshare default(shared) proc_bind(close)
            self%cmat_exp = CMPLX_ZERO
            self%rho_exp  = 0.0
            !$omp end parallel workshare
        endif
    end subroutine reset_exp

    ! Multiply matrices by a scalar
    ! the workshare pragma is faster than a parallel do
    subroutine apply_weight( self, w )
        class(reconstructor), intent(inout) :: self
        real,                 intent(in)    :: w
        if(allocated(self%cmat_exp) .and. allocated(self%rho_exp) )then
            !$omp parallel workshare default(shared) proc_bind(close)
            self%cmat_exp = w * self%cmat_exp
            self%rho_exp  = w * self%rho_exp
            !$omp end parallel workshare
        endif
    end subroutine apply_weight

    subroutine set_sh_lim(self, sh_lim)
        class(reconstructor), intent(inout) :: self !< this instance
        integer,              intent(in)    :: sh_lim
        if( sh_lim < 1 )then
            self%sh_lim = self%nyq
        else
            self%sh_lim = min(sh_lim,self%nyq)
        endif
    end subroutine set_sh_lim

    subroutine pad_with_zeros( self, vol_prev, rho_prev )
        class(reconstructor),      intent(inout) :: self !< this instance
        class(image),               intent(in)   :: vol_prev
        real(kind=c_float_complex), intent(in)   :: rho_prev(:,:,:)
        call self%pad_mats(self%rho, vol_prev, rho_prev)
    end subroutine pad_with_zeros

    ! GETTERS

    !> get the kbinterpol window
    function get_kbwin( self ) result( wf )
        class(reconstructor), intent(inout) :: self !< this instance
        type(kbinterpol) :: wf                      !< return kbintpol window
        wf = kbinterpol(KBWINSZ, KBALPHA)
    end function get_kbwin

    ! I/O
    !>Write reconstructed image
    subroutine write_rho( self, kernam )
        class(reconstructor), intent(in) :: self   !< this instance
        class(string),        intent(in) :: kernam !< kernel name
        integer :: filnum, ierr
        call del_file(kernam)
        call fopen(filnum, kernam, status='NEW', action='WRITE', access='STREAM', iostat=ierr)
        call fileiochk( 'simple_reconstructor ; write rho '//kernam%to_char(), ierr)
        write(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileiochk('read_rho; simple_reconstructor writing '//kernam%to_char(), ierr)
        call fclose(filnum)
    end subroutine write_rho

    !> Read sampling density matrix
    subroutine read_rho( self, kernam )
        class(reconstructor), intent(inout) :: self !< this instance
        class(string),        intent(in)    :: kernam !< kernel name
        integer :: filnum, ierr
        call fopen(filnum, file=kernam, status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('read_rho; simple_reconstructor opening '//kernam%to_char(), ierr)
        read(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileiochk('simple_reconstructor::read_rho; simple_reconstructor reading '//kernam%to_char(), ierr)
        call fclose(filnum)
    end subroutine read_rho

    !> Serial read of sampling density matrix given an open pipe
    subroutine read_raw_rho( self, fhandle )
        class(reconstructor), intent(inout) :: self
        integer,              intent(in)    :: fhandle
        integer :: ierr
        read(fhandle, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileiochk('simple_reconstructor::raw read_rho; simple_reconstructor reading from pipe', ierr)
    end subroutine read_raw_rho

    ! CONVOLUTION INTERPOLATION

    !> insert Fourier plane into NATIVE (unpadded) expanded Fourier volume,
    !> sampling the PADDED plane by STRIDING 
    subroutine insert_plane_strided( self, se, o, fpl, pwght )
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(reconstructor), intent(inout) :: self
        class(sym),           intent(inout) :: se
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl
        real,                 intent(in)    :: pwght
        type(ori) :: o_sym
        complex   :: comp
        real      :: rotmats(se%get_nsym(),3,3), w(self%wdim,self%wdim,self%wdim), loc(3), ctfval
        integer   :: win(2,3), i, h, k, l, nsym, isym, iwinsz, sh, stride, fpllims_pd(3,2)
        integer   :: fpllims(3,2), fplnyq, hp, kp, pf
        real      :: pf2
        if( pwght < TINY ) return
        ! window size
        iwinsz = ceiling(KBWINSZ - 0.5)
        ! stride along h dimension for interpolation: all threads are at least
        ! wdim pixels away from each other to avoid race conditions
        stride = self%wdim ! this is striding for OpenMP, unrelated to the Fourier striding
        ! setup rotation matrices
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 ) then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        ! Plane limits/nyq are for the PADDED plane (input)
        fpllims_pd = fpl%frlims
        ! Native (unpadded) iteration limits so that hp=h*pf and kp=k*pf are in-bounds
        fpllims_pd = fpl%frlims
        pf         = STRIDE_GRID_PAD_FAC
        pf2        = real(pf*pf)
        fpllims    = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf)
        ! KB interpolation / insertion
        !$omp parallel default(shared) private(h,k,l,sh,comp,ctfval,w,win,loc,hp,kp) proc_bind(close)
        do isym = 1, nsym
            do l = 0, stride-1
                !$omp do schedule(static)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    hp = h * pf
                    do k = fpllims(2,1), fpllims(2,2)
                        kp = k * pf
                        sh = nint(sqrt(real(h*h + k*k)))
                        if (sh > self%nyq) cycle
                        ! non-uniform sampling location on the ORIGINAL (native) lattice
                        loc = matmul(real([h,k,0]), rotmats(isym,:,:))
                        ! window on ORIGINAL (native) lattice
                        win(1,:) = nint(loc)
                        win(2,:) = win(1,:) + iwinsz
                        win(1,:) = win(1,:) - iwinsz
                        ! no need to update outside the non-redundant Friedel limits consistent with compress_exp
                        if( win(2,1) < self%lims(1,1) ) cycle
                        ! Fourier component & CTF from PADDED plane at STRIDED indices
                        comp   = pwght * pf2 * fplane_get_cmplx(fpl, hp, kp)
                        ! CTF values are calculated analytically, no FFTW/padding scaling to account for
                        ctfval = pwght * fplane_get_ctfsq(fpl, hp, kp)
                        ! KB weights evaluated in ORIGINAL coordinates / geometry
                        call self%kbwin%apod_mat_3d(loc, iwinsz, self%wdim, w)
                        ! expanded matrices update (NATIVE volume)
                        self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) = &
                            self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + comp*w
                        self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) = &
                            self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + ctfval*w
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call o_sym%kill
    end subroutine insert_plane_strided

    subroutine sampl_dens_correct( self )
        class(reconstructor), intent(inout) :: self
        integer :: h,k,m,phys(3),sh
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,sh) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh   = nint(sqrt(real(h*h + k*k + m*m)))
                    phys = self%comp_addr_phys(h, k, m )
                    if( sh > self%sh_lim )then
                        ! outside Nyqvist, zero
                        call self%set_cmat_at(phys(1),phys(2),phys(3), CMPLX_ZERO)
                    else
                        call self%div_cmat_at(phys(1),phys(2),phys(3), self%rho(phys(1),phys(2),phys(3)))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine sampl_dens_correct

    subroutine compress_exp( self )
        class(reconstructor), intent(inout) :: self
        integer :: phys(3), h, k, m
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            THROW_HARD('expanded complex or rho matrices do not exist; compress_exp')
        endif
        ! Fourier components & rho matrices compression
        call self%reset        
        !$omp parallel do collapse(3) private(h,k,m,phys) schedule(static) default(shared) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    if(abs(self%cmat_exp(h,k,m)) < TINY) cycle
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(self%ldim_img(2),0,k < 0)
                        phys(3) = m + 1 + MERGE(self%ldim_img(3),0,m < 0)
                        call self%set_cmat_at(phys(1),phys(2),phys(3), self%cmat_exp(h,k,m))
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(self%ldim_img(2),0,-k < 0)
                        phys(3) = -m + 1 + MERGE(self%ldim_img(3),0,-m < 0)
                        call self%set_cmat_at(phys(1),phys(2),phys(3), conjg(self%cmat_exp(h,k,m)))
                    endif
                    self%rho(phys(1),phys(2),phys(3)) = self%rho_exp(h,k,m)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compress_exp

    subroutine expand_exp( self )
        class(reconstructor), intent(inout) :: self
        integer :: phys(3), h, k, m, logi(3)
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            THROW_HARD('expanded complex or rho matrices do not exist; expand_exp')
        endif
        call self%reset_exp
        ! Fourier components & rho matrices expansion
        !$omp parallel do collapse(3) private(h,k,m,phys,logi) schedule(static) default(shared) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    logi = [h,k,m]
                    phys = self%comp_addr_phys(h,k,m)
                    ! this should be safe even if there isn't a 1-to-1 correspondence
                    ! btw logi and phys since we are accessing shared data.
                    self%cmat_exp(h,k,m) = self%get_fcomp(logi, phys)
                    self%rho_exp(h,k,m)  = self%rho(phys(1),phys(2),phys(3))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine expand_exp

    subroutine write_rho_as_mrc( self, fname )
        class(reconstructor), intent(inout) :: self
        class(string),        intent(in)    :: fname
        type(image) :: img
        integer :: c,phys(3),h,k,m
        call img%new([self%rho_shape(2),self%rho_shape(2),self%rho_shape(2)], 1.0)
        c = self%rho_shape(2)/2+1
        !$omp parallel do collapse(3) private(h,k,m,phys) schedule(static) default(shared) proc_bind(close)
        do h = 0,self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(self%ldim_img(2),0,k < 0)
                        phys(3) = m + 1 + MERGE(self%ldim_img(3),0,m < 0)
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(self%ldim_img(2),0,-k < 0)
                        phys(3) = -m + 1 + MERGE(self%ldim_img(3),0,-m < 0)
                    endif
                    call img%set([1+h,k+c,m+c], self%rho(phys(1),phys(2),phys(3)))
                end do
            end do
        end do
        !$omp end parallel do
        call img%write(fname)
        call img%kill
    end subroutine write_rho_as_mrc

    subroutine write_absfc_as_mrc( self, fname )
        class(reconstructor), intent(inout) :: self
        class(string),        intent(in)    :: fname
        type(image) :: img
        integer :: c,phys(3),h,k,m
        call img%new([self%rho_shape(2),self%rho_shape(2),self%rho_shape(2)], 1.0)
        c = self%rho_shape(2)/2+1
        !$omp parallel do collapse(3) private(h,k,m,phys) schedule(static) default(shared) proc_bind(close)
        do h = 0,self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(self%ldim_img(2),0,k < 0)
                        phys(3) = m + 1 + MERGE(self%ldim_img(3),0,m < 0)
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(self%ldim_img(2),0,-k < 0)
                        phys(3) = -m + 1 + MERGE(self%ldim_img(3),0,-m < 0)
                    endif
                    call img%set([1+h,k+c,m+c], abs(self%get_cmat_at(phys(1),phys(2),phys(3))))
                end do
            end do
        end do
        !$omp end parallel do
        call img%write(fname)
        call img%kill
    end subroutine write_absfc_as_mrc

    ! SUMMATION

    !> for summing reconstructors generated by parallel execution
    subroutine sum_reduce( self, self_in )
        class(reconstructor), intent(inout) :: self    !< this instance
        class(reconstructor), intent(in)    :: self_in !< other instance
        call self%sum_reduce_mats(self_in, self%rho, self_in%rho)
    end subroutine sum_reduce

    subroutine add_invtausq2rho( self, fsc )
        class(reconstructor),  intent(inout) :: self !< this instance
        real,                  intent(in)    :: fsc(:)
        real,     allocatable :: sig2(:), tau2(:), ssnr(:)
        integer,  allocatable :: cnt(:)
        real(dp), allocatable :: rsum(:)
        real              :: fudge, cc, scale, pad_factor, invtau2
        integer           :: h, k, m, sh, phys(3), sz, reslim_ind
        ! logical           :: l_combined
        ! l_combined = trim(params_glob%combine_eo).eq.'yes'
        sz = size(fsc)
        allocate(ssnr(0:sz), rsum(0:sz), cnt(0:sz), tau2(0:sz), sig2(0:sz))
        rsum = 0.d0
        cnt  = 0
        ssnr = 0.0
        tau2 = 0.0
        sig2 = 0.0
        scale = real(params_glob%box_crop) / real(params_glob%box_croppd)
        pad_factor = 1.0 / scale**3
        fudge = params_glob%tau
        ! SSNR
        do k = 1,sz
            cc = max(0.001,fsc(k))
            ! if( l_combined )then
                ! update to filtering scheme since e/o were identical during alignment
                ! cc = sqrt(2.*cc / (cc+1.))
            ! endif
            cc      = min(0.999,cc)
            ssnr(k) = fudge * cc / (1.-cc)
        enddo
        ! Noise
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,sh) proc_bind(close) reduction(+:cnt,rsum)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(scale * sqrt(real(h*h + k*k + m*m)))
                    if( sh > sz ) cycle
                    phys     = self%comp_addr_phys(h, k, m)
                    cnt(sh)  = cnt(sh) + 1
                    rsum(sh) = rsum(sh) + real(self%rho(phys(1),phys(2),phys(3)),dp)
                enddo
            enddo
        enddo
        !$omp end parallel do
        rsum = rsum * pad_factor
        where( rsum > 1.d-10 )
            sig2 = real(real(cnt,dp) / rsum)
        else where
            sig2 = 0.0
        end where
        ! Signal
        tau2 = ssnr * sig2
        ! add Tau2 inverse to denominator
        ! because signal assumed infinite at very low resolution there is no addition
        reslim_ind = max(6, calc_fourier_index(params_glob%hp, params_glob%box_crop, params_glob%smpd_crop))
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,sh,invtau2) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(scale*sqrt(real(h*h + k*k + m*m)))
                    if( (sh < reslim_ind) .or. (sh > sz) ) cycle
                    phys = self%comp_addr_phys(h, k, m)
                    if( tau2(sh) > TINY)then
                        invtau2 = 1.0/(pad_factor*fudge*tau2(sh))
                    else
                        invtau2 = min(1.e3, 1.e3 * self%rho(phys(1),phys(2),phys(3)))
                    endif
                    self%rho(phys(1),phys(2),phys(3)) = self%rho(phys(1),phys(2),phys(3)) + invtau2
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine add_invtausq2rho

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine dealloc_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        if( allocated(self%rho_exp)  ) deallocate(self%rho_exp)
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
