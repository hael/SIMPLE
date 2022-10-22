! 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
!$ use omp_lib
include 'simple_lib.f08'
use simple_sym,        only: sym
use simple_kbinterpol, only: kbinterpol
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
    real                        :: winsz          = RECWINSZ    !< window half-width
    real                        :: alpha          = KBALPHA     !< oversampling ratio
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
    logical                     :: linear_interp  = .false.     !< Reconstruction interpolation false=>kb|true=>trilinear
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
    ! GETTER
    procedure          :: get_kbwin
    ! I/O
    procedure          :: write_rho, write_rho_as_mrc
    procedure          :: read_rho
    ! CONVOLUTION INTERPOLATION
    procedure          :: insert_plane
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

integer(timer_int_kind) :: trec

contains

    ! CONSTRUCTORS

    subroutine alloc_rho( self, spproj, expand )
        use simple_sp_project, only: sp_project
        class(reconstructor), intent(inout) :: self   !< this instance
        class(sp_project),    intent(inout) :: spproj !< project description
        logical, optional,    intent(in)    :: expand !< expand flag
        integer :: dim
        logical :: l_expand
        l_expand = .true.
        if(.not. self%exists() ) THROW_HARD('construct image before allocating rho; alloc_rho')
        if(      self%is_2d()  ) THROW_HARD('only for volumes; alloc_rho')
        if( present(expand) )l_expand = expand
        call self%dealloc_rho
        l_expand = .true.
        if( present(expand) ) l_expand = expand
        self%ldim_img       =  self%get_ldim()
        self%nyq            =  self%get_lfny(1)
        self%sh_lim         =  self%nyq
        self%winsz          =  params_glob%winsz
        self%alpha          =  params_glob%alpha
        self%ctfflag        =  spproj%get_ctfflag_type(params_glob%oritype)
        self%phaseplate     =  spproj%has_phaseplate(params_glob%oritype)
        self%kbwin          =  kbinterpol(self%winsz,self%alpha)
        self%wdim           =  self%kbwin%get_wdim()
        self%lims           =  self%loop_lims(2)
        self%cyc_lims       =  self%loop_lims(3)
        self%shconst_rec    =  self%get_shconst()
        self%linear_interp  = trim(params_glob%interpfun) == 'linear'
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
            dim = maxval(abs(self%lims)) + ceiling(self%winsz)
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
    ! The shared memory used in a parallel section should be initialised
    ! with a (redundant) parallel section, because of how pages are organised.
    ! Memory otherwise becomes associated with the single thread used for
    ! allocation, causing load imbalance. This will reduce cache misses.
    subroutine reset( self )
        class(reconstructor), intent(inout) :: self !< this instance
        integer :: i, j, k
        call self%set_ft(.true.)
        !$omp parallel do collapse(3) default(shared) schedule(static) private(i,j,k) proc_bind(close)
        do i=1,self%rho_shape(1)
            do j=1,self%rho_shape(2)
                do k=1,self%rho_shape(3)
                    call self%set_cmat_at([i,j,k], cmplx(0.,0.))
                    self%rho(i,j,k) = 0.
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine reset

    ! resets the reconstructor expanded matrices before reconstruction
    ! The shared memory used in a parallel section should be initialised
    ! with a (redundant) parallel section, because of how pages are organised.
    ! Memory otherwise becomes associated with the single thread used for
    ! allocation, causing load imbalance. This will reduce cache misses.
    subroutine reset_exp( self )
        class(reconstructor), intent(inout) :: self !< this instance
        integer :: h, k, l
        if(allocated(self%cmat_exp) .and. allocated(self%rho_exp) )then
            !$omp parallel do collapse(3) default(shared) schedule(static) private(h,k,l) proc_bind(close)
            do h=self%ldim_exp(1,1),self%ldim_exp(1,2)
                do k=self%ldim_exp(2,1),self%ldim_exp(2,2)
                    do l=self%ldim_exp(3,1),self%ldim_exp(3,2)
                        self%cmat_exp(h,k,l) = cmplx(0.,0.)
                        self%rho_exp(h,k,l)  = 0.
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine reset_exp

    ! the same trick is applied here (see above) since this is after (single-threaded) read
    subroutine apply_weight( self, w )
        class(reconstructor), intent(inout) :: self
        real,                 intent(in)    :: w
        integer :: h, k
        if(allocated(self%cmat_exp) .and. allocated(self%rho_exp) )then
            !$omp parallel do collapse(2) default(shared) schedule(static) private(h,k) proc_bind(close)
            do h=self%ldim_exp(1,1),self%ldim_exp(1,2)
                do k=self%ldim_exp(2,1),self%ldim_exp(2,2)
                    self%cmat_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2)) = &
                        w * self%cmat_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))
                    self%rho_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))  = &
                        w * self%rho_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))
                end do
            end do
            !$omp end parallel do
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

    ! GETTERS

    !> get the kbinterpol window
    function get_kbwin( self ) result( wf )
        class(reconstructor), intent(inout) :: self !< this instance
        type(kbinterpol) :: wf                      !< return kbintpol window
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
        call fileiochk( 'simple_reconstructor ; write rho '//trim(kernam), ierr)
        write(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileiochk('read_rho; simple_reconstructor writing '//trim(kernam), ierr)
        call fclose(filnum)
    end subroutine write_rho

    !> Read sampling density matrix
    subroutine read_rho( self, kernam )
        class(reconstructor), intent(inout) :: self !< this instance
        character(len=*),     intent(in)    :: kernam !< kernel name
        integer :: filnum, ierr
        call fopen(filnum, file=trim(kernam), status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('read_rho; simple_reconstructor opening '//trim(kernam), ierr)
        read(filnum, pos=1, iostat=ierr) self%rho
        if( ierr .ne. 0 ) &
            call fileiochk('simple_reconstructor::read_rho; simple_reconstructor reading '&
            &// trim(kernam), ierr)
        call fclose(filnum)
    end subroutine read_rho

    ! CONVOLUTION INTERPOLATION

    !> insert Fourier plane, single orientation
    subroutine insert_plane( self, se, o, fpl, pwght )
        use simple_fplane, only: fplane
        class(reconstructor), intent(inout) :: self    !< instance
        class(sym),           intent(inout) :: se      !< symmetry elements
        class(ori),           intent(inout) :: o       !< orientation
        class(fplane),        intent(in)    :: fpl     !< Fourier plane
        real,                 intent(in)    :: pwght   !< external particle weight (affects both fplane and rho)
        type(ori) :: o_sym
        complex   :: comp, oshift
        real      :: rotmats(se%get_nsym(),3,3), w(self%wdim,self%wdim,self%wdim)
        real      :: vec(3), loc(3), odists(3), dists(3), shconst_here(2), scale, arg, ctfval
        real      :: w000, w001, w010, w011, w100, w101, w110, w111
        integer   :: i, h, k, nsym, isym, iwinsz, sh, win(2,3), floc(3), cloc(3)
        if( pwght < TINY )return
        ! window size
        iwinsz = ceiling(self%winsz - 0.5)
        ! setup rotation matrices
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym=2,nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        ! scale & memoize for origin shifting
        scale        = real(self%ldim_img(1)) / real(fpl%ldim(1))
        shconst_here = -o%get_2Dshift() * fpl%shconst(1:2)
        if( self%linear_interp )then
            !$omp parallel default(shared) proc_bind(close)&
            !$omp private(h,k,sh,comp,arg,oshift,ctfval,vec,loc,dists,odists,floc,cloc,w000,w001,w010,w011,w100,w101,w110,w111)
            do isym=1,nsym
                !$omp do collapse(2) schedule(static)
                do h=fpl%frlims(1,1),fpl%frlims(1,2)
                    do k=fpl%frlims(2,1),fpl%frlims(2,2)
                        sh = nint(sqrt(real(h*h + k*k)))
                        if( sh > fpl%nyq ) cycle
                        vec  = real([h,k,0])
                        ! non-uniform sampling location
                        loc  = scale * matmul(vec, rotmats(isym,:,:))
                        ! no need to update outside the non-redundant Friedel limits consistent with compress_exp
                        floc = floor(loc)
                        cloc = floc + 1
                        if( cloc(1) < self%lims(1,1) )cycle
                        ! shift
                        arg    = dot_product(shconst_here, vec(1:2))
                        oshift = cmplx(cos(arg), sin(arg))
                        ! Fourier component x particle weight x shift & CTF
                        comp   = (pwght * fpl%cmplx_plane(h,k)) * oshift
                        ctfval =  pwght * fpl%ctfsq_plane(h,k)
                        ! interpolation Fcs
                        dists  = loc - real(floc)
                        odists = 1.0 - dists
                        w000 = product(odists)
                        w001 = odists(1) * odists(2) *  dists(3)
                        w010 = odists(1) *  dists(2) * odists(3)
                        w011 = odists(1) *  dists(2) *  dists(3)
                        w100 =  dists(1) * odists(2) * odists(3)
                        w101 =  dists(1) * odists(2) *  dists(3)
                        w110 =  dists(1) *  dists(2) * odists(3)
                        w111 = product(dists)
                        self%cmat_exp(floc(1), floc(2), floc(3)) = self%cmat_exp(floc(1), floc(2), floc(3)) + w000 * comp
                        self%cmat_exp(floc(1), floc(2), cloc(3)) = self%cmat_exp(floc(1), floc(2), cloc(3)) + w001 * comp
                        self%cmat_exp(floc(1), cloc(2), floc(3)) = self%cmat_exp(floc(1), cloc(2), floc(3)) + w010 * comp
                        self%cmat_exp(floc(1), cloc(2), cloc(3)) = self%cmat_exp(floc(1), cloc(2), cloc(3)) + w011 * comp
                        self%cmat_exp(cloc(1), floc(2), floc(3)) = self%cmat_exp(cloc(1), floc(2), floc(3)) + w100 * comp
                        self%cmat_exp(cloc(1), floc(2), cloc(3)) = self%cmat_exp(cloc(1), floc(2), cloc(3)) + w101 * comp
                        self%cmat_exp(cloc(1), cloc(2), floc(3)) = self%cmat_exp(cloc(1), cloc(2), floc(3)) + w110 * comp
                        self%cmat_exp(cloc(1), cloc(2), cloc(3)) = self%cmat_exp(cloc(1), cloc(2), cloc(3)) + w111 * comp
                        ! interpolation ctf^2
                        self%rho_exp(floc(1), floc(2), floc(3))  = self%rho_exp(floc(1), floc(2), floc(3))  + w000 * ctfval
                        self%rho_exp(floc(1), floc(2), cloc(3))  = self%rho_exp(floc(1), floc(2), cloc(3))  + w001 * ctfval
                        self%rho_exp(floc(1), cloc(2), floc(3))  = self%rho_exp(floc(1), cloc(2), floc(3))  + w010 * ctfval
                        self%rho_exp(floc(1), cloc(2), cloc(3))  = self%rho_exp(floc(1), cloc(2), cloc(3))  + w011 * ctfval
                        self%rho_exp(cloc(1), floc(2), floc(3))  = self%rho_exp(cloc(1), floc(2), floc(3))  + w100 * ctfval
                        self%rho_exp(cloc(1), floc(2), cloc(3))  = self%rho_exp(cloc(1), floc(2), cloc(3))  + w101 * ctfval
                        self%rho_exp(cloc(1), cloc(2), floc(3))  = self%rho_exp(cloc(1), cloc(2), floc(3))  + w110 * ctfval
                        self%rho_exp(cloc(1), cloc(2), cloc(3))  = self%rho_exp(cloc(1), cloc(2), cloc(3))  + w111 * ctfval
                    end do
                end do
                !$omp end do
            end do
            !$omp end parallel
        else
            ! KB interpolation
            !$omp parallel default(shared) private(i,h,k,sh,comp,arg,oshift,ctfval,w,win,vec,loc,dists)&
            !$omp proc_bind(close)
            do isym=1,nsym
                !$omp do collapse(2) schedule(static)
                do h=fpl%frlims(1,1),fpl%frlims(1,2)
                    do k=fpl%frlims(2,1),fpl%frlims(2,2)
                        sh = nint(sqrt(real(h*h + k*k)))
                        if( sh > fpl%nyq ) cycle
                        vec  = real([h,k,0])
                        ! non-uniform sampling location
                        loc  = scale * matmul(vec, rotmats(isym,:,:))
                        ! window
                        win(1,:) = nint(loc)
                        win(2,:) = win(1,:) + iwinsz
                        win(1,:) = win(1,:) - iwinsz
                        ! no need to update outside the non-redundant Friedel limits consistent with compress_exp
                        if( win(2,1) < self%lims(1,1) )cycle
                        ! Fourier component & CTF
                        comp   = fpl%cmplx_plane(h,k)
                        ctfval = fpl%ctfsq_plane(h,k)
                        ! shift
                        arg    = dot_product(shconst_here, vec(1:2))
                        oshift = cmplx(cos(arg), sin(arg))
                        ! (weighted) kernel & CTF values
                        w = 1.
                        do i=1,self%wdim
                            dists    = real(win(1,:) + i - 1) - loc
                            w(i,:,:) = w(i,:,:) * self%kbwin%apod(dists(1))
                            w(:,i,:) = w(:,i,:) * self%kbwin%apod(dists(2))
                            w(:,:,i) = w(:,:,i) * self%kbwin%apod(dists(3))
                        enddo
                        w = w / sum(w)
                        w = w * pwght
                        ! expanded matrices update
                        self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                            &self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + (comp*w)*oshift
                        self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                            &self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + ctfval*w
                    end do
                end do
                !$omp end do
            end do
            !$omp end parallel
        endif
        call o_sym%kill
    end subroutine insert_plane

    !>  is for uneven distribution of orientations correction
    !>  from Pipe & Menon 1999
    subroutine sampl_dens_correct( self, do_gridcorr )
        use simple_gridding, only: mul_w_instr
        class(reconstructor),    intent(inout) :: self
        logical, optional,       intent(in)    :: do_gridcorr
        complex(kind=c_float_complex), pointer :: cmatW(:,:,:)    =>null()
        complex(kind=c_float_complex), pointer :: cmatWprev(:,:,:)=>null()
        complex, parameter :: one   = cmplx(1.,0.)
        complex, parameter :: zero  = cmplx(0.,0.)
        type(kbinterpol)   :: kbwin
        type(image)        :: W_img, Wprev_img
        real, allocatable  :: antialw(:)
        real               :: winsz, val_prev, val, invrho, rsh_sq
        integer            :: h,k,m, phys(3), iter, sh, cmat_shape(3), i,j,l
        logical            :: l_gridcorr, l_lastiter
        logical, parameter :: skip_pipemenon = .false.
        logical, parameter :: do_hann_window = .true.
        ! kernel
        winsz   = max(1., 2.*self%kbwin%get_winsz())
        kbwin   = kbinterpol(winsz, self%alpha)
        antialw = self%hannw()
        l_gridcorr = .true.
        if( present(do_gridcorr) ) l_gridcorr = do_gridcorr
        l_gridcorr = l_gridcorr .and. (GRIDCORR_MAXITS > 0)
        if( l_gridcorr )then
            cmat_shape = self%get_array_shape()
            call W_img%new(self%ldim_img, self%get_smpd())
            call Wprev_img%new(self%ldim_img, self%get_smpd())
            call W_img%set_ft(.true.)
            call Wprev_img%set_ft(.true.)
            call W_img%get_cmat_ptr(cmatW)
            call Wprev_img%get_cmat_ptr(cmatWprev)
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(i,j,l) proc_bind(close)
            do l = 1,cmat_shape(3)
                do j = 1,cmat_shape(2)
                    do i = 1,cmat_shape(1)
                        ! init
                        cmatWprev(i,j,l) = one
                        ! W <- W * rho
                        cmatW(i,j,l) = cmplx(self%rho(i,j,l),0.)
                    end do
                end do
            end do
            !$omp end parallel do
            if( .not. skip_pipemenon )then
                do iter = 1, GRIDCORR_MAXITS
                    l_lastiter = (iter == GRIDCORR_MAXITS)
                    ! W <- (W / rho) x kernel
                    call W_img%ifft()
                    call mul_w_instr(W_img, params_glob%interpfun, kbwin=kbwin)
                    call W_img%fft()
                    !$omp parallel do default(shared) private(i,j,l,val,val_prev) proc_bind(close)&
                    !$omp collapse(3) schedule(static)
                    do l = 1,cmat_shape(3)
                        do j = 1,cmat_shape(2)
                            do i = 1,cmat_shape(1)
                                ! W <- Wprev / ((W / rho) x kernel)
                                val      = mycabs(cmatW(i,j,l))
                                if( val > 1.0e38 )then
                                    cmatW(i,j,l) = zero
                                else
                                    val_prev     = real(cmatWprev(i,j,l))
                                    cmatW(i,j,l) = cmplx(min(val_prev/val, 1.e20),0.)
                                endif
                                if( l_lastiter )then
                                    cycle
                                else
                                    ! W <- W * rho
                                    cmatWprev(i,j,l) = cmatW(i,j,l)
                                    cmatW(i,j,l)     = self%rho(i,j,l)*cmatW(i,j,l)
                                endif
                            end do
                        end do
                    end do
                    !$omp end parallel do
                enddo
            end if
            nullify(cmatW)
            nullify(cmatWprev)
            call Wprev_img%kill
            ! Fourier comps / rho
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys,rsh_sq,sh,invrho) proc_bind(close)
            do h = self%lims(1,1),self%lims(1,2)
                do k = self%lims(2,1),self%lims(2,2)
                    do m = self%lims(3,1),self%lims(3,2)
                        rsh_sq = real(h*h + k*k + m*m)
                        sh     = nint(sqrt(rsh_sq))
                        phys   = W_img%comp_addr_phys(h, k, m)
                        if( sh > self%sh_lim )then
                            ! outside Nyqvist, zero
                            call self%set_cmat_at(phys(1),phys(2),phys(3), zero)
                        else
                            if( skip_pipemenon )then
                                invrho = 1. / (1.e-2+self%rho(phys(1),phys(2),phys(3)))
                            else
                                invrho = real(W_img%get_cmat_at(phys(1), phys(2), phys(3)))
                            end if
                            if( do_hann_window )then
                                invrho = invrho * antialw(max(1,abs(h)))*antialw(max(1,abs(k)))*antialw(max(1,abs(m)))
                            endif
                            call self%mul_cmat_at(phys(1),phys(2),phys(3),invrho)
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            ! cleanup
            call W_img%kill
        else
            ! division by rho
            !$omp parallel do collapse(3) default(shared) schedule(static)&
            !$omp private(h,k,m,phys,sh) proc_bind(close)
            do h = self%lims(1,1),self%lims(1,2)
                do k = self%lims(2,1),self%lims(2,2)
                    do m = self%lims(3,1),self%lims(3,2)
                        sh   = nint(sqrt(real(h*h + k*k + m*m)))
                        phys = self%comp_addr_phys(h, k, m )
                        if( sh > self%sh_lim )then
                            ! outside Nyqvist, zero
                            call self%set_cmat_at(phys(1),phys(2),phys(3), zero)
                        else
                            call self%div_cmat_at(phys, self%rho(phys(1),phys(2),phys(3)))
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine sampl_dens_correct

    subroutine compress_exp( self )
        class(reconstructor), intent(inout) :: self
        integer :: phys(3), h, k, m
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            THROW_HARD('expanded complex or rho matrices do not exist; compress_exp')
        endif
        call self%reset
        ! Fourier components & rho matrices compression
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
        character(len=*), intent(in) :: fname
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

    ! SUMMATION

    !> for summing reconstructors generated by parallel execution
    subroutine sum_reduce( self, self_in )
        class(reconstructor), intent(inout) :: self    !< this instance
        class(reconstructor), intent(in)    :: self_in !< other instance
        complex(kind=c_float_complex), pointer :: ptr_self(:,:,:)    => null()
        complex(kind=c_float_complex), pointer :: ptr_self_in(:,:,:) => null()
        call self%get_cmat_ptr(ptr_self)
        call self_in%get_cmat_ptr(ptr_self_in)
        !$omp parallel workshare proc_bind(close)
        ptr_self = ptr_self + ptr_self_in
        self%rho = self%rho + self_in%rho
        !$omp end parallel workshare
    end subroutine sum_reduce

    subroutine add_invtausq2rho( self, fsc)
        use simple_estimate_ssnr, only: fsc2optlp_sub
        class(reconstructor),  intent(inout) :: self !< this instance
        real,     allocatable, intent(in)    :: fsc(:)
        real,     allocatable :: sig2(:), tau2(:), ssnr(:)
        integer,  allocatable :: cnt(:)
        real(dp), allocatable :: rsum(:)
        real, parameter   :: fudge = 1.0
        real              :: cc, scale, pad_factor, invtau2
        integer           :: h, k, m, sh, phys(3), sz, reslim_ind
        logical           :: l_combined
        l_combined = trim(params_glob%combine_eo).eq.'yes'
        sz = size(fsc)
        allocate(ssnr(0:sz), rsum(0:sz), cnt(0:sz), tau2(0:sz), sig2(0:sz))
        rsum = 0.d0
        cnt  = 0
        ssnr = 0.0
        tau2 = 0.0
        sig2 = 0.0
        scale = real(params_glob%box) / real(params_glob%boxpd)
        pad_factor = 1.0 / scale**3
        ! SSNR
        do k = 1,sz
            cc = max(0.001,fsc(k))
            if( l_combined )then
                ! update to filtering scheme since e/o were identical during alignment
                cc = sqrt(2.*cc / (cc+1.))
            endif
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
        reslim_ind = max(6, calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
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
            self%linear_interp = .false.
        endif
    end subroutine dealloc_rho

end module simple_reconstructor
