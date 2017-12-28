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
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable        :: cmat_exp(:,:,:)              !< Fourier components of expanded reconstructor
    real,    allocatable        :: rho_exp(:,:,:)               !< sampling+CTF**2 density of expanded reconstructor
    real,    allocatable        :: ctf_sqSpatFreq(:,:)          !< CTF squared reciprocal pixels
    real,    allocatable        :: ctf_ang(:,:)                 !< CTF effective defocus
    integer, allocatable        :: ind_map(:,:,:)               !< logical to physical index mapping (2D)
    real                        :: winsz          = RECWINSZ    !< window half-width
    real                        :: alpha          = KBALPHA     !< oversampling ratio
    real                        :: shconst_rec(3) = 0.          !< memoized constants for origin shifting 
    integer                     :: wdim           = 0           !< dim of interpolation matrix
    integer                     :: nyq            = 0           !< Nyqvist Fourier index
    integer                     :: ldim_img(3)    = 0           !< logical dimension of the original image
    integer                     :: ldim_exp(3,2)  = 0           !< logical dimension of the expanded complex matrix
    integer                     :: lims(3,2)      = 0           !< Friedel limits
    integer                     :: rho_shape(3)   = 0           !< shape of sampling density matrix
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
    procedure          :: write_rho
    procedure          :: read_rho
    ! CONVOLUTION INTERPOLATION
    procedure, private :: insert_fplane_1
    procedure, private :: insert_fplane_2
    generic            :: insert_fplane => insert_fplane_1, insert_fplane_2
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    procedure          :: expand_exp
    ! SUMMATION
    procedure          :: sum
    ! RECONSTRUCTION
    procedure          :: rec
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

real, parameter :: SHTHRESH=0.0001
#include "simple_local_flags.inc"

contains

    ! CONSTRUCTORS

    subroutine alloc_rho( self, p, expand )
        class(reconstructor), intent(inout) :: self   !< this instance
        class(params),        intent(in)    :: p      !< parameters object
        logical, optional,    intent(in)    :: expand !< expand flag
        real    :: inv1, inv2
        integer :: dim, h, k, sh
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
        self%rho_shape(1)    = fdim(self%ldim_img(1))
        self%rho_shape(2:3)  = self%ldim_img(2:3)
        ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
        self%kp = fftwf_alloc_real(int(product(self%rho_shape),c_size_t))
        ! Set up the rho array which will point at the allocated memory
        call c_f_pointer(self%kp,self%rho,self%rho_shape)
        self%rho_allocated = .true.
        if( l_expand )then
            ! setup expanded matrices
            dim  = maxval(abs(self%lims)) + ceiling(KBWINSZ)
            self%ldim_exp(1,:) = [self%lims(1,1)-self%wdim, dim]
            self%ldim_exp(2,:) = [-dim, dim]
            self%ldim_exp(3,:) = [-dim, dim]
            allocate(self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),self%ldim_exp(2,1):self%ldim_exp(2,2),&
                &self%ldim_exp(3,1):self%ldim_exp(3,2)), source=cmplx(0.,0.), stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor cmat_exp")
            allocate(self%rho_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),self%ldim_exp(2,1):self%ldim_exp(2,2),&
                &self%ldim_exp(3,1):self%ldim_exp(3,2)), source=0., stat=alloc_stat)
            allocchk("In: alloc_rho; simple_reconstructor rho_exp")
        end if
        ! build CTF related matrices
        if( self%ctf%flag .ne. CTFFLAG_NO)then
            allocate(self%ctf_ang(self%cyc_lims(1,1):self%cyc_lims(1,2), self%cyc_lims(2,1):self%cyc_lims(2,2)),&
            &self%ctf_sqSpatFreq(self%cyc_lims(1,1):self%cyc_lims(1,2), self%cyc_lims(2,1):self%cyc_lims(2,2)), source=0.)
            !$omp parallel do collapse(2) default(shared) schedule(static) private(h,k,sh,inv1,inv2) proc_bind(close)
            do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
                do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
                    sh = nint(hyp(real(h),real(k)))
                    if( sh > self%nyq + 1 )cycle
                    ! evaluate the transfer function
                    inv1 = real(h)*(1./real(self%ldim_img(1)))
                    inv2 = real(k)*(1./real(self%ldim_img(2)))
                    self%ctf_sqSpatFreq(h,k) = inv1 * inv1 + inv2 * inv2
                    self%ctf_ang(h,k) = atan2(real(k), real(h))
                enddo
            enddo
            !$omp end parallel do
        endif
        ! generate index map
        allocate( self%ind_map(self%cyc_lims(1,1):self%cyc_lims(1,2),self%cyc_lims(2,1):self%cyc_lims(2,2),3), source=0)
        call self%get_2Dphys_ind_mapping(self%cyc_lims(1:2,:), self%ind_map)
        !
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
                    self%cmat_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2)) = w * self%cmat_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))
                    self%rho_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))  = w * self%rho_exp(h,k,self%ldim_exp(3,1):self%ldim_exp(3,2))
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine apply_weight

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

    ! CONVOLUTION INTERPOLATION

    !> insert Fourier plane, single orientation
    ! subroutine insert_fplane_1( self, se, o, fpl, pwght )
    !     use simple_sym, only: sym
    !     class(reconstructor), intent(inout) :: self  !< instance
    !     class(sym),           intent(inout) :: se    !< symmetry elements
    !     class(ori),           intent(inout) :: o     !< orientation
    !     class(image),         intent(inout) :: fpl   !< Fourier plane
    !     real,                 intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
    !     real, allocatable :: rotmats(:,:,:)
    !     type(ori) :: o_sym
    !     type(ctf) :: tfun
    !     integer   :: logi(3), win(3,2), sh, i, h, k, nsym, isym
    !     complex   :: comp, oshift
    !     real      :: vec(3), loc(3), shconst_here(2), dfx, dfy, angast, phshift
    !     real      :: w(self%wdim,self%wdim,self%wdim), arg, tval, tvalsq
    !     ! setup CTF
    !     if( self%ctf%flag /= CTFFLAG_NO )then
    !         ! make CTF object & get CTF info
    !         tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
    !         dfx  = o%get('dfx')
    !         if( self%tfastig )then ! astigmatic CTF model
    !             dfy    = o%get('dfy')
    !             angast = o%get('angast')
    !         else                   ! non-astigmatic CTF model
    !             dfy    = dfx
    !             angast = 0.
    !         endif
    !         call tfun%init(dfx, dfy, angast)
    !         ! additional phase shift from the Volta
    !         phshift = 0.
    !         if( self%phaseplate ) phshift = o%get('phshift')
    !     endif
    !     ! setup rotation matrices
    !     nsym = se%get_nsym()
    !     allocate(rotmats(nsym,3,3), source=0.0)
    !     rotmats(1,:,:) = o%get_mat()
    !     if( nsym > 1 )then
    !         do isym=2,nsym
    !             o_sym = se%apply(o, isym)
    !             rotmats(isym,:,:) = o_sym%get_mat()
    !         end do
    !     endif
    !     ! memoize for origin shifting
    !     shconst_here = (-o%get_2Dshift()) * self%shconst_rec(1:2)
    !     ! the parallellisation must run over one plane @ the time to avoid race conditions
    !     ! but by starting the parallel section here we reduce thread creation O/H
    !     ! and lower the serial slack, while preserving a low memory footprint. The speeduop 
    !     ! (on 5000 images of betagal) is modest (10%) but significant
    !     !$omp parallel default(shared) private(i,h,k,sh,comp,arg,oshift,logi,tval,tvalsq,w,win,vec,loc) proc_bind(close)
    !     do isym=1,nsym
    !         !$omp do collapse(2) schedule(static)
    !         do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
    !             do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
    !                 sh = nint(hyp(real(h),real(k)))
    !                 if( sh > self%nyq + 1 ) cycle
    !                 logi = [h,k,0]
    !                 vec  = real(logi)
    !                 ! non-uniform sampling location
    !                 loc  = matmul(vec, rotmats(isym,:,:))
    !                 ! window
    !                 call sqwin_3d(loc(1), loc(2), loc(3), self%winsz, win)
    !                 ! no need to update outside the non-redundant Friedel limits
    !                 ! consistent with compress_exp
    !                 if( win(1,2) < self%lims(1,1) )cycle
    !                 ! Fourier component
    !                 comp   = fpl%get_fcomp(logi, self%ind_map(h,k,:))
    !                 ! shift
    !                 arg    = dot_product(shconst_here, vec(1:2))
    !                 oshift = cmplx(cos(arg), sin(arg))
    !                 ! transfer function
    !                 if( self%ctf%flag /= CTFFLAG_NO )then
    !                     ! CTF and CTF**2 values
    !                     if( self%phaseplate )then
    !                         tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k), phshift)
    !                     else
    !                         tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k))
    !                     endif
    !                     tvalsq = tval * tval
    !                     if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
    !                 else
    !                     tval   = 1.
    !                     tvalsq = tval
    !                 endif
    !                 ! (weighted) kernel values
    !                 w = pwght
    !                 do i=1,self%wdim
    !                     w(i,:,:) = w(i,:,:) * self%kbwin%apod(real(win(1,1) + i - 1) - loc(1))
    !                     w(:,i,:) = w(:,i,:) * self%kbwin%apod(real(win(2,1) + i - 1) - loc(2))
    !                     w(:,:,i) = w(:,:,i) * self%kbwin%apod(real(win(3,1) + i - 1) - loc(3))
    !                 enddo
    !                 ! expanded matrices update
    !                 ! CTF and w modulates the component before origin shift
    !                 self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
    !                 &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + (comp*tval*w)*oshift
    !                 self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
    !                 &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + tvalsq*w
    !             end do
    !         end do
    !         !$omp end do nowait
    !     end do
    !     !$omp end parallel 
    ! end subroutine insert_fplane_1

    ! !> insert Fourier plane, distribution of orientations (with weights)
    ! subroutine insert_fplane_2( self, se, os, fpl, pwght, state )
    !     use simple_sym, only: sym
    !     class(reconstructor), intent(inout) :: self  !< instance
    !     class(sym),           intent(inout) :: se    !< symmetry elements
    !     class(oris),          intent(inout) :: os    !< orientations
    !     class(image),         intent(inout) :: fpl   !< Fourier plane
    !     real,                 intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
    !     integer, optional,    intent(in)    :: state !< state to reconstruct
    !     type(ori) :: o_sym, o
    !     type(ctf) :: tfun
    !     integer   :: logi(3), win(3,2), sh, i, h, k, nsym, isym, iori, noris, sstate, states(os%get_noris())
    !     complex   :: comp, oshift
    !     real      :: vec(3), loc(3), shifts(os%get_noris(),2), ows(os%get_noris()), dfx, dfy, angast, phshift
    !     real      :: w(self%wdim,self%wdim,self%wdim), arg, tval, tvalsq, rotmats(os%get_noris(),se%get_nsym(),3,3)
    !     ! take care of optional state flag
    !     sstate = 1
    !     if( present(state) ) sstate = state
    !     ! setup CTF
    !     if( self%ctf%flag /= CTFFLAG_NO )then
    !         ! make CTF object & get CTF info
    !         o    = os%get_ori(1)
    !         tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
    !         dfx  = o%get('dfx')
    !         if( self%tfastig )then ! astigmatic CTF model
    !             dfy    = o%get('dfy')
    !             angast = o%get('angast')
    !         else                   ! non-astigmatic CTF model
    !             dfy    = dfx
    !             angast = 0.
    !         endif
    !         call tfun%init(dfx, dfy, angast)
    !         ! additional phase shift from the Volta
    !         phshift = 0.
    !         if( self%phaseplate ) phshift = o%get('phshift')
    !     endif
    !     ! setup orientation weights/states/rotation matrices/shifts
    !     nsym  = se%get_nsym()
    !     noris = os%get_noris()
    !     do iori=1,noris
    !         o            = os%get_ori(iori)
    !         ows(iori)    = pwght * o%get('ow')
    !         states(iori) = nint(o%get('state'))
    !         if( ows(iori) < TINY ) cycle
    !         rotmats(iori,1,:,:) = o%get_mat()
    !         if( nsym > 1 )then
    !             do isym=2,nsym
    !                 o_sym = se%apply(o, isym)
    !                 rotmats(iori,isym,:,:) = o_sym%get_mat()
    !             end do
    !         endif
    !         shifts(iori,:) = (-o%get_2Dshift()) * self%shconst_rec(1:2)
    !     enddo
    !     ! the parallellisation must run over one plane @ the time to avoid race conditions
    !     ! but by starting the parallel section here we reduce thread creation O/H
    !     ! and lower the serial slack, while preserving a low memory footprint
    !     !$omp parallel default(shared) private(i,h,k,sh,comp,arg,oshift,logi,tval,tvalsq,w,win,vec,loc) proc_bind(close)
    !     do isym=1,nsym
    !         do iori=1,noris
    !             if( ows(iori) < TINY .or. states(iori) /= sstate ) cycle
    !             !$omp do collapse(2) schedule(static)
    !             do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
    !                 do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
    !                     sh = nint(hyp(real(h),real(k)))
    !                     if( sh > self%nyq + 1 ) cycle
    !                     logi = [h,k,0]
    !                     vec  = real(logi)
    !                     ! non-uniform sampling location
    !                     loc  = matmul(vec, rotmats(iori,isym,:,:))
    !                     ! window
    !                     call sqwin_3d(loc(1), loc(2), loc(3), self%winsz, win)
    !                     ! no need to update outside the non-redundant Friedel limits
    !                     ! consistent with compress_exp
    !                     if( win(1,2) < self%lims(1,1) )cycle
    !                     ! Fourier component
    !                     comp   = fpl%get_fcomp(logi, self%ind_map(h,k,:))
    !                     ! shift
    !                     arg    = dot_product(shifts(iori,:), vec(1:2))
    !                     oshift = cmplx(cos(arg), sin(arg))
    !                     ! transfer function
    !                     if( self%ctf%flag /= CTFFLAG_NO )then
    !                         ! CTF and CTF**2 values
    !                         if( self%phaseplate )then
    !                             tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k), phshift)
    !                         else
    !                             tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k))
    !                         endif
    !                         tvalsq = tval * tval
    !                         if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
    !                     else
    !                         tval   = 1.
    !                         tvalsq = tval
    !                     endif
    !                     ! (weighted) kernel values
    !                     w = ows(iori)
    !                     do i=1,self%wdim
    !                         w(i,:,:) = w(i,:,:) * self%kbwin%apod(real(win(1,1) + i - 1) - loc(1))
    !                         w(:,i,:) = w(:,i,:) * self%kbwin%apod(real(win(2,1) + i - 1) - loc(2))
    !                         w(:,:,i) = w(:,:,i) * self%kbwin%apod(real(win(3,1) + i - 1) - loc(3))
    !                     enddo
    !                     ! expanded matrices update
    !                     ! CTF and w modulates the component before origin shift
    !                     self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
    !                     &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + (comp*tval*w)*oshift
    !                     self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
    !                     &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + tvalsq*w
    !                 end do
    !             end do
    !             !$omp end do nowait
    !         end do
    !     end do
    !     !$omp end parallel 
    ! end subroutine insert_fplane_2

    !> insert Fourier plane, single orientation
    subroutine insert_fplane_1( self, se, o, fpl, pwght )
        use simple_sym, only: sym
        class(reconstructor), intent(inout) :: self  !< instance
        class(sym),           intent(inout) :: se    !< symmetry elements
        class(ori),           intent(inout) :: o     !< orientation
        class(image),         intent(inout) :: fpl   !< Fourier plane
        real,                 intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
        real, allocatable :: rotmats(:,:,:)
        type(ori) :: o_sym
        type(ctf) :: tfun
        integer   :: logi(3), phys(3), sh, i, h, k, nsym, isym, iwinsz
        integer   :: win(2,3) ! window boundary array in fortran contiguous format
        complex   :: comp, oshift
        real      :: w(self%wdim,self%wdim,self%wdim), vec(3), loc(3), dists(3), shconst_here(2)
        real      :: arg, dfx, dfy, angast, phshift, tval, tvalsq
        ! window size
        iwinsz = ceiling(self%winsz - 0.5)
        ! setup CTF
        if( self%ctf%flag /= CTFFLAG_NO )then
            ! make CTF object & get CTF info
            tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            dfx  = o%get('dfx')
            if( self%tfastig )then ! astigmatic CTF model
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else                   ! non-astigmatic CTF model
                dfy    = dfx
                angast = 0.
            endif
            call tfun%init(dfx, dfy, angast)
            ! additional phase shift from the Volta
            phshift = 0.
            if( self%phaseplate ) phshift = o%get('phshift')
        endif
        ! setup rotation matrices
        nsym = se%get_nsym()
        allocate(rotmats(nsym,3,3), source=0.0)
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym=2,nsym
                o_sym = se%apply(o, isym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        ! memoize for origin shifting
        shconst_here = (-o%get_2Dshift()) * self%shconst_rec(1:2)
        ! the parallellisation must run over one plane @ the time to avoid race conditions
        ! but by starting the parallel section here we reduce thread creation O/H
        ! and lower the serial slack, while preserving a low memory footprint. The speeduop 
        ! (on 5000 images of betagal) is modest (10%) but significant
        !$omp parallel default(shared) private(i,h,k,sh,comp,arg,oshift,logi,tval,tvalsq,w,win,vec,loc,dists,phys) proc_bind(close)
        do isym=1,nsym
            !$omp do collapse(2) schedule(static)
            do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
                do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
                    sh = nint(hyp(real(h),real(k)))
                    if( sh > self%nyq + 1 ) cycle
                    logi = [h,k,0]
                    vec  = real(logi)
                    ! non-uniform sampling location
                    loc  = matmul(vec, rotmats(isym,:,:))
                    ! window
                    win(1,:) = nint(loc)
                    win(2,:) = win(1,:) + iwinsz
                    win(1,:) = win(1,:) - iwinsz
                    ! no need to update outside the non-redundant Friedel limits
                    ! consistent with compress_exp
                    if( win(2,1) < self%lims(1,1) )cycle
                    ! Fourier component
                    phys = self%ind_map(h,k,:)
                    comp = fpl%get_fcomp(logi, phys)
                    ! shift
                    arg    = dot_product(shconst_here, vec(1:2))
                    oshift = cmplx(cos(arg), sin(arg))
                    ! transfer function
                    if( self%ctf%flag /= CTFFLAG_NO )then
                        ! CTF and CTF**2 values
                        if( self%phaseplate )then
                            tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k), phshift)
                        else
                            tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k))
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
                        dists    = real(win(1,:) + i - 1) - loc
                        w(i,:,:) = w(i,:,:) * self%kbwin%apod(dists(1))
                        w(:,i,:) = w(:,i,:) * self%kbwin%apod(dists(2))
                        w(:,:,i) = w(:,:,i) * self%kbwin%apod(dists(3))
                    enddo
                    ! expanded matrices update
                    ! CTF and w modulates the component before origin shift
                    self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                    &self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + (comp*tval*w)*oshift
                    self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                    &self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + tvalsq*w
                end do
            end do
            !$omp end do nowait
        end do
        !$omp end parallel
    end subroutine insert_fplane_1

    !> insert Fourier plane, distribution of orientations (with weights)
    subroutine insert_fplane_2( self, se, os, fpl, pwght, state )
        use simple_sym, only: sym
        class(reconstructor), intent(inout) :: self  !< instance
        class(sym),           intent(inout) :: se    !< symmetry elements
        class(oris),          intent(inout) :: os    !< orientations
        class(image),         intent(inout) :: fpl   !< Fourier plane
        real,                 intent(in)    :: pwght !< external particle weight (affects both fplane and rho)
        integer, optional,    intent(in)    :: state !< state to reconstruct
        type(ori) :: o_sym, o
        type(ctf) :: tfun
        complex   :: comp, oshift
        integer   :: logi(3), sh, i, h, k, nsym, isym, iori, noris, sstate, states(os%get_noris()), iwinsz
        integer   :: win(2,3) ! window boundary array in fortran contiguous format
        real      :: vec(3), loc(3), shifts(os%get_noris(),2), ows(os%get_noris()), dfx, dfy, angast, phshift
        real      :: w(self%wdim,self%wdim,self%wdim), arg, tval, tvalsq, rotmats(os%get_noris(),se%get_nsym(),3,3)
        ! take care of optional state flag
        sstate = 1
        if( present(state) ) sstate = state
        ! window size
        iwinsz = ceiling(self%winsz - 0.5)
        ! setup CTF
        if( self%ctf%flag /= CTFFLAG_NO )then
            ! make CTF object & get CTF info
            o    = os%get_ori(1)
            tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            dfx  = o%get('dfx')
            if( self%tfastig )then ! astigmatic CTF model
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else                   ! non-astigmatic CTF model
                dfy    = dfx
                angast = 0.
            endif
            call tfun%init(dfx, dfy, angast)
            ! additional phase shift from the Volta
            phshift = 0.
            if( self%phaseplate ) phshift = o%get('phshift')
        endif
        ! setup orientation weights/states/rotation matrices/shifts
        nsym  = se%get_nsym()
        noris = os%get_noris()
        do iori=1,noris
            o            = os%get_ori(iori)
            ows(iori)    = pwght * o%get('ow')
            states(iori) = nint(o%get('state'))
            if( ows(iori) < TINY ) cycle
            rotmats(iori,1,:,:) = o%get_mat()
            if( nsym > 1 )then
                do isym=2,nsym
                    o_sym = se%apply(o, isym)
                    rotmats(iori,isym,:,:) = o_sym%get_mat()
                end do
            endif
            shifts(iori,:) = (-o%get_2Dshift()) * self%shconst_rec(1:2)
        enddo
        ! the parallellisation must run over one plane @ the time to avoid race conditions
        ! but by starting the parallel section here we reduce thread creation O/H
        ! and lower the serial slack, while preserving a low memory footprint
        !$omp parallel default(shared) private(i,h,k,sh,comp,arg,oshift,logi,tval,tvalsq,w,win,vec,loc) proc_bind(close)
        do isym=1,nsym
            do iori=1,noris
                if( ows(iori) < TINY .or. states(iori) /= sstate ) cycle
                !$omp do collapse(2) schedule(static)
                do h=self%cyc_lims(1,1),self%cyc_lims(1,2)
                    do k=self%cyc_lims(2,1),self%cyc_lims(2,2)
                        sh = nint(hyp(real(h),real(k)))
                        if( sh > self%nyq + 1 ) cycle
                        logi = [h,k,0]
                        vec  = real(logi)
                        ! non-uniform sampling location
                        loc  = matmul(vec, rotmats(iori,isym,:,:))
                        ! window
                        win(1,:) = nint(loc)
                        win(2,:) = win(1,:) + iwinsz
                        win(1,:) = win(1,:) - iwinsz
                        ! no need to update outside the non-redundant Friedel limits
                        ! consistent with compress_exp
                        if( win(2,1) < self%lims(1,1) )cycle
                        ! Fourier component
                        comp   = fpl%get_fcomp(logi, self%ind_map(h,k,:))
                        ! shift
                        arg    = dot_product(shifts(iori,:), vec(1:2))
                        oshift = cmplx(cos(arg), sin(arg))
                        ! transfer function
                        if( self%ctf%flag /= CTFFLAG_NO )then
                            ! CTF and CTF**2 values
                            if( self%phaseplate )then
                                tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k), phshift)
                            else
                                tval = tfun%eval(self%ctf_sqSpatFreq(h,k), self%ctf_ang(h,k))
                            endif
                            tvalsq = tval * tval
                            if( self%ctf%flag == CTFFLAG_FLIP ) tval = abs(tval)
                        else
                            tval   = 1.
                            tvalsq = tval
                        endif
                        ! (weighted) kernel values
                        w = ows(iori)
                        do i=1,self%wdim
                            w(i,:,:) = w(i,:,:) * self%kbwin%apod(real(win(1,1) + i - 1) - loc(1))
                            w(:,i,:) = w(:,i,:) * self%kbwin%apod(real(win(1,2) + i - 1) - loc(2))
                            w(:,:,i) = w(:,:,i) * self%kbwin%apod(real(win(1,3) + i - 1) - loc(3))
                        enddo
                        ! expanded matrices update
                        ! CTF and w modulates the component before origin shift
                        self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                        &self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + (comp*tval*w)*oshift
                        self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) =&
                        &self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)) + tvalsq*w
                    end do
                end do
                !$omp end do nowait
            end do
        end do
        !$omp end parallel 
    end subroutine insert_fplane_2

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
        ! redundant parallel initialisation because the First Touch policy of OpenMP 
        ! distributes the pages over the memory of the system to allow better cache
        ! utilisation
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    phys  = W_img%comp_addr_phys([h,k,m])
                    call W_img%set_cmat_at(phys,cmplx(0., 0.))
                    call Wprev_img%set_cmat_at(phys,cmplx(0., 0.))
                end do
            end do
        end do
        !$omp end parallel do
        ! kernel
        kbwin = kbinterpol(winsz, self%alpha)
        ! weights init to 1.
        W_img = one
        do iter = 1, maxits_here
            Wprev_img = W_img 
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
        integer :: phys(3), h, k, m
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            stop 'expanded complex or rho matrices do not exist; simple_reconstructor::compress_exp'
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

    subroutine expand_exp( self )
        class(reconstructor), intent(inout) :: self
        integer :: phys(3), h, k, m, logi(3)
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            stop 'expanded complex or rho matrices do not exist; simple_reconstructor::expand_exp'
        endif
        call self%reset_exp
        ! Fourier components & rho matrices expansion
        !$omp parallel do collapse(3) private(h,k,m,phys,logi) schedule(static) default(shared) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    logi = [h,k,m]
                    phys = self%comp_addr_phys([h,k,m])
                    ! this should be safe even if there isn't a 1-to-1 correspondence
                    ! btw logi and phys since we are accessing shared data.
                    self%cmat_exp(h,k,m) = self%get_fcomp(logi, phys)
                    self%rho_exp(h,k,m)  = self%rho(phys(1),phys(2),phys(3))
                end do
            end do
        end do
    end subroutine expand_exp

    ! SUMMATION

    !> for summing reconstructors generated by parallel execution
    subroutine sum( self, self_in )
         class(reconstructor), intent(inout) :: self    !< this instance
         class(reconstructor), intent(in)    :: self_in !< other instance
         call self%add_workshare(self_in, self%rho, self_in%rho)
    end subroutine sum

    ! RECONSTRUCTION

    !> reconstruction routine
    subroutine rec( self, p, o, se, state, part )
        use simple_prep4cgrid, only: prep4cgrid
        class(reconstructor), intent(inout) :: self      !< this object
        class(params),        intent(in)    :: p         !< parameters
        class(oris),          intent(inout) :: o         !< orientations
        class(sym),           intent(inout) :: se        !< symmetry element
        integer,              intent(in)    :: state     !< state to reconstruct
        integer, optional,    intent(in)    :: part      !< partition (4 parallel rec)
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
                type(ori) :: orientation
                integer   :: state, state_balance, ind
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
                    call self%insert_fplane(se, orientation, img_pad, pwght=pw)
                    deallocate(stkname)
                endif
            end subroutine rec_dens

    end subroutine rec

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
        if(allocated(self%ctf_ang)        ) deallocate(self%ctf_ang)
        if(allocated(self%ctf_sqSpatFreq) ) deallocate(self%ctf_sqSpatFreq)
        if(allocated(self%ind_map)        ) deallocate(self%ind_map)
    end subroutine dealloc_rho

end module simple_reconstructor
