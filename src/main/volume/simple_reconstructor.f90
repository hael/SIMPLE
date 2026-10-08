!@descr: 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
use simple_core_module_api
use simple_fftw3
use simple_image,      only: image
use simple_parameters, only: parameters
use simple_gridding,   only: kb_stencil_inv_envelope_1d, deapodize3D_inplace
implicit none

public :: reconstructor, gridding_half_restore, exp_samples, insert_planes_multi
private
#include "simple_local_flags.inc"

!> Dense outputs produced while restoring one gridding half. The base image is
!! deliberately undeapodized because FSC and cFAR are estimated on it;
!! final is the deapodized map handed to downstream consumers.
type :: gridding_half_restore
    type(image) :: base
    type(image) :: final
  contains
    procedure :: new                => new_gridding_half_restore
    procedure :: prepare_final      => prepare_gridding_half_final
    procedure :: finalize_from_base => finalize_gridding_half_restore
    procedure :: kill               => kill_gridding_half_restore
end type gridding_half_restore

!> The KB sample geometry of points on an expanded Fourier lattice (window corners, normalized separable
!! weights, Friedel flags), built once and gathered from any number of volumes on the lattice; a plane
!! set (new_plane) also holds the plane positions its samples are written back to (put_plane)
type :: exp_samples
    private
    integer :: n    = 0
    integer :: wdim = 0
    integer, allocatable :: win(:,:)    !< (3,n) window lower corners
    real,    allocatable :: w(:,:,:)    !< (wdim,3,n) normalized separable KB weights
    logical, allocatable :: l_conj(:)   !< (n) the point was reflected into h >= 0
    integer, allocatable :: hk(:,:)     !< (2,n) plane positions (new_plane)
    complex, allocatable :: vbuf(:)     !< gathered values of project_fplane, kept with the sample set
    logical :: l_plane = .false.
  contains
    procedure :: new       => exp_samples_new
    procedure :: new_plane => exp_samples_new_plane
    procedure :: gather    => exp_samples_gather
    procedure :: put_plane => exp_samples_put_plane
    procedure :: get_n     => exp_samples_get_n
    procedure :: kill      => exp_samples_kill
end type exp_samples

type, extends(image) :: reconstructor
    private
    class(parameters),  pointer :: p_ptr => null()              !< pointer to parameters object
    type(kbinterpol)            :: kbwin                        !< window function object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable, public :: cmat_exp(:,:,:)             !< Fourier components of expanded reconstructor, only made public for the sake of GPU implementation
    real,    allocatable, public :: rho_exp(:,:,:)              !< sampling+CTF**2 density of expanded reconstructor, only made public for the sake of GPU implementation
    real,    allocatable         :: invenv1d(:)                  !< inverse KB-stencil envelope on the native lattice
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
    logical                     :: rho_allocated  = .false.     !< existence of rho matrix
  contains
    ! CONSTRUCTORS
    procedure          :: new_accumulator
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    procedure          :: reset_exp
    procedure          :: apply_weight
    procedure          :: apply_weight_sums
    procedure          :: set_sh_lim
    procedure          :: pad_with_zeros
    ! GETTERS
    procedure          :: get_kbwin
    procedure          :: get_rho_copy
    ! I/O
    procedure          :: write_raw_accum
    procedure          :: write_rho, write_rho_as_mrc, write_absfc_as_mrc
    procedure          :: read_rho, read_raw_rho
    ! CONVOLUTION INTERPOLATION
    procedure          :: insert_plane_oversamp
    procedure          :: sampl_dens_correct
    procedure          :: deapodize => deapodize_self
    procedure          :: deapodize_volume
    procedure          :: restore_base
    procedure          :: restore_final
    procedure          :: floor_rho_shellwise
    procedure          :: compress_exp
    procedure          :: expand_exp
    procedure          :: project_fplane
    ! SUMMATION
    procedure          :: sum_reduce
    procedure          :: add_invtausq2rho
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
    procedure          :: kill => kill_reconstructor
end type reconstructor

contains

    ! CONSTRUCTORS

    subroutine new_gridding_half_restore( self, backend, wthreads )
        class(gridding_half_restore), intent(inout) :: self
        class(reconstructor),         intent(in)    :: backend
        logical, optional,            intent(in)    :: wthreads
        integer :: ldim(3)
        call self%kill
        if( .not. backend%exists() ) THROW_HARD('construct gridding backend before restoration result')
        ldim = backend%get_ldim()
        call self%base%new(ldim, backend%get_smpd(), wthreads)
    end subroutine new_gridding_half_restore

    subroutine prepare_gridding_half_final( self, backend, wthreads )
        class(gridding_half_restore), intent(inout) :: self
        class(reconstructor),         intent(in)    :: backend
        logical, optional,            intent(in)    :: wthreads
        if( .not. backend%exists() ) THROW_HARD('construct gridding backend before final restoration')
        call self%final%new(backend%get_ldim(), backend%get_smpd(), wthreads)
    end subroutine prepare_gridding_half_final

    !> Construct a complete one-half gridding accumulator on the active crop
    !! grid. This is the reconstruction-backend lifecycle entry point; callers
    !! no longer need to coordinate the inherited image allocation with rho.
    subroutine new_accumulator( self, params, spproj, expand, wthreads )
        use simple_sp_project, only: sp_project
        class(reconstructor),      intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        class(sp_project),         intent(inout) :: spproj
        logical,                   intent(in)    :: expand
        logical, optional,         intent(in)    :: wthreads
        integer :: ldim(3)
        call self%kill
        ldim = [params%box_crop, params%box_crop, params%box_crop]
        call self%image%new(ldim, params%smpd_crop, wthreads)
        call self%alloc_rho(params, spproj, expand)
        if( expand ) call self%reset_exp
    end subroutine new_accumulator

    subroutine alloc_rho( self, params, spproj, expand )
        use simple_sp_project, only: sp_project
        class(reconstructor),      intent(inout) :: self           !< this instance
        class(parameters), target, intent(inout) :: params         !< parameters object
        class(sp_project),         intent(inout) :: spproj         !< project description
        logical,                   intent(in)    :: expand         !< expand flag
        integer :: dim
        logical :: l_expand
        if(.not. self%exists() ) THROW_HARD('construct image before allocating rho; alloc_rho')
        if(      self%is_2d()  ) THROW_HARD('only for volumes; alloc_rho')
        call self%dealloc_rho
        self%p_ptr          => params
        l_expand            = expand
        self%ldim_img       = self%get_ldim()
        self%nyq            = self%get_lfny(1)
        self%sh_lim         = self%nyq
        self%ctfflag        = spproj%get_ctfflag_type(self%p_ptr %oritype)
        self%kbwin          = kbinterpol(KBWINSZ,KBALPHA)
        self%wdim           = self%kbwin%get_wdim()
        self%lims           = self%loop_lims(2)
        self%cyc_lims       = self%loop_lims(3)
        self%shconst_rec    = self%get_shconst()
        call kb_stencil_inv_envelope_1d(self%ldim_img(1), self%invenv1d)
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

    ! Multiply the non-expanded accumulator (Fourier sums & sampling density) by
    ! a scalar; used to decay the persistent trailing-reconstruction chain
    subroutine apply_weight_sums( self, w )
        class(reconstructor), intent(inout) :: self
        real,                 intent(in)    :: w
        call self%scale_mats(self%rho, w)
    end subroutine apply_weight_sums

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

    !> copy of the sampling density; test/diagnostic boundary only
    subroutine get_rho_copy( self, rho )
        class(reconstructor), intent(in)  :: self
        real, allocatable,    intent(out) :: rho(:,:,:)
        if( .not. self%rho_allocated ) THROW_HARD('gridding accumulator density is not allocated')
        allocate(rho(self%rho_shape(1),self%rho_shape(2),self%rho_shape(3)), source=self%rho)
    end subroutine get_rho_copy

    ! I/O

    !> Write one raw gridding accumulator without imposing state/half naming
    !! policy. The caller owns the explicit numerator and density filenames.
    subroutine write_raw_accum( self, numerator_fname, density_fname )
        class(reconstructor), intent(inout) :: self
        class(string),        intent(in)    :: numerator_fname, density_fname
        if( .not. self%exists() ) THROW_HARD('gridding accumulator image is not constructed')
        if( .not. self%rho_allocated ) THROW_HARD('gridding accumulator density is not allocated')
        call self%write(numerator_fname, del_if_exists=.true.)
        call self%write_rho(density_fname)
    end subroutine write_raw_accum

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

    !> Divide a native-grid real-space map by this gridding operator's
    !! discrete KB deposition envelope.
    subroutine deapodize_self( self )
        class(reconstructor), intent(inout) :: self
        if( .not. allocated(self%invenv1d) ) THROW_HARD('gridding deapodization envelope is not built')
        call deapodize3D_inplace(self, self%invenv1d)
    end subroutine deapodize_self

    subroutine deapodize_volume( self, vol )
        class(reconstructor), intent(in)    :: self
        class(image),         intent(inout) :: vol
        if( .not. allocated(self%invenv1d) ) THROW_HARD('gridding deapodization envelope is not built')
        call deapodize3D_inplace(vol, self%invenv1d)
    end subroutine deapodize_volume

    !> Promote the undeapodized FSC base representation to the deapodized final-map
    !! representation without changing the base image.
    subroutine finalize_gridding_half_restore( self, backend )
        class(gridding_half_restore), intent(inout) :: self
        class(reconstructor),         intent(in)    :: backend
        if( .not. self%base%exists() ) THROW_HARD('gridding half restore has no base image')
        call self%final%copy(self%base)
        call backend%deapodize_volume(self%final)
    end subroutine finalize_gridding_half_restore

    !> Restore the dense, undeapodized native-grid map on which the gridding
    !! FSC is estimated. With preserve_numerator=.true., sampling-density
    !! correction runs on a temporary image and the raw Fourier numerator is
    !! restored for a later FSC-prior replay. Rho is deliberately not copied:
    !! attaching a prior to it remains an in-place backend operation.
    subroutine restore_base( self, vol, preserve_numerator )
        class(reconstructor), intent(inout) :: self
        class(image),         intent(inout) :: vol
        logical, optional,    intent(in)    :: preserve_numerator
        complex, allocatable :: cmat(:,:,:)
        logical :: l_preserve
        call validate_restore_source(self)
        l_preserve = .false.
        if( present(preserve_numerator) ) l_preserve = preserve_numerator
        if( l_preserve )then
            cmat = self%get_cmat()
            call self%sampl_dens_correct
            vol = self
            call self%set_cmat(cmat)
            deallocate(cmat)
            call vol%ifft()
            call vol%clip_inplace(self%ldim_img)
        else
            call validate_restore_target(self, vol)
            call self%sampl_dens_correct
            call self%ifft()
            call self%clip(vol)
        endif
    end subroutine restore_base

    !> Restore the dense, deapodized native-grid map shipped to downstream
    !! consumers. Numerator preservation supports the ML replay path after an
    !! FSC-derived prior has been attached to rho.
    subroutine restore_final( self, vol, preserve_numerator )
        class(reconstructor), intent(inout) :: self
        class(image),         intent(inout) :: vol
        logical, optional,    intent(in)    :: preserve_numerator
        complex, allocatable :: cmat(:,:,:)
        logical :: l_preserve
        call validate_restore_source(self)
        l_preserve = .false.
        if( present(preserve_numerator) ) l_preserve = preserve_numerator
        call validate_restore_target(self, vol)
        if( l_preserve ) cmat = self%get_cmat()
        call self%sampl_dens_correct
        call self%ifft()
        if( l_preserve )then
            call vol%zero_and_unflag_ft
            call self%clip(vol)
            call self%deapodize_volume(vol)
            call self%set_cmat(cmat)
            deallocate(cmat)
        else
            call self%deapodize
            call self%clip(vol)
        endif
    end subroutine restore_final

    subroutine validate_restore_source( self )
        class(reconstructor), intent(in) :: self
        if( .not. self%exists() ) THROW_HARD('gridding restoration source is not constructed')
        if( .not. self%rho_allocated ) THROW_HARD('gridding restoration density is not allocated')
    end subroutine validate_restore_source

    subroutine validate_restore_target( self, vol )
        class(reconstructor), intent(in) :: self
        class(image),         intent(in) :: vol
        if( .not. vol%exists() ) THROW_HARD('construct gridding restoration target before use')
        if( any(vol%get_ldim() /= self%ldim_img) ) THROW_HARD('gridding restoration target has incompatible dimensions')
    end subroutine validate_restore_target

    !> KB gridding of one prepared plane into the expanded accumulators (h-strided OpenMP race
    !! avoidance, separable-KB stencil). The optional weight w scales the data and the density terms
    !! alike (a fractional state weight); without it the plane enters with weight one.
    subroutine insert_plane_oversamp( self, se, o, fpl, w )
        use simple_math, only: ceil_div, floor_div
        class(reconstructor), intent(inout) :: self
        class(sym),           intent(inout) :: se
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl
        real, optional,       intent(in)    :: w
        type(ori) :: o_sym
        complex   :: comp, cmplx_raw
        real      :: rotmats(se%get_nsym(),3,3), loc(3), hrow(3), ctfval
        real      :: wx(self%wdim), wy(self%wdim), wz(self%wdim), ww, ctfsq_raw
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2, 3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3, 2)
        integer   :: fpllims(3, 2), pf, ix, iy, iz, hx, ky, mz
        integer   :: nyq_disk, h_sq, k_max_h, k_lo, k_hi
        real      :: source_scale, eps_norm, inv_wdim, dens_w
        logical   :: l_weighted
        ! window size
        iwinsz = ceiling(KBWINSZ - 0.5)
        ! After rotation, source h-lines separated by only wdim can still
        ! update the same 3-D interpolation cell.  sqrt(3)*wdim guarantees
        ! separation by a full window width along at least one target axis.
        stride = ceiling(sqrt(3.0) * real(self%wdim))
        ! setup rotation matrices
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 ) then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        ! Native iteration limits: the planes are stored on the native lattice
        fpllims_pd      = fpl%frlims
        pf           = OSMPL_PAD_FAC
        source_scale = real(pf*pf)
        l_weighted   = present(w)
        dens_w       = 1.0
        if( l_weighted )then
            source_scale = source_scale * w
            dens_w       = w
        endif
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf)
        eps_norm     = epsilon(1.0)
        inv_wdim     = 1.0 / real(self%wdim)
        ! integer disk gate: bit-equivalent to nint(sqrt(h*h+k*k)) > nyq
        nyq_disk = self%nyq * (self%nyq + 1)
        ! KB interpolation / insertion
        !$omp parallel default(shared) private(h,k,l,h_sq,k_max_h,k_lo,k_hi,comp,cmplx_raw,&
        !$omp& ctfsq_raw,ctfval,wx,wy,wz,ww,win,loc,hrow,r11,r12,r13,r21,r22,r23,&
        !$omp& isym,ix,iy,iz,hx,ky,mz) proc_bind(close)
        do isym = 1, nsym
            r11 = rotmats(isym,1,1); r12 = rotmats(isym,1,2); r13 = rotmats(isym,1,3)
            r21 = rotmats(isym,2,1); r22 = rotmats(isym,2,2); r23 = rotmats(isym,2,3)
            do l = 0, stride-1
                !$omp do schedule(static,1)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    h_sq = h*h
                    if( h_sq > nyq_disk ) cycle
                    k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                    k_lo    = max(fpllims(2,1), -k_max_h)
                    k_hi    = min(fpllims(2,2),  k_max_h)
                    hrow(1) = real(h) * r11
                    hrow(2) = real(h) * r12
                    hrow(3) = real(h) * r13
                    do k = k_lo, k_hi
                        ! gen_fplane4rec stores only k<=0, on the native lattice; Friedel symmetry for k>0.
                        if( k <= 0 )then
                            cmplx_raw = fpl%cmplx_plane(h,k)
                            ctfsq_raw = fpl%ctfsq_plane(h,k)
                        else
                            cmplx_raw = conjg(fpl%cmplx_plane(-h,-k))
                            ctfsq_raw = fpl%ctfsq_plane(-h,-k)
                        endif
                        if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. &
                            ctfsq_raw <= TINY ) cycle
                        ! The expanded reconstruction volume is indexed by native
                        ! Fourier-grid cells. Keep sample location and KB weights in
                        ! that same coordinate system.
                        loc(1) = hrow(1) + real(k) * r21
                        loc(2) = hrow(2) + real(k) * r22
                        loc(3) = hrow(3) + real(k) * r23
                        win(1,:) = nint(loc)
                        win(2,:) = win(1,:) + iwinsz
                        win(1,:) = win(1,:) - iwinsz
                        ! no need to update outside the non-redundant Friedel limits consistent with compress_exp
                        if( win(2,1) < self%lims(1,1) ) cycle
                        comp   = source_scale * cmplx_raw
                        ! CTF values are calculated analytically, no FFTW/padding scaling to account for
                        ctfval = ctfsq_raw
                        if( l_weighted ) ctfval = dens_w * ctfval
                        call kb_apod_vecs_3d_fast(loc, wx, wy, wz)
                        do iz = 1, self%wdim
                            mz = win(1,3) + iz - 1
                            do iy = 1, self%wdim
                                ky = win(1,2) + iy - 1
                                do ix = 1, self%wdim
                                    hx = win(1,1) + ix - 1
                                    ww = wx(ix) * (wy(iy) * wz(iz))
                                    self%cmat_exp(hx,ky,mz) = self%cmat_exp(hx,ky,mz) + comp * ww
                                    self%rho_exp( hx,ky,mz) = self%rho_exp( hx,ky,mz) + ctfval * ww
                                end do
                            end do
                        end do
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call o_sym%kill

    contains

        subroutine kb_apod_vecs_3d_fast( loc, wx, wy, wz )
            real, intent(in)  :: loc(3)
            real, intent(out) :: wx(:), wy(:), wz(:)
            integer :: i, win_lo(3)
            real    :: base(3), ww(3), sx, sy, sz
            win_lo = nint(loc) - iwinsz
            base   = real(win_lo) - loc
            do i = 1, self%wdim
                ww    = self%kbwin%apod_fast(base + real(i-1))
                wx(i) = ww(1)
                wy(i) = ww(2)
                wz(i) = ww(3)
            end do
            sx = sum(wx)
            sy = sum(wy)
            sz = sum(wz)
            if( abs(sx) > eps_norm )then
                wx = wx * (1.0 / sx)
            else
                wx = inv_wdim
            endif
            if( abs(sy) > eps_norm )then
                wy = wy * (1.0 / sy)
            else
                wy = inv_wdim
            endif
            if( abs(sz) > eps_norm )then
                wz = wz * (1.0 / sz)
            else
                wz = inv_wdim
            endif
        end subroutine kb_apod_vecs_3d_fast

    end subroutine insert_plane_oversamp

    !>  Floor rho at its shell mean / frac (default 1000, RELION-style; clamped >= 1) before
    !!  sampl_dens_correct, which divides unfloored. Needed for kernel-weighted (flex) state maps,
    !!  whose rho is small and noisy in low-occupancy regions.
    subroutine floor_rho_shellwise( self, frac )
        class(reconstructor), intent(inout) :: self
        real, optional,       intent(in)    :: frac
        real(dp), allocatable :: shsum(:)
        integer,  allocatable :: shcnt(:)
        integer  :: h, k, m, phys(3), sh, nsh
        real     :: ffrac, floor_here
        ffrac = 1000.0
        if( present(frac) ) ffrac = max(1.0, frac)
        nsh = self%sh_lim
        allocate(shsum(0:nsh), source=0.d0)
        allocate(shcnt(0:nsh), source=0)
        ! pass 1: spherical average of rho per shell
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > nsh ) cycle
                    phys = self%comp_addr_phys(h, k, m)
                    shsum(sh) = shsum(sh) + real(self%rho(phys(1),phys(2),phys(3)), dp)
                    shcnt(sh) = shcnt(sh) + 1
                end do
            end do
        end do
        ! pass 2: floor each voxel at its shell mean / frac
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,phys,sh,floor_here) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > nsh ) cycle
                    if( shcnt(sh) < 1 ) cycle
                    phys = self%comp_addr_phys(h, k, m)
                    floor_here = real(shsum(sh)/real(shcnt(sh),dp)) / ffrac
                    if( self%rho(phys(1),phys(2),phys(3)) < floor_here ) &
                        &self%rho(phys(1),phys(2),phys(3)) = floor_here
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(shsum, shcnt)
    end subroutine floor_rho_shellwise

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
        real(dp), allocatable :: rho_shell_sum(:), rho_shell_mean(:)
        integer,  allocatable :: rho_shell_cnt(:)
        complex :: comp_here
        real    :: rho_here, rho_floor
        integer :: phys(3), h, k, m, sh
        if(.not. allocated(self%cmat_exp) .or. .not.allocated(self%rho_exp))then
            THROW_HARD('expanded complex or rho matrices do not exist; compress_exp')
        endif
        ! Fourier components & rho matrices compression
        call self%reset        
        !$omp parallel do collapse(3) private(h,k,m,phys,sh,rho_here,rho_floor,comp_here)&
        !$omp& schedule(static) default(shared) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    comp_here = self%cmat_exp(h,k,m)
                    rho_here  = self%rho_exp(h,k,m)
                    if( abs(comp_here) < TINY .and. rho_here <= TINY ) cycle
                    if (h > 0) then
                        phys(1) = h + 1
                        phys(2) = k + 1 + MERGE(self%ldim_img(2),0,k < 0)
                        phys(3) = m + 1 + MERGE(self%ldim_img(3),0,m < 0)
                        if( abs(comp_here) >= TINY ) call self%set_cmat_at(phys(1),phys(2),phys(3), comp_here)
                    else
                        phys(1) = -h + 1
                        phys(2) = -k + 1 + MERGE(self%ldim_img(2),0,-k < 0)
                        phys(3) = -m + 1 + MERGE(self%ldim_img(3),0,-m < 0)
                        if( abs(comp_here) >= TINY ) call self%set_cmat_at(phys(1),phys(2),phys(3), conjg(comp_here))
                    endif
                    self%rho(phys(1),phys(2),phys(3)) = rho_here
                end do
            end do
        end do
        !$omp end parallel do
        if( allocated(rho_shell_sum) ) deallocate(rho_shell_sum, rho_shell_mean, rho_shell_cnt)
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
        do m = self%lims(3,1),self%lims(3,2)
            do k = self%lims(2,1),self%lims(2,2)
                do h = self%lims(1,1),self%lims(1,2)
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

    ! ---------------------------------------------------------------- multi-volume gather

    !> The KB sample geometry of n native Fourier locations on an expanded lattice: per point the window's
    !! lower corner, the normalized separable weights and whether the point was Friedel-reflected into
    !! h >= 0. Built once, it reads any number of expanded volumes on the lattice (gather). The arrays keep
    !! their capacity across calls, so a sample set reused particle by particle does not reallocate.
    subroutine exp_samples_new( self, rec, locs )
        class(exp_samples),   intent(inout) :: self
        class(reconstructor), intent(in)    :: rec
        real,                 intent(in)    :: locs(:,:)   !< (3,n) native Fourier locations
        integer :: j, n
        n = size(locs,2)
        call alloc_exp_samples(self, n, rec%wdim)
        do j = 1, n
            call exp_samples_set(self, rec, j, locs(:,j))
        end do
    end subroutine exp_samples_new

    !> sample j at the native location loc: Friedel reflection into h >= 0, window lower corner and
    !! normalized separable KB weights
    subroutine exp_samples_set( self, rec, j, loc_in )
        class(exp_samples),   intent(inout) :: self
        class(reconstructor), intent(in)    :: rec
        integer,              intent(in)    :: j
        real,                 intent(in)    :: loc_in(3)
        real    :: loc(3), base(3), ww3(3), s
        integer :: i, ax
        loc = loc_in
        self%l_conj(j) = loc(1) < 0.
        if( self%l_conj(j) ) loc = -loc
        self%win(:,j) = nint(loc) - ceiling(KBWINSZ - 0.5)
        base = real(self%win(:,j)) - loc
        do i = 1, self%wdim
            ww3 = rec%kbwin%apod(base + real(i-1))
            self%w(i,:,j) = ww3
        end do
        do ax = 1, 3
            s = sum(self%w(:,ax,j))
            if( abs(s) > epsilon(1.0) )then
                self%w(:,ax,j) = self%w(:,ax,j) * (1.0 / s)
            else
                self%w(:,ax,j) = 1.0 / real(self%wdim)
            endif
        end do
    end subroutine exp_samples_set

    !> The central-section samples of the plane fpl_ref at orientation o, on the native lattice, k <= 0,
    !! inside the disc of the volume's and the plane's Nyquist (the planes of gen_fplane4rec)
    subroutine exp_samples_new_plane( self, rec, o, fpl_ref )
        use simple_math, only: ceil_div, floor_div
        class(exp_samples),   intent(inout) :: self
        class(reconstructor), intent(in)    :: rec
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl_ref
        real    :: rotmat(3,3), hrow(3), loc(3)
        integer :: fpllims(2,2), h, k, h_sq, k_max_h, k_lo, k_hi, nyq_eff, nyq_disk, ns, ipass
        rotmat       = o%get_mat()
        fpllims(1,1) = ceil_div (fpl_ref%frlims(1,1), OSMPL_PAD_FAC)
        fpllims(1,2) = floor_div(fpl_ref%frlims(1,2), OSMPL_PAD_FAC)
        fpllims(2,1) = ceil_div (fpl_ref%frlims(2,1), OSMPL_PAD_FAC)
        fpllims(2,2) = floor_div(fpl_ref%frlims(2,2), OSMPL_PAD_FAC)
        nyq_eff = rec%nyq
        if( fpl_ref%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl_ref%nyq / OSMPL_PAD_FAC))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        ! pass 1 counts the disc samples, pass 2 fills the (grow-only) sample arrays in place: a sample
        ! set kept across particles allocates nothing per particle
        do ipass = 1, 2
            ns = 0
            do h = fpllims(1,1), fpllims(1,2)
                h_sq = h*h
                if( h_sq > nyq_disk ) cycle
                k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                k_lo    = max(fpllims(2,1), -k_max_h)
                k_hi    = min(0, min(fpllims(2,2), k_max_h))
                if( ipass == 1 )then
                    ns = ns + max(0, k_hi - k_lo + 1)
                    cycle
                endif
                hrow(1) = real(h) * rotmat(1,1)
                hrow(2) = real(h) * rotmat(1,2)
                hrow(3) = real(h) * rotmat(1,3)
                do k = k_lo, k_hi
                    ns = ns + 1
                    loc(1) = hrow(1) + real(k) * rotmat(2,1)
                    loc(2) = hrow(2) + real(k) * rotmat(2,2)
                    loc(3) = hrow(3) + real(k) * rotmat(2,3)
                    call exp_samples_set(self, rec, ns, loc)
                    self%hk(:,ns) = [h, k]
                end do
            end do
            if( ipass == 1 )then
                call alloc_exp_samples(self, ns, rec%wdim)
                if( allocated(self%hk) )then
                    if( size(self%hk,2) < ns ) deallocate(self%hk)
                endif
                if( .not. allocated(self%hk) ) allocate(self%hk(2,max(1,ns)))
            endif
        end do
        self%l_plane = .true.
    end subroutine exp_samples_new_plane

    !> vals(j) is the volume at point j, conjugated where the point was reflected; a point whose window
    !! leaves this volume's expanded lattice reads zero
    subroutine exp_samples_gather( self, rec, vals )
        class(exp_samples),   intent(in)  :: self
        class(reconstructor), intent(in)  :: rec
        complex,              intent(out) :: vals(:)
        call gather_samples(self%n, self%wdim, self%win, self%w, self%l_conj, rec, vals)
    end subroutine exp_samples_gather

    !> exp_samples_gather into the set's own grow-only buffer (project_fplane)
    subroutine exp_samples_gather_vbuf( self, rec )
        class(exp_samples),   intent(inout) :: self
        class(reconstructor), intent(in)    :: rec
        if( allocated(self%vbuf) )then
            if( size(self%vbuf) < self%n ) deallocate(self%vbuf)
        endif
        if( .not. allocated(self%vbuf) ) allocate(self%vbuf(max(1,self%n)))
        call gather_samples(self%n, self%wdim, self%win, self%w, self%l_conj, rec, self%vbuf)
    end subroutine exp_samples_gather_vbuf

    !> the KB-window gather of n samples from rec%cmat_exp; a sample whose window leaves the lattice reads zero
    subroutine gather_samples( n, wdim, win, w, l_conj, rec, vals )
        integer,              intent(in)  :: n, wdim
        integer,              intent(in)  :: win(3,n)
        real,                 intent(in)  :: w(wdim,3,n)
        logical,              intent(in)  :: l_conj(n)
        class(reconstructor), intent(in)  :: rec
        complex,              intent(out) :: vals(:)
        complex :: val
        real    :: wyz
        integer :: lb(3), ub(3), j, ix, iy, iz, hx, ky, mz
        lb = lbound(rec%cmat_exp)
        ub = ubound(rec%cmat_exp)
        do j = 1, n
            if( any(win(:,j) < lb) .or. any(win(:,j) + wdim - 1 > ub) )then
                vals(j) = CMPLX_ZERO
                cycle
            endif
            val = CMPLX_ZERO
            do iz = 1, wdim
                mz = win(3,j) + iz - 1
                do iy = 1, wdim
                    ky  = win(2,j) + iy - 1
                    wyz = w(iy,2,j) * w(iz,3,j)
                    do ix = 1, wdim
                        hx  = win(1,j) + ix - 1
                        val = val + rec%cmat_exp(hx,ky,mz) * (w(ix,1,j) * wyz)
                    end do
                end do
            end do
            if( l_conj(j) ) val = conjg(val)
            vals(j) = val
        end do
    end subroutine gather_samples

    !> Write the gathered plane samples vals into fpl_out at the plane positions of new_plane. fpl_out takes
    !! the reference plane's geometry and carries no ctfsq or transfer plane; it is zeroed only when that
    !! geometry (bounds, Nyquist) changes, so the out-of-disc remainder stays zero. With apply_ctf_amp the
    !! samples are multiplied by the reference plane's forward transfer, or sqrt(ctf^2) without one.
    subroutine exp_samples_put_plane( self, vals, fpl_ref, fpl_out, apply_ctf_amp )
        class(exp_samples), intent(in)    :: self
        complex,            intent(in)    :: vals(:)
        class(fplane_type), intent(in)    :: fpl_ref
        type(fplane_type),  intent(inout) :: fpl_out
        logical,            intent(in)    :: apply_ctf_amp
        logical :: l_realloc
        integer :: j, h, k
        if( .not. self%l_plane ) THROW_HARD('sample set has no plane positions; exp_samples_put_plane')
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_ref%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_ref%cmplx_plane)) .or. fpl_out%nyq /= fpl_ref%nyq
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_ref%cmplx_plane,1):ubound(fpl_ref%cmplx_plane,1), &
                &lbound(fpl_ref%cmplx_plane,2):ubound(fpl_ref%cmplx_plane,2)), source=CMPLX_ZERO)
        endif
        if( allocated(fpl_out%ctfsq_plane)    ) deallocate(fpl_out%ctfsq_plane)
        if( allocated(fpl_out%transfer_plane) ) deallocate(fpl_out%transfer_plane)
        fpl_out%frlims  = fpl_ref%frlims
        fpl_out%shconst = fpl_ref%shconst
        fpl_out%nyq     = fpl_ref%nyq
        do j = 1, self%n
            h = self%hk(1,j)
            k = self%hk(2,j)
            if( apply_ctf_amp )then
                if( allocated(fpl_ref%transfer_plane) )then
                    fpl_out%cmplx_plane(h,k) = fpl_ref%transfer_plane(h,k) * vals(j)
                else
                    fpl_out%cmplx_plane(h,k) = sqrt(max(0., fpl_ref%ctfsq_plane(h,k))) * vals(j)
                endif
            else
                fpl_out%cmplx_plane(h,k) = vals(j)
            endif
        end do
    end subroutine exp_samples_put_plane

    integer pure function exp_samples_get_n( self ) result( n )
        class(exp_samples), intent(in) :: self
        n = self%n
    end function exp_samples_get_n

    subroutine alloc_exp_samples( self, n, wdim )
        class(exp_samples), intent(inout) :: self
        integer,            intent(in)    :: n, wdim
        if( allocated(self%win) )then
            if( size(self%win,2) < n .or. size(self%w,1) /= wdim ) deallocate(self%win, self%w, self%l_conj)
        endif
        if( .not. allocated(self%win) ) allocate(self%win(3,max(1,n)), self%w(wdim,3,max(1,n)), self%l_conj(max(1,n)))
        self%n       = n
        self%wdim    = wdim
        self%l_plane = .false.
    end subroutine alloc_exp_samples

    subroutine exp_samples_kill( self )
        class(exp_samples), intent(inout) :: self
        if( allocated(self%win)    ) deallocate(self%win)
        if( allocated(self%w)      ) deallocate(self%w)
        if( allocated(self%l_conj) ) deallocate(self%l_conj)
        if( allocated(self%hk)     ) deallocate(self%hk)
        if( allocated(self%vbuf)   ) deallocate(self%vbuf)
        self%n       = 0
        self%wdim    = 0
        self%l_plane = .false.
    end subroutine exp_samples_kill

    !> Project this volume (cmat_exp current) into the Cartesian Fourier-plane storage of gen_fplane4rec:
    !! the central-section samples of fpl_ref (exp_samples new_plane) gathered and written into fpl_out
    !! (put_plane: no ctfsq or transfer plane; zeroed only when the geometry changes). samples is optional
    !! work storage a caller reuses particle by particle.
    subroutine project_fplane( self, o, fpl_ref, fpl_out, apply_ctf_amp, samples )
        class(reconstructor),        intent(in)    :: self
        class(ori),                  intent(inout) :: o
        class(fplane_type),          intent(in)    :: fpl_ref
        type(fplane_type),           intent(inout) :: fpl_out
        logical,           optional, intent(in)    :: apply_ctf_amp
        type(exp_samples), optional, intent(inout) :: samples
        type(exp_samples) :: work
        logical :: l_apply_ctf_amp
        if( .not. allocated(self%cmat_exp) ) THROW_HARD('expanded matrix does not exist; reconstructor project_fplane')
        if( .not. allocated(fpl_ref%cmplx_plane) .or. .not. allocated(fpl_ref%ctfsq_plane) )then
            THROW_HARD('reference Fourier plane does not exist; reconstructor project_fplane')
        endif
        l_apply_ctf_amp = .false.
        if( present(apply_ctf_amp) ) l_apply_ctf_amp = apply_ctf_amp
        ! callers keep cmat_exp current: expanding here would refresh the whole 3D lattice per particle.
        ! A caller-held sample set (one per thread) is reused across particles and allocates nothing
        if( present(samples) )then
            call samples%new_plane(self, o, fpl_ref)
            call exp_samples_gather_vbuf(samples, self)
            call samples%put_plane(samples%vbuf, fpl_ref, fpl_out, l_apply_ctf_amp)
        else
            call work%new_plane(self, o, fpl_ref)
            call exp_samples_gather_vbuf(work, self)
            call work%put_plane(work%vbuf, fpl_ref, fpl_out, l_apply_ctf_amp)
            call work%kill
        endif
    end subroutine project_fplane

    ! ---------------------------------------------------------------- multi-target insertion

    !> Batch KB insertion of nrecords observation-model planes (gen_fplane4rec observation_model) into
    !! K expanded numerators: target q receives conj(T)y of plane i scaled by data_w(q,i). The density is
    !! one array rho(P,lattice) beside the targets: P=1 shared (|T|^2), P=K diagonal (dens_w(q,q,i)|T|^2),
    !! P=K(K+1)/2 the packed upper triangle (dens_w(q,r,i)|T|^2 at (r(r-1))/2+q). Per-record geometry is
    !! derived serially (se%apply is not guaranteed thread-safe); the h-line sweep is threaded with a
    !! stride that keeps concurrent lines off any one interpolation cell. insert_plane_oversamp stays the
    !! single-target path of refinement.
    subroutine insert_planes_multi( recs, rho, se, orientations, fpls, data_w, dens_w, valid, nrecords )
        use simple_math, only: ceil_div, floor_div
        type(reconstructor), intent(inout) :: recs(:)
        real,                intent(inout) :: rho(:,:,:,:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: data_w(:,:), dens_w(:,:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(ori) :: o_sym
        type(kbinterpol) :: kbwin
        complex   :: comp_base, cmplx_raw
        real, allocatable :: rotmats(:,:,:,:), data_w_sp(:,:), dens_w_packed(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        integer, parameter :: WDIM = 2 * ceiling(KBWINSZ - 0.5) + 1
        real      :: loc(3), hrow(3), ctfsq_raw, ww, r2(3), pf2, eps_norm, inv_wdim
        real      :: wx(WDIM), wy(WDIM), wz(WDIM)
        integer   :: win(3,2), h, k, l, nsym, isym, iwinsz, stride
        integer   :: ix, iy, iz, hx, ky, mz, q, r, i, ncomp, ipair
        integer   :: h_sq, k_max_h, k_lo, k_hi, ih, ik, im, nyq_eff
        integer   :: exp_lb(3), exp_ub(3), exp_shape(3), npairs
        logical   :: l_shared, l_diagonal
        ncomp = size(recs)
        if( ncomp <= 0 .or. nrecords <= 0 ) return
        npairs     = (ncomp * (ncomp + 1)) / 2
        l_diagonal = size(rho,1) == ncomp
        l_shared   = size(rho,1) == 1 .and. .not. l_diagonal
        if( size(orientations) < nrecords .or. size(fpls) < nrecords .or. size(valid) < nrecords )then
            THROW_HARD('record array smaller than batch; insert_planes_multi')
        endif
        if( size(data_w,1) < ncomp .or. size(data_w,2) < nrecords .or. size(dens_w,1) < ncomp .or. &
            &size(dens_w,2) < ncomp .or. size(dens_w,3) < nrecords )then
            THROW_HARD('weight array smaller than batch; insert_planes_multi')
        endif
        if( .not. allocated(recs(1)%cmat_exp) ) THROW_HARD('expanded matrix does not exist; insert_planes_multi')
        exp_lb    = lbound(recs(1)%cmat_exp)
        exp_ub    = ubound(recs(1)%cmat_exp)
        exp_shape = shape(recs(1)%cmat_exp)
        if( (.not. l_shared .and. .not. l_diagonal .and. size(rho,1) < npairs) .or. size(rho,2) < exp_shape(1) .or. &
            &size(rho,3) < exp_shape(2) .or. size(rho,4) < exp_shape(3) )then
            THROW_HARD('density array shape mismatch; insert_planes_multi')
        endif
        if( recs(1)%wdim /= WDIM ) THROW_HARD('unexpected KB window dimension; insert_planes_multi')
        kbwin    = recs(1)%kbwin
        nsym     = se%get_nsym()
        iwinsz   = ceiling(KBWINSZ - 0.5)
        ! source h-lines in one colour must map to non-overlapping windows for every rotation:
        ! sqrt(3)*WDIM guarantees a full window width along at least one target axis
        stride   = ceiling(sqrt(3.0) * real(WDIM))
        pf2      = real(OSMPL_PAD_FAC*OSMPL_PAD_FAC)
        eps_norm = epsilon(1.0)
        inv_wdim = 1.0 / real(WDIM)
        allocate(rotmats(3,3,nsym,nrecords), data_w_sp(ncomp,nrecords), source=0.)
        if( .not. l_shared .and. .not. l_diagonal ) allocate(dens_w_packed(npairs,nrecords), source=0.)
        allocate(fpllims(2,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            if( .not. allocated(fpls(i)%transfer_plane) )then
                THROW_HARD('forward transfer plane does not exist; insert_planes_multi')
            endif
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims(1,1,i) = ceil_div (fpls(i)%frlims(1,1), OSMPL_PAD_FAC)
            fpllims(1,2,i) = floor_div(fpls(i)%frlims(1,2), OSMPL_PAD_FAC)
            fpllims(2,1,i) = ceil_div (fpls(i)%frlims(2,1), OSMPL_PAD_FAC)
            fpllims(2,2,i) = floor_div(fpls(i)%frlims(2,2), OSMPL_PAD_FAC)
            nyq_eff = recs(1)%nyq
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq / OSMPL_PAD_FAC))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            do q = 1, ncomp
                data_w_sp(q,i) = real(data_w(q,i))
            end do
            if( .not. l_shared .and. .not. l_diagonal )then
                do r = 1, ncomp
                    do q = 1, r
                        dens_w_packed((r*(r-1))/2 + q,i) = real(dens_w(q,r,i))
                    end do
                end do
            endif
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,h,k,l,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,ctfsq_raw,&
        !$omp& comp_base,wx,wy,wz,ww,win,loc,hrow,r2,isym,ix,iy,iz,hx,ky,mz,ih,ik,im,q,ipair) proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            do isym = 1, nsym
                r2 = rotmats(2,:,isym,i)
                do l = 0, stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i) - h_sq)))
                        k_lo = max(fpllims(2,1,i), -k_max_h)
                        k_hi = min(fpllims(2,2,i),  k_max_h)
                        hrow = real(h) * rotmats(1,:,isym,i)
                        do k = k_lo, k_hi
                            loc = hrow + real(k) * r2
                            ! native lattice, k<=0 stored; Friedel symmetry for k>0
                            if( k <= 0 )then
                                cmplx_raw = conjg(fpls(i)%transfer_plane(h,k)) * fpls(i)%cmplx_plane(h,k)
                                ctfsq_raw = fpls(i)%ctfsq_plane(h,k)
                            else
                                cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-h,-k)) * fpls(i)%cmplx_plane(-h,-k))
                                ctfsq_raw = fpls(i)%ctfsq_plane(-h,-k)
                            endif
                            if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. ctfsq_raw <= TINY ) cycle
                            win(:,1) = nint(loc) - iwinsz
                            win(:,2) = nint(loc) + iwinsz
                            ! cmat_exp stores h>=0 as the independent Friedel half
                            if( win(1,2) < 0 ) cycle
                            if( any(win(:,1) < exp_lb) .or. any(win(:,2) > exp_ub) ) cycle
                            comp_base = pf2 * cmplx_raw
                            call kb_weights(loc, win(:,1), wx, wy, wz)
                            do iz = 1, WDIM
                                mz = win(3,1) + iz - 1
                                im = mz - exp_lb(3) + 1
                                do iy = 1, WDIM
                                    ky = win(2,1) + iy - 1
                                    ik = ky - exp_lb(2) + 1
                                    do ix = 1, WDIM
                                        hx = win(1,1) + ix - 1
                                        ih = hx - exp_lb(1) + 1
                                        ww = wx(ix) * (wy(iy) * wz(iz))
                                        do q = 1, ncomp
                                            recs(q)%cmat_exp(hx,ky,mz) = recs(q)%cmat_exp(hx,ky,mz) + &
                                                &(data_w_sp(q,i) * comp_base) * ww
                                        end do
                                        if( l_shared )then
                                            rho(1,ih,ik,im) = rho(1,ih,ik,im) + ctfsq_raw * ww
                                        else if( l_diagonal )then
                                            do q = 1, ncomp
                                                rho(q,ih,ik,im) = rho(q,ih,ik,im) + real(dens_w(q,q,i)) * ctfsq_raw * ww
                                            end do
                                        else
                                            !$omp simd
                                            do ipair = 1, npairs
                                                rho(ipair,ih,ik,im) = rho(ipair,ih,ik,im) + &
                                                    &(dens_w_packed(ipair,i) * ctfsq_raw) * ww
                                            end do
                                        endif
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !$omp end do
                end do
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, data_w_sp, fpllims, nyq_disks)
        if( allocated(dens_w_packed) ) deallocate(dens_w_packed)

    contains

        !> normalized separable KB weights at loc from the window's lower corner (apod_fast, as the
        !! refinement's insert_plane_oversamp)
        subroutine kb_weights( loc, win_lo, wx, wy, wz )
            real,    intent(in)  :: loc(3)
            integer, intent(in)  :: win_lo(3)
            real,    intent(out) :: wx(:), wy(:), wz(:)
            integer :: j
            real    :: base(3), ww3(3), sx, sy, sz
            base = real(win_lo) - loc
            do j = 1, WDIM
                ww3   = kbwin%apod_fast(base + real(j-1))
                wx(j) = ww3(1)
                wy(j) = ww3(2)
                wz(j) = ww3(3)
            end do
            sx = sum(wx); sy = sum(wy); sz = sum(wz)
            if( abs(sx) > eps_norm )then; wx = wx/sx; else; wx = inv_wdim; endif
            if( abs(sy) > eps_norm )then; wy = wy/sy; else; wy = inv_wdim; endif
            if( abs(sz) > eps_norm )then; wz = wz/sz; else; wz = inv_wdim; endif
        end subroutine kb_weights

    end subroutine insert_planes_multi

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
        real              :: fudge, cc, invtau2
        integer           :: h, k, m, sh, phys(3), sz, reslim_ind
        sz = size(fsc)
        allocate(ssnr(0:sz), rsum(0:sz), cnt(0:sz), tau2(0:sz), sig2(0:sz))
        rsum  = 0.d0
        cnt   = 0
        ssnr  = 0.0
        tau2  = 0.0
        sig2  = 0.0
        fudge = self%p_ptr %tau
        ! SSNR
        do k = 1,sz
            cc      = max(0.001,fsc(k))
            cc      = min(0.999,cc)
            ssnr(k) = cc / (1.-cc)
        enddo
        ! Noise
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,sh) proc_bind(close) reduction(+:cnt,rsum)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > sz ) cycle
                    phys     = self%comp_addr_phys(h, k, m)
                    cnt(sh)  = cnt(sh) + 1
                    rsum(sh) = rsum(sh) + real(self%rho(phys(1),phys(2),phys(3)),dp)
                enddo
            enddo
        enddo
        !$omp end parallel do
        where( rsum > 1.d-10 )
            sig2 = real(real(cnt,dp) / rsum)
        else where
            sig2 = 0.0
        end where
        ! Signal
        tau2 = ssnr * sig2
        ! add Tau2 inverse to denominator
        ! because signal assumed infinite at very low resolution there is no addition
        reslim_ind = max(6, calc_fourier_index(self%p_ptr %hp, self%p_ptr %box_crop, self%p_ptr %smpd_crop))
        !$omp parallel do collapse(3) default(shared) schedule(static)&
        !$omp private(h,k,m,phys,sh,invtau2) proc_bind(close)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                do m = self%lims(3,1),self%lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( (sh < reslim_ind) .or. (sh > sz) ) cycle
                    phys = self%comp_addr_phys(h, k, m)
                    if( tau2(sh) > TINY)then
                        invtau2 = 1.0/(fudge*tau2(sh))
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

    subroutine kill_gridding_half_restore( self )
        class(gridding_half_restore), intent(inout) :: self
        call self%base%kill
        call self%final%kill
    end subroutine kill_gridding_half_restore

    !> Complete destructor for the extended image. A bare reconstructor%kill
    !! must release the inherited image, rho, and expanded accumulator storage.
    subroutine kill_reconstructor( self )
        class(reconstructor), intent(inout) :: self
        call self%dealloc_rho
        self%p_ptr => null()
        self%shconst_rec = 0.
        self%wdim        = 0
        self%nyq         = 0
        self%sh_lim      = 0
        self%ldim_img    = 0
        self%ldim_exp    = 0
        self%lims        = 0
        self%rho_shape   = 0
        self%cyc_lims    = 0
        call self%image%kill
    end subroutine kill_reconstructor

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
        if( allocated(self%invenv1d) ) deallocate(self%invenv1d)
        if( self%rho_allocated )then
            call fftwf_free(self%kp)
            self%rho => null()
            self%rho_allocated = .false.
        endif
    end subroutine dealloc_rho

end module simple_reconstructor
