!@descr: 3D reconstruction from projections using convolution interpolation (gridding)
module simple_reconstructor
use simple_core_module_api
use simple_fftw3
use simple_image,      only: image
use simple_parameters, only: parameters
implicit none

public :: reconstructor
private
#include "simple_local_flags.inc"

type, extends(image) :: reconstructor
    private
    class(parameters),  pointer :: p_ptr => null()              !< pointer to parameters object
    type(kbinterpol)            :: kbwin                        !< window function object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex, allocatable, public :: cmat_exp(:,:,:)             !< Fourier components of expanded reconstructor, only made public for the sake of GPU implementation
    real,    allocatable, public :: rho_exp(:,:,:)              !< sampling+CTF**2 density of expanded reconstructor, only made public for the sake of GPU implementation
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
    procedure          :: insert_plane_oversamp
    procedure          :: insert_plane_oversamp_opt
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    procedure          :: expand_exp
    procedure          :: project_polar
    procedure          :: project_fplane
    ! SUMMATION
    procedure          :: sum_reduce
    procedure          :: add_invtausq2rho
    procedure          :: add_conical_invtausq2rho
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

contains

    ! CONSTRUCTORS

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
        self%phaseplate     = spproj%has_phaseplate(self%p_ptr %oritype)
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

    !> Experimental variant of insert_plane_oversamp for benchmarking.
    !! It preserves the original h-strided OpenMP race-avoidance scheme, but
    !! scalarizes rotation, inlines fplane access, and updates the two expanded
    !! matrices in one explicit separable-KB stencil pass.
    subroutine insert_plane_oversamp( self, se, o, fpl, compact_source )
        use simple_math, only: ceil_div, floor_div
        class(reconstructor), intent(inout) :: self
        class(sym),           intent(inout) :: se
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl
        logical, optional,    intent(in)    :: compact_source
        type(ori) :: o_sym
        complex   :: comp, cmplx_raw
        real      :: rotmats(se%get_nsym(),3,3), loc(3), hrow(3), ctfval
        real      :: wx(self%wdim), wy(self%wdim), wz(self%wdim), ww, ctfsq_raw
        real      :: r11, r12, r13, r21, r22, r23
        integer   :: win(2, 3), h, k, l, nsym, isym, iwinsz, stride, fpllims_pd(3, 2)
        integer   :: fpllims(3, 2), hp, kp, pf, ix, iy, iz, hx, ky, mz
        integer   :: nyq_disk, h_sq, k_max_h, k_lo, k_hi
        real      :: source_scale, eps_norm, inv_wdim
        logical   :: l_compact_source
        ! window size
        iwinsz = ceiling(KBWINSZ - 0.5)
        ! stride along h dimension for interpolation: all threads are at least
        ! wdim pixels away from each other to avoid race conditions
        stride = self%wdim
        ! setup rotation matrices
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 ) then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            end do
        endif
        ! Native (unpadded) iteration limits so that hp=h*pf and kp=k*pf are in-bounds
        fpllims_pd      = fpl%frlims
        l_compact_source = .false.
        if( present(compact_source) ) l_compact_source = compact_source
        if( l_compact_source )then
            ! The source is already a native-grid 2D KB numerator/CTF^2 sum.
            ! Its padded-FFT amplitude scaling was applied during 2D assembly.
            pf           = 1
            source_scale = 1.0
        else
            pf           = OSMPL_PAD_FAC
            source_scale = real(pf*pf)
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
        !$omp& ctfsq_raw,ctfval,wx,wy,wz,ww,win,loc,hrow,hp,kp,r11,r12,r13,r21,r22,r23,&
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
                    hp      = h * pf
                    hrow(1) = real(h) * r11
                    hrow(2) = real(h) * r12
                    hrow(3) = real(h) * r13
                    do k = k_lo, k_hi
                        kp = k * pf
                        ! gen_fplane4rec stores only k<=0; use Friedel symmetry for kp>0.
                        if( kp <= 0 )then
                            cmplx_raw = fpl%cmplx_plane(hp,kp)
                            ctfsq_raw = fpl%ctfsq_plane(hp,kp)
                        else
                            cmplx_raw = conjg(fpl%cmplx_plane(-hp,-kp))
                            ctfsq_raw = fpl%ctfsq_plane(-hp,-kp)
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

    subroutine insert_plane_oversamp_opt( self, se, o, fpl )
        use simple_math,       only: ceil_div, floor_div
        use simple_kbinterpol, only: apod_kb15_a2
        class(reconstructor), intent(inout) :: self
        class(sym),           intent(inout) :: se
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl
        integer, parameter :: WDIM   = 3
        integer, parameter :: STRIDE = WDIM
        type(ori) :: o_sym
        real      :: rotmats(3,3,se%get_nsym()), pf2
        integer   :: fpllims(3, 2), fpllims_pd(3, 2)
        integer   :: clb3D(3), cdim3D(3), clb2D(2), cdim2D(2)
        integer   :: nyq_disk, jsym, nsym, iwinsz
        if( self%wdim /= WDIM ) then
            THROW_HARD('insert_plane_oversamp_opt only implemented for wdim=3!')
        endif
        ! window size
        iwinsz = ceiling(KBWINSZ - 0.5)
        ! Setup rotation matrices
        nsym = se%get_nsym()
        rotmats(:,:,1) = o%get_mat()
        if( nsym > 1 ) then
            do jsym = 2, nsym
                call se%apply(o, jsym, o_sym)
                rotmats(:,:,jsym) = o_sym%get_mat()
            end do
        endif
        call o_sym%kill
        ! Native iteration limits so that hp=h*pf and kp=k*pf are in-bounds
        fpllims_pd   = fpl%frlims
        pf2          = real(OSMPL_PAD_FAC**2)
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), OSMPL_PAD_FAC)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), OSMPL_PAD_FAC)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), OSMPL_PAD_FAC)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), OSMPL_PAD_FAC)
        ! bit-equivalent to nint(sqrt(h*h+k*k)) > nyq
        nyq_disk = self%nyq * (self%nyq + 1)
        ! 3D arrays boundaries
        clb3D  = lbound(self%cmat_exp)
        cdim3D = ubound(self%cmat_exp) - clb3D + 1
        ! 2D planes boundaries
        clb2D  = lbound(fpl%cmplx_plane)
        cdim2D = ubound(fpl%cmplx_plane) - clb2D + 1
        ! KB interpolation / insertion with KB interpolation
        call kernel(self%cmat_exp, self%rho_exp, fpl%cmplx_plane, fpl%ctfsq_plane)
        contains

            subroutine kernel( cmatexp, rhoexp, fcomp_plane, ctfsq_plane )
                complex(sp),      intent(inout) :: cmatexp(cdim3D(1),cdim3D(2),cdim3D(3))
                real(sp),         intent(inout) :: rhoexp(cdim3D(1),cdim3D(2),cdim3D(3))
                complex(sp),      intent(in)    :: fcomp_plane(cdim2D(1), cdim2D(2))
                real(sp),         intent(in)    :: ctfsq_plane(cdim2D(1), cdim2D(2))
                complex(sp) :: comp
                real    :: wx(WDIM), wy(WDIM), wz(WDIM), base(3), loc(3)
                real    :: r21, r22, r23, sx, sy, sz, comp_scale, ctfsq, wyz
                integer :: h,k,l, h_sq, k_max_h, k_lo,k_hi, hp,kp,hpb,kpb, iy,iz, ky,mz, i
                integer :: win(3, 2), isym
                comp_scale = pf2
                !$omp parallel default(shared) private(h,k,l,h_sq,k_max_h,k_lo,k_hi,comp,&
                !$omp& ctfsq,wx,wy,wz,sx,sy,sz,i,win,loc,r21,r22,r23,isym,iy,iz,ky,mz,wyz,&
                !$omp& base,hp,kp,hpb,kpb) proc_bind(close)
                do isym = 1, nsym
                    r21 = rotmats(2,1,isym); r22 = rotmats(2,2,isym); r23 = rotmats(2,3,isym)
                    do l = 0, STRIDE-1
                        !$omp do schedule(static,1)
                        do h = fpllims(1,1)+l, fpllims(1,2), STRIDE
                            h_sq = h*h
                            if( h_sq > nyq_disk ) cycle
                            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
                            k_lo    = max(fpllims(2,1), -k_max_h)
                            k_hi    = min(fpllims(2,2),  k_max_h)
                            loc     = real(h) * rotmats(1,1:3,isym)
                            loc     = loc + real(k_lo-1) * [r21, r22, r23]
                            ! padded h coordinate
                            hp = h * OSMPL_PAD_FAC
                            do k = k_lo, k_hi
                                ! rotation
                                loc(1) = loc(1) + r21
                                loc(2) = loc(2) + r22
                                loc(3) = loc(3) + r23
                                win(:,1) = nint(loc)
                                win(:,2) = win(:,1) + iwinsz
                                win(:,1) = win(:,1) - iwinsz
                                ! no need to update outside the non-redundant Friedel limits consistent with compress_exp
                                if( win(1,2) < self%lims(1,1) ) cycle
                                ! padded coordinate
                                kp = k * OSMPL_PAD_FAC
                                ! gen_fplane4rec stores only k<=0; use Friedel symmetry for kp>0.
                                if( kp <= 0 )then
                                    hpb   = hp - clb2D(1) + 1
                                    kpb   = kp - clb2D(2) + 1
                                    comp  = fcomp_plane(hpb,kpb)
                                    ctfsq = ctfsq_plane(hpb,kpb)
                                else
                                    hpb   = -hp - clb2D(1) + 1
                                    kpb   = -kp - clb2D(2) + 1
                                    comp  = conjg(fcomp_plane(hpb,kpb))
                                    ctfsq =       ctfsq_plane(hpb,kpb)
                                endif
                                if( abs(real(comp)) + abs(aimag(comp)) <= TINY .and. ctfsq <= TINY ) cycle
                                ! FFTW padding scaling
                                comp  = comp_scale * comp
                                ! precompute and normalize weights
                                base = real(win(:,1)) - loc
                                sx = 0.0; sy = 0.0; sz = 0.0
                                do i = 1, WDIM
                                    wx(i) = apod_kb15_a2(base(1) + real(i-1))
                                    wy(i) = apod_kb15_a2(base(2) + real(i-1))
                                    wz(i) = apod_kb15_a2(base(3) + real(i-1))
                                    sx = sx + wx(i)
                                    sy = sy + wy(i)
                                    sz = sz + wz(i)
                                end do
                                wx = wx / sx; wy = wy / sy; wz = wz / sz
                                ! Arrays updates
                                win(:,1) = win(:,1) - clb3D ! adjust for 1-based indexing
                                do iz = 1, WDIM
                                    mz = win(3,1) + iz
                                    do iy = 1, WDIM
                                        ky  = win(2,1) + iy
                                        wyz = wy(iy) * wz(iz)
                                        cmatexp(win(1,1)+1:win(1,1)+WDIM, ky, mz) = &
                                            &cmatexp(win(1,1)+1:win(1,1)+WDIM, ky, mz) + (wyz*comp) * wx(:WDIM)
                                        rhoexp(win(1,1)+1:win(1,1)+WDIM, ky, mz)  = &
                                            &rhoexp(win(1,1)+1:win(1,1)+WDIM, ky, mz)  + (wyz*ctfsq) * wx(:WDIM)
                                    end do
                                end do
                            end do
                        end do
                        !$omp end do
                    end do
                end do
                !$omp end parallel
            end subroutine kernel

    end subroutine insert_plane_oversamp_opt

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

    subroutine project_polar( self, eulspace, nspace, kfromto, polar_x, polar_y, pfts_state, ctf2_state )
        class(reconstructor), intent(inout) :: self
        class(oris),          intent(inout) :: eulspace
        integer,              intent(in)    :: nspace, kfromto(2)
        real(sp),             intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(dp),          intent(out)   :: pfts_state(:,:,:)
        real(dp),             intent(out)   :: ctf2_state(:,:,:)
        integer :: iproj, irot, k, kloc, kfrom, kto, pftsz, noris
        real    :: e_rotmat(3,3), loc(3), px, py
        kfrom = kfromto(1)
        kto   = kfromto(2)
        pftsz = size(polar_x,1)
        noris = eulspace%get_noris()
        if( .not. allocated(self%cmat_exp) .or. .not. allocated(self%rho_exp) )then
            THROW_HARD('expanded matrices do not exist; reconstructor project_polar')
        endif
        if( kfrom < 1 .or. kfrom > kto ) THROW_HARD('invalid kfromto; reconstructor project_polar')
        if( noris < nspace )             THROW_HARD('eulspace smaller than nspace; reconstructor project_polar')
        if( size(polar_x,2) /= (kto-kfrom+1) ) THROW_HARD('polar_x k-span mismatch; reconstructor project_polar')
        if( size(polar_y,1) /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )then
            THROW_HARD('polar coordinate shape mismatch; reconstructor project_polar')
        endif
        if( size(pfts_state,1) /= pftsz .or. size(pfts_state,2) /= size(polar_x,2) .or. &
            &size(pfts_state,3) /= nspace )then
            THROW_HARD('pfts_state shape mismatch; reconstructor project_polar')
        endif
        if( any(shape(ctf2_state) /= shape(pfts_state)) )then
            THROW_HARD('ctf2_state shape mismatch; reconstructor project_polar')
        endif
        call self%expand_exp
        !$omp parallel do default(shared) private(iproj,e_rotmat,irot,k,kloc,loc,px,py) schedule(static) proc_bind(close)
        do iproj = 1, nspace
            e_rotmat = eulspace%get_mat(iproj)
            do k = kfrom, kto
                kloc = k - kfrom + 1
                do irot = 1, pftsz
                    px     = real(polar_x(irot,kloc))
                    py     = real(polar_y(irot,kloc))
                    loc(1) = px*e_rotmat(1,1) + py*e_rotmat(2,1)
                    loc(2) = px*e_rotmat(1,2) + py*e_rotmat(2,2)
                    loc(3) = px*e_rotmat(1,3) + py*e_rotmat(2,3)
                    pfts_state(irot,kloc,iproj) = cmplx(interp_cmat_exp(self, loc), kind=dp)
                    ctf2_state(irot,kloc,iproj) = real(max(0., interp_rho_exp(self, loc)), kind=dp)
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine project_polar

    !> Project this volume into the same Cartesian Fourier-plane storage used
    !! by gen_fplane4rec.  If apply_ctf_amp is true, the model is multiplied by
    !! the stored shift/CTF/sigma transfer plane when available, falling back to
    !! sqrt(ctf^2) only for legacy fplanes that do not carry the transfer.
    subroutine project_fplane( self, o, fpl_ref, fpl_out, apply_ctf_amp )
        use simple_math, only: ceil_div, floor_div
        class(reconstructor), intent(inout) :: self
        class(ori),           intent(inout) :: o
        class(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),    intent(inout) :: fpl_out
        logical, optional,    intent(in)    :: apply_ctf_amp
        real    :: rotmat(3,3), loc(3), hrow(3), ctfamp
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        logical :: l_apply_ctf_amp, l_realloc
        if( .not. allocated(self%cmat_exp) .or. .not. allocated(self%rho_exp) )then
            THROW_HARD('expanded matrices do not exist; reconstructor project_fplane')
        endif
        if( .not. allocated(fpl_ref%cmplx_plane) .or. .not. allocated(fpl_ref%ctfsq_plane) )then
            THROW_HARD('reference Fourier plane does not exist; reconstructor project_fplane')
        endif
        l_apply_ctf_amp = .false.
        if( present(apply_ctf_amp) ) l_apply_ctf_amp = apply_ctf_amp
        fpl_out%frlims  = fpl_ref%frlims
        fpl_out%shconst = fpl_ref%shconst
        fpl_out%nyq     = fpl_ref%nyq
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_ref%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_ref%cmplx_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_ref%cmplx_plane,1):ubound(fpl_ref%cmplx_plane,1), &
                &lbound(fpl_ref%cmplx_plane,2):ubound(fpl_ref%cmplx_plane,2)))
        endif
        fpl_out%cmplx_plane = CMPLX_ZERO
        l_realloc = .not. allocated(fpl_out%ctfsq_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%ctfsq_plane) /= lbound(fpl_ref%ctfsq_plane)) .or. &
                &any(ubound(fpl_out%ctfsq_plane) /= ubound(fpl_ref%ctfsq_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%ctfsq_plane) ) deallocate(fpl_out%ctfsq_plane)
            allocate(fpl_out%ctfsq_plane(lbound(fpl_ref%ctfsq_plane,1):ubound(fpl_ref%ctfsq_plane,1), &
                &lbound(fpl_ref%ctfsq_plane,2):ubound(fpl_ref%ctfsq_plane,2)))
        endif
        fpl_out%ctfsq_plane = fpl_ref%ctfsq_plane
        if( allocated(fpl_ref%transfer_plane) )then
            l_realloc = .not. allocated(fpl_out%transfer_plane)
            if( .not. l_realloc )then
                l_realloc = any(lbound(fpl_out%transfer_plane) /= lbound(fpl_ref%transfer_plane)) .or. &
                    &any(ubound(fpl_out%transfer_plane) /= ubound(fpl_ref%transfer_plane))
            endif
            if( l_realloc )then
                if( allocated(fpl_out%transfer_plane) ) deallocate(fpl_out%transfer_plane)
                allocate(fpl_out%transfer_plane(lbound(fpl_ref%transfer_plane,1):ubound(fpl_ref%transfer_plane,1), &
                    &lbound(fpl_ref%transfer_plane,2):ubound(fpl_ref%transfer_plane,2)))
            endif
            fpl_out%transfer_plane = fpl_ref%transfer_plane
        else
            if( allocated(fpl_out%transfer_plane) ) deallocate(fpl_out%transfer_plane)
        endif
        ! Callers using this per-particle projector must keep cmat_exp current.
        ! Expanding here would refresh the full 3D Fourier lattice for every
        ! particle projection.
        rotmat      = o%get_mat()
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl_ref%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = self%nyq
        if( fpl_ref%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl_ref%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            hrow(1) = real(h) * rotmat(1,1)
            hrow(2) = real(h) * rotmat(1,2)
            hrow(3) = real(h) * rotmat(1,3)
            do k = k_lo, k_hi
                kp     = k * pf
                loc(1) = hrow(1) + real(k) * rotmat(2,1)
                loc(2) = hrow(2) + real(k) * rotmat(2,2)
                loc(3) = hrow(3) + real(k) * rotmat(2,3)
                fpl_out%cmplx_plane(hp,kp) = interp_cmat_exp(self, loc)
                if( l_apply_ctf_amp )then
                    if( allocated(fpl_ref%transfer_plane) )then
                        fpl_out%cmplx_plane(hp,kp) = fpl_ref%transfer_plane(hp,kp) * fpl_out%cmplx_plane(hp,kp)
                    else
                        ctfamp = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp)))
                        fpl_out%cmplx_plane(hp,kp) = ctfamp * fpl_out%cmplx_plane(hp,kp)
                    endif
                endif
            end do
        end do
    end subroutine project_fplane

    pure function interp_cmat_exp( self, loc ) result( comp )
        class(reconstructor), intent(in) :: self
        real,                 intent(in) :: loc(3)
        complex :: comp
        real    :: w(1:self%wdim,1:self%wdim,1:self%wdim)
        integer :: iwinsz, win(2,3)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz
        call self%kbwin%apod_mat_3d(loc, iwinsz, self%wdim, w)
        comp = sum(w * self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)))
    end function interp_cmat_exp

    pure real function interp_rho_exp( self, loc ) result( rho )
        class(reconstructor), intent(in) :: self
        real,                 intent(in) :: loc(3)
        real    :: w(1:self%wdim,1:self%wdim,1:self%wdim)
        integer :: iwinsz, win(2,3)
        iwinsz   = ceiling(KBWINSZ - 0.5)
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + iwinsz
        win(1,:) = win(1,:) - iwinsz
        call self%kbwin%apod_mat_3d(loc, iwinsz, self%wdim, w)
        rho = sum(w * self%rho_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)))
    end function interp_rho_exp

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

    subroutine add_conical_invtausq2rho( self, cones_obj )
        use simple_fsc, only: fsc_area_score_result
        class(reconstructor),        intent(inout) :: self !< this instance
        type(fsc_area_score_result), intent(in)    :: cones_obj
        real,     allocatable :: tau2(:,:)
        real(dp), allocatable :: cctfsq_sums(:,:)
        real(dp) :: ckcl, norm, dot, sig2
        real     :: cos_half, fudge, cc, tau2_thr, invtau2, snr, maxdot
        integer  :: phys(3), logi(3), maxdir
        integer  :: reslim_ind, nfreqs, ndirs, idir, freqlim, ifreq, h,i,j,k,l, sh
        freqlim  = self%get_filtsz()
        if( freqlim /= cones_obj%nfreqs ) then
            THROW_HARD('invalid freqlim; reconstructor add_conical_invtausq2rho')
        endif
        cos_half  = cos(deg2rad(real(cones_obj%cone_half_angle_deg, kind=dp)))
        nfreqs    = cones_obj%nfreqs
        ndirs     = cones_obj%ndirs
        fudge     = self%p_ptr %tau
        ! Conical spectral CTF2 sums
        allocate(cctfsq_sums(nfreqs,ndirs),source=0.d0)
        !$omp parallel default(shared) private(idir,l,k,ckcl,h,norm,sh,dot,phys) proc_bind(close)
        !$omp do schedule(guided) collapse(3) reduction(+:cctfsq_sums)
        do idir = 1, ndirs
            do l = self%lims(3,1), self%lims(3,2)
                do k = self%lims(2,1), self%lims(2,2)
                    ckcl = real(k,kind=dp) * cones_obj%dirs(2,idir) + real(l,kind=dp) * cones_obj%dirs(3,idir)
                    do h = self%lims(1,1), self%lims(1,2)
                        norm = sqrt(real(h*h + k*k + l*l, kind=dp))
                        if( norm <= 0.0d0 ) cycle
                        sh = nint(norm)
                        if( (sh == 0) .or. (sh > freqlim) ) cycle
                        dot = abs((real(h,kind=dp)*cones_obj%dirs(1,idir) + ckcl) / norm)
                        if( dot < cos_half ) cycle
                        phys = self%comp_addr_phys(h,k,l)
                        cctfsq_sums(sh,idir) = cctfsq_sums(sh,idir) + real(self%rho(phys(1),phys(2),phys(3)),dp)
                    end do
                enddo
            end do
        end do
        !$omp end do
        !$omp end parallel
        ! Conical spectral SNR & Tau2
        allocate(tau2(nfreqs,ndirs), source=0.0)
        do idir = 1, ndirs
            do ifreq = 1, nfreqs
                ! SSNR
                cc  = min(0.999, max(0.001, cones_obj%cfsc(ifreq,idir)))
                snr = cc / (1.-cc)
                ! Tau2
                sig2 = 0.0
                if( cctfsq_sums(ifreq,idir) > 1.d-10 )then
                    sig2 = real(real(cones_obj%counts(ifreq,idir),dp) / cctfsq_sums(ifreq,idir))
                end if
                tau2(ifreq,idir) = snr * sig2
            end do
        enddo
        ! add Tau2 inverse to denominator
        ! because signal assumed infinite at very low resolution there is no addition
        reslim_ind = max(6, calc_fourier_index(self%p_ptr %hp, self%p_ptr %box_crop, self%p_ptr %smpd_crop))
        !$omp parallel do default(shared) schedule(guided) proc_bind(close) collapse(3)&
        !$omp& private(i,j,k,logi,norm,sh,tau2_thr,idir,maxdir,dot,maxdot,invtau2)
        do k = 1, self%rho_shape(3)
            do j = 1, self%rho_shape(2)
                do i = 1, self%rho_shape(1)
                    logi = self%comp_addr_logi(i,j,k)
                    norm = sqrt(real(sum(logi*logi), kind=dp))
                    if( norm <= 0.0d0 ) cycle
                    sh = nint(norm)
                    if( (sh < reslim_ind) .or. (sh > freqlim) ) cycle
                    ! Nearest cone
                    maxdir   = 0
                    maxdot   = -1.
                    do idir = 1, ndirs
                        dot = abs((real(logi(1),kind=dp) * cones_obj%dirs(1,idir) + &
                            &real(logi(2),kind=dp) * cones_obj%dirs(2,idir) + &
                            &real(logi(3),kind=dp) * cones_obj%dirs(3,idir)) / norm)
                        if( dot > maxdot ) then
                            maxdot = dot
                            maxdir = idir
                        endif
                    enddo
                    tau2_thr = tau2(sh,maxdir)
                    ! CTF2 + 1/Tau2
                    if( tau2_thr > TINY)then
                        invtau2 = 1.0 / (fudge * tau2_thr)
                    else
                        invtau2 = min(1.e3, 1.e3 * self%rho(i,j,k))
                    endif
                    self%rho(i,j,k) = self%rho(i,j,k) + invtau2
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine add_conical_invtausq2rho

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
