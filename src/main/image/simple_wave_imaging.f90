!@descr: the 2D image formation of one particle under a single-scatter wave model, its linearisation and adjoint
! (doc/implementation_notes/planned/ghostbuster_operators.md, sections 3.2, 3.3 and 3.5); the lens and detection of
! a given exit wave (multislice); the detector step of simulate_scattering
module simple_wave_imaging
use simple_defs,       only: PI
use simple_type_defs,  only: ctfparams, ctfvars, CTFFLAG_NO, CTFFLAG_FLIP
use simple_error,      only: simple_exception
use simple_ctf,        only: ctf
use simple_math_ctf,   only: electron_interaction_constant
use simple_wave_field, only: wave_field, wave_freq
use simple_image,      only: image
use simple_rnd,        only: poisson_rnd, gasdev
implicit none
private
#include "simple_local_flags.inc"

public :: wave_imaging, WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE, detect_intensity

integer, parameter :: WAVE_PROJECTION = 1   !< exp(c u) of the central section
integer, parameter :: WAVE_FIRST_BORN = 2   !< 1 + c u
integer, parameter :: WAVE_RYTOV      = 3   !< exp(c u) of the paraboloid section
integer, parameter :: WAVE_MSLICE = 4   !< the exit wave is given (simple_wave_multislice): lens and detection

!> One thread's image formation. The input s(h,k) is the native-disk spectrum (-R <= h,k <= R,
!! h^2+k^2 <= R^2, R = nat/2) of the specimen's projected potential, in volts times voxels and in the units
!! of an unnormalised native DFT (the projector's gathered section times nat^3 and the shift phase). On the
!! fine grid of n = osf*nat points per side:
!!
!!   a     = c bwd(embed(P s)),          c = i sigma a_px (sqrt(1-f^2) + i f), f = fraca
!!   psi   = exp(a) (projection, Rytov) or 1 + a (first Born)
!!   psi_T = bwd(T P^-1 fwd(psi)),       T = exp(i eval_phase) exp(-B |k|^2/4), Kirkland's exp(-i chi),
!!                                       T(0) = 1 (the phase shift acts on the scattered waves)
!!
!! s refers to the specimen's midplane (the paraboloid is sampled about the box centre), and so does the
!! nominal defocus. Rytov exponentiates the scattered field where the wave leaves the specimen, at the exit
!! face, exit_depth (nat/2 for a curved model, the specimen filling the box) downstream:
!! P = exp(-i pi lambda d |k|^2) carries the field there and P^-1 carries the wave back to the midplane for
!! the lens (Kirkland's propagator; specter's Rytov does the same). For first Born the two cancel exactly,
!! and without curvature the model has no depth, so P = 1 there and Rytov is projection.
!!   ihat  = crop(fwd(|psi_T|^2))        the native DFT of the relative intensity (vacuum 1)
!!
!! forward stores psi and psi_T; linearize maps ds to d ihat at the stored point, adjoint is its transpose
!! under the real inner product Re<a,b> = sum Re(conj(a) b). The real-space origin of the fine grid is its
!! first point (no checkerboard): every step is shift-equivariant, so ihat keeps the origin convention of s.
type :: wave_imaging
    private
    integer              :: nat   = 0
    integer              :: R     = 0
    integer              :: n     = 0
    integer              :: osf   = 0
    integer              :: model = 0
    complex              :: c     = cmplx(0.,0.)  !< the scattering coefficient of the current particle
    real                 :: wl    = 0.             !< the current particle's wavelength in native pixels
    type(wave_field)     :: fpsi                  !< psi, the exit wave on the fine grid
    type(wave_field)     :: fpt                   !< psi_T, the wave after the lens
    type(wave_field)     :: fw                    !< work field
    complex, allocatable :: tlens(:,:)            !< T P^-1 on the fine grid, FFT order
    complex, allocatable :: pexit(:,:)            !< P on the native disk (Rytov), 1 otherwise
    real,    allocatable :: fsq(:,:)              !< |k|^2 of the fine grid in native cycles per pixel
    real,    allocatable :: fang(:,:)             !< azimuth atan2(k,h) of the fine grid
    logical              :: l_optics = .false.
    logical              :: l_point  = .false.   !< forward has stored a point
  contains
    procedure :: new
    procedure :: set_optics
    procedure :: forward
    procedure :: linearize
    procedure :: adjoint
    procedure :: get_exit_wave
    procedure :: forward_wave
    procedure :: linearize_wave
    procedure :: adjoint_wave
    procedure :: get_n
    procedure :: get_c
    procedure :: get_wl
    procedure, private :: lens_forward
    procedure, private :: lens_linearize
    procedure, private :: lens_adjoint
    procedure :: kill
end type wave_imaging

contains

    !> nat the native box (even), osf the oversampling of the fine grid, model one of WAVE_*
    subroutine new( self, nat, osf, model )
        class(wave_imaging), intent(inout) :: self
        integer,             intent(in)    :: nat, osf, model
        integer :: i, j, hf, kf
        call self%kill
        if( nat < 2 .or. mod(nat,2) /= 0 ) THROW_HARD('the native box must be even; wave_imaging%new')
        if( osf < 1 ) THROW_HARD('the oversampling must be at least 1; wave_imaging%new')
        select case(model)
            case(WAVE_PROJECTION, WAVE_FIRST_BORN, WAVE_RYTOV, WAVE_MSLICE)
            case DEFAULT
                THROW_HARD('unknown wave model; wave_imaging%new')
        end select
        self%nat   = nat
        self%R     = nat / 2
        self%osf   = osf
        self%n     = osf * nat
        self%model = model
        call self%fpsi%new(self%n)
        call self%fpt%new(self%n)
        call self%fw%new(self%n)
        allocate(self%tlens(self%n,self%n), source=cmplx(1.,0.))
        allocate(self%pexit(-self%R:self%R,-self%R:self%R), source=cmplx(1.,0.))
        ! the particle-independent geometry of the lens, once: the fine grid's frequencies h/nat
        allocate(self%fsq(self%n,self%n), self%fang(self%n,self%n))
        do j = 1, self%n
            kf = wave_freq(j, self%n)
            do i = 1, self%n
                hf = wave_freq(i, self%n)
                self%fsq(i,j)  = real(hf*hf + kf*kf) / real(self%nat)**2
                self%fang(i,j) = atan2(real(kf), real(hf))
            end do
        end do
    end subroutine new

    !> The particle's lens and scattering coefficient: ctfparms at the native sampling (smpd), bfac the
    !! B-factor (A^2) of the lens envelope, exit_depth (native pixels) the distance from the midplane to
    !! where Rytov exponentiates (0: at the midplane). ctfflag no is an aberration-free lens; a
    !! phase-flipped particle has no wave model
    subroutine set_optics( self, ctfparms, bfac, exit_depth )
        class(wave_imaging), intent(inout) :: self
        type(ctfparams),     intent(in)    :: ctfparms
        real,                intent(in)    :: bfac, exit_depth
        type(ctf)     :: tfun
        type(ctfvars) :: vars
        real    :: sig, beta, ph, env, dprop
        integer :: i, j, h, k
        if( self%n == 0 ) THROW_HARD('new has not been called; wave_imaging%set_optics')
        if( ctfparms%ctfflag == CTFFLAG_FLIP ) THROW_HARD('phase-flipped particles have no wave model; set_optics')
        if( ctfparms%fraca < 0. .or. ctfparms%fraca >= 1. ) THROW_HARD('fraca must be in [0,1); set_optics')
        sig    = electron_interaction_constant(ctfparms%kv)
        beta   = sqrt(1. - ctfparms%fraca**2)
        self%c = cmplx(0., sig * ctfparms%smpd) * cmplx(beta, ctfparms%fraca)
        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        vars = tfun%get_ctfvars(ctfparms%phshift)
        self%wl = vars%wl
        ! Rytov's exit face exit_depth pixels downstream of the midplane, P = exp(-i pi wl d s^2), wl and s
        ! in pixels; P = 1 for the other models
        dprop = 0.
        if( self%model == WAVE_RYTOV ) dprop = exit_depth
        do k = -self%R, self%R
            do h = -self%R, self%R
                ph = -PI * vars%wl * dprop * real(h*h + k*k) / real(self%nat)**2
                self%pexit(h,k) = cmplx(cos(ph), sin(ph))
            end do
        end do
        if( ctfparms%ctfflag == CTFFLAG_NO )then
            self%tlens = cmplx(1.,0.)
        else
            do j = 1, self%n
                do i = 1, self%n
                    ph  = tfun%eval_phase(self%fsq(i,j), self%fang(i,j), vars%phshift)
                    env = exp(-0.25 * bfac * self%fsq(i,j) / ctfparms%smpd**2)
                    self%tlens(i,j) = env * cmplx(cos(ph), sin(ph))
                end do
            end do
            ! the unscattered beam keeps its phase: a phase plate (phshift) shifts the scattered waves
            ! relative to it, and a common phase of all would vanish from the intensity
            self%tlens(1,1) = cmplx(1.,0.)
        endif
        ! P^-1 back to the midplane, on the fine grid
        if( dprop /= 0. )then
            do j = 1, self%n
                do i = 1, self%n
                    ph = PI * vars%wl * dprop * self%fsq(i,j)
                    self%tlens(i,j) = self%tlens(i,j) * cmplx(cos(ph), sin(ph))
                end do
            end do
        endif
        self%l_optics = .true.
        self%l_point  = .false.
    end subroutine set_optics

    !> ihat, the native DFT of the relative intensity, of the spectrum s; stores psi and psi_T
    subroutine forward( self, s, ihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: s(-self%R:self%R,-self%R:self%R)
        complex,             intent(out)   :: ihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:)
        if( .not. self%l_optics ) THROW_HARD('set_optics has not been called; wave_imaging%forward')
        if( self%model == WAVE_MSLICE ) THROW_HARD('multislice gives the exit wave: forward_wave; wave_imaging%forward')
        psi => self%fpsi%ptr()
        call self%fpsi%embed(s * self%pexit, self%R)
        call self%fpsi%bwd
        psi = self%c * psi
        select case(self%model)
            case(WAVE_PROJECTION, WAVE_RYTOV)
                psi = exp(psi)
            case(WAVE_FIRST_BORN)
                psi = cmplx(1.,0.) + psi
        end select
        call self%lens_forward(ihat)
    end subroutine forward

    !> ihat of the exit wave psi_in (fine grid, real space, origin at the first point, referred to the
    !! midplane); stores it as the point
    subroutine forward_wave( self, psi_in, ihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: psi_in(self%n,self%n)
        complex,             intent(out)   :: ihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:)
        if( .not. self%l_optics ) THROW_HARD('set_optics has not been called; wave_imaging%forward_wave')
        if( self%model /= WAVE_MSLICE ) THROW_HARD('the exit wave is given only to the multislice model; forward_wave')
        psi => self%fpsi%ptr()
        psi =  psi_in
        call self%lens_forward(ihat)
    end subroutine forward_wave

    !> the lens and the detector from the exit wave in fpsi: psi_T = bwd(T fwd(psi)), ihat = crop(fwd |psi_T|^2)
    subroutine lens_forward( self, ihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(out)   :: ihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:), pt(:,:), w(:,:)
        psi => self%fpsi%ptr()
        pt  => self%fpt%ptr()
        w   => self%fw%ptr()
        pt = psi
        call self%fpt%fwd
        pt = pt * self%tlens
        call self%fpt%bwd
        w = cmplx(real(pt)**2 + aimag(pt)**2, 0.)
        call self%fw%fwd
        call self%fw%crop(ihat, self%R)
        self%l_point = .true.
    end subroutine lens_forward

    !> d ihat = J ds at the point forward stored
    subroutine linearize( self, ds, dihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: ds(-self%R:self%R,-self%R:self%R)
        complex,             intent(out)   :: dihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:), w(:,:)
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_imaging%linearize')
        if( self%model == WAVE_MSLICE ) THROW_HARD('multislice gives the exit wave: linearize_wave')
        psi => self%fpsi%ptr()
        w   => self%fw%ptr()
        call self%fw%embed(ds * self%pexit, self%R)
        call self%fw%bwd
        w = self%c * w                              ! da
        if( self%model /= WAVE_FIRST_BORN ) w = psi * w   ! d psi = psi da
        call self%lens_linearize(dihat)
    end subroutine linearize

    !> d ihat of a change dpsi of the exit wave at the point stored
    subroutine linearize_wave( self, dpsi, dihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: dpsi(self%n,self%n)
        complex,             intent(out)   :: dihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: w(:,:)
        if( .not. self%l_point ) THROW_HARD('forward_wave has not been called; wave_imaging%linearize_wave')
        if( self%model /= WAVE_MSLICE ) THROW_HARD('the exit wave is given only to the multislice model; linearize_wave')
        w => self%fw%ptr()
        w =  dpsi
        call self%lens_linearize(dihat)
    end subroutine linearize_wave

    !> the lens and detector linearised, from d psi in fw
    subroutine lens_linearize( self, dihat )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(out)   :: dihat(-self%R:self%R,-self%R:self%R)
        complex, pointer :: pt(:,:), w(:,:)
        pt => self%fpt%ptr()
        w  => self%fw%ptr()
        call self%fw%fwd
        w = w * self%tlens
        call self%fw%bwd                            ! d psi_T
        w = cmplx(2. * real(conjg(pt) * w), 0.)     ! dI = 2 Re(conj(psi_T) d psi_T)
        call self%fw%fwd
        call self%fw%crop(dihat, self%R)
    end subroutine lens_linearize

    !> ds_adj = J^T q at the point forward stored, under Re<a,b>
    subroutine adjoint( self, q, ds_adj )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: q(-self%R:self%R,-self%R:self%R)
        complex,             intent(out)   :: ds_adj(-self%R:self%R,-self%R:self%R)
        complex, pointer :: psi(:,:), w(:,:)
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_imaging%adjoint')
        if( self%model == WAVE_MSLICE ) THROW_HARD('multislice gives the exit wave: adjoint_wave')
        psi => self%fpsi%ptr()
        w   => self%fw%ptr()
        call self%lens_adjoint(q)
        if( self%model /= WAVE_FIRST_BORN ) w = conjg(psi) * w
        ! (c bwd embed)^T = embed_t bwd^H conj(c) = embed_t fwd conj(c) / n^2
        w = (conjg(self%c) / real(self%n)**2) * w
        call self%fw%fwd
        call self%fw%embed_t(ds_adj, self%R)
        ds_adj = ds_adj * conjg(self%pexit)
    end subroutine adjoint

    !> dpsi_adj, the transpose of linearize_wave applied to q (under Re<a,b>), at the point stored
    subroutine adjoint_wave( self, q, dpsi_adj )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: q(-self%R:self%R,-self%R:self%R)
        complex,             intent(out)   :: dpsi_adj(self%n,self%n)
        complex, pointer :: w(:,:)
        if( .not. self%l_point ) THROW_HARD('forward_wave has not been called; wave_imaging%adjoint_wave')
        if( self%model /= WAVE_MSLICE ) THROW_HARD('the exit wave is given only to the multislice model; adjoint_wave')
        w => self%fw%ptr()
        call self%lens_adjoint(q)
        dpsi_adj = w
    end subroutine adjoint_wave

    !> the transpose of the lens and detector: fw receives the exit wave's cotangent of q
    subroutine lens_adjoint( self, q )
        class(wave_imaging), intent(inout) :: self
        complex,             intent(in)    :: q(-self%R:self%R,-self%R:self%R)
        complex, pointer :: pt(:,:), w(:,:)
        pt => self%fpt%ptr()
        w  => self%fw%ptr()
        ! (crop o fwd)^T = fwd^H crop_t = n^2 bwd crop_t, onto the real field dI
        call self%fw%crop_t(q, self%R)
        call self%fw%bwd
        ! d psi_T* = 2 Re(.) psi_T
        w = (2. * real(self%n)**2 * real(w)) * pt
        ! (bwd T fwd)^T = bwd conj(T) fwd
        call self%fw%fwd
        w = w * conjg(self%tlens)
        call self%fw%bwd
    end subroutine lens_adjoint

    !> the exit wave psi on the fine grid (real space, origin at the first point), of the point forward stored
    subroutine get_exit_wave( self, psi_out )
        class(wave_imaging), intent(in)  :: self
        complex,             intent(out) :: psi_out(self%n,self%n)
        complex, pointer :: psi(:,:)
        if( .not. self%l_point ) THROW_HARD('forward has not been called; wave_imaging%get_exit_wave')
        psi => self%fpsi%ptr()
        psi_out = psi
    end subroutine get_exit_wave

    pure integer function get_n( self )
        class(wave_imaging), intent(in) :: self
        get_n = self%n
    end function get_n

    !> the scattering coefficient c = i sigma a (sqrt(1-fraca^2) + i fraca) of the particle set_optics set
    pure complex function get_c( self )
        class(wave_imaging), intent(in) :: self
        get_c = self%c
    end function get_c

    !> the particle's wavelength in native pixels
    pure real function get_wl( self )
        class(wave_imaging), intent(in) :: self
        get_wl = self%wl
    end function get_wl

    !> idempotent
    subroutine kill( self )
        class(wave_imaging), intent(inout) :: self
        call self%fpsi%kill
        call self%fpt%kill
        call self%fw%kill
        if( allocated(self%tlens) ) deallocate(self%tlens)
        if( allocated(self%pexit) ) deallocate(self%pexit)
        if( allocated(self%fsq)   ) deallocate(self%fsq)
        if( allocated(self%fang)  ) deallocate(self%fang)
        self%nat = 0; self%R = 0; self%n = 0; self%osf = 0; self%model = 0
        self%c = cmplx(0.,0.)
        self%wl = 0.
        self%l_optics = .false.
        self%l_point  = .false.
    end subroutine kill

    !> The detector of simulate_scattering on a relative-intensity image (vacuum 1): with total_dose > 0 the
    !! intensity becomes expected electrons per pixel (total_dose smpd^2), Poisson-sampled and normalised back
    !! by the same dose; noise_sdev adds Gaussian detector noise after that. The result stays nonnegative
    subroutine detect_intensity( img, total_dose, noise_sdev )
        class(image), intent(inout) :: img
        real,         intent(in)    :: total_dose, noise_sdev
        real, allocatable :: intensity(:,:,:)
        real    :: dose_per_pixel, expected
        integer :: i, j, ldim(3)
        if( total_dose < 0. .or. noise_sdev < 0. ) THROW_HARD('total_dose and noise_sdev must be nonnegative; detect_intensity')
        ldim = img%get_ldim()
        intensity = img%get_rmat()
        if( total_dose > 0.0 )then
            dose_per_pixel = total_dose * img%get_smpd()**2
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    expected = max(0.0, intensity(i,j,1) * dose_per_pixel)
                    intensity(i,j,1) = real(poisson_rnd(expected)) / dose_per_pixel
                enddo
            enddo
        endif
        if( noise_sdev > 0.0 )then
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    intensity(i,j,1) = max(0.0, intensity(i,j,1) + noise_sdev * gasdev())
                enddo
            enddo
        endif
        call img%set_rmat(intensity, .false.)
    end subroutine detect_intensity

end module simple_wave_imaging
