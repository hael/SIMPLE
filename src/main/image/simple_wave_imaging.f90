!@descr: the 2D image formation of one particle under a single-scatter wave model, its linearisation and adjoint;
! the lens and detection of a given exit wave (multislice); the detector step of simulate_scattering and
! simulate_particles (the chain: doc/policies/3D/wave_model_policy.md, section 3)
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

!> One thread's image formation from a native-disk spectrum s (volts times voxels, unnormalised native DFT):
!! exit wave, lens, intensity, crop (wave_model_policy.md, section 3); forward stores psi and psi_T, linearize and
!! adjoint act at that point under Re<a,b>. The fine grid's real-space origin is its first point
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
    real,    allocatable :: c2ang(:,:), s2ang(:,:) !< cos and sin of twice the azimuth atan2(k,h) of the fine grid
    real,    allocatable :: env(:,:)              !< the B-factor envelope of env_bfac at env_smpd
    complex, allocatable :: pinv(:,:)             !< P^-1 on the fine grid, of pexit_wl and pexit_d
    real                 :: env_bfac = 0., env_smpd = 0.
    real                 :: pexit_wl = 0., pexit_d = 0.   !< the wavelength and depth of pexit
    logical              :: l_env = .false., l_pexit = .false.
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
        call self%fw%new(self%n, with_real=.true.)
        allocate(self%tlens(self%n,self%n), source=cmplx(1.,0.))
        allocate(self%pexit(-self%R:self%R,-self%R:self%R), source=cmplx(1.,0.))
        ! the particle-independent geometry of the lens, once: the fine grid's frequencies h/nat
        allocate(self%fsq(self%n,self%n), self%c2ang(self%n,self%n), self%s2ang(self%n,self%n), self%env(self%n,self%n))
        allocate(self%pinv(self%n,self%n))
        do j = 1, self%n
            kf = wave_freq(j, self%n)
            do i = 1, self%n
                hf = wave_freq(i, self%n)
                self%fsq(i,j)   = real(hf*hf + kf*kf) / real(self%nat)**2
                self%c2ang(i,j) = cos(2. * atan2(real(kf), real(hf)))
                self%s2ang(i,j) = sin(2. * atan2(real(kf), real(hf)))
            end do
        end do
    end subroutine new

    !> The particle's lens and scattering coefficient: ctfparms at the native sampling (smpd), bfac the
    !! B-factor (A^2) of the lens envelope, exit_depth (native pixels) the distance from the midplane to
    !! where Rytov exponentiates (0: at the midplane). ctfflag no is an aberration-free lens
    subroutine set_optics( self, ctfparms, bfac, exit_depth )
        class(wave_imaging), intent(inout) :: self
        type(ctfparams),     intent(in)    :: ctfparms
        real,                intent(in)    :: bfac, exit_depth
        type(ctf)     :: tfun
        type(ctfvars) :: vars
        real    :: sig, beta, ph, dprop, dsum, ddiff, ca, sa, pwl, hwl2cs
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
        ! in pixels; P = 1 for the other models. It and the envelope change with the wavelength and B only
        dprop = 0.
        if( self%model == WAVE_RYTOV ) dprop = exit_depth
        if( .not. self%l_pexit .or. vars%wl /= self%pexit_wl .or. dprop /= self%pexit_d )then
            do k = -self%R, self%R
                do h = -self%R, self%R
                    ph = -PI * vars%wl * dprop * real(h*h + k*k) / real(self%nat)**2
                    self%pexit(h,k) = cmplx(cos(ph), sin(ph))
                end do
            end do
            ! and P^-1 back to the midplane, on the fine grid
            do j = 1, self%n
                do i = 1, self%n
                    ph = PI * vars%wl * dprop * self%fsq(i,j)
                    self%pinv(i,j) = cmplx(cos(ph), sin(ph))
                end do
            end do
            self%pexit_wl = vars%wl
            self%pexit_d  = dprop
            self%l_pexit  = .true.
        endif
        if( .not. self%l_env .or. bfac /= self%env_bfac .or. ctfparms%smpd /= self%env_smpd )then
            self%env      = exp(-0.25 * bfac * self%fsq / ctfparms%smpd**2)
            self%env_bfac = bfac
            self%env_smpd = ctfparms%smpd
            self%l_env    = .true.
        endif
        ! T on the fine grid: the aberration phase of ctf%eval_phase with the astigmatism's angle terms expanded;
        ! P^-1 a factor of its own, so that P and P^-1 cancel to rounding (first Born is Rytov linearised)
        if( ctfparms%ctfflag == CTFFLAG_NO )then
            self%tlens = cmplx(1.,0.)
        else
            dsum   = vars%dfx + vars%dfy
            ddiff  = vars%dfx - vars%dfy
            ca     = cos(2. * vars%angast)
            sa     = sin(2. * vars%angast)
            pwl    = PI * vars%wl
            hwl2cs = 0.5 * vars%wl * vars%wl * vars%cs
            do j = 1, self%n
                do i = 1, self%n
                    ph = pwl * self%fsq(i,j) * (0.5 * (dsum + (self%c2ang(i,j) * ca + self%s2ang(i,j) * sa) * ddiff) &
                        &- hwl2cs * self%fsq(i,j)) + vars%phshift
                    self%tlens(i,j) = self%env(i,j) * cmplx(cos(ph), sin(ph))
                end do
            end do
            ! the unscattered beam keeps its phase: a phase plate (phshift) shifts the scattered waves
            ! relative to it, and a common phase of all would vanish from the intensity
            self%tlens(1,1) = cmplx(1.,0.)
        endif
        if( dprop /= 0. ) self%tlens = self%tlens * self%pinv
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
        complex, pointer :: psi(:,:), pt(:,:)
        real,    pointer :: r(:,:)
        psi => self%fpsi%ptr()
        pt  => self%fpt%ptr()
        r   => self%fw%rptr()
        pt = psi
        call self%fpt%fwd
        pt = pt * self%tlens
        call self%fpt%bwd
        r = real(pt)**2 + aimag(pt)**2
        call self%fw%fwd_real_crop(ihat, self%R)
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
        real,    pointer :: r(:,:)
        pt => self%fpt%ptr()
        w  => self%fw%ptr()
        r  => self%fw%rptr()
        call self%fw%fwd
        w = w * self%tlens
        call self%fw%bwd                            ! d psi_T
        r = 2. * real(conjg(pt) * w)                ! dI = 2 Re(conj(psi_T) d psi_T)
        call self%fw%fwd_real_crop(dihat, self%R)
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
        real,    pointer :: r(:,:)
        pt => self%fpt%ptr()
        w  => self%fw%ptr()
        r  => self%fw%rptr()
        ! (crop o fwd)^T = fwd^H crop_t = n^2 bwd crop_t, onto the real field dI: r = n^2 Re(bwd crop_t q)
        call self%fw%crop_t_real(q, self%R)
        ! d psi_T* = 2 Re(.) psi_T
        w = (2. * r) * pt
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
        if( allocated(self%c2ang) ) deallocate(self%c2ang)
        if( allocated(self%s2ang) ) deallocate(self%s2ang)
        if( allocated(self%env)   ) deallocate(self%env)
        if( allocated(self%pinv)  ) deallocate(self%pinv)
        self%env_bfac = 0.; self%env_smpd = 0.; self%pexit_wl = 0.; self%pexit_d = 0.
        self%l_env = .false.; self%l_pexit = .false.
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
