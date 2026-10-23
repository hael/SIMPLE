!@descr: Ghostbuster-compatible first-Born and multislice electron-scattering forward models
module simple_scattering_simulator
use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_float, c_float_complex, c_ptr, c_size_t, c_null_ptr
use simple_core_module_api, only: PI, gasdev, ran3
use simple_fftw3, only: fftwf_alloc_complex, fftwf_free, fftwf_plan_dft_2d, fftwf_execute_dft, &
    &fftwf_destroy_plan, FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
use simple_image, only: image
implicit none
private

public :: scattering_simulator, electron_interaction_constant

complex(c_float_complex), parameter :: IMAG_UNIT = cmplx(0.0_c_float, 1.0_c_float, kind=c_float)

type :: scattering_simulator
    private
    integer                           :: box = 0
    real                              :: smpd = 0.0
    real                              :: wavelength = 0.0
    real                              :: sigma_e = 0.0
    real                              :: absorption = 0.1
    type(c_ptr)                       :: spatial_ptr = c_null_ptr
    type(c_ptr)                       :: fourier_ptr = c_null_ptr
    type(c_ptr)                       :: plan_fwd = c_null_ptr
    type(c_ptr)                       :: plan_bwd = c_null_ptr
    complex(c_float_complex), pointer :: spatial(:,:) => null()
    complex(c_float_complex), pointer :: fourier(:,:) => null()
contains
    procedure :: new      => scattering_new
    procedure :: simulate => scattering_simulate
    procedure :: kill     => scattering_kill
end type scattering_simulator

contains

    subroutine scattering_new( self, box, smpd, kv, absorption )
        class(scattering_simulator), intent(inout) :: self
        integer,                     intent(in)    :: box
        real,                        intent(in)    :: smpd, kv, absorption
        integer(c_size_t) :: ncomplex
        if( box < 2 ) error stop 'scattering_simulator: box must be at least two pixels'
        if( smpd <= 0.0 ) error stop 'scattering_simulator: smpd must be positive'
        if( kv <= 0.0 ) error stop 'scattering_simulator: kv must be positive'
        if( absorption < 0.0 .or. absorption >= 1.0 ) &
            &error stop 'scattering_simulator: absorption must be in [0,1)'
        call self%kill
        self%box        = box
        self%smpd       = smpd
        self%wavelength = electron_wavelength(kv)
        self%sigma_e    = electron_interaction_constant(kv)
        self%absorption = absorption
        ncomplex = int(box, c_size_t) * int(box, c_size_t)
        self%spatial_ptr = fftwf_alloc_complex(ncomplex)
        self%fourier_ptr = fftwf_alloc_complex(ncomplex)
        if( .not. c_associated(self%spatial_ptr) .or. .not. c_associated(self%fourier_ptr) ) &
            &error stop 'scattering_simulator: FFTW allocation failed'
        call c_f_pointer(self%spatial_ptr, self%spatial, [box, box])
        call c_f_pointer(self%fourier_ptr, self%fourier, [box, box])
        self%spatial = cmplx(0.0_c_float, 0.0_c_float, kind=c_float)
        self%fourier = cmplx(0.0_c_float, 0.0_c_float, kind=c_float)
        !$omp critical (fftw_planner)
        self%plan_fwd = fftwf_plan_dft_2d(box, box, self%spatial, self%fourier, FFTW_FORWARD, FFTW_ESTIMATE)
        self%plan_bwd = fftwf_plan_dft_2d(box, box, self%fourier, self%spatial, FFTW_BACKWARD, FFTW_ESTIMATE)
        !$omp end critical (fftw_planner)
        if( .not. c_associated(self%plan_fwd) .or. .not. c_associated(self%plan_bwd) ) &
            &error stop 'scattering_simulator: FFTW plan creation failed'
    end subroutine scattering_new

    subroutine scattering_simulate( self, potential, detector, model, dfx, dfy, angast, cs, bfac, total_dose, noise_sdev )
        class(scattering_simulator), intent(inout) :: self
        class(image),                intent(in)    :: potential
        class(image),                intent(inout) :: detector
        character(len=*),            intent(in)    :: model
        real,                        intent(in)    :: dfx, dfy, angast, cs, bfac, total_dose, noise_sdev
        real(c_float), pointer                :: v(:,:,:)
        complex(c_float_complex), allocatable :: wave(:,:)
        real, allocatable                     :: intensity(:,:,:)
        integer                               :: i, j, iz, ldim(3)
        real                                  :: coeff_real, distance, dose_per_pixel, expected
        complex(c_float_complex)              :: potential_coeff
        if( self%box == 0 ) error stop 'scattering_simulator: new must be called before simulate'
        ldim = potential%get_ldim()
        if( ldim(1) /= self%box .or. ldim(2) /= self%box ) &
            &error stop 'scattering_simulator: potential x/y dimensions differ from simulator box'
        if( abs(potential%get_smpd() - self%smpd) > 10.0 * epsilon(self%smpd) ) &
            &error stop 'scattering_simulator: potential sampling differs from simulator sampling'
        if( cs < 0.0 .or. bfac < 0.0 .or. total_dose < 0.0 .or. noise_sdev < 0.0 ) &
            &error stop 'scattering_simulator: cs, bfac, total_dose, and noise_sdev must be nonnegative'
        call potential%get_rmat_ptr(v)
        allocate(wave(self%box,self%box))
        coeff_real      = self%sigma_e * self%smpd
        potential_coeff = cmplx(sqrt(1.0 - self%absorption**2), self%absorption, kind=c_float)
        select case(trim(model))
        case('first_born')
            wave = cmplx(1.0_c_float, 0.0_c_float, kind=c_float)
            do iz = 1, ldim(3)
                self%spatial = cmplx(coeff_real * v(:,:,iz), 0.0, kind=c_float) * potential_coeff
                call fftwf_execute_dft(self%plan_fwd, self%spatial, self%fourier)
                distance = real(ldim(3) - iz) * self%smpd
                call apply_fresnel(self, distance)
                call inverse_fft(self)
                wave = wave + IMAG_UNIT * self%spatial
            enddo
        case('multislice')
            wave = cmplx(1.0_c_float, 0.0_c_float, kind=c_float)
            do iz = 1, ldim(3)
                wave = wave * exp(IMAG_UNIT * coeff_real * potential_coeff * v(:,:,iz))
                if( iz < ldim(3) )then
                    self%spatial = wave
                    call fftwf_execute_dft(self%plan_fwd, self%spatial, self%fourier)
                    call apply_fresnel(self, self%smpd)
                    call inverse_fft(self)
                    wave = self%spatial
                endif
            enddo
        case default
            error stop 'scattering_simulator: model must be first_born or multislice'
        end select
        self%spatial = wave
        call fftwf_execute_dft(self%plan_fwd, self%spatial, self%fourier)
        call apply_objective_lens(self, dfx, dfy, angast, cs, bfac)
        call inverse_fft(self)
        allocate(intensity(self%box,self%box,1))
        intensity(:,:,1) = real(abs(self%spatial)**2)
        if( total_dose > 0.0 )then
            dose_per_pixel = total_dose * self%smpd**2
            do j = 1, self%box
                do i = 1, self%box
                    expected = max(0.0, intensity(i,j,1) * dose_per_pixel)
                    intensity(i,j,1) = real(poisson_sample(expected)) / dose_per_pixel
                enddo
            enddo
        endif
        if( noise_sdev > 0.0 )then
            do j = 1, self%box
                do i = 1, self%box
                    intensity(i,j,1) = max(0.0, intensity(i,j,1) + noise_sdev * gasdev())
                enddo
            enddo
        endif
        call detector%set_rmat(intensity, .false.)
        deallocate(wave, intensity)
    end subroutine scattering_simulate

    subroutine apply_fresnel( self, distance )
        class(scattering_simulator), intent(inout) :: self
        real,                        intent(in)    :: distance
        integer :: i, j
        real    :: kx, ky, ksq, phase
        do j = 1, self%box
            ky = fourier_frequency(j, self%box, self%smpd)
            do i = 1, self%box
                kx = fourier_frequency(i, self%box, self%smpd)
                ksq = kx*kx + ky*ky
                phase = PI * distance * self%wavelength * ksq
                self%fourier(i,j) = self%fourier(i,j) * exp(IMAG_UNIT * phase)
            enddo
        enddo
    end subroutine apply_fresnel

    subroutine apply_objective_lens( self, dfx, dfy, angast, cs, bfac )
        class(scattering_simulator), intent(inout) :: self
        real,                        intent(in)    :: dfx, dfy, angast, cs, bfac
        integer :: i, j
        real    :: kx, ky, ksq, theta, defocus_a, cs_a, chi, envelope, angast_rad
        cs_a = cs * 1.0e7
        angast_rad = angast * PI / 180.0
        do j = 1, self%box
            ky = fourier_frequency(j, self%box, self%smpd)
            do i = 1, self%box
                kx = fourier_frequency(i, self%box, self%smpd)
                ksq = kx*kx + ky*ky
                theta = atan2(ky, kx)
                defocus_a = 1.0e4 * 0.5 * (dfx + dfy + (dfx - dfy) * cos(2.0 * (theta - angast_rad)))
                chi = PI * (-defocus_a * self%wavelength * ksq + &
                    &0.5 * cs_a * self%wavelength**3 * ksq**2)
                envelope = exp(-0.25 * bfac * ksq)
                self%fourier(i,j) = self%fourier(i,j) * envelope * exp(IMAG_UNIT * chi)
            enddo
        enddo
    end subroutine apply_objective_lens

    subroutine inverse_fft( self )
        class(scattering_simulator), intent(inout) :: self
        call fftwf_execute_dft(self%plan_bwd, self%fourier, self%spatial)
        self%spatial = self%spatial / real(self%box * self%box, c_float)
    end subroutine inverse_fft

    subroutine scattering_kill( self )
        class(scattering_simulator), intent(inout) :: self
        if( c_associated(self%plan_fwd) )then
            !$omp critical (fftw_planner)
            call fftwf_destroy_plan(self%plan_fwd)
            !$omp end critical (fftw_planner)
        endif
        if( c_associated(self%plan_bwd) )then
            !$omp critical (fftw_planner)
            call fftwf_destroy_plan(self%plan_bwd)
            !$omp end critical (fftw_planner)
        endif
        if( c_associated(self%spatial_ptr) ) call fftwf_free(self%spatial_ptr)
        if( c_associated(self%fourier_ptr) ) call fftwf_free(self%fourier_ptr)
        nullify(self%spatial, self%fourier)
        self%spatial_ptr = c_null_ptr
        self%fourier_ptr = c_null_ptr
        self%plan_fwd = c_null_ptr
        self%plan_bwd = c_null_ptr
        self%box = 0
    end subroutine scattering_kill

    pure real function electron_wavelength( kv )
        real, intent(in) :: kv
        electron_wavelength = 12.26 / sqrt(1000.0 * kv + 0.9784 * (1000.0 * kv)**2 / 1.0e6)
    end function electron_wavelength

    pure real function electron_interaction_constant( kv )
        real, intent(in) :: kv
        real(kind=8), parameter :: PLANCK = 6.62607015d-34
        real(kind=8), parameter :: ELECTRON_MASS = 9.1093837139d-31
        real(kind=8), parameter :: ELEMENTARY_CHARGE = 1.602176634d-19
        real(kind=8), parameter :: LIGHT_SPEED = 2.99792458d8
        real(kind=8) :: voltage, wavelength_m, relativistic_mass
        voltage = 1000.0d0 * real(kv, kind=8)
        wavelength_m = 1.0d-10 * real(electron_wavelength(kv), kind=8)
        relativistic_mass = ELECTRON_MASS * &
            &(1.0d0 + ELEMENTARY_CHARGE * voltage / (ELECTRON_MASS * LIGHT_SPEED**2))
        electron_interaction_constant = real(2.0d0 * real(PI,kind=8) * relativistic_mass * &
            &ELEMENTARY_CHARGE * wavelength_m / PLANCK**2 * 1.0d-10)
    end function electron_interaction_constant

    pure real function fourier_frequency( index, box, smpd )
        integer, intent(in) :: index, box
        real,    intent(in) :: smpd
        integer :: signed_index
        signed_index = index - 1
        if( signed_index > box / 2 ) signed_index = signed_index - box
        fourier_frequency = real(signed_index) / (real(box) * smpd)
    end function fourier_frequency

    integer function poisson_sample( mean )
        real, intent(in) :: mean
        real    :: a, b, inv_alpha, kreal, log_mean, product, sqrt_mean, u, us, v, vr
        integer :: k
        if( mean <= 0.0 )then
            poisson_sample = 0
        else if( mean < 30.0 )then
            product = 1.0
            k = 0
            do
                k = k + 1
                product = product * max(ran3(), tiny(1.0))
                if( product <= exp(-mean) ) exit
            enddo
            poisson_sample = k - 1
        else
            sqrt_mean = sqrt(mean)
            log_mean  = log(mean)
            b         = 0.931 + 2.53 * sqrt_mean
            a         = -0.059 + 0.02483 * b
            inv_alpha = 1.1239 + 1.1328 / (b - 3.4)
            vr        = 0.9277 - 3.6224 / (b - 2.0)
            do
                u  = ran3() - 0.5
                v  = ran3()
                us = 0.5 - abs(u)
                if( us <= tiny(us) ) cycle
                kreal = floor((2.0 * a / us + b) * u + mean + 0.43)
                k  = int(kreal)
                if( us >= 0.07 .and. v <= vr ) exit
                if( k < 0 .or. (us < 0.013 .and. v > us) ) cycle
                if( log(v * inv_alpha / (a / (us * us) + b)) <= &
                    &-mean + kreal * log_mean - log_gamma(kreal + 1.0) ) exit
            enddo
            poisson_sample = k
        endif
    end function poisson_sample

end module simple_scattering_simulator
