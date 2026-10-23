!@descr: analytic tests of the Ghostbuster-compatible scattering forward models
module simple_scattering_simulator_tester
use simple_image,                only: image
use simple_scattering_simulator, only: scattering_simulator, electron_interaction_constant
use simple_test_utils,           only: assert_real, assert_true
implicit none
private

public :: run_all_scattering_simulator_tests

contains

    subroutine run_all_scattering_simulator_tests
        integer, parameter :: BOX = 8
        real,    parameter :: SMPD = 1.5, KV = 200.0, POTENTIAL_VALUE = 0.25
        class(scattering_simulator), allocatable :: simulator
        type(image) :: potential, detector
        real, allocatable :: volume(:,:,:), intensity(:,:,:)
        real :: expected, phase, sigma_e
        sigma_e = electron_interaction_constant(KV)
        call assert_real(7.29e-4, sigma_e, 2.0e-6, '200 kV interaction constant matches the reference value')
        call potential%new([BOX,BOX,BOX], SMPD)
        call detector%new([BOX,BOX,1], SMPD)
        allocate(volume(BOX,BOX,BOX))
        allocate(simulator)
        call simulator%new(BOX, SMPD, KV, 0.0)

        volume = 0.0
        call potential%set_rmat(volume, .false.)
        call simulator%simulate(potential, detector, 'first_born', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(maxval(abs(intensity - 1.0)) < 2.0e-5, 'zero potential gives unit first-Born intensity')
        deallocate(intensity)
        call simulator%simulate(potential, detector, 'multislice', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(maxval(abs(intensity - 1.0)) < 2.0e-5, 'zero potential gives unit multislice intensity')
        deallocate(intensity)

        volume = POTENTIAL_VALUE
        call potential%set_rmat(volume, .false.)
        phase = sigma_e * SMPD * real(BOX) * POTENTIAL_VALUE
        expected = 1.0 + phase**2
        call simulator%simulate(potential, detector, 'first_born', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(maxval(abs(intensity - expected)) < 3.0e-5, &
            &'constant potential gives the analytic first-Born interference intensity')
        deallocate(intensity)
        call simulator%simulate(potential, detector, 'multislice', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(maxval(abs(intensity - 1.0)) < 3.0e-5, &
            &'constant real potential is a pure phase object under multislice')
        deallocate(intensity)
        call simulator%kill

        call simulator%new(BOX, SMPD, KV, 0.1)
        expected = exp(-2.0 * sigma_e * SMPD * real(BOX) * 0.1 * POTENTIAL_VALUE)
        call simulator%simulate(potential, detector, 'multislice', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(maxval(abs(intensity - expected)) < 3.0e-5, &
            &'absorptive multislice potential gives analytic exponential attenuation')
        deallocate(intensity)
        call simulator%simulate(potential, detector, 'multislice', 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0)
        intensity = detector%get_rmat()
        call assert_true(all(intensity >= 0.0), 'Poisson-sampled detector intensity is nonnegative')
        call assert_true(maxval(abs(intensity * (20.0 * SMPD**2) - &
            &real(nint(intensity * (20.0 * SMPD**2))))) < 2.0e-5, &
            &'dose-normalized detector values correspond to integral electron counts')
        deallocate(intensity)

        call simulator%kill
        deallocate(simulator, volume)
        call detector%kill
        call potential%kill
    end subroutine run_all_scattering_simulator_tests

end module simple_scattering_simulator_tester
