!@descr: analytic tests for three-dimensional shape descriptors
module simple_volume_shape_tester
use simple_test_utils, only: assert_real
use simple_image,      only: image
implicit none
private
public :: run_all_volume_shape_tests

integer, parameter :: BOX = 32, CENTRE = BOX / 2 + 1, HALF_WIDTH = 6
real,    parameter :: TOL = 1.e-4

contains

    subroutine run_all_volume_shape_tests()
        type(image) :: vol
        real, allocatable :: density(:,:,:)
        real :: q
        call vol%new([BOX, BOX, BOX], 1.0, wthreads=.false.)
        allocate(density(BOX, BOX, BOX), source=0.0)
        q = real(HALF_WIDTH * (HALF_WIDTH + 1)) / 3.0

        density(CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH, &
            &CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH, CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH) = 1.0
        call vol%set_rmat(density, .false.)
        call verify_shape(vol, 'cube', 0.0, 0.0, 0.0, 0.0, 3.0 * q)

        density = 0.0
        density(CENTRE, CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH, &
            &CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH) = 1.0
        call vol%set_rmat(density, .false.)
        call verify_shape(vol, 'plate', 1.0, 0.25, 0.25, 0.5, 2.0 * q)

        density = 0.0
        density(CENTRE, CENTRE, CENTRE-HALF_WIDTH:CENTRE+HALF_WIDTH) = 1.0
        call vol%set_rmat(density, .false.)
        call verify_shape(vol, 'rod', 1.0, 1.0, 1.0, 0.0, q)

        density = 0.0
        call vol%set_rmat(density, .false.)
        call verify_shape(vol, 'empty volume', 0.0, 0.0, 0.0, 0.0, 0.0)

        call vol%kill
        deallocate(density)
    end subroutine run_all_volume_shape_tests

    subroutine verify_shape(vol, label, expected_ecc, expected_aniso, expected_asph, expected_acyl, expected_rg_sq)
        type(image),      intent(in) :: vol
        character(len=*), intent(in) :: label
        real, intent(in) :: expected_ecc, expected_aniso, expected_asph, expected_acyl, expected_rg_sq
        real :: eccentricity, anisotropy, asphericity, acylindricity, rg_sq
        call vol%calc_3D_shape_descriptors(real(BOX), eccentricity, anisotropy, asphericity, acylindricity, rg_sq)
        call assert_real(expected_ecc,   eccentricity,  TOL, 'volume shape: '//label//' eccentricity')
        call assert_real(expected_aniso, anisotropy,    TOL, 'volume shape: '//label//' anisotropy')
        call assert_real(expected_asph,  asphericity,   TOL, 'volume shape: '//label//' asphericity')
        call assert_real(expected_acyl,  acylindricity, TOL, 'volume shape: '//label//' acylindricity')
        call assert_real(expected_rg_sq, rg_sq,         TOL, 'volume shape: '//label//' radius of gyration squared')
        call assert_real(anisotropy, asphericity**2 + 0.75 * acylindricity**2, TOL, &
            &'volume shape: '//label//' descriptor identity')
    end subroutine verify_shape

end module simple_volume_shape_tester
