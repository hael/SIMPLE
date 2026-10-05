!@descr: analytic tests for three-dimensional shape descriptors
module simple_volume_shape_tester
use simple_test_utils, only: assert_real, assert_int, enter_fixture, leave_fixture, tests_failed
use simple_string,     only: string
use simple_image,      only: image
use simple_image_bin,  only: image_bin
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
        call test_vol_shape_descr_components()
    end subroutine run_all_volume_shape_tests

    !> the components vol_shape_descr counts: one blob is one object; two blobs of a size are two;
    !! a blob outside the mask does not count; a speck below the size fraction does not count
    subroutine test_vol_shape_descr_components()
        integer, parameter :: NB = 48, C = NB / 2 + 1, HW = 4, SHIFT = 14
        real,    parameter :: SMPD_VOL = 2.0, LP = 6.0, FRAC = 0.1
        type(image)       :: vol
        type(image_bin)   :: bin
        type(string)      :: cwd_saved, root
        real, allocatable :: density(:,:,:)
        integer :: nccs, nfail0
        nfail0 = tests_failed
        call enter_fixture('vol_shape_descr', cwd_saved, root)
        call vol%new([NB, NB, NB], SMPD_VOL, wthreads=.false.)
        allocate(density(NB, NB, NB), source=0.0)
        ! one blob at the centre
        density(C-HW:C+HW, C-HW:C+HW, C-HW:C+HW) = 1.0
        call vol%set_rmat(density, .false.)
        call bin%vol_shape_descr(vol, LP, 23.0, nccs, min_frac=FRAC, tag='_one')
        call assert_int(1, nccs, 'vol_shape_descr: one blob is one object')
        call bin%kill_bimg
        ! a second blob of the same size, inside the mask
        density(C+SHIFT-HW:C+SHIFT+HW, C-HW:C+HW, C-HW:C+HW) = 1.0
        call vol%set_rmat(density, .false.)
        call bin%vol_shape_descr(vol, LP, 23.0, nccs, min_frac=FRAC, tag='_two')
        call assert_int(2, nccs, 'vol_shape_descr: two blobs of a size are two objects')
        call bin%kill_bimg
        ! the same with a mask that leaves the second blob out
        call bin%vol_shape_descr(vol, LP, 8.0, nccs, min_frac=FRAC, tag='_masked')
        call assert_int(1, nccs, 'vol_shape_descr: a blob outside the mask does not count')
        call bin%kill_bimg
        ! the central blob with a speck
        density = 0.0
        density(C-HW:C+HW, C-HW:C+HW, C-HW:C+HW) = 1.0
        density(C-SHIFT:C-SHIFT+1, C:C+1, C:C+1) = 1.0
        call vol%set_rmat(density, .false.)
        call bin%vol_shape_descr(vol, LP, 23.0, nccs, min_frac=FRAC, tag='_speck')
        call assert_int(1, nccs, 'vol_shape_descr: a speck below the size fraction does not count')
        call bin%kill_bimg
        call vol%kill
        deallocate(density)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_vol_shape_descr_components

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
