!@descr: self-contained MRC volume validation tests
module simple_mrc_validate_tester
use simple_test_utils, only: assert_double, assert_int, assert_real, assert_true
use simple_defs,      only: dp
use simple_string,    only: string
use simple_syslib,    only: del_file, file_exists
use simple_imghead,   only: find_img_smpd, find_ldim_nptcls
use simple_image,     only: image
implicit none
private
public :: run_all_mrc_validate_tests

integer, parameter :: LDIM(3) = [12, 10, 8]
real,    parameter :: SMPD = 1.75
character(len=*), parameter :: INPUT_FILE  = 'tmp_mrc_validate_input.mrc'
character(len=*), parameter :: OUTPUT_FILE = 'tmp_mrc_validate_output.mrc'

contains

    subroutine run_all_mrc_validate_tests()
        type(image) :: source, loaded, roundtrip
        real, allocatable :: reference(:,:,:), actual(:,:,:)
        real(dp) :: expected_sum, expected_sumsq
        integer  :: x, y, z

        call cleanup_files
        call source%new(LDIM, SMPD, wthreads=.false.)
        expected_sum   = 0._dp
        expected_sumsq = 0._dp
        do z = 1, LDIM(3)
            do y = 1, LDIM(2)
                do x = 1, LDIM(1)
                    call source%set([x, y, z], voxel_value(x, y, z))
                    expected_sum   = expected_sum   + real(voxel_value(x, y, z), dp)
                    expected_sumsq = expected_sumsq + real(voxel_value(x, y, z), dp)**2
                end do
            end do
        end do
        reference = source%get_rmat()
        call source%write(string(INPUT_FILE))
        call source%kill

        call assert_true(file_exists(string(INPUT_FILE)), 'mrc_validate creates its source volume internally')
        if( .not. file_exists(string(INPUT_FILE)) ) return
        call verify_header(INPUT_FILE, 'source')

        call loaded%new(LDIM, SMPD, wthreads=.false.)
        call loaded%read(string(INPUT_FILE))
        actual = loaded%get_rmat()
        call verify_density(actual, reference, expected_sum, expected_sumsq, 'source')
        call loaded%write(string(OUTPUT_FILE))
        call loaded%kill

        call assert_true(file_exists(string(OUTPUT_FILE)), 'mrc_validate writes the round-trip volume')
        if( .not. file_exists(string(OUTPUT_FILE)) )then
            call cleanup_files
            return
        endif
        call verify_header(OUTPUT_FILE, 'round-trip')
        call roundtrip%new(LDIM, SMPD, wthreads=.false.)
        call roundtrip%read(string(OUTPUT_FILE))
        actual = roundtrip%get_rmat()
        call verify_density(actual, reference, expected_sum, expected_sumsq, 'round-trip')
        call roundtrip%kill
        call cleanup_files
    end subroutine run_all_mrc_validate_tests

    subroutine verify_header(filename, label)
        character(len=*), intent(in) :: filename, label
        integer :: ldim_read(3), nsections
        real :: smpd_read
        call find_ldim_nptcls(string(filename), ldim_read, nsections)
        smpd_read = find_img_smpd(string(filename))
        call assert_true(all(ldim_read == LDIM), 'mrc_validate '//label//' dimensions are 12 x 10 x 8')
        call assert_int(LDIM(3), nsections, 'mrc_validate '//label//' contains eight sections')
        call assert_real(SMPD, smpd_read, 1.e-6, 'mrc_validate '//label//' sampling distance is 1.75 A')
    end subroutine verify_header

    subroutine verify_density(actual, reference, expected_sum, expected_sumsq, label)
        real,             intent(in) :: actual(:,:,:), reference(:,:,:)
        real(dp),         intent(in) :: expected_sum, expected_sumsq
        character(len=*), intent(in) :: label
        real(dp) :: actual_sum, actual_sumsq
        actual_sum   = sum(real(actual, dp))
        actual_sumsq = sum(real(actual, dp)**2)
        call assert_int(product(LDIM), size(actual), 'mrc_validate '//label//' contains 960 voxels')
        call assert_real(0., maxval(abs(actual - reference)), 0., &
            &'mrc_validate '//label//' maximum voxel error is zero')
        call assert_double(expected_sum, actual_sum, 'mrc_validate '//label//' voxel sum is exact')
        call assert_double(expected_sumsq, actual_sumsq, 'mrc_validate '//label//' squared norm is exact')
    end subroutine verify_density

    pure real function voxel_value(x, y, z) result(value)
        integer, intent(in) :: x, y, z
        value = real(x + 10 * y + 100 * z - 500)
    end function voxel_value

    subroutine cleanup_files()
        call del_file(INPUT_FILE)
        call del_file(OUTPUT_FILE)
    end subroutine cleanup_files

end module simple_mrc_validate_tester
