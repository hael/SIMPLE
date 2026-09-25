!@descr: unit tests for binary persistence of simple_motion_model
module simple_motion_model_tester
use, intrinsic :: iso_fortran_env, only: int8, int32, int64
use simple_core_module_api
use simple_cmdline,      only: cmdline
use simple_motion_model, only: motion_model
use simple_parameters,   only: parameters
use simple_test_utils
implicit none
private
#include "simple_local_flags.inc"

public :: run_all_motion_model_tests

contains

    subroutine run_all_motion_model_tests()
        write(*,'(A)') '**** running all motion model tests ****'
        call test_binary_roundtrip_with_optional_arrays()
        call test_binary_roundtrip_without_optional_arrays()
    end subroutine run_all_motion_model_tests

    subroutine test_binary_roundtrip_with_optional_arrays()
        write(*,'(A)') 'test_binary_roundtrip_with_optional_arrays'
        call test_binary_roundtrip(&
            &string('tmp_motion_model_patched_input.bin'),&
            &string('tmp_motion_model_patched_output.bin'),&
            &string('tmp_motion_model_patched_output.star'), .true.)
    end subroutine test_binary_roundtrip_with_optional_arrays

    subroutine test_binary_roundtrip_without_optional_arrays()
        write(*,'(A)') 'test_binary_roundtrip_without_optional_arrays'
        call test_binary_roundtrip(&
            &string('tmp_motion_model_minimal_input.bin'),&
            &string('tmp_motion_model_minimal_output.bin'),&
            &string('tmp_motion_model_minimal_output.star'), .false.)
    end subroutine test_binary_roundtrip_without_optional_arrays

    subroutine test_binary_roundtrip( input_bin, output_bin, output_star, with_optional_arrays )
        type(string), intent(in) :: input_bin, output_bin, output_star
        logical,      intent(in) :: with_optional_arrays
        type(motion_model)       :: model
        type(parameters), target :: params
        type(cmdline)            :: cline
        call cline%set('prg', 'motion_correct')
        call cline%set('mkdir', 'no')
        call cline%set('smpd', 1.5)
        call cline%set('kv', 300.)
        call cline%set('total_dose', 40.)
        call cline%set('fraction_dose_target', 1.5)
        call cline%set('scale_movies', 0.5)
        call params%new(cline)
        call del_file(input_bin)
        call del_file(output_bin)
        call del_file(output_star)
        call write_reference_binary(input_bin, with_optional_arrays)
        call model%read(input_bin, params)
        call model%write(output_star, output_bin, .true.)
        call assert_true(file_exists(output_bin), 'motion model roundtrip writes a binary model')
        call assert_true(binary_files_equal(input_bin, output_bin),&
            &'motion model read/write preserves the independently generated binary payload')
        call model%kill()
        call del_file(input_bin)
        call del_file(output_bin)
        call del_file(output_star)
        call cline%kill()
    end subroutine test_binary_roundtrip

    subroutine write_reference_binary( fname, with_optional_arrays )
        type(string), intent(in) :: fname
        logical,      intent(in) :: with_optional_arrays
        integer(int8) :: flag
        integer       :: funit, ios, i, noutliers
        integer       :: file_version, model_version
        integer       :: ldim_movie(2), ldim(2), ldim_patch(2)
        integer       :: nframes, total_nframes, npatch, nx_patch, ny_patch
        integer       :: fixed_frame, eer_fraction
        integer, allocatable :: patch_bounds(:,:,:,:), outlier_coords(:,:)
        real          :: smpd_movie, smpd, binning
        real          :: voltage, dose_per_frame, target_dose_per_frame
        real          :: total_dose, accumulated_dose, rmsd_fit(2)
        real, allocatable :: drift_x(:), drift_y(:), frameweights(:)
        real, allocatable :: patch_coords(:,:,:), local_x(:,:,:), local_y(:,:,:)
        real(dp)      :: coeffs_x(18), coeffs_y(18)
        character(len=:), allocatable :: gain
        file_version         = 0
        model_version        = 0
        smpd_movie           = 1.5
        ldim_movie           = [9, 7]
        smpd                 = 1.5
        ldim                 = [9, 7]
        binning              = 1.0
        nframes              = 3
        total_nframes        = 4
        voltage              = 300.0
        dose_per_frame       = 1.0
        target_dose_per_frame = 2.5
        total_dose           = 10.0
        accumulated_dose     = 3.0
        rmsd_fit             = [0.125, 0.25]
        fixed_frame          = 2
        eer_fraction         = 0
        ldim_patch           = [4, 4]
        allocate(drift_x(nframes), drift_y(nframes), frameweights(nframes))
        drift_x              = [0.0, 0.25, -0.5]
        drift_y              = [0.0, -0.75, 0.125]
        frameweights         = [0.2, 0.3, 0.5]
        coeffs_x             = [(real(i,dp) / 100.0_dp, i=1,18)]
        coeffs_y             = -coeffs_x
        if( with_optional_arrays )then
            gain       = 'tmp_motion_model_gain.mrc'
            nx_patch   = 3
            ny_patch   = 2
            npatch     = nx_patch * ny_patch
            noutliers  = 2
            allocate(patch_bounds(nx_patch,ny_patch,2,2))
            allocate(patch_coords(nx_patch,ny_patch,2))
            allocate(local_x(nframes,nx_patch,ny_patch), local_y(nframes,nx_patch,ny_patch))
            allocate(outlier_coords(2,noutliers))
            patch_bounds  = reshape([(i, i=1,size(patch_bounds))], shape(patch_bounds))
            patch_coords  = reshape([(real(i) / 10.0, i=1,size(patch_coords))], shape(patch_coords))
            local_x       = reshape([(real(i) / 20.0, i=1,size(local_x))], shape(local_x))
            local_y       = -local_x
            outlier_coords = reshape([1, 2, 7, 8], shape(outlier_coords))
        else
            gain       = ''
            nx_patch   = 0
            ny_patch   = 0
            npatch     = 0
            noutliers  = 0
        endif
        open(newunit=funit, file=fname%to_char(), access='stream', form='unformatted',&
            &status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot create motion model test fixture')
        write(funit) file_version, model_version
        call write_reference_string(funit, 'tmp_motion_model_movie.mrc')
        call write_reference_string(funit, gain)
        write(funit) smpd_movie, ldim_movie, smpd, ldim, binning
        write(funit) nframes, total_nframes
        write(funit) voltage, dose_per_frame, target_dose_per_frame, total_dose, accumulated_dose
        write(funit) rmsd_fit, fixed_frame
        flag = 1_int8
        write(funit) flag
        flag = 0_int8
        write(funit) flag
        write(funit) eer_fraction
        write(funit) drift_x, drift_y, frameweights
        write(funit) npatch, nx_patch, ny_patch, ldim_patch
        if( with_optional_arrays )then
            write(funit) patch_bounds, patch_coords
            write(funit) coeffs_x, coeffs_y
            write(funit) local_x, local_y
        endif
        write(funit) noutliers
        if( with_optional_arrays ) write(funit) outlier_coords
        close(funit)
        deallocate(drift_x, drift_y, frameweights)
        if( allocated(patch_bounds) ) deallocate(patch_bounds)
        if( allocated(patch_coords) ) deallocate(patch_coords)
        if( allocated(local_x) ) deallocate(local_x)
        if( allocated(local_y) ) deallocate(local_y)
        if( allocated(outlier_coords) ) deallocate(outlier_coords)
    end subroutine write_reference_binary

    subroutine write_reference_string( funit, text )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: text
        integer(int32) :: n
        n = len(text)
        write(funit) n
        if( n > 0 ) write(funit) text
    end subroutine write_reference_string

    logical function binary_files_equal( fname1, fname2 )
        type(string), intent(in) :: fname1, fname2
        integer(int64) :: nbytes1, nbytes2
        integer     :: unit1, unit2, ios
        integer(int8), allocatable :: bytes1(:), bytes2(:)
        binary_files_equal = .false.
        inquire(file=fname1%to_char(), size=nbytes1, iostat=ios)
        if( ios /= 0 ) return
        inquire(file=fname2%to_char(), size=nbytes2, iostat=ios)
        if( ios /= 0 .or. nbytes1 /= nbytes2 ) return
        allocate(bytes1(nbytes1), bytes2(nbytes2))
        open(newunit=unit1, file=fname1%to_char(), access='stream', form='unformatted',&
            &status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        open(newunit=unit2, file=fname2%to_char(), access='stream', form='unformatted',&
            &status='old', action='read', iostat=ios)
        if( ios /= 0 )then
            close(unit1)
            return
        endif
        read(unit1) bytes1
        read(unit2) bytes2
        close(unit1)
        close(unit2)
        binary_files_equal = all(bytes1 == bytes2)
        deallocate(bytes1, bytes2)
    end function binary_files_equal

end module simple_motion_model_tester
