!@descr: unit tests for binary persistence and polynomial refitting of simple_motion_model
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

character(len=2), parameter :: FLIP_MODES(5) = ['no', 'x ', 'y ', 'xy', 'yx']

contains

    subroutine run_all_motion_model_tests()
        write(*,'(A)') '**** running all motion model tests ****'
        call test_flipgain_lifecycle()
        call test_binary_roundtrip_with_optional_arrays()
        call test_binary_roundtrip_with_rejected_patch()
        call test_binary_roundtrip_without_optional_arrays()
        call test_refit_polynomial_known_coefficients()
    end subroutine run_all_motion_model_tests

    subroutine test_flipgain_lifecycle()
        use simple_image, only: image
        character(len=2), parameter :: EXPECTED_MODES(5) = ['NO', 'X ', 'Y ', 'XY', 'YX']
        class(motion_model), allocatable :: model
        class(parameters), allocatable, target :: params
        class(image), pointer :: frames(:)
        type(string) :: movie, gain
        integer :: imode
        write(*,'(A)') 'test_flipgain_lifecycle'
        allocate(model, params)
        allocate(frames(1))
        call frames(1)%new([6,4,1], 1.5, wthreads=.false.)
        movie = 'tmp_motion_model_flipgain_movie.mrc'
        gain  = 'tmp_motion_model_flipgain_gain.mrc'
        call frames(1)%write(movie, del_if_exists=.true.)
        params%fraction_dose_target = 1.5
        params%total_dose = 1.0
        params%eer_upsampling = 1
        call assert_char('no', model%flipgain, 'motion model gain flipping defaults to no')
        do imode = 1,size(FLIP_MODES)
            params%flipgain = FLIP_MODES(imode)
            call model%new(params, movie, [6,4], 1.5, frames, 1, 1, 300., 1., 0, gain=gain)
            call assert_char(EXPECTED_MODES(imode), model%flipgain,&
                &'constructor stores the gain flipping mode in uppercase')
            params%flipgain = FLIP_MODES(mod(imode,size(FLIP_MODES))+1)
            call assert_char(EXPECTED_MODES(imode), model%flipgain,&
                &'gain flipping is independent of later parameter changes')
            call model%kill()
            call assert_char('no', model%flipgain, 'killing an active model resets gain flipping')
        enddo
        model%flipgain = 'xy'
        call model%kill()
        call assert_char('no', model%flipgain, 'killing an inactive model also resets gain flipping')
        call frames(1)%kill()
        deallocate(frames)
        call del_file(movie)
    end subroutine test_flipgain_lifecycle

    subroutine test_refit_polynomial_known_coefficients()
        integer, parameter :: NFRAMES_TEST = 7, NGRID = 3, REF_FRAMES(2) = [1,4]
        ! Three samples per spatial axis and six nonzero times give full rank.
        ! Dyadic coordinates/coefficients keep the fixture exact in single precision.
        real(dp), parameter :: GRID(NGRID) = [-0.5_dp, 0.0_dp, 0.5_dp]
        real(dp), parameter :: COEFF_TOL = 1.e-10_dp
        real, parameter :: RMSD_TOL = 32.0 * epsilon(1.0) ! single-precision evaluation/accumulation
        real(dp), parameter :: RESIDUAL_SCALE = 1.0_dp / 65536.0_dp
        type(motion_model) :: model
        real(dp) :: coeffs(3,6,2), expected(3,6,2), spatial(6), tau, delta, value(2)
        real(dp) :: residual(NFRAMES_TEST), expected_rmsd
        real :: local_x(NFRAMES_TEST,NGRID,NGRID), local_y(NFRAMES_TEST,NGRID,NGRID)
        integer :: i, j, iframe, itest, ref_frame, term, axis
        character(len=96) :: message
        write(*,'(A)') 'test_refit_polynomial_known_coefficients'
        ! Every temporal/spatial coefficient is nonzero; X and Y are independent.
        coeffs(:,:,1) = reshape([ &
            & 3, -2,  1,   5,  3, -2,  -4,  2,  3, &
            & 7, -3,  2,  -5,  4, -1,   6, -2, -3], [3,6]) / 1024.0_dp
        coeffs(:,:,2) = reshape([ &
            &-2,  5, -3,   4, -1,  2,   7, -4,  1, &
            &-6,  2,  3,   3, -5,  2,  -1,  4, -2], [3,6]) / 1024.0_dp
        model%nframes = NFRAMES_TEST
        model%ldim = [129,97] ! non-square dimensions catch an X/Y normalization mix-up
        model%nx_patch = NGRID
        model%ny_patch = NGRID
        model%npatch = NGRID*NGRID
        model%fixed_frame = 1
        allocate(model%patch_coords(NGRID,NGRID,2))
        ! This in-memory fixture needs neither movie images nor parameter parsing.
        model%exists = .true.
        do j = 1,NGRID
            do i = 1,NGRID
                model%patch_coords(i,j,:) = real(1.0_dp + (GRID([i,j])+0.5_dp)*real(model%ldim-1,dp))
                spatial = [1.0_dp, GRID(i), GRID(i)**2, GRID(j), GRID(j)**2, GRID(i)*GRID(j)]
                do iframe = 1,NFRAMES_TEST
                    tau = real(iframe-1,dp)
                    ! Independent Horner evaluation of the separable cubic field;
                    ! do not use patch_poly/apply_patch_poly to generate the truth.
                    do axis = 1,2
                        value(axis) = sum(spatial * &
                            &((coeffs(3,:,axis)*tau + coeffs(2,:,axis))*tau + coeffs(1,:,axis))) * tau
                    enddo
                    local_x(iframe,i,j) = real(value(1))
                    local_y(iframe,i,j) = real(value(2))
                enddo
            enddo
        enddo
        call model%set_local_offsets(local_x, local_y)
        do itest = 1,size(REF_FRAMES)
            ref_frame = REF_FRAMES(itest)
            delta = real(ref_frame-1,dp)
            ! For q(t)=a*t+b*t^2+c*t^3, q(t+d)-q(d) has coefficients
            ! [a+2*b*d+3*c*d^2, b+3*c*d, c]. Apply this to each spatial monomial.
            expected(1,:,:) = coeffs(1,:,:) + 2.0_dp*delta*coeffs(2,:,:) + 3.0_dp*delta**2*coeffs(3,:,:)
            expected(2,:,:) = coeffs(2,:,:) + 3.0_dp*delta*coeffs(3,:,:)
            expected(3,:,:) = coeffs(3,:,:)
            call model%refit_polynomial(ref_frame)
            do term = 1,6
                write(message,'(A,I0,A,I0)') 'refit recovers X cubic coefficients, spatial term ',term,', reference ',ref_frame
                call assert_true(all(abs(model%model_coeffs_x(3*term-2:3*term)-expected(:,term,1)) <= COEFF_TOL), message)
                write(message,'(A,I0,A,I0)') 'refit recovers Y cubic coefficients, spatial term ',term,', reference ',ref_frame
                call assert_true(all(abs(model%model_coeffs_y(3*term-2:3*term)-expected(:,term,2)) <= COEFF_TOL), message)
            enddo
            call assert_true(all(abs(model%rmsd_fit) <= RMSD_TOL), 'an exact cubic motion field has zero refit RMSD')
            call assert_true(all(model%local_offsets_x == local_x) .and. all(model%local_offsets_y == local_y), &
                &'refitting preserves the measured local offsets')
            call assert_int(1, model%fixed_frame, 'refitting leaves stored reference-frame metadata unchanged')
        enddo
        ! Add a known residual orthogonal to t, t^2 and t^3 on t=-3,...,3.
        ! sum(t^4)=196 and sum(t^6)=1588, so r=49*t^4-397*t^2 is orthogonal
        ! to t^2; symmetry handles odd powers. r(0)=0 preserves re-referencing.
        do iframe = 1,NFRAMES_TEST
            tau = real(iframe-REF_FRAMES(2),dp)
            residual(iframe) = (49.0_dp*tau**4 - 397.0_dp*tau**2) * RESIDUAL_SCALE
            model%local_offsets_x(iframe,:,:) = local_x(iframe,:,:) + real(residual(iframe))
            model%local_offsets_y(iframe,:,:) = local_y(iframe,:,:) - real(2.0_dp*residual(iframe))
        enddo
        call model%refit_polynomial(REF_FRAMES(2))
        call assert_true(all(abs(model%model_coeffs_x-reshape(expected(:,:,1),[18])) <= COEFF_TOL), &
            &'orthogonal residuals leave the fitted X coefficients unchanged')
        call assert_true(all(abs(model%model_coeffs_y-reshape(expected(:,:,2),[18])) <= COEFF_TOL), &
            &'orthogonal residuals leave the fitted Y coefficients unchanged')
        expected_rmsd = sqrt(2.0_dp*(348.0_dp**2+804.0_dp**2+396.0_dp**2)/7.0_dp) * RESIDUAL_SCALE
        call assert_real(real(expected_rmsd), model%rmsd_fit(1), RMSD_TOL, 'X fit RMSD matches the known residual')
        call assert_real(real(2.0_dp*expected_rmsd), model%rmsd_fit(2), RMSD_TOL, 'Y fit RMSD matches the known residual')
        call model%kill()
    end subroutine test_refit_polynomial_known_coefficients

    subroutine test_binary_roundtrip_with_optional_arrays()
        write(*,'(A)') 'test_binary_roundtrip_with_optional_arrays'
        call test_binary_roundtrip(&
            &string('tmp_motion_model_patched_input.bin'),&
            &string('tmp_motion_model_patched_output.bin'),&
            &string('tmp_motion_model_patched_output.star'), .true., .true.)
    end subroutine test_binary_roundtrip_with_optional_arrays

    subroutine test_binary_roundtrip_with_rejected_patch()
        write(*,'(A)') 'test_binary_roundtrip_with_rejected_patch'
        call test_binary_roundtrip(&
            &string('tmp_motion_model_rejected_patch_input.bin'),&
            &string('tmp_motion_model_rejected_patch_output.bin'),&
            &string('tmp_motion_model_rejected_patch_output.star'), .true., .false.)
    end subroutine test_binary_roundtrip_with_rejected_patch

    subroutine test_binary_roundtrip_without_optional_arrays()
        write(*,'(A)') 'test_binary_roundtrip_without_optional_arrays'
        call test_binary_roundtrip(&
            &string('tmp_motion_model_minimal_input.bin'),&
            &string('tmp_motion_model_minimal_output.bin'),&
            &string('tmp_motion_model_minimal_output.star'), .false., .false.)
    end subroutine test_binary_roundtrip_without_optional_arrays

    subroutine test_binary_roundtrip( input_bin, output_bin, output_star, with_optional_arrays, patch_accepted )
        type(string), intent(in) :: input_bin, output_bin, output_star
        logical,      intent(in) :: with_optional_arrays, patch_accepted
        class(motion_model), allocatable :: model
        class(parameters), allocatable, target :: params
        class(cmdline), allocatable :: cline
        integer :: imode
        allocate(model, params, cline)
        ! a program outside every UI table: parameters is only the container
        ! here, and a registered program (motion_correct requires a project)
        ! made the result depend on whether an earlier suite built the UI
        ! (the SIMPLE_UNIT_ORDER=reverse run of test=units failed)
        call cline%set('prg', 'motion_model_tester')
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
        do imode = 1,size(FLIP_MODES)
            ! Deliberately disagree with the file: read must use stored metadata.
            params%flipgain = FLIP_MODES(mod(imode,size(FLIP_MODES))+1)
            call write_reference_binary(input_bin, with_optional_arrays, patch_accepted, FLIP_MODES(imode))
            call model%read(input_bin, params)
            call assert_char(FLIP_MODES(imode), model%flipgain, 'binary gain flipping overrides runtime parameters')
            call assert_true(model%patch_accepted .eqv. patch_accepted,&
                &'motion model read preserves patch acceptance')
            call model%write(output_star, output_bin, patch_accepted)
            call assert_true(file_exists(output_bin), 'motion model roundtrip writes a binary model')
            call assert_true(binary_files_equal(input_bin, output_bin),&
                &'motion model read/write preserves the independently generated binary payload')
        enddo
        call model%kill()
        call del_file(input_bin)
        call del_file(output_bin)
        call del_file(output_star)
        call cline%kill()
    end subroutine test_binary_roundtrip

    subroutine write_reference_binary( fname, with_optional_arrays, patch_accepted, flipgain )
        type(string), intent(in) :: fname
        logical,      intent(in) :: with_optional_arrays, patch_accepted
        character(len=*), intent(in) :: flipgain
        integer(int8) :: flag
        integer       :: funit, ios, i, noutliers
        integer       :: file_version, model_version
        integer       :: ldim_movie(2), ldim(2), ldim_patch(2)
        integer       :: nframes, total_nframes, npatch, nx_patch, ny_patch
        integer       :: fixed_frame, eer_fraction, eer_upsampling
        integer, allocatable :: patch_bounds(:,:,:,:), outlier_coords(:,:)
        real          :: smpd_movie, smpd, binning
        real          :: voltage, dose_per_frame, target_dose_per_frame
        real          :: total_dose, accumulated_dose, rmsd_fit(2)
        real, allocatable :: drift_x(:), drift_y(:), frameweights(:)
        real, allocatable :: patch_coords(:,:,:), local_x(:,:,:), local_y(:,:,:)
        real(dp)      :: coeffs_x(18), coeffs_y(18)
        character(len=:), allocatable :: gain
        file_version         = 1
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
        eer_upsampling       = 1
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
        call write_reference_string(funit, trim(flipgain))
        write(funit) smpd_movie, ldim_movie, smpd, ldim, binning
        write(funit) nframes, total_nframes
        write(funit) voltage, dose_per_frame, target_dose_per_frame, total_dose, accumulated_dose
        write(funit) rmsd_fit, fixed_frame
        flag = 1_int8
        write(funit) flag
        flag = 0_int8
        write(funit) flag
        flag = merge(1_int8, 0_int8, patch_accepted)
        write(funit) flag
        write(funit) eer_fraction, eer_upsampling
        write(funit) drift_x, drift_y
        write(funit) frameweights
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
