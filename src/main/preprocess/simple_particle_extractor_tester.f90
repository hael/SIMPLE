!@descr: unit tests for binary-model fallback and initialization in simple_particle_extractor
module simple_particle_extractor_tester
use simple_defs,               only: logfhandle
use simple_ori,                only: ori
use simple_parameters,         only: parameters
use simple_particle_extractor, only: ptcl_extractor
use simple_syslib,             only: file_exists
use simple_test_utils,         only: assert_true, assert_false, assert_int
implicit none
private
public :: run_all_particle_extractor_tests

contains

    subroutine run_all_particle_extractor_tests()
        integer :: saved_log, scratch_log
        write(*,'(A)') '**** running all particle extractor tests ****'
        saved_log = logfhandle
        open(newunit=scratch_log, status='scratch', action='write')
        logfhandle = scratch_log
        call test_missing_motion_model()
        call test_fallback_releases_movie_state()
        call test_init_mic_selects_micrograph_mode()
        call test_binary_model_constructors()
        logfhandle = saved_log
        close(scratch_log)
    end subroutine run_all_particle_extractor_tests

    subroutine test_missing_motion_model()
        character(len=*), parameter :: MISSING_MODEL = 'tmp_particle_extractor_missing.mmodel'
        class(ptcl_extractor), allocatable :: extractor
        class(parameters), allocatable, target :: params
        type(ori) :: omic
        write(*,'(A)') 'test_missing_motion_model'
        allocate(extractor, params)
        call omic%new(is_ptcl=.false.)
        params%box = 8
        params%pcontrast = 'black'
        call extractor%init_mov(omic, params)
        call assert_fallback(extractor, 8, .true.)
        ! A legacy sidecar entry must not revive the retired metadata reader.
        call omic%set('mc_starfile', 'tmp_particle_extractor_unused.star')
        call omic%set('mcmodel', '')
        params%box = 12
        params%pcontrast = 'white'
        call extractor%init_mov(omic, params)
        call assert_fallback(extractor, 12, .false.)
        call assert_false(file_exists(MISSING_MODEL), 'the missing-model fixture path does not exist')
        if( .not.file_exists(MISSING_MODEL) )then
            call omic%set('mcmodel', MISSING_MODEL)
            call extractor%init_mov(omic, params)
            call assert_fallback(extractor, 12, .false.)
        endif
        call extractor%kill()
        call omic%kill()
    end subroutine test_missing_motion_model

    subroutine test_fallback_releases_movie_state()
        class(ptcl_extractor), allocatable :: extractor
        class(parameters), allocatable, target :: params
        type(ori) :: omic
        write(*,'(A)') 'test_fallback_releases_movie_state'
        allocate(extractor, params)
        call omic%new(is_ptcl=.false.)
        params%box = 8
        params%pcontrast = 'black'
        ! Seed retained frame metadata to exercise reuse across micrographs.
        call extractor%set_test_movie_state(2)
        call assert_int(2, extractor%get_nframes(), 'the retained-state fixture has two frames')
        call assert_true(extractor%does_exist(), 'the retained-state fixture marks the extractor active')
        call assert_true(extractor%from_mov(), 'the retained-state fixture selects movie mode')
        call assert_true(extractor%from_model(), 'the retained-state fixture has model metadata')
        call assert_true(extractor%has_frames(), 'the retained-state fixture allocates frames')
        call assert_true(extractor%has_isoshifts(), 'the retained-state fixture allocates movie shifts')
        call assert_true(extractor%has_weights(), 'the retained-state fixture allocates frame weights')
        call assert_true(extractor%has_hotpix_coords(), 'the retained-state fixture allocates defect coordinates')
        call assert_true(extractor%has_model_frameweights(), 'the retained-state fixture allocates model frame weights')
        call extractor%init_mov(omic, params)
        call assert_fallback(extractor, 8, .true.)
        call assert_false(extractor%has_isoshifts(), 'fallback releases the previous movie shifts')
        call assert_false(extractor%has_weights(), 'fallback releases the previous frame weights')
        call assert_false(extractor%has_hotpix_coords(), 'fallback releases the previous defect coordinates')
        call assert_false(extractor%has_model_frameweights(), 'fallback releases the previous binary model')
        call extractor%kill()
        call omic%kill()
    end subroutine test_fallback_releases_movie_state

    subroutine test_init_mic_selects_micrograph_mode()
        class(ptcl_extractor), allocatable :: extractor
        logical, allocatable :: mask(:,:,:)
        write(*,'(A)') 'test_init_mic_selects_micrograph_mode'
        allocate(extractor)
        call extractor%get_particle_mask(mask)
        call assert_false(allocated(mask), 'an uninitialized extractor returns an unallocated mask')
        call extractor%init_mic(8, .false.)
        call assert_false(extractor%from_mov(), 'direct micrograph initialization selects micrograph extraction')
        call assert_false(extractor%get_l_neg(), 'direct micrograph initialization preserves contrast')
        call extractor%get_particle_mask(mask)
        call assert_true(allocated(mask), 'direct micrograph initialization returns an allocated mask')
        if( allocated(mask) )then
            mask = .false.
            call extractor%get_particle_mask(mask)
            call assert_true(any(mask), 'changing a returned mask does not modify the private normalization mask')
        endif
        call extractor%kill()
        call extractor%get_particle_mask(mask)
        call assert_false(allocated(mask), 'a killed extractor clears a previously allocated output mask')
        call assert_false(extractor%does_exist(), 'killing the extractor clears its active state')
    end subroutine test_init_mic_selects_micrograph_mode

    subroutine test_binary_model_constructors()
        use simple_defs,              only: dp, nthr_glob
        use simple_image,             only: image
        use simple_motion_model,      only: motion_model
        use simple_particle_extractor, only: dealloc_particles_frames
        use simple_string,            only: string
        use simple_syslib,            only: del_file
        use simple_test_utils,        only: assert_real
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        integer, parameter :: NX = 48, NY = 40, NFR = 2, BOX_TEST = 8
        integer, parameter :: DRIFT(2,NFR) = reshape([0,0,2,-1], [2,NFR])
        real, parameter :: PIXEL_SIZE = 1.5, DOSE = 6., FRAME_WEIGHTS(NFR) = [0.25,0.75]
        ! Allow single-precision I/O, FFTs and background normalization roundoff.
        real, parameter :: PIXEL_TOL = 1024. * epsilon(1.)
        class(ptcl_extractor), allocatable :: extractor, mic_extractor
        class(parameters), allocatable, target :: params
        class(motion_model), allocatable :: model
        class(image), pointer :: frames(:)
        type(image), allocatable :: actual(:), expected(:)
        type(image) :: gain_image, truth_image
        type(ori) :: omic
        type(string) :: movie, gain, model_file, star_file
        integer, allocatable :: pinds(:)
        real :: pixels(NX,NY,1), gain_pixels(NX,NY,1), vmin, vmax, vmean, vsdev
        real(dp) :: omega(2), frequency(2), critical_dose(2), dw(2,NFR), weights(NFR), phase(2), pi
        integer :: coords(2,1), shift(2), i, j, iframe, imode, ref_frame, saved_nthr
        !$ integer :: saved_omp
        character(len=96) :: message
        write(*,'(A)') 'test_binary_model_constructors'
        saved_nthr = nthr_glob
        !$ saved_omp = omp_get_max_threads()
        nthr_glob = 1
        !$ call omp_set_num_threads(1)
        allocate(extractor, mic_extractor, params, model, frames(NFR), actual(1), expected(1), pinds(1))
        movie = 'tmp_particle_extractor_movie.mrc'
        gain = 'tmp_particle_extractor_gain.mrc'
        model_file = 'tmp_particle_extractor_model.bin'
        star_file = 'tmp_particle_extractor_model.star'
        pi = acos(-1._dp)
        omega = 2._dp*pi*[4._dp/real(NX,dp), 7._dp/real(NY,dp)]
        frequency = omega / (2._dp*pi*real(PIXEL_SIZE,dp))
        do i = 1,NX
            gain_pixels(i,:,1) = 1. + real(i*i)/1024.
        enddo
        call gain_image%new([NX,NY,1], PIXEL_SIZE, wthreads=.false.)
        call gain_image%set_rmat(gain_pixels, .false.)
        call gain_image%write(gain, del_if_exists=.true.)
        do iframe = 1,NFR
            do j = 1,NY
                do i = 1,NX
                    pixels(i,j,1) = real(sin(omega(1)*real(i-1,dp) + 0.3_dp*iframe) + &
                        &0.7_dp*cos(omega(2)*real(j-1,dp) - 0.2_dp*iframe)) / gain_pixels(NX-i+1,j,1)
                enddo
            enddo
            call frames(iframe)%new([NX,NY,1], PIXEL_SIZE, wthreads=.false.)
            call frames(iframe)%set_rmat(pixels, .false.)
            call frames(iframe)%write(movie, iframe, del_if_exists=(iframe == 1))
        enddo
        params%box = BOX_TEST
        params%pcontrast = 'white'
        params%tof = NFR
        params%flipgain = 'x'
        params%fraction_dose_target = DOSE
        params%total_dose = real(NFR)*DOSE
        params%eer_upsampling = 1
        call model%new(params, movie, [NX,NY], PIXEL_SIZE, frames, NFR, 2, 300., DOSE, 0, gain=gain)
        call model%set_drift_offsets(real(DRIFT(1,:)), real(DRIFT(2,:)))
        call model%set_frameweights(FRAME_WEIGHTS)
        ! The constructor, not the saved DW flag or live flip parameter, selects preparation.
        model%dw = .false.
        call model%write(star_file, model_file, .false.)
        call del_file(star_file)
        params%flipgain = 'no'
        call omic%new(is_ptcl=.false.)
        call omic%set('mcmodel', model_file)
        call assert_false(file_exists(star_file), 'valid-model extraction needs no STAR sidecar')
        coords(:,1) = [20,16]
        pinds = 1
        call actual(1)%new([BOX_TEST,BOX_TEST,1], PIXEL_SIZE, wthreads=.false.)
        call expected(1)%new([BOX_TEST,BOX_TEST,1], PIXEL_SIZE, wthreads=.false.)
        call truth_image%new([NX,NY,1], PIXEL_SIZE, wthreads=.false.)
        call mic_extractor%init_mic(BOX_TEST, .false.)
        do imode = 1,2
            if( imode == 1 )then
                call extractor%init_model(omic, params)
                ref_frame = 2
                weights = 1._dp
                dw = 1._dp
            else
                call extractor%init_mov(omic, params)
                ref_frame = 1
                weights = real(FRAME_WEIGHTS,dp)
                ! Two sinusoidal Fourier modes: normalized dose weights in closed form.
                critical_dose = 0.245_dp*frequency**(-1.665_dp) + 2.81_dp
                dw(:,1) = 1._dp / sqrt(1._dp + exp(-real(DOSE,dp)/critical_dose))
                dw(:,2) = dw(:,1)*exp(-real(DOSE,dp)/(2._dp*critical_dose))
            endif
            call assert_true(extractor%does_exist(), 'both binary-model constructors complete initialization')
            call assert_true(extractor%from_model(), 'both constructors retain binary-model metadata')
            call assert_true(extractor%from_mov() .eqv. (imode == 2), 'constructors retain distinct extraction modes')
            call assert_int(NFR, extractor%get_nframes(), 'both constructors copy the model frame count')
            pixels = 0.
            do iframe = 1,NFR
                shift = DRIFT(:,iframe) - DRIFT(:,ref_frame)
                do j = 1,NY
                    do i = 1,NX
                        phase = omega*real([i-1,j-1]-shift,dp) + [0.3_dp,-0.2_dp]*iframe
                        pixels(i,j,1) = pixels(i,j,1) + real(weights(iframe) * &
                            &(dw(1,iframe)*sin(phase(1)) + 0.7_dp*dw(2,iframe)*cos(phase(2))))
                    enddo
                enddo
            enddo
            call truth_image%set_rmat(pixels, .false.)
            ! Share only the downstream normalization, not movie preparation or motion evaluation.
            call mic_extractor%extract_particles_from_mic(truth_image, pinds, coords, expected, vmin, vmax, vmean, vsdev)
            if( imode == 1 )then
                call extractor%extract_all_particles_frames(coords, [1,NFR], NFR)
                call extractor%generate_cumulative_particles_stack(1, actual, [1,NFR], vmin, vmax, vmean, vsdev)
                call dealloc_particles_frames()
            else
                call extractor%extract_particles(pinds, coords, actual, vmin, vmax, vmean, vsdev)
            endif
            write(message,'(A,I0)') 'gain, drift reference and weighting match analytic particles for constructor ',imode
            call assert_real(0., maxval(abs(actual(1)%get_rmat()-expected(1)%get_rmat())), PIXEL_TOL, trim(message))
        enddo
        call extractor%kill()
        call mic_extractor%kill()
        call model%kill()
        call omic%kill()
        call gain_image%kill()
        call truth_image%kill()
        call actual(1)%kill()
        call expected(1)%kill()
        do iframe = 1,NFR
            call frames(iframe)%kill()
        enddo
        deallocate(frames)
        call del_file(movie)
        call del_file(gain)
        call del_file(model_file)
        nthr_glob = saved_nthr
        !$ call omp_set_num_threads(saved_omp)
    end subroutine test_binary_model_constructors

    subroutine assert_fallback( extractor, box, neg )
        class(ptcl_extractor), intent(in) :: extractor
        integer,               intent(in) :: box
        logical,               intent(in) :: neg
        logical, allocatable :: mask(:,:,:)
        call assert_true(extractor%does_exist(), 'fallback leaves an initialized extractor')
        call assert_false(extractor%from_mov(), 'missing binary metadata selects micrograph extraction')
        call assert_false(extractor%from_model(), 'missing binary metadata leaves no active motion model')
        call assert_int(0, extractor%get_nframes(), 'fallback retains no movie frames')
        call assert_false(extractor%has_frames(), 'fallback does not allocate movie frames')
        call assert_int(box, extractor%get_box(), 'fallback uses the requested particle box')
        call assert_true(extractor%get_l_neg() .eqv. neg, 'fallback uses the requested particle contrast')
        call extractor%get_particle_mask(mask)
        call assert_true(allocated(mask), 'fallback initializes the normalization mask')
        if( allocated(mask) )then
            call assert_true(all(shape(mask) == [box,box,1]), 'fallback mask matches the particle box')
            call assert_true(any(mask), 'fallback mask contains foreground pixels')
            call assert_true(any(.not.mask), 'fallback mask leaves background pixels')
        endif
    end subroutine assert_fallback

end module simple_particle_extractor_tester
