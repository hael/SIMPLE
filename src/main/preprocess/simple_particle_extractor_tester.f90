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
        call test_init_mic_selects_micrograph_mode()
        call test_binary_model_constructors()
        call test_particle_frame_cache()
        call test_movie_outlier_curation()
        logfhandle = saved_log
        close(scratch_log)
    end subroutine run_all_particle_extractor_tests

    subroutine test_movie_outlier_curation()
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        class(ptcl_extractor), allocatable :: extractor
        logical :: gain_only, combined, flat_gain_only, flat_combined, no_gain, non_eer, clean_gain
        !$ integer :: saved_omp
        write(*,'(A)') 'test_movie_outlier_curation'
        !$ saved_omp = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        allocate(extractor)
        call extractor%test_movie_outlier_curation(gain_only, combined, flat_gain_only, flat_combined,&
            &no_gain, non_eer, clean_gain)
        call assert_true(gain_only, 'movie curation repairs EER gain defects with no saved outliers')
        call assert_true(combined, 'movie curation repairs statistical and gain defects including saved overlap')
        call assert_true(flat_gain_only, 'gain-only repair survives constant frame-sum statistics')
        call assert_true(flat_combined, 'statistical curation still repairs gain defects with a constant frame sum')
        call assert_true(no_gain, 'no-gain movie curation leaves frames unchanged with no saved outliers')
        call assert_true(non_eer, 'non-EER movie curation leaves frames unchanged with no saved outliers')
        call assert_true(clean_gain, 'clean EER gain does not enable statistical re-detection without saved outliers')
        deallocate(extractor)
        !$ call omp_set_num_threads(saved_omp)
    end subroutine test_movie_outlier_curation

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
        call assert_micrograph_state(extractor, 8, .true.)
        ! A legacy sidecar entry must not revive the retired metadata reader.
        call omic%set('mc_starfile', 'tmp_particle_extractor_unused.star')
        call omic%set('mcmodel', '')
        params%box = 12
        params%pcontrast = 'white'
        call extractor%init_mov(omic, params)
        call assert_micrograph_state(extractor, 12, .false.)
        call assert_false(file_exists(MISSING_MODEL), 'the missing-model fixture path does not exist')
        if( .not.file_exists(MISSING_MODEL) )then
            call omic%set('mcmodel', MISSING_MODEL)
            call extractor%init_mov(omic, params)
            call assert_micrograph_state(extractor, 12, .false.)
        endif
        call extractor%kill()
        call omic%kill()
    end subroutine test_missing_motion_model

    subroutine test_init_mic_selects_micrograph_mode()
        class(ptcl_extractor), allocatable :: extractor
        logical, allocatable :: mask(:,:,:)
        write(*,'(A)') 'test_init_mic_selects_micrograph_mode'
        allocate(extractor)
        call extractor%get_particle_mask(mask)
        call assert_false(allocated(mask), 'an uninitialized extractor returns an unallocated mask')
        call extractor%init_mic(8, .false.)
        call assert_micrograph_state(extractor, 8, .false.)
        call extractor%get_particle_mask(mask)
        call assert_true(allocated(mask), 'direct micrograph initialization returns an allocated mask')
        if( allocated(mask) )then
            mask = .false.
            call extractor%get_particle_mask(mask)
            call assert_true(any(mask), 'changing a returned mask does not modify the private normalization mask')
        endif
        call extractor%init_mic(12, .true.)
        call assert_micrograph_state(extractor, 12, .true.)
        call extractor%init_mic(6, .false.)
        call assert_micrograph_state(extractor, 6, .false.)
        call extractor%kill()
        call extractor%kill()
        call extractor%get_particle_mask(mask)
        call assert_false(allocated(mask), 'a killed extractor clears a previously allocated output mask')
        call assert_false(extractor%does_exist(), 'repeated kill leaves the extractor inactive')
        call assert_int(0, extractor%get_nframes(), 'repeated kill leaves no movie frames')
        call assert_int(0, extractor%get_box(), 'repeated kill clears the particle box')
        call extractor%init_mic(8, .true.)
        call assert_micrograph_state(extractor, 8, .true.)
        call extractor%kill()
    end subroutine test_init_mic_selects_micrograph_mode

    subroutine test_binary_model_constructors()
        use simple_defs,              only: dp
        use simple_image,             only: image
        use simple_motion_model,      only: motion_model
        use simple_string,            only: string
        use simple_syslib,            only: del_file
        use simple_test_utils,        only: assert_real
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        integer, parameter :: NX = 48, NY = 40, NFR = 3, BOX_TEST = 8, REF_FRAME = 2
        integer, parameter :: DRIFT(2,NFR) = reshape([0,0,2,-1,-1,2], [2,NFR])
        real, parameter :: PIXEL_SIZE = 1.5, DOSE = 6., FRAME_WEIGHTS(NFR) = [0.25,0.5,0.25]
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
        integer :: coords(2,1), shift(2), i, j, iframe, imode
        !$ integer :: saved_omp
        character(len=96) :: message
        write(*,'(A)') 'test_binary_model_constructors'
        !$ saved_omp = omp_get_max_threads()
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
        call model%new(params, movie, [NX,NY], PIXEL_SIZE, frames, NFR, REF_FRAME, 300., DOSE, 0, gain=gain)
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
                weights = 1._dp
                dw = 1._dp
            else
                call extractor%init_mov(omic, params)
                weights = real(FRAME_WEIGHTS,dp)
                ! Two sinusoidal Fourier modes: normalized dose weights in closed form.
                critical_dose = 0.245_dp*frequency**(-1.665_dp) + 2.81_dp
                do iframe = 1,NFR
                    dw(:,iframe) = exp(-real(iframe-1,dp)*real(DOSE,dp)/(2._dp*critical_dose))
                enddo
                do i = 1,size(dw,dim=1)
                    dw(i,:) = dw(i,:) / sqrt(sum(dw(i,:)**2))
                enddo
            endif
            call assert_true(extractor%does_exist(), 'both binary-model constructors complete initialization')
            call assert_false(extractor%has_particle_frames([1,NFR], 1), 'both constructors invalidate cached particles')
            call assert_true(extractor%from_model(), 'both constructors retain binary-model metadata')
            call assert_true(extractor%from_mov() .eqv. (imode == 2), 'constructors retain distinct extraction modes')
            call assert_int(NFR, extractor%get_nframes(), 'both constructors copy the model frame count')
            call assert_dose_weighting_report(extractor, imode == 2, .false.)
            pixels = 0.
            do iframe = 1,NFR
                shift = DRIFT(:,iframe) - DRIFT(:,REF_FRAME)
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
                call extractor%extract_all_particles_frames(coords, [1,NFR])
                call extractor%generate_cumulative_particles_stack(1, actual, [1,NFR], vmin, vmax, vmean, vsdev)
            else
                call extractor%extract_particles(pinds, coords, actual, vmin, vmax, vmean, vsdev)
            endif
            write(message,'(A,I0)') 'gain, drift reference and weighting match analytic particles for constructor ',imode
            call assert_real(0., maxval(abs(actual(1)%get_rmat()-expected(1)%get_rmat())), PIXEL_TOL, trim(message))
        enddo
        model%dw = .true.
        call model%write(star_file, model_file, .false.)
        call del_file(star_file)
        do imode = 1,2
            if( imode == 1 )then
                call extractor%init_model(omic, params)
            else
                call extractor%init_mov(omic, params)
            endif
            call assert_dose_weighting_report(extractor, imode == 2, .true.)
            call extractor%extract_all_particles_frames(coords, [1,NFR])
            call extractor%init_mic(BOX_TEST, imode == 2)
            call assert_micrograph_state(extractor, BOX_TEST, imode == 2)
            call mic_extractor%init_mic(BOX_TEST, imode == 2)
            call mic_extractor%extract_particles_from_mic(truth_image, pinds, coords, expected, vmin, vmax, vmean, vsdev)
            call extractor%extract_particles_from_mic(truth_image, pinds, coords, actual, vmin, vmax, vmean, vsdev)
            write(message,'(A,I0)') 'reused micrograph extractor matches fresh initialization after constructor ',imode
            call assert_real(0., maxval(abs(actual(1)%get_rmat()-expected(1)%get_rmat())), PIXEL_TOL, trim(message))
        enddo
        call extractor%kill()
        call assert_dose_weighting_report(extractor, .false.)
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
        !$ call omp_set_num_threads(saved_omp)
    end subroutine test_binary_model_constructors

    subroutine test_particle_frame_cache()
        use simple_image,        only: image
        use simple_motion_model, only: motion_model
        use simple_string,       only: string
        use simple_syslib,       only: del_file
        use simple_test_utils,   only: assert_real
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        integer, parameter :: NX = 48, NY = 40, NFR = 3, REF_FRAME = 2, NPTCLS = 2
        integer, parameter :: BOXES(3) = [8,10,10]
        real, parameter :: SAMPLINGS(3) = [1.5,1.5,2.0], PIXEL_TOL = 1024. * epsilon(1.)
        integer, parameter :: COORDS(2,NPTCLS) = reshape([10,10,24,19], [2,NPTCLS])
        class(ptcl_extractor), allocatable :: extractor, other, mic_extractor
        class(parameters), allocatable, target :: params
        class(motion_model), allocatable :: model
        class(image), pointer :: frames(:)
        type(image), allocatable :: actual(:), expected(:)
        type(image) :: truth_image
        type(ori) :: omic, missing_model
        type(string) :: movie, model_file, star_file
        real :: movie_pixels(NX,NY,1,NFR), vmin, vmax, vmean, vsdev
        integer :: i, j, iframe, icase
        !$ integer :: saved_omp
        write(*,'(A)') 'test_particle_frame_cache'
        !$ saved_omp = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        allocate(extractor, other, mic_extractor, params, model, frames(NFR), actual(NPTCLS), expected(NPTCLS))
        movie = 'tmp_particle_frame_cache.mrc'
        model_file = 'tmp_particle_frame_cache.mmodel'
        star_file = 'tmp_particle_frame_cache.star'
        call omic%new(is_ptcl=.false.)
        call omic%set('mcmodel', model_file)
        call missing_model%new(is_ptcl=.false.)
        params%pcontrast = 'white'
        params%tof = NFR
        params%flipgain = 'no'
        params%fraction_dose_target = 0.
        params%total_dose = 0.
        params%eer_upsampling = 1
        do iframe = 1,NFR
            do j = 1,NY
                do i = 1,NX
                    movie_pixels(i,j,1,iframe) = sin(0.3*real(i)+0.2*real(iframe)) + &
                        &0.7*cos(0.5*real(j)-0.4*real(iframe))
                enddo
            enddo
        enddo
        call assert_false(extractor%has_particle_frames([1,NFR], 1), 'a new extractor has no particle-frame cache')
        ! Change box size, then sampling alone, while preserving frame bounds and particle count.
        do icase = 1,size(BOXES)
            params%box = BOXES(icase)
            do iframe = 1,NFR
                call frames(iframe)%new([NX,NY,1], SAMPLINGS(icase), wthreads=.false.)
                call frames(iframe)%set_rmat(movie_pixels(:,:,:,iframe), .false.)
                call frames(iframe)%write(movie, iframe, del_if_exists=(iframe == 1))
            enddo
            call model%new(params, movie, [NX,NY], SAMPLINGS(icase), frames, NFR, REF_FRAME, 300., 0., 0)
            call model%write(star_file, model_file, .false.)
            call del_file(star_file)
            call extractor%init_model(omic, params)
            call assert_false(extractor%has_particle_frames([1,NFR], 1), 'model initialization clears the previous cache')
            do i = 1,NPTCLS
                call actual(i)%new([params%box,params%box,1], SAMPLINGS(icase), wthreads=.false.)
                call expected(i)%new([params%box,params%box,1], SAMPLINGS(icase), wthreads=.false.)
            enddo
            call truth_image%new([NX,NY,1], SAMPLINGS(icase), wthreads=.false.)
            call mic_extractor%init_mic(params%box, .false.)
            call extractor%extract_all_particles_frames(COORDS, [1,NFR])
            call assert_true(extractor%has_particle_frames([1,NFR], NPTCLS), 'extraction populates the requested cache')
            call check_cached_sum([1,NFR], NPTCLS)
            if( icase == 1 )then
                ! Matching cache shapes must not let one extractor overwrite or destroy another's pixels.
                call other%init_model(omic, params)
                call other%extract_all_particles_frames(COORDS+3, [1,NFR])
                call check_cached_sum([1,NFR], NPTCLS)
                call other%kill()
                call other%kill()
                call assert_false(other%has_particle_frames([1,NFR], 1), 'repeated kill leaves no cache')
                call check_cached_sum([1,NFR], NPTCLS)
                call extractor%extract_all_particles_frames(COORDS(:,1:1), [2,NFR])
                call assert_true(extractor%has_particle_frames([2,NFR], 1), 'a partial cache keeps its actual lower bound')
                call assert_true(extractor%has_particle_frames([2,2], 1), 'a single cached frame is a valid interval')
                call assert_false(extractor%has_particle_frames([1,NFR], 1), 'movie bounds do not imply cached bounds')
                call assert_false(extractor%has_particle_frames([2,NFR+1], 1), 'frames above the cached interval are invalid')
                call assert_false(extractor%has_particle_frames([3,2], 1), 'reversed cached bounds are invalid')
                call assert_false(extractor%has_particle_frames([2,NFR], 0), 'an empty particle request is invalid')
                call assert_false(extractor%has_particle_frames([2,NFR], NPTCLS), 'a smaller extraction drops old particles')
                call check_cached_sum([2,NFR], 1)
                call check_cached_sum([2,2], 1)
                call extractor%extract_all_particles_frames(COORDS(:,1:1), [2,2])
                call assert_false(extractor%has_particle_frames([2,NFR], 1), 'a movie frame above the cache is invalid')
                call check_cached_sum([2,2], 1)
                call extractor%extract_all_particles_frames(COORDS, [1,NFR])
            endif
        enddo
        call extractor%init_mic(params%box, .false.)
        call assert_micrograph_state(extractor, params%box, .false.)
        call extractor%init_model(omic, params)
        call extractor%extract_all_particles_frames(COORDS, [1,NFR])
        call extractor%init_mov(missing_model, params)
        call assert_false(extractor%has_particle_frames([1,NFR], 1), 'movie fallback clears the cache')
        call extractor%kill()
        call mic_extractor%kill()
        call model%kill()
        call omic%kill()
        call missing_model%kill()
        call truth_image%kill()
        do i = 1,NPTCLS
            call actual(i)%kill()
            call expected(i)%kill()
        enddo
        do iframe = 1,NFR
            call frames(iframe)%kill()
        enddo
        deallocate(frames)
        call del_file(movie)
        call del_file(model_file)
        !$ call omp_set_num_threads(saved_omp)
      contains
        subroutine check_cached_sum(ffromto, n)
            integer, intent(in) :: ffromto(2), n
            integer, allocatable :: pinds(:)
            integer :: iptcl
            ! Zero motion: sum the raw fixture pixels, sharing only downstream normalization.
            call truth_image%set_rmat(sum(movie_pixels(:,:,:,ffromto(1):ffromto(2)), dim=4), .false.)
            pinds = [(iptcl, iptcl=1,n)]
            call mic_extractor%extract_particles_from_mic(truth_image, pinds, COORDS, expected, vmin, vmax, vmean, vsdev)
            call extractor%generate_cumulative_particles_stack(n, actual, ffromto, vmin, vmax, vmean, vsdev)
            do iptcl = 1,n
                call assert_real(0., maxval(abs(actual(iptcl)%get_rmat()-expected(iptcl)%get_rmat())), PIXEL_TOL, &
                    &'cached particle sums match the current extractor and requested interval')
            enddo
        end subroutine check_cached_sum
    end subroutine test_particle_frame_cache

    subroutine assert_dose_weighting_report( extractor, applied, model_dw )
        class(ptcl_extractor), intent(in) :: extractor
        logical,              intent(in) :: applied
        logical, optional,    intent(in) :: model_dw
        integer :: saved_log, scratch_log, ios, nmodel, napplied
        character(len=256) :: line
        character(len=32) :: label
        logical :: value
        saved_log = logfhandle
        open(newunit=scratch_log, status='scratch', action='readwrite')
        logfhandle = scratch_log
        call extractor%display()
        logfhandle = saved_log
        rewind(scratch_log)
        nmodel = 0
        napplied = 0
        do
            read(scratch_log,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, 'doseweighting_') == 0 ) cycle
            read(line,*,iostat=ios) label, value
            call assert_int(0, ios, 'dose-weighting diagnostics contain a logical value')
            if( ios /= 0 ) cycle
            select case(trim(label))
                case('doseweighting_applied')
                    napplied = napplied + 1
                    call assert_true(value .eqv. applied, 'display reports the weighting applied by the constructor')
                case('doseweighting_model')
                    nmodel = nmodel + 1
                    if( present(model_dw) )then
                        call assert_true(value .eqv. model_dw, 'display preserves the saved model weighting flag')
                    endif
            end select
        enddo
        close(scratch_log)
        call assert_int(1, napplied, 'display reports applied dose weighting exactly once')
        call assert_int(merge(1,0,present(model_dw)), nmodel, 'display reports model weighting only for an active model')
    end subroutine assert_dose_weighting_report

    subroutine assert_micrograph_state( extractor, box, neg )
        class(ptcl_extractor), intent(in) :: extractor
        integer,               intent(in) :: box
        logical,               intent(in) :: neg
        logical, allocatable :: mask(:,:,:)
        call assert_dose_weighting_report(extractor, .false.)
        call assert_true(extractor%does_exist(), 'micrograph initialization leaves an active extractor')
        call assert_false(extractor%from_mov(), 'micrograph initialization selects micrograph extraction')
        call assert_false(extractor%from_model(), 'micrograph initialization leaves no active motion model')
        call assert_int(0, extractor%get_nframes(), 'micrograph initialization retains no movie frames')
        call assert_false(extractor%has_frames(), 'micrograph initialization releases movie frames')
        call assert_false(extractor%has_isoshifts(), 'micrograph initialization releases movie shifts')
        call assert_false(extractor%has_weights(), 'micrograph initialization releases frame weights')
        call assert_false(extractor%has_hotpix_coords(), 'micrograph initialization releases defect coordinates')
        call assert_false(extractor%has_model_frameweights(), 'micrograph initialization releases model frame weights')
        call assert_false(extractor%has_particle_frames([1,1], 1), 'micrograph initialization releases the particle cache')
        call assert_int(box, extractor%get_box(), 'micrograph initialization uses the requested particle box')
        call assert_true(extractor%get_l_neg() .eqv. neg, 'micrograph initialization uses the requested particle contrast')
        call extractor%get_particle_mask(mask)
        call assert_true(allocated(mask), 'micrograph initialization creates the normalization mask')
        if( allocated(mask) )then
            call assert_true(all(shape(mask) == [box,box,1]), 'micrograph normalization mask matches the particle box')
            call assert_true(any(mask), 'micrograph normalization mask contains foreground pixels')
            call assert_true(any(.not.mask), 'micrograph normalization mask leaves background pixels')
        endif
    end subroutine assert_micrograph_state

end module simple_particle_extractor_tester
