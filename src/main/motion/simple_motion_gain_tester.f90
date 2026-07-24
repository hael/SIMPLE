!@descr: unit tests for motion gain helper and analyzer
module simple_motion_gain_tester
use simple_core_module_api
use simple_image,                        only: image
use simple_motion_gain_analysis,         only: gain_flip_analyzer
use simple_motion_gain_analysis_helpers, only: read_movies_and_sum_frames
use simple_test_utils
implicit none
private
#include "simple_local_flags.inc"

public :: run_all_motion_gain_tests

contains

    subroutine run_all_motion_gain_tests()
        write(*,'(A)') '**** running all motion gain tests ****'
        call test_read_movies_and_sum_frames_counts()
        call test_gain_flip_analyzer_batch_updates()
    end subroutine run_all_motion_gain_tests

    subroutine test_read_movies_and_sum_frames_counts()
        type(string) :: movie1, movie2
        type(string) :: movie_fnames(2)
        type(image)  :: sum_img
        integer      :: ldim(3), n_movies, total_frames
        real         :: smpd, val

        ldim = [4,4,1]
        smpd = 1.5
        movie1 = 'motion_gain_test_movie1.mrcs'
        movie2 = 'motion_gain_test_movie2.mrcs'

        call create_movie_stack(movie1, ldim, smpd, [1.0, 2.0, 3.0])
        call create_movie_stack(movie2, ldim, smpd, [4.0, 5.0])

        movie_fnames(1) = movie1
        movie_fnames(2) = movie2

        call read_movies_and_sum_frames(movie_fnames, smpd, sum_img, n_movies, total_frames)

        call assert_int(2, n_movies, 'read_movies_and_sum_frames should count movies')
        call assert_int(5, total_frames, 'read_movies_and_sum_frames should count frames')

        val = sum_img%get_rmat_at(1,1,1)
        call assert_real(15.0, val, 1.0e-6, 'sum image pixel should equal frame-value sum')

        call sum_img%kill()
        call del_file(movie1)
        call del_file(movie2)
    end subroutine test_read_movies_and_sum_frames_counts

    subroutine test_gain_flip_analyzer_batch_updates()
        type(gain_flip_analyzer) :: analyzer
        type(image)              :: part_sum
        type(string)             :: gainref
        integer                  :: ldim(3)
        real                     :: smpd
        logical                  :: ran_analysis

        ldim    = [4,4,1]
        smpd    = 1.5
        gainref = 'motion_gain_test_gainref.mrc'

        call create_gainref(gainref, ldim, smpd)
        call create_constant_image(part_sum, ldim, smpd, 2.0)

        call analyzer%new(gainref, smpd)

        call analyzer%analyze_if_due(part_sum, 10, 10, ran_analysis)
        call assert_true(ran_analysis, 'analysis should run when first threshold is reached')
        call assert_int(1, analyzer%analysis_count, 'analysis_count should increase on first due batch')

        call analyzer%analyze_if_due(part_sum, 5, 5, ran_analysis)
        call assert_true(ran_analysis, 'analysis should run on configured step interval')
        call assert_int(2, analyzer%analysis_count, 'analysis_count should increase on second due batch')

        call analyzer%analyze_if_due(part_sum, 1, 1, ran_analysis)
        call assert_false(ran_analysis, 'analysis should skip when movie step is not due')
        call assert_int(16, analyzer%total_frames, 'analyzer should accumulate total_frames from batches')
        call assert_int(16, analyzer%total_movies, 'analyzer should accumulate total_movies from batches')
        call assert_true(file_exists(string('average_all_movies.mrc')), 'analyzer should write averaged image output')
        call assert_true(file_exists(string('best_gainref_by_corr.mrc')), 'analyzer should write best-gain output')

        call analyzer%kill()
        call part_sum%kill()
        call del_file(gainref)
        call del_file(string('average_all_movies.mrc'))
        call del_file(string('best_gainref_by_corr.mrc'))
    end subroutine test_gain_flip_analyzer_batch_updates

    subroutine create_movie_stack(fname, ldim, smpd, frame_values)
        type(string), intent(in) :: fname
        integer,      intent(in) :: ldim(3)
        real,         intent(in) :: smpd
        real,         intent(in) :: frame_values(:)
        type(image)              :: frame
        real, allocatable        :: rmat(:,:,:)
        integer                  :: iframe

        call frame%new(ldim, smpd, wthreads=.false.)
        allocate(rmat(ldim(1),ldim(2),1))

        do iframe=1,size(frame_values)
            rmat = frame_values(iframe)
            call frame%set_rmat(rmat, .false.)
            call frame%write(fname, iframe, del_if_exists=(iframe == 1))
        enddo

        deallocate(rmat)
        call frame%kill()
    end subroutine create_movie_stack

    subroutine create_gainref(fname, ldim, smpd)
        type(string), intent(in) :: fname
        integer,      intent(in) :: ldim(3)
        real,         intent(in) :: smpd
        type(image)              :: img
        real, allocatable        :: rmat(:,:,:)
        integer                  :: i, j

        call img%new(ldim, smpd, wthreads=.false.)
        allocate(rmat(ldim(1),ldim(2),1))
        do j=1,ldim(2)
            do i=1,ldim(1)
                rmat(i,j,1) = real(i + 10 * j)
            enddo
        enddo
        call img%set_rmat(rmat, .false.)
        call img%write(fname, del_if_exists=.true.)

        deallocate(rmat)
        call img%kill()
    end subroutine create_gainref

    subroutine create_constant_image(img, ldim, smpd, val)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: ldim(3)
        real,        intent(in)    :: smpd, val
        real, allocatable          :: rmat(:,:,:)

        call img%new(ldim, smpd, wthreads=.false.)
        allocate(rmat(ldim(1),ldim(2),1))
        rmat = val
        call img%set_rmat(rmat, .false.)
        deallocate(rmat)
    end subroutine create_constant_image

end module simple_motion_gain_tester
