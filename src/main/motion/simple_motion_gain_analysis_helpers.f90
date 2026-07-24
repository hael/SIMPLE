module simple_motion_gain_analysis_helpers
use simple_core_module_api
use simple_image, only: image
implicit none
private
#include "simple_local_flags.inc"

public :: read_movies_and_sum_frames

contains

    subroutine read_movies_and_sum_frames(movie_fnames, smpd, sum_img, n_movies, total_frames)
        class(string), intent(in)    :: movie_fnames(:)
        real,          intent(in)    :: smpd
        type(image),   intent(inout) :: sum_img
        integer,       intent(out)   :: n_movies
        integer,       intent(out)   :: total_frames
        type(image) :: frame
        integer     :: ldim_ref(3), ldim_cur(3), nframes, imov, iframe, ifoo
        logical     :: have_sum

        n_movies     = 0
        total_frames = 0
        have_sum     = .false.

        if( size(movie_fnames) < 1 )then
            THROW_HARD('No movies provided to read_movies_and_sum_frames')
        endif

        do imov=1,size(movie_fnames)
            if( .not. file_exists(movie_fnames(imov)) )then
                THROW_HARD('Movie file does not exist: '//movie_fnames(imov)%to_char())
            endif

            select case(fname2format(movie_fnames(imov)))
            case('K')
                THROW_HARD('EER movies are not supported by this test: '//movie_fnames(imov)%to_char())
            case DEFAULT
                call find_ldim_nptcls(movie_fnames(imov), ldim_cur, nframes)
            end select

            if( nframes < 1 )then
                THROW_HARD('No frames in movie stack: '//movie_fnames(imov)%to_char())
            endif

            ldim_cur(3) = 1
            if( .not. have_sum )then
                ldim_ref = ldim_cur
                call sum_img%new(ldim_ref, smpd, wthreads=.false.)
                call sum_img%zero()
                call frame%new(ldim_ref, smpd, wthreads=.false.)
                have_sum = .true.
            else if( ldim_cur(1) /= ldim_ref(1) .or. ldim_cur(2) /= ldim_ref(2) )then
                THROW_HARD('Movie dimensions differ from first movie dimensions: '//movie_fnames(imov)%to_char())
            endif

            do iframe=1,nframes
                call frame%read(movie_fnames(imov), iframe)
                call sum_img%add_workshare(frame)
                total_frames = total_frames + 1
            enddo
            n_movies = n_movies + 1
        enddo

        if( have_sum ) call frame%kill()
    end subroutine read_movies_and_sum_frames

end module simple_motion_gain_analysis_helpers
