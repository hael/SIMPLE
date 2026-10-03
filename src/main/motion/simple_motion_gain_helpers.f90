!@descr: gain-reference helpers: summed movie frames, normalised inverse average intensity and gain-reference previews
module simple_motion_gain_helpers
use simple_core_module_api
use simple_image,           only: image
use simple_eer_factory,     only: eer_decoder

implicit none
private
#include "simple_local_flags.inc"

integer, parameter  :: EER_THUMB_UPSAMPLING = 1

public :: read_movies_and_sum_frames
public :: add_movies_to_gain_sum
public :: normalized_inverse_average_intensity
public :: write_gain_from_sum
public :: gainref_to_jpg

contains

    subroutine read_movies_and_sum_frames(movie_fnames, smpd, sum_img, n_movies, total_frames)
        class(string), intent(in)    :: movie_fnames(:)
        real,          intent(in)    :: smpd
        type(image),   intent(inout) :: sum_img
        integer,       intent(out)   :: n_movies
        integer,       intent(out)   :: total_frames
        type(image)                  :: frame, eer_frame(1)
        type(eer_decoder)            :: eer
        integer                      :: ldim_ref(3), ldim_cur(3), nframes, imov, iframe
        logical                      :: have_sum

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
                ! Decode raw EER events at native 4K sampling before thumbnail downscaling.
                call eer%new(movie_fnames(imov), smpd, EER_THUMB_UPSAMPLING)
                nframes = eer%get_nframes()
                ldim_cur = eer%get_ldim()
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
                call eer_frame(1)%new(ldim_ref, smpd, wthreads=.false.)
                have_sum = .true.
            else if( ldim_cur(1) /= ldim_ref(1) .or. ldim_cur(2) /= ldim_ref(2) )then
                THROW_HARD('Movie dimensions differ from first movie dimensions: '//movie_fnames(imov)%to_char())
            endif

            select case(fname2format(movie_fnames(imov)))
            case('K')
                call eer%decode(eer_frame, nframes)
                call sum_img%add_workshare(eer_frame(1))
                total_frames = total_frames + 1
                call eer%kill()
            case DEFAULT 
                do iframe=1,nframes
                    call frame%read(movie_fnames(imov), iframe)
                    call sum_img%add_workshare(frame)
                    total_frames = total_frames + 1
                enddo
            end select
            n_movies = n_movies + 1
        enddo

        if( have_sum ) then
            call frame%kill()
            call eer_frame(1)%kill()
        endif
    end subroutine read_movies_and_sum_frames

    !> Sums the frames of one batch of movies into the running sum gain_sum, which the
    !> first call creates; nmovies and nframes accumulate across calls.
    subroutine add_movies_to_gain_sum(movie_fnames, smpd, gain_sum, nmovies, nframes)
        class(string), intent(in)    :: movie_fnames(:)
        real,          intent(in)    :: smpd
        type(image),   intent(inout) :: gain_sum
        integer,       intent(inout) :: nmovies, nframes
        type(image) :: batch_sum
        integer     :: batch_movies, batch_frames, ldim_sum(3), ldim_batch(3)

        call read_movies_and_sum_frames(movie_fnames, smpd, batch_sum, batch_movies, batch_frames)
        if( gain_sum%exists() )then
            ldim_sum   = gain_sum%get_ldim()
            ldim_batch = batch_sum%get_ldim()
            if( any(ldim_sum /= ldim_batch) )then
                THROW_HARD('Movie batch dimensions differ from the running gain sum; add_movies_to_gain_sum')
            endif
            call gain_sum%add_workshare(batch_sum)
        else
            call gain_sum%copy(batch_sum)
        endif
        nmovies = nmovies + batch_movies
        nframes = nframes + batch_frames
        call batch_sum%kill()
    end subroutine add_movies_to_gain_sum

    subroutine normalized_inverse_average_intensity(sum_img, nframes, inv_avg_img, avg_value)
        class(image), intent(in)    :: sum_img
        integer,      intent(in)    :: nframes
        type(image),  intent(out)   :: inv_avg_img
        real,         intent(out)   :: avg_value
        real, parameter :: INV_AVG_EPS = 1.0e-12
        integer :: ldim(3)
        real    :: npix
        real, allocatable :: rmat(:,:,:)

        if( nframes < 1 )then
            THROW_HARD('nframes must be >= 1; normalized_inverse_average_intensity')
        endif

        ldim = sum_img%get_ldim()
        if( any(ldim < 1) )then
            THROW_HARD('sum_img has invalid dimensions; normalized_inverse_average_intensity')
        endif

        npix = real(ldim(1) * ldim(2) * ldim(3))
        rmat = sum_img%get_rmat()
        avg_value = sum(rmat) / (npix * real(nframes))

        if( abs(avg_value) <= INV_AVG_EPS )then
            THROW_HARD('Average intensity is near zero; cannot compute stable inverse scaling')
        endif

        rmat = rmat / real(nframes)
        where( abs(rmat) > INV_AVG_EPS )
            rmat = avg_value / rmat
        elsewhere
            rmat = 0.
        end where
      
        call inv_avg_img%new(ldim, sum_img%get_smpd(), wthreads=.false.)
        call inv_avg_img%set_rmat(rmat, .false.)
        
        deallocate(rmat)

    end subroutine normalized_inverse_average_intensity

    !> Writes the gain reference estimated from gain_sum, the sum of nframes movie frames,
    !> to fname: the normalised inverse average intensity.
    subroutine write_gain_from_sum(gain_sum, nframes, fname)
        class(image),  intent(in) :: gain_sum
        integer,       intent(in) :: nframes
        class(string), intent(in) :: fname
        type(image) :: gain_img
        real        :: avg_value

        call normalized_inverse_average_intensity(gain_sum, nframes, gain_img, avg_value)
        call gain_img%write(fname, del_if_exists=.true.)
        call gain_img%kill()
    end subroutine write_gain_from_sum

    !> Reads a gain reference from disk and writes a normalized jpeg preview
    !> resized to the standard GUI micrograph thumbnail size.
    subroutine gainref_to_jpg(gainref_fname, jpg_fname, quality)
        class(string),     intent(in) :: gainref_fname
        class(string),     intent(in) :: jpg_fname
        integer, optional, intent(in) :: quality
        type(image)                   :: gain_img, thumb_img
        integer                       :: ldim(3), ldim_thumb(3), ifoo
        real                          :: smpd, scale

        if( .not. file_exists(gainref_fname) )then
            THROW_HARD('Gain reference file does not exist: '//gainref_fname%to_char()//'; gainref_to_jpg')
        endif

        call find_ldim_nptcls(gainref_fname, ldim, ifoo)
        ldim(3) = 1
        smpd    = find_img_smpd(gainref_fname)
        call gain_img%new(ldim, smpd, wthreads=.false.)
        call gain_img%read(gainref_fname)

        ! .gain/EER-TIFF references are stored bottom-up relative to the .mrc convention
        if( fname2format(gainref_fname) == 'J' .or. fname2format(gainref_fname) == 'L' )then
            call gain_img%flip('Y')
        endif

        ! resize to the standard GUI micrograph thumbnail size
        scale           = real(GUI_PSPECSZ)/real(maxval(ldim(1:2)))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_thumb(3)   = 1
        call thumb_img%new(ldim_thumb, smpd)
        call gain_img%fft()
        call gain_img%clip(thumb_img)
        call thumb_img%ifft()
        call gain_img%kill()

        call thumb_img%write_jpg(jpg_fname, norm=.true., quality=quality)
        call thumb_img%kill()
    end subroutine gainref_to_jpg

end module simple_motion_gain_helpers