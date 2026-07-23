!@descr: utility functions for motion correction
module simple_motion_correct_utils
use simple_core_module_api
use simple_image,       only: image
use simple_eer_factory, only: eer_decoder
implicit none

public :: correct_gain, flip_gain, calc_eer_fraction, extract_outliers
private
#include "simple_local_flags.inc"

contains

    ! gain correction, calculate image sum and identify outliers
    subroutine correct_gain( frames_here, gainref_fname, gainimg, eerdecoder, frames_range )
        type(image),     allocatable, intent(inout) :: frames_here(:)
        class(string),                intent(in)    :: gainref_fname
        class(image),                 intent(inout) :: gainimg
        class(eer_decoder), optional, intent(in)    :: eerdecoder
        integer,            optional, intent(in)    :: frames_range(2)
        integer :: ldim_here(3), ldim_gain(3), iframe, ifoo, nframes_here, from, to
        write(logfhandle,'(a)') '>>> PERFORMING GAIN CORRECTION'
        nframes_here = size(frames_here)
        if( present(frames_range) )then
            from = frames_range(1)
            to   = frames_range(2)
        else
            from = 1
            to   = nframes_here
        endif
        if( to > nframes_here ) THROW_HARD('Invalid frame range 1; correct_gain')
        if( (from < 1) .or. (from > to) ) THROW_HARD('Invalid frame range 2; correct_gain')
        if( present(eerdecoder) )then
            call eerdecoder%prep_gainref(gainref_fname, gainimg)
        else
            if( fname2format(gainref_fname)=='L' )then
                THROW_HARD('''.gain'' files only for use with EER movies! correct_gain')
            endif
            ldim_here = frames_here(from)%get_ldim()
            call find_ldim_nptcls(gainref_fname,ldim_gain,ifoo)
            if( ldim_gain(1).ne.ldim_here(1) .or. ldim_gain(2).ne.ldim_here(2) )then
                THROW_HARD('Inconsistent dimensions between movie frames & gain reference! correct_gain')
            endif
            call gainimg%new(ldim_gain, frames_here(from)%get_smpd())
            call gainimg%read(gainref_fname)
        endif
        !$omp parallel do schedule(static) default(shared) private(iframe) proc_bind(close)
        do iframe = from,to
            call frames_here(iframe)%mul(gainimg)
        enddo
        !$omp end parallel do
    end subroutine correct_gain

    ! flips, writes gain reference & update command line
    subroutine flip_gain( cline, gainref_fname, mode )
        use simple_cmdline, only: cmdline
        class(cmdline),   intent(inout) :: cline
        class(string),    intent(inout) :: gainref_fname
        character(len=*),  intent(in)   :: mode
        type(image)  :: gain
        type(string) :: new_fname, ext
        real    :: smpd
        integer :: ldim(3), n
        if( .not.cline%defined('gainref') ) return
        if( .not.cline%defined('flipgain') ) return
        select case(trim(uppercase(mode)))
        case('NO')
            return
        case('X','Y','XY','YX')
            ! supported
        case DEFAULT
            THROW_HARD('UNSUPPORTED GAIN REFERENCE FLIPPING MODE: '//trim(mode))
        end select
        if( .not.file_exists(gainref_fname) )then
            THROW_HARD('Could not find gain reference: '//gainref_fname%to_char())
        endif
        call find_ldim_nptcls(gainref_fname, ldim, n)
        smpd = find_img_smpd(gainref_fname)
        ldim(3) = 1
        call gain%new(ldim, smpd, wthreads=.false.)
        call gain%read(gainref_fname)
        call gain%flip(mode)
        gainref_fname = basename(gainref_fname)
        ext           = fname2ext(gainref_fname)
        gainref_fname = get_fbody(gainref_fname, ext, separator=.true.)
        gainref_fname = './'//gainref_fname%to_char()//'_flip'//trim(uppercase(mode))//'.mrc'
        call gain%write(gainref_fname)
        new_fname     = simple_abspath(gainref_fname)
        gainref_fname = new_fname%to_char()
        call cline%set('gainref', gainref_fname%to_char())
        call gain%kill
        call new_fname%kill
        call ext%kill
    end subroutine flip_gain

    !>  Utility to calculate the number fractions, # eer frames per fraction and adjusted total_dose
    !   while minimizing the number of leftover frames given total dose, total # of eer frames & a dose target
    subroutine calc_eer_fraction(n_eer_frames, fraction_dose_target, tot_dose, nfractions, eerfraction)
        integer, intent(in)    :: n_eer_frames
        real,    intent(in)    :: fraction_dose_target
        real,    intent(inout) :: tot_dose
        integer, intent(out)   :: nfractions, eerfraction
        real    :: dose_per_eer_frame
        integer :: ffloor, fceil, nffloor, nfceil
        ffloor  = floor(real(n_eer_frames) / (tot_dose / fraction_dose_target))
        fceil   = ffloor+1
        nffloor = floor(real(n_eer_frames) / real(ffloor))
        nfceil  = floor(real(n_eer_frames) / real(fceil))
        if( nffloor*ffloor > nfceil*fceil )then
            nfractions  = nffloor
            eerfraction = ffloor
        else
            nfractions  = nfceil
            eerfraction = fceil
        endif
        ! adjusting total dose
        dose_per_eer_frame = tot_dose / real(n_eer_frames)
        tot_dose           = dose_per_eer_frame * real(eerfraction) * real(nfractions)
    end subroutine calc_eer_fraction

    subroutine extract_outliers(movie_fname, smpd, low_outliers_out, high_outliers_out)
        class(string), intent(in)  :: movie_fname
        real,          intent(in)  :: smpd
        logical, allocatable, optional, intent(out) :: low_outliers_out(:,:), high_outliers_out(:,:)
        type(image)                :: movie_sum, frame
        type(string)               :: sum_fname, ext
        real,          allocatable :: sum_rmat(:,:,:)
        integer                    :: ldim(3), nframes, iframe
        integer                    :: n_low, n_high
        real, parameter            :: NSIGMAS = 4.
        real, parameter            :: LOW_ZERO_EPS = 1.0e-6
        real                       :: ave, sdev, var, lthresh, uthresh
        logical                    :: err
        logical, allocatable       :: low_outliers(:,:), high_outliers(:,:)
        select case(fname2format(movie_fname))
        case('K')
            write(logfhandle,'(a)') '>>> EER MOVIES NOT YET SUPPORTED FOR DEAD PIXEL EXTRACTION'
            return
        case DEFAULT
            call find_ldim_nptcls(movie_fname, ldim, nframes)
        end select
        if( nframes < 1 )then
            THROW_HARD('No frames in movie stack: '//movie_fname%to_char())
        endif
        ldim(3) = 1

        write(logfhandle,'(a)') '>>> SUMMING MOVIE FRAMES FOR DEAD PIXEL EXTRACTION'
        call movie_sum%new(ldim, smpd, wthreads=.false.)
        call frame%new(ldim, smpd, wthreads=.false.)

        call movie_sum%read(movie_fname, 1)
        do iframe=2,nframes
            call frame%read(movie_fname, iframe)
            call movie_sum%add_workshare(frame)
        enddo

        sum_rmat = movie_sum%get_rmat()
        call moment(sum_rmat(:,:,1), ave, sdev, var, err)
        if( err )then
            THROW_HARD('Failed to compute statistics for movie sum: '//movie_fname%to_char())
        endif
        if( sdev < TINY )then
            write(logfhandle,'(a)') '>>> Movie sum has near-zero variance; no outliers detected'
            n_low  = 0
            n_high = 0
            allocate(low_outliers(ldim(1), ldim(2)), source=.false.)
            allocate(high_outliers(ldim(1), ldim(2)), source=.false.)
        else
            lthresh = ave - NSIGMAS * sdev
            uthresh = ave + NSIGMAS * sdev
            allocate(low_outliers(ldim(1), ldim(2)), source=.false.)
            allocate(high_outliers(ldim(1), ldim(2)), source=.false.)
            !$omp workshare
            where(abs(sum_rmat(:,:,1)) <= LOW_ZERO_EPS)
                low_outliers = .true.
            elsewhere
                low_outliers = .false.
            end where
            where(sum_rmat(:,:,1) > uthresh)
                high_outliers = .true.
            elsewhere
                high_outliers = .false.
            end where
            !$omp end workshare
            n_low  = count(low_outliers)
            n_high = count(high_outliers)
            write(logfhandle,'(a,1x,i0,1x,a,1x,i0,1x,a,f12.6,1x,a,f12.3)') &
                '>>> LOW/HIGH OUTLIERS:', n_low, '/', n_high, ' LOW|x|<=', LOW_ZERO_EPS, ' HIGH>x:', uthresh
        endif

        if( present(low_outliers_out) ) then
            if( allocated(low_outliers_out) ) deallocate(low_outliers_out)
            allocate(low_outliers_out(ldim(1), ldim(2)), source=low_outliers)
        endif
        if( present(high_outliers_out) ) then
            if( allocated(high_outliers_out) ) deallocate(high_outliers_out)
            allocate(high_outliers_out(ldim(1), ldim(2)), source=high_outliers)
        endif

        
        
        ! ext       = fname2ext(movie_fname)
        ! sum_fname = get_fbody(basename(movie_fname), ext, separator=.true.)
        ! sum_fname = sum_fname%to_char()//'_sum.mrc'
        ! call movie_sum%write(sum_fname)
        ! write(logfhandle,'(a)') '>>> Wrote movie-frame sum to: '//sum_fname%to_char()

        if( allocated(low_outliers) )  deallocate(low_outliers)
        if( allocated(high_outliers) ) deallocate(high_outliers)
        if( allocated(sum_rmat) ) deallocate(sum_rmat)
        call movie_sum%kill()
        call frame%kill()
        call sum_fname%kill()
        call ext%kill()
    end subroutine extract_outliers

end module simple_motion_correct_utils
