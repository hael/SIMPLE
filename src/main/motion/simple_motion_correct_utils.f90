!@descr: utility functions for motion correction
module simple_motion_correct_utils
use simple_core_module_api
use simple_image,       only: image
use simple_eer_factory, only: eer_decoder
implicit none

public :: correct_gain, flip_gain, calc_eer_fraction, pix2polycoords, patch_poly, apply_patch_poly
public :: eer_scale_movie_convention
private
#include "simple_local_flags.inc"

contains

    !> Return decoded EER sampling and dimensions from physical input values.
    ! Mode 1 leaves both unchanged; mode 2 doubles X/Y and halves sampling.
    subroutine eer_scale_movie_convention( smpd, ldim, eer_upsampling, smpd_out, ldim_out )
        real,    intent(in)  :: smpd
        integer, intent(in)  :: ldim(2), eer_upsampling
        real,    intent(out) :: smpd_out
        integer, intent(out) :: ldim_out(2)
        if( any(ldim < 1) .or. smpd <= 0. ) THROW_HARD('Invalid movie dimensions or sampling!')
        select case(eer_upsampling)
        case(1)
            ldim_out = ldim
            smpd_out = smpd
        case(2)
            ldim_out = 2 * ldim
            smpd_out = smpd / 2.
        case DEFAULT
            THROW_HARD('Unsupported EER up-sampling: '//int2str(eer_upsampling))
        end select
    end subroutine eer_scale_movie_convention

    !> Converts a pixel coordinate to the normalized stage coordinate
    pure real(dp) function pix2polycoords( pixel_coord, image_dim )
        real(dp), intent(in) :: pixel_coord
        integer,  intent(in) :: image_dim
        pix2polycoords = (pixel_coord-1.d0) / real(image_dim-1,dp) - 0.5d0
    end function pix2polycoords

    !> Polynomial basis for the 18-parameter local-motion model
    pure function patch_poly( p, n ) result( res )
        real(dp), intent(in) :: p(:)
        integer,  intent(in) :: n
        real(dp) :: res(n)
        real(dp) :: x, y, t
        x = p(1)
        y = p(2)
        t = p(3)
        res(    1) = t
        res(    2) = t**2
        res(    3) = t**3
        res( 4: 6) = x * res( 1: 3)  ! x   * {t,t^2,t^3}
        res( 7: 9) = x * res( 4: 6)  ! x^2 * {t,t^2,t^3}
        res(10:12) = y * res( 1: 3)  ! y   * {t,t^2,t^3}
        res(13:15) = y * res(10:12)  ! y^2 * {t,t^2,t^3}
        res(16:18) = y * res( 4: 6)  ! x*y * {t,t^2,t^3}
    end function patch_poly

    !> Evaluates the 18-parameter local-motion polynomial
    pure function apply_patch_poly( c, x, y, t ) result( res_sp )
        real(dp), intent(in) :: c(18), x, y, t
        real(sp) :: res_sp
        real(dp) :: res
        real(dp) :: x2, y2, xy, t2, t3
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        res =       c( 1) * t      + c( 2) * t2      + c( 3) * t3
        res = res + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
        res = res + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
        res = res + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
        res = res + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
        res = res + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
        res_sp = real(res)
    end function apply_patch_poly

    ! gain correction, calculate image sum and identify outliers
    subroutine correct_gain( frames_here, gainref_fname, gainimg, eerdecoder, frames_range, flipgain )
        type(image),     allocatable, intent(inout) :: frames_here(:)
        class(string),                intent(in)    :: gainref_fname
        class(image),                 intent(inout) :: gainimg
        class(eer_decoder), optional, intent(in)    :: eerdecoder
        integer,            optional, intent(in)    :: frames_range(2)
        character(len=*),   optional, intent(in)    :: flipgain
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
        ! Orient the loaded/prepared gain before applying it to any frame. Keep
        ! gainimg in this orientation so EER defect coordinates use the same axes.
        if( present(flipgain) )then
            write(logfhandle,'(a,a)') '>>> GAIN FLIPPING: ', trim(flipgain)
            if( trim(uppercase(flipgain)) /= 'NO' ) call gainimg%flip(flipgain)
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

end module simple_motion_correct_utils
