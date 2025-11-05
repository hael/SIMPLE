module simple_motion_correct_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,       only: image
use simple_eer_factory, only: eer_decoder
implicit none

public :: correct_gain, flip_gain, apply_dose_weighing, micrograph_interp, calc_eer_fraction
private
#include "simple_local_flags.inc"

contains

    ! gain correction, calculate image sum and identify outliers
    subroutine correct_gain( frames_here, gainref_fname, gainimg, eerdecoder, frames_range )
        type(image),     allocatable, intent(inout) :: frames_here(:)
        character(len=*),             intent(in)    :: gainref_fname
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
        character(len=*), intent(inout) :: gainref_fname, mode
        type(image)                   :: gain
        character(len=:), allocatable :: new_fname
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
            THROW_HARD('Could not find gain reference: '//trim(gainref_fname))
        endif
        call find_ldim_nptcls(gainref_fname, ldim, n, smpd=smpd)
        ldim(3) = 1
        call gain%new(ldim, smpd, wthreads=.false.)
        call gain%read(gainref_fname)
        call gain%flip(mode)
        gainref_fname = basename(gainref_fname)
        gainref_fname = get_fbody(gainref_fname, trim(fname2ext(gainref_fname)), separator=.true.)
        gainref_fname = './'//trim(gainref_fname)//'_flip'//trim(uppercase(mode))//'.mrc'
        call gain%write(gainref_fname)
        new_fname     = simple_abspath(gainref_fname)
        gainref_fname = trim(new_fname)
        call cline%set('gainref', gainref_fname)
        call gain%kill
        deallocate(new_fname)
    end subroutine flip_gain

    ! Following Grant & Grigorieff; eLife 2015;4:e06980
    ! Frames assumed in fourier space
    subroutine apply_dose_weighing( nframes, frames, frange, total_dose, voltage )
        integer,      intent(in)    :: nframes
        class(image), intent(inout) :: frames(nframes)
        integer,      intent(in)    :: frange(2)
        real,         intent(in)    :: total_dose   ! in e-/A2
        real,         intent(in)    :: voltage      ! in kV
        real, parameter   :: A=0.245, B=-1.665, C=2.81
        real    :: qs(frange(1):frange(2)), acc_doses(nframes)
        real    :: spaFreqk, dose_per_frame, twoNe, smpd, spafreq, limhsq,limksq
        integer :: nrflims(3,2), ldim(3), hphys,kphys, i, h,k
        if( .not.frames(1)%is_ft() ) THROW_HARD('Frames should be in in the Fourier domain')
        nrflims = frames(frange(1))%loop_lims(2)
        smpd    = frames(frange(1))%get_smpd()
        qs      = 0.
        ldim    = frames(frange(1))%get_ldim()
        limhsq  = (real(ldim(1))*smpd)**2.
        limksq  = (real(ldim(2))*smpd)**2.
        ! Accumulated frame dosess
        dose_per_frame = total_dose / real(nframes)
        acc_doses      = (/(real(i)*dose_per_frame, i=1,nframes)/)
        if( is_equal(voltage,200.) )then
            acc_doses = acc_doses / 0.8
        else if( is_equal(voltage,100.) )then
            acc_doses = acc_doses / 0.64
        endif
        ! dose normalization
        !$omp parallel do private(h,k,spafreq,spafreqk,twone,kphys,hphys,i,qs)&
        !$omp default(shared) schedule(static) proc_bind(close)
        do k = nrflims(2,1),nrflims(2,2)
            kphys    = k + 1 + merge(ldim(2),0,k<0)
            spaFreqk = real(k*k)/limksq
            do h = nrflims(1,1),nrflims(1,2)
                hphys   = h + 1
                spaFreq = sqrt( real(h*h)/limhsq + spaFreqk )
                twoNe   = 2.*(A*spaFreq**B + C)
                qs = exp(-acc_doses(frange(1):frange(2))/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                do i = frange(1),frange(2)
                    call frames(i)%mul_cmat_at([hphys,kphys,1], qs(i))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_dose_weighing

    !>  Per frame real space polynomial interpolation
    subroutine micrograph_interp(interp_fixed_frame, fixed_frame, nframes, frames,&
        &weights, poly_coeffs_dim, poly_coeffs, frame_output)
        integer,               intent(in)    :: interp_fixed_frame, fixed_frame, nframes, poly_coeffs_dim
        type(image),           intent(inout) :: frames(nframes)
        real,                  intent(in)    :: weights(nframes)
        real(dp),              intent(in)    :: poly_coeffs(poly_coeffs_dim,2)
        type(image),           intent(inout) :: frame_output
        real, pointer :: rmatin(:,:,:), rmatout(:,:,:)
        real(dp)      :: t,ti, dt,dt2,dt3, x,x2,y,y2,xy, A1,A2, B1x,B1x2,B1xy,B2x,B2x2,B2xy
        integer       :: ldim(3), i, j, iframe
        real          :: w, pixx,pixy
        ldim = frames(1)%get_ldim()
        call frame_output%zero_and_unflag_ft
        call frame_output%get_rmat_ptr(rmatout)
        ti = real(interp_fixed_frame-fixed_frame, dp)
        do iframe = 1,nframes
            call frames(iframe)%get_rmat_ptr(rmatin)
            w = weights(iframe)
            t = real(iframe-fixed_frame, dp)
            dt  = ti-t
            dt2 = ti*ti - t*t
            dt3 = ti*ti*ti - t*t*t
            B1x  = sum(poly_coeffs(4:6,1)   * [dt,dt2,dt3])
            B1x2 = sum(poly_coeffs(7:9,1)   * [dt,dt2,dt3])
            B1xy = sum(poly_coeffs(16:18,1) * [dt,dt2,dt3])
            B2x  = sum(poly_coeffs(4:6,2)   * [dt,dt2,dt3])
            B2x2 = sum(poly_coeffs(7:9,2)   * [dt,dt2,dt3])
            B2xy = sum(poly_coeffs(16:18,2) * [dt,dt2,dt3])
            !$omp parallel do default(shared) private(i,j,x,x2,y,y2,xy,A1,A2,pixx,pixy)&
            !$omp proc_bind(close) schedule(static)
            do j = 1, ldim(2)
                y  = real(j-1,dp) / real(ldim(2)-1,dp) - 0.5d0
                y2 = y*y
                A1 =           sum(poly_coeffs(1:3,1)   * [dt,dt2,dt3])
                A1 = A1 + y  * sum(poly_coeffs(10:12,1) * [dt,dt2,dt3])
                A1 = A1 + y2 * sum(poly_coeffs(13:15,1) * [dt,dt2,dt3])
                A2 =           sum(poly_coeffs(1:3,2)   * [dt,dt2,dt3])
                A2 = A2 + y  * sum(poly_coeffs(10:12,2) * [dt,dt2,dt3])
                A2 = A2 + y2 * sum(poly_coeffs(13:15,2) * [dt,dt2,dt3])
                do i = 1, ldim(1)
                    x  = real(i-1,dp) / real(ldim(1)-1,dp) - 0.5d0
                    x2 = x*x
                    xy = x*y
                    pixx = real(i) + real(A1 + B1x*x + B1x2*x2 + B1xy*xy)
                    pixy = real(j) + real(A2 + B2x*x + B2x2*x2 + B2xy*xy)
                    rmatout(i,j,1) = rmatout(i,j,1) + w*interp_bilin(pixx,pixy)
                end do
            end do
            !$omp end parallel do
        enddo
        nullify(rmatin,rmatout)
        contains

        pure real function interp_bilin( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: t, u
            logical  :: outside
            outside = .false.
            x1_h = floor(xval)
            x2_h = x1_h + 1
            if( x1_h<1 .or. x2_h<1 )then
                x1_h    = 1
                outside = .true.
            endif
            if( x1_h>ldim(1) .or. x2_h>ldim(1) )then
                x1_h    = ldim(1)
                outside = .true.
            endif
            y1_h = floor(yval)
            y2_h = y1_h + 1
            if( y1_h<1 .or. y2_h<1 )then
                y1_h    = 1
                outside = .true.
            endif
            if( y1_h>ldim(2) .or. y2_h>ldim(2) )then
                y1_h    = ldim(2)
                outside = .true.
            endif
            if( outside )then
                interp_bilin = rmatin(x1_h, y1_h, 1)
                return
            endif
            t  = xval - real(x1_h)
            u  = yval - real(y1_h)
            interp_bilin =  (1. - t) * (1. - u) * rmatin(x1_h, y1_h, 1) + &
                &      t  * (1. - u) * rmatin(x2_h, y1_h, 1) + &
                &      t  *       u  * rmatin(x2_h, y2_h, 1) + &
                &(1. - t) *       u  * rmatin(x1_h, y2_h, 1)
        end function interp_bilin

        pure real function interp_nn( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x, y
            x = min(max(nint(xval),1),ldim(1))
            y = min(max(nint(yval),1),ldim(2))
            interp_nn = rmatin(x,y,1)
        end function interp_nn

    end subroutine micrograph_interp

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
