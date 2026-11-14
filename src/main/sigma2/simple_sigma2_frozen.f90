!@descr: frozen log-Savitzky-Golay bootstrap spectra for 2D matching
module simple_sigma2_frozen
use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_global, SIGMA2_PROV_PSPEC, SIGMA2_STATE_COMMITTED
use simple_syslib, only: simple_atomic_replace
implicit none
private
public :: prepare_frozen_sigma2, read_frozen_sigma2, frozen_sigma2_band_stats

character(len=16), parameter :: MAGIC = 'SIMPLE_S2FRZ001'
integer(int64), parameter :: HEADER_BYTES = 56_int64

contains

    !> Unweighted shell means of log(measured/anchor), the low-pass edge and counts above the anchor.
    subroutine frozen_sigma2_band_stats( committed, band, mean_log, edge_log, nabove, status, message )
        character(len=*), intent(in) :: committed
        integer, intent(in) :: band(2)
        real(real64), intent(out) :: mean_log(2), edge_log(2)
        integer, intent(out) :: nabove(2), status
        character(len=*), intent(out) :: message
        real(real32), allocatable :: measured(:,:), anchor(:,:)
        real(real64), allocatable :: log_ratio(:)
        integer :: half
        mean_log = 0.d0; edge_log = 0.d0; nabove = 0
        call read_frozen_sigma2(committed, anchor, status, message)
        if( status /= 0 ) return
        call sigma2_state_read_global(committed, measured, status, message)
        if( status /= 0 ) return
        status = 1
        message = 'invalid matching band for frozen sigma2 readout'
        if( any(shape(measured) /= shape(anchor)) ) return
        if( band(1) < lbound(anchor,1) .or. band(2) > ubound(anchor,1) .or. band(2) < band(1) ) return
        message = 'measured sigma2 must be finite and positive for the anchor readout'
        if( .not. all(ieee_is_finite(measured(band(1):band(2),:))) ) return
        if( any(measured(band(1):band(2),:) <= 0.0_real32) ) return
        do half = 1, 2
            log_ratio = log(real(measured(band(1):band(2),half),real64)) - &
                &log(real(anchor(band(1):band(2),half),real64))
            mean_log(half) = sum(log_ratio)/real(size(log_ratio),real64)
            edge_log(half) = log_ratio(size(log_ratio))
            nabove(half) = count(log_ratio > 0.d0)
        enddo
        status = 0
        message = ''
    end subroutine frozen_sigma2_band_stats

    !> The owner saves a bootstrap once; later stages and restarts verify the same companion.
    subroutine prepare_frozen_sigma2( committed, status, message )
        character(len=*), intent(in) :: committed
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        real(real32), allocatable :: spectra(:,:), smoothed(:,:)
        character(len=:), allocatable :: fname, candidate
        integer :: unit, half
        logical :: exists
        call sigma2_state_read_header(committed, header, status, message)
        if( status /= 0 ) return
        fname = trim(committed)//'.frozen'
        inquire(file=fname, exist=exists)
        if( exists .and. header%provenance /= SIGMA2_PROV_PSPEC )then
            call read_frozen_sigma2(committed, spectra, status, message)
            return
        endif
        if( header%provenance /= SIGMA2_PROV_PSPEC .or. header%state /= SIGMA2_STATE_COMMITTED )then
            status = 1
            message = 'sigma2_frozen needs its saved bootstrap; start a fresh run to create it'
            return
        endif
        call sigma2_state_read_global(committed, spectra, status, message)
        if( status /= 0 ) return
        allocate(smoothed(size(spectra,1),2))
        do half = 1, 2
            call smooth_log_spectrum(spectra(:,half), smoothed(:,half), status, message)
            if( status /= 0 ) return
        enddo
        candidate = fname//'.next'
        open(newunit=unit, file=candidate, access='stream', form='unformatted', &
            &status='replace', action='write', iostat=status)
        if( status /= 0 )then
            message = 'cannot create frozen sigma2 companion'; return
        endif
        write(unit, iostat=status) MAGIC, header%kfrom, header%kto, header%box, header%nptcls, &
            &header%smpd, header%layout_digest, header%generation, smoothed
        if( status /= 0 )then
            close(unit, status='delete')
            message = 'cannot write frozen sigma2 companion'; return
        endif
        close(unit, iostat=status)
        if( status /= 0 )then
            message = 'cannot close frozen sigma2 companion'; return
        endif
        call simple_atomic_replace(candidate, fname, status)
        message = ''
        if( status /= 0 ) message = 'cannot publish frozen sigma2 companion'
    end subroutine prepare_frozen_sigma2

    !> Read-only consumer: check provenance metadata and length before allocating the native-grid curves.
    subroutine read_frozen_sigma2( committed, spectra, status, message )
        character(len=*), intent(in) :: committed
        real(real32), allocatable, intent(out) :: spectra(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        character(len=16) :: file_magic
        integer(int32) :: kfrom, kto, box, nptcls
        integer(int64) :: digest, generation, nbytes
        real(real64) :: smpd
        integer :: unit
        call sigma2_state_read_header(committed, header, status, message)
        if( status /= 0 ) return
        open(newunit=unit, file=trim(committed)//'.frozen', access='stream', form='unformatted', &
            &status='old', action='read', iostat=status)
        if( status /= 0 )then
            message = 'missing frozen sigma2 bootstrap companion'; return
        endif
        read(unit, iostat=status) file_magic, kfrom, kto, box, nptcls, smpd, digest, generation
        if( status /= 0 )then
            close(unit)
            message = 'truncated frozen sigma2 header'; return
        endif
        if( .not. ieee_is_finite(smpd) )then
            close(unit)
            status = 1
            message = 'nonfinite sampling in frozen sigma2 companion'; return
        endif
        if( file_magic /= MAGIC .or. kfrom /= header%kfrom .or. kto /= header%kto .or. &
            &box /= header%box .or. nptcls /= header%nptcls .or. digest /= header%layout_digest .or. &
            &smpd /= header%smpd .or. generation < 1_int64 .or. generation > header%generation .or. &
            &header%state /= SIGMA2_STATE_COMMITTED )then
            close(unit)
            status = 1
            message = 'frozen sigma2 companion does not match the canonical grid or particle layout'; return
        endif
        inquire(unit=unit, size=nbytes, iostat=status)
        if( status == 0 )then
            if( nbytes /= HEADER_BYTES + 8_int64*int(kto-kfrom+1,int64) ) status = 1
        endif
        if( status /= 0 )then
            close(unit)
            message = 'frozen sigma2 companion has the wrong size'; return
        endif
        allocate(spectra(kfrom:kto,2))
        read(unit, iostat=status) spectra
        close(unit)
        if( status == 0 )then
            if( .not. all(ieee_is_finite(spectra)) ) status = 1
        endif
        if( status == 0 )then
            if( any(spectra <= 0.0_real32) ) status = 1
        endif
        message = ''
        if( status /= 0 ) message = 'frozen sigma2 spectrum must be finite and positive'
    end subroutine read_frozen_sigma2

    !> Quadratic least squares on seven consecutive shells, evaluated at each shell including the endpoints.
    subroutine smooth_log_spectrum( input, output, status, message )
        real(real32), intent(in) :: input(:)
        real(real32), intent(out) :: output(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(real64), parameter :: X(7) = [-3.d0,-2.d0,-1.d0,0.d0,1.d0,2.d0,3.d0]
        real(real64) :: y(7), z, value
        integer :: k, first
        status = 1
        message = 'frozen sigma2 requires at least seven finite positive bootstrap shells'
        if( size(input) < 7 .or. size(output) /= size(input) ) return
        if( .not. all(ieee_is_finite(input)) ) return
        if( any(input <= 0.0_real32) ) return
        do k = 1, size(input)
            first = max(1, min(k-3, size(input)-6))
            y = log(real(input(first:first+6),real64))
            z = real(k-first-3,real64)
            value = sum(y)/7.d0 + dot_product(X,y)*z/28.d0 + &
                &dot_product(X*X-4.d0,y)*(z*z-4.d0)/84.d0
            if( value < log(real(tiny(1.0_real32),real64)) .or. &
                &value > log(real(huge(1.0_real32),real64)) )then
                message = 'smoothed bootstrap sigma2 is outside the representable positive range'; return
            endif
            output(k) = real(exp(value),real32)
        enddo
        status = 0
        message = ''
    end subroutine smooth_log_spectrum

end module simple_sigma2_frozen
