!@descr: the power-spectrum bootstrap of the canonical sigma2 store in per-half sums: the partition files the
!  calc_pspec workers write and the master reads, and the seed spectra the master forms from them
! A worker sums, per even/odd half, the scaled and sanitised power spectra of the sampled particles of its
! partition (calc_pspec_strategy); the master adds the partitions' sums and forms one seed spectrum per half
! (calc_pspec_assemble). The normalisation is the bootstrap's: each half's sum over its number of active
! particles, sampled or not, minus the average image's power over two, then non-positive values repaired.
! The average image is the mean of the images the workers summed (sampled, with a valid spectrum), so each
! partition file also carries how many images its Fourier sum holds.
module simple_sigma2_pspec_sums
use, intrinsic :: iso_fortran_env, only: int64, real32, real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,         only: DTINY, logfhandle
use simple_error,        only: simple_exception
use simple_string_utils, only: int2str
implicit none
private
#include "simple_local_flags.inc"

public :: pspec_part_sums_fname, write_pspec_part_sums, read_pspec_part_sums, bootstrap_seed_spectra, sum_valid_records

character(len=16), parameter :: PSPEC_SUMS_MAGIC = 'SIMPLE_PSPECSUM2'

contains

    !> the per-half sums file of partition ipart
    function pspec_part_sums_fname( ipart ) result( fname )
        integer, intent(in) :: ipart
        character(len=:), allocatable :: fname
        fname = 'init_pspec_part'//int2str(ipart)//'.dat'
    end function pspec_part_sums_fname

    !> writes the per-half sums sums(shells,2) of the partition covering particles fromp:top, and the number
    !! nimgs of images its Fourier sum holds
    subroutine write_pspec_part_sums( fname, fromp, top, nimgs, sums )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: fromp, top, nimgs
        real(real64),     intent(in) :: sums(:,:)
        integer :: funit, io_stat
        if( size(sums,2) /= 2 ) THROW_HARD('one sum per even/odd half; write_pspec_part_sums')
        if( nimgs < 0 .or. nimgs > max(0,top-fromp+1) ) THROW_HARD('image count outside the partition; write_pspec_part_sums')
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='replace', &
            &action='write', iostat=io_stat)
        if( io_stat /= 0 ) THROW_HARD('cannot open the bootstrap sums of a partition for writing')
        write(funit, iostat=io_stat) PSPEC_SUMS_MAGIC, fromp, top, size(sums,1), nimgs, sums
        close(funit)
        if( io_stat /= 0 ) THROW_HARD('cannot write the bootstrap sums of a partition')
    end subroutine write_pspec_part_sums

    !> reads a partition's per-half sums over nshell_expected shells, the particle range it covers and the
    !! number of images its Fourier sum holds; status /= 0 with a message for a missing, foreign, truncated
    !! or oversized file, one summed over another number of shells, or an image count outside the range.
    !! Shell count and file size are checked before anything is allocated.
    subroutine read_pspec_part_sums( fname, nshell_expected, fromp, top, nimgs, sums, status, message )
        character(len=*),          intent(in)  :: fname
        integer,                   intent(in)  :: nshell_expected
        integer,                   intent(out) :: fromp, top, nimgs
        real(real64), allocatable, intent(out) :: sums(:,:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        character(len=16) :: magic
        integer(int64)    :: file_bytes, expected_bytes
        integer :: funit, nshell
        fromp = 0
        top   = -1
        nimgs = 0
        if( nshell_expected < 1 ) THROW_HARD('read_pspec_part_sums needs a positive shell count')
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', &
            &action='read', iostat=status)
        if( status /= 0 )then
            message = 'cannot open the bootstrap sums of a partition: '//trim(fname); return
        endif
        read(funit, iostat=status) magic, fromp, top, nshell, nimgs
        if( status == 0 .and. magic /= PSPEC_SUMS_MAGIC ) status = 1
        if( status /= 0 )then
            close(funit)
            message = 'not a bootstrap sums file of a partition: '//trim(fname); return
        endif
        if( nshell /= nshell_expected )then
            close(funit)
            status  = 1
            message = 'bootstrap sums of a partition over another number of shells: '//trim(fname); return
        endif
        if( nimgs < 0 .or. nimgs > max(0,top-fromp+1) )then
            close(funit)
            status  = 1
            message = 'bootstrap sums of a partition with an image count outside its range: '//trim(fname); return
        endif
        ! the header (magic, fromp, top, nshell, nimgs) and the sums, nothing more
        expected_bytes = int(len(magic),int64) + 4_int64*int(storage_size(nshell)/8,int64) + &
            &2_int64*int(nshell,int64)*int(storage_size(1.0_real64)/8,int64)
        inquire(unit=funit, size=file_bytes, iostat=status)
        if( status /= 0 .or. file_bytes /= expected_bytes )then
            close(funit)
            status  = 1
            message = 'bootstrap sums of a partition of the wrong size: '//trim(fname); return
        endif
        allocate(sums(nshell,2))
        read(funit, iostat=status) sums
        close(funit)
        if( status /= 0 )then
            message = 'truncated bootstrap sums of a partition: '//trim(fname); return
        endif
        message = ''
    end subroutine read_pspec_part_sums

    !> The seed spectrum of each half (seeds(shells,2)) from the sums over all partitions: the sum over the
    !! half's number of active particles nactive (sampled or not), minus pspec_ave (the average image's
    !! power over two), then repaired (remove_negative_sigmas); a half without active particles gets 1.
    subroutine bootstrap_seed_spectra( sums, nactive, pspec_ave, seeds )
        real(real64), intent(in)  :: sums(:,:)
        integer,      intent(in)  :: nactive(2)
        real,         intent(in)  :: pspec_ave(:)
        real(real64), intent(out) :: seeds(:,:)
        integer :: half
        if( size(sums,2) /= 2 .or. size(seeds,2) /= 2 .or. size(seeds,1) /= size(sums,1) .or. &
            &size(pspec_ave) /= size(sums,1) ) THROW_HARD('inconsistent shells; bootstrap_seed_spectra')
        do half = 1, 2
            if( nactive(half) > 0 )then
                seeds(:,half) = sums(:,half) / real(nactive(half),real64)
                seeds(:,half) = seeds(:,half) - real(pspec_ave,real64)
                call remove_negative_sigmas(seeds(:,half), half, nactive(half))
            else
                seeds(:,half) = 1.d0
            endif
        enddo
    end subroutine bootstrap_seed_spectra

    !> Adds the records block(shells, rows first_row:) of the active particles whose every shell is finite and
    !! positive to their half's sums and counts, in row order, as a reduction of the store sums them; the
    !! active records it skips are counted in ninvalid
    subroutine sum_valid_records( block, first_row, active, eo, sums, counts, ninvalid )
        real(real32),   intent(in)    :: block(:,:)
        integer,        intent(in)    :: first_row
        logical,        intent(in)    :: active(:)
        integer,        intent(in)    :: eo(:)
        real(real64),   intent(inout) :: sums(:,:)
        integer(int64), intent(inout) :: counts(2), ninvalid
        integer :: i, row, half
        do i = 1, size(block,2)
            row = first_row + i - 1
            if( .not. active(row) ) cycle
            if( .not.(all(ieee_is_finite(block(:,i))) .and. all(block(:,i) > 0.0_real32)) )then
                ninvalid = ninvalid + 1_int64
                cycle
            endif
            if( eo(row) < 0 .or. eo(row) > 1 ) THROW_HARD('active particle has invalid even/odd assignment; sum_valid_records')
            half = eo(row) + 1
            sums(:,half) = sums(:,half) + real(block(:,i),real64)
            counts(half) = counts(half) + 1_int64
        enddo
    end subroutine sum_valid_records

    !> non-finite values become 1; a spectrum with nothing above DTINY becomes 1 throughout; any other
    !! value at or below DTINY takes the previous shell's value when that is positive, else the next
    !! positive value, else 1
    subroutine remove_negative_sigmas( spec, half, weight )
        real(real64), intent(inout) :: spec(:)
        integer,      intent(in)    :: half, weight
        logical :: is_positive, fixed_from_prev
        integer :: nn, idx
        where( .not. ieee_is_finite(spec) ) spec = 1.d0
        if( .not. any(spec > DTINY) )then
            write(logfhandle,*) '>>> WARNING: calc_pspec_assemble floored empty global sigma2 spectrum to 1.0; ', &
                'eo/weight: ', half, weight
            spec = 1.d0
            return
        endif
        do idx = 1, size(spec)
            if( spec(idx) <= DTINY )then
                ! first try the previous value
                fixed_from_prev = .false.
                if( idx - 1 >= 1 )then
                    if( spec(idx-1) > DTINY )then
                        spec(idx) = spec(idx-1)
                        fixed_from_prev = .true.
                    end if
                end if
                if( .not. fixed_from_prev )then
                    is_positive = .false.
                    nn          = idx
                    do while (.not. is_positive)
                        nn = nn + 1
                        if( nn > size(spec) )then
                            spec(idx) = 1.d0
                            exit
                        end if
                        if( spec(nn) > DTINY )then
                            is_positive = .true.
                            spec(idx) = spec(nn)
                        end if
                    end do
                end if
            end if
        end do
    end subroutine remove_negative_sigmas

end module simple_sigma2_pspec_sums
