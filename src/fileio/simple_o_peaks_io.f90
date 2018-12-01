module simple_o_peaks_io
use simple_defs
use simple_fileio
use simple_error, only: simple_exception
use simple_oris,  only: oris
implicit none

public :: open_o_peaks_io, close_o_peaks_io, write_empty_o_peaks_file, write_o_peak, read_o_peak
private
#include "simple_local_flags.inc"

! type for one peak orientation
type o_peak_ori
    real    :: eul(3)   = 0.
    real    :: shift(2) = 0.
    real    :: corr     = 0.
    real    :: ow       = 0.
    integer :: iptcl    = 0
    integer :: proj     = 0
    integer :: state    = 0
end type o_peak_ori

integer :: recsz = 0, funit = 0
logical :: l_open

contains

    subroutine open_o_peaks_io( fname )
        character(len=*), intent(in) :: fname
        type(o_peak_ori) :: o_peak_record(NPEAKS2REFINE)
        integer :: tmpunit, io_stat
        if( .not. l_open )then
            inquire(iolength=recsz) o_peak_record
            call fopen(tmpunit, trim(fname), access='DIRECT', action='READWRITE',&
                &status='UNKNOWN', form='UNFORMATTED', recl=recsz, iostat=io_stat)
            call fileiochk('open_o_peaks_io '//trim(fname), io_stat)
            funit  = tmpunit
            l_open = .true.
        endif
    end subroutine open_o_peaks_io

    subroutine close_o_peaks_io
        if( l_open )then
            call fclose(funit)
            funit  = 0
            l_open = .false.
        endif
    end subroutine close_o_peaks_io

    subroutine write_empty_o_peaks_file( fname, fromto )
        character(len=*), intent(in)    :: fname
        integer,          intent(in)    :: fromto(2)
        type(o_peak_ori) :: o_peak_record(NPEAKS2REFINE)
        integer :: recind, iptcl
        call open_o_peaks_io( fname )
        do iptcl=fromto(1),fromto(2)
            ! index stuff
            recind = iptcl - fromto(1) + 1
            if( recind < 1 )then
                write(logfhandle,*) 'iptcl:  ', iptcl
                write(logfhandle,*) 'fromto: ', fromto
                write(logfhandle,*) 'recind: ', recind
                THROW_HARD('nonsensical record index; write_o_peak')
            endif
            ! write record
            write(funit,rec=recind) o_peak_record
        end do
        call close_o_peaks_io
    end subroutine write_empty_o_peaks_file

    subroutine write_o_peak( o_peak, fromto, iptcl )
        integer,     intent(in)    :: fromto(2), iptcl
        class(oris), intent(inout) :: o_peak
        type(o_peak_ori) :: o_peak_record(NPEAKS2REFINE)
        integer :: recind, ipeak, npeaks
        real    :: ow
        if( l_open )then
            ! index stuff
            recind = iptcl - fromto(1) + 1
            if( recind < 1 )then
                write(logfhandle,*) 'iptcl:  ', iptcl
                write(logfhandle,*) 'fromto: ', fromto
                write(logfhandle,*) 'recind: ', recind
                THROW_HARD('nonsensical record index; write_o_peak')
            endif
            npeaks = o_peak%get_noris()
            if( npeaks == 0 )then
                THROW_WARN('npeaks == 0 for peak set; write_o_peak')
                return
            endif
            if( npeaks > NPEAKS2REFINE )then
                THROW_HARD('npeaks too large to conform with recl; write_o_peak')
            endif
            ! transfer record data
            do ipeak=1,npeaks
                ow = o_peak%get(ipeak,'ow')
                if( ow > TINY )then
                    o_peak_record(ipeak)%eul   = o_peak%get_euler(ipeak)
                    o_peak_record(ipeak)%shift = o_peak%get_2Dshift(ipeak)
                    o_peak_record(ipeak)%corr  = o_peak%get(ipeak, 'corr')
                    o_peak_record(ipeak)%ow    = ow
                    o_peak_record(ipeak)%iptcl = iptcl
                    o_peak_record(ipeak)%proj  = nint(o_peak%get(ipeak, 'proj'))
                    o_peak_record(ipeak)%state = o_peak%get_state(ipeak)
                endif
            end do
            ! write record
            write(funit,rec=recind) o_peak_record
        else
            THROW_HARD('file not open; write_o_peak')
        endif
    end subroutine write_o_peak

    subroutine read_o_peak( o_peak, fromto, iptcl, n_nozero )
        integer,     intent(in)    :: fromto(2), iptcl
        class(oris), intent(inout) :: o_peak
        integer,     intent(out)   :: n_nozero
        type(o_peak_ori) :: o_peak_record(NPEAKS2REFINE)
        integer :: recind, ipeak, npeaks, noris
        if( l_open )then
            ! index stuff
            recind = iptcl - fromto(1) + 1
            if( recind < 1 )then
                write(logfhandle,*) 'iptcl:  ', iptcl
                write(logfhandle,*) 'fromto: ', fromto
                write(logfhandle,*) 'recind: ', recind
                THROW_HARD('nonsensical record index; write_o_peak')
            endif
            ! read record
            read(funit,rec=recind) o_peak_record
            ! count # peaks
            npeaks = 0
            do ipeak=1,NPEAKS2REFINE
                if( o_peak_record(ipeak)%ow > TINY ) npeaks = npeaks + 1
            end do
            ! transfer record to o_peak
            noris = o_peak%get_noris()
            if( noris < npeaks ) call o_peak%new(npeaks)
            n_nozero = 0
            do ipeak=1,NPEAKS2REFINE
                if( o_peak_record(ipeak)%ow > TINY )then
                    n_nozero = n_nozero + 1
                    call o_peak%set_euler(n_nozero,    o_peak_record(ipeak)%eul)
                    call o_peak%set_shift(n_nozero,    o_peak_record(ipeak)%shift)
                    call o_peak%set(n_nozero, 'corr',  o_peak_record(ipeak)%corr)
                    call o_peak%set(n_nozero, 'ow',    o_peak_record(ipeak)%ow)
                    call o_peak%set(n_nozero, 'iptcl', real(o_peak_record(ipeak)%iptcl))
                    call o_peak%set(n_nozero, 'proj',  real(o_peak_record(ipeak)%proj))
                    call o_peak%set(n_nozero, 'state', real(o_peak_record(ipeak)%state))
                endif
            end do
        endif
    end subroutine read_o_peak

end module simple_o_peaks_io
