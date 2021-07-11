module simple_ori_io
use simple_defs
use simple_fileio
use simple_error, only: simple_exception
use simple_oris,  only: oris
implicit none

public :: open_ori_io, close_ori_io, write_empty_ori_file, write_oris, read_oris
public :: get_ori_filesz
private
#include "simple_local_flags.inc"

type ori_mini
    real    :: eul(3)    = 0.
    real    :: shift(2)  = 0.
    real    :: corr      = 0.
    real    :: specscore = 0.
    integer :: iptcl     = 0
    integer :: class     = 0
    integer :: proj      = 0
    integer :: state     = 0
end type ori_mini

integer :: recsz = 0, funit = 0
logical :: l_open

contains

    subroutine open_ori_io( fname )
        character(len=*), intent(in) :: fname
        type(ori_mini) :: ori_record
        integer :: tmpunit, io_stat
        if( .not. l_open )then
            inquire(iolength=recsz) ori_record
            call fopen(tmpunit, trim(fname), access='DIRECT', action='READWRITE',&
                &status='UNKNOWN', form='UNFORMATTED', recl=recsz, iostat=io_stat)
            call fileiochk('open_ori_io '//trim(fname), io_stat)
            funit  = tmpunit
            l_open = .true.
        endif
    end subroutine open_ori_io

    subroutine close_ori_io
        if( l_open )then
            call fclose(funit)
            funit  = 0
            l_open = .false.
        endif
    end subroutine close_ori_io

    subroutine write_empty_ori_file( fname, fromto )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: fromto(2)
        type(ori_mini) :: ori_record
        integer        :: recind, iptcl
        call open_ori_io( fname )
        do iptcl=fromto(1),fromto(2)
            ! indexing
            recind = iptcl - fromto(1) + 1
            if( recind < 1 )then
                write(logfhandle,*) 'iptcl:  ', iptcl
                write(logfhandle,*) 'fromto: ', fromto
                write(logfhandle,*) 'recind: ', recind
                THROW_HARD('nonsensical record index; write_empty_ori_file')
            endif
            ! write record
            write(funit,rec=recind) ori_record
        end do
        call close_ori_io
    end subroutine write_empty_ori_file

    subroutine write_oris( os, fromto )
        class(oris), intent(inout) :: os
        integer,     intent(in)    :: fromto(2)
        type(ori_mini) :: ori_record
        integer        :: iptcl, recind
        if( l_open )then
            do iptcl=fromto(1),fromto(2)
                ! indexing
                recind = iptcl - fromto(1) + 1
                if( recind < 1 )then
                    write(logfhandle,*) 'iptcl:  ', iptcl
                    write(logfhandle,*) 'fromto: ', fromto
                    write(logfhandle,*) 'recind: ', recind
                    THROW_HARD('nonsensical record index; write_ori')
                endif
                ! transfer record data
                ori_record%eul       = os%get_euler(iptcl)
                ori_record%shift     = os%get_2Dshift(iptcl)
                ori_record%corr      = os%get(iptcl, 'corr')
                ori_record%specscore = os%get(iptcl, 'specscore')
                ori_record%iptcl     = iptcl
                ori_record%class     = nint(os%get(iptcl, 'class'))
                ori_record%proj      = nint(os%get(iptcl, 'proj'))
                ori_record%state     = os%get_state(iptcl)
                ! write record
                write(funit,rec=recind) ori_record
            end do
        else
            THROW_HARD('file not open; write_ori')
        endif
    end subroutine write_oris

    subroutine read_oris( os, fromto )
        class(oris), intent(inout) :: os
        integer,     intent(in)    :: fromto(2)
        type(ori_mini) :: ori_record
        integer        :: iptcl, recind
        if( l_open )then
            do iptcl=fromto(1),fromto(2)
                ! indexing
                recind = iptcl - fromto(1) + 1
                if( recind < 1 )then
                    write(logfhandle,*) 'iptcl:  ', iptcl
                    write(logfhandle,*) 'fromto: ', fromto
                    write(logfhandle,*) 'recind: ', recind
                    THROW_HARD('nonsensical record index; read_ori')
                endif
                ! read record
                read(funit,rec=recind) ori_record
                ! transfer record data
                call os%set_euler(iptcl, ori_record%eul)
                call os%set_shift(iptcl, ori_record%shift)
                call os%set(iptcl, 'corr',      ori_record%corr)
                call os%set(iptcl, 'specscore', ori_record%specscore)
                call os%set(iptcl, 'iptcl',     real(ori_record%iptcl))
                call os%set(iptcl, 'class',     real(ori_record%class))
                call os%set(iptcl, 'proj',      real(ori_record%proj))
                call os%set(iptcl, 'state',     real(ori_record%state))
            end do
        endif
    end subroutine read_oris

    integer function get_ori_filesz( fname )
        character(len=*), intent(in) :: fname
        type(ori_mini) :: ori_record
        integer :: recsz, ori_recsz
        get_ori_filesz = 0
        inquire(iolength=ori_recsz) ori_record
        inquire(file=fname, size=recsz)
        if( mod(recsz,ori_recsz) /= 0 ) THROW_HARD('ori-file has inconsistent size: '//trim(fname))
        get_ori_filesz = recsz / ori_recsz
    end function get_ori_filesz

end module simple_ori_io
