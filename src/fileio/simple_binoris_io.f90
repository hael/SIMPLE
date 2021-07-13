module simple_binoris_io
use simple_defs
use simple_strings
use simple_fileio
use simple_error
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
use simple_syslib,     only: file_exists
implicit none

public :: binread_oritab, binread_ctfparams_state_eo, binread_nlines, binwrite_oritab
private
#include "simple_local_flags.inc"

contains

    subroutine binread_oritab( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        if( .not. file_exists(fname) )then
            THROW_HARD('file: '//trim(fname)//' does not exist in cwd')
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read(fname)
            case('T')
                call a%read(fname, fromto=fromto)
            case DEFAULT
                THROW_HARD('format: '//fname2format(fname)//' unsupported')
        end select
    end subroutine binread_oritab

    subroutine binread_ctfparams_state_eo( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in), optional    :: fromto(2)
        if( .not. file_exists(fname) )then
            THROW_HARD('file: '//trim(fname)//' does not exist in cwd')
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read_ctfparams_state_eo(fname)
            case('T')
                call a%read_ctfparams_state_eo(fname)
            case DEFAULT
                THROW_HARD('format of file: '//trim(fname)//' unsupported')
        end select
    end subroutine binread_ctfparams_state_eo

    function binread_nlines( fname ) result( nl )
        use simple_binoris,    only: binoris
        use simple_parameters, only: params_glob
        character(len=*), intent(in) :: fname
        integer       :: nl
        type(binoris) :: bos
        nl = -1
        if( .not. file_exists(fname) )then
            THROW_HARD('file: '//trim(fname)//' does not exist in cwd')
        endif
        select case(fname2format(fname))
        case('O')
            call bos%open(fname)
            nl = bos%get_n_records(params_glob%spproj_iseg)
            call bos%close
        case('T')
            nl = nlines(fname)
        case DEFAULT
            THROW_HARD('format of file: '//trim(fname)//' unsupported')
        end select
    end function binread_nlines

    subroutine binwrite_oritab( fname, spproj, a, fromto, isegment )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        integer(kind(ENUM_ORISEG)), optional, intent(in) :: isegment
        select case(fname2format(fname))
            case('O')
                call spproj%write(fname, fromto, isegment)
            case('T')
                call a%write(fname, fromto)
            case DEFAULT
                THROW_HARD('format of file: '//trim(fname)//' unsupported')
        end select
    end subroutine binwrite_oritab

end module simple_binoris_io
