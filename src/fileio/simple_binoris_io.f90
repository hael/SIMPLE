
module simple_binoris_io
use simple_defs
use simple_ori,        only: ori
use simple_oris,       only: oris
use simple_fileio,     only: file_exists, nlines, fileiochk, fname2format
use simple_strings,    only: str_has_substr
use simple_sp_project, only: sp_project
implicit none

contains

    subroutine binread_oritab( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_1'
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read(fname)
            case('T')
                call a%read(fname, fromto=fromto)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_oritab'
        end select
    end subroutine binread_oritab

    subroutine binread_ctfparams_state_eo( fname, spproj, a, fromto )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in), optional    :: fromto(2)
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_ctfparams_and_state'
        endif
        select case(fname2format(fname))
            case('O')
                call spproj%read_ctfparams_state_eo(fname)
            case('T')
                call a%read_ctfparams_state_eo(fname)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_ctfparams_state_eo'
        end select
    end subroutine binread_ctfparams_state_eo

    function binread_nlines( p, fname ) result( nl )
        use simple_params,  only: params
        use simple_binoris, only: binoris
        class(params),    intent(in) :: p
        character(len=*), intent(in) :: fname
        integer       :: nl
        type(binoris) :: bos
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_nlines'
        endif
        select case(fname2format(fname))
            case('O')
                call bos%open(fname)
                nl = bos%get_n_records(p%spproj_a_seg)
                call bos%close
            case('T')
                nl = nlines(fname)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binread_nlines'
        end select
    end function binread_nlines

    subroutine binwrite_oritab( fname, spproj, a, fromto, isegment )
        character(len=*),  intent(in)    :: fname
        class(sp_project), intent(inout) :: spproj
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        integer, optional, intent(in)    :: isegment
        select case(fname2format(fname))
            case('O')
                call spproj%write(fname, fromto, isegment)
            case('T')
                call a%write(fname, fromto)
            case DEFAULT
                write(*,*) 'file: ', trim(fname)
                stop 'format unsupported; simple_binoris_io :: binwrite_oritab'
        end select
    end subroutine binwrite_oritab

end module simple_binoris_io
