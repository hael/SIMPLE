
module simple_binoris_io
use simple_defs
use simple_ori,     only: ori
use simple_oris,    only: oris
use simple_fileio,  only: file_exists, nlines, fileio_errmsg
use simple_strings, only: str_has_substr
implicit none

interface binread_oritab
    module procedure binread_oritab_1
end interface

interface binwrite_oritab
    module procedure binwrite_oritab_1
end interface

contains

    subroutine binread_oritab_1( fname, a, fromto, nst )
        character(len=*),  intent(in)    :: fname
        class(oris),       intent(inout) :: a
        integer,           intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
        integer :: irec
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_oritab_1'
        endif
        call a%read(fname, fromto=fromto, nst=nst)
    end subroutine binread_oritab_1

    subroutine binread_ctfparams_state_eo( fname, a, fromto )
        character(len=*), intent(in)    :: fname
        class(oris),      intent(inout) :: a
        integer,          intent(in)    :: fromto(2)
        integer :: irec
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_ctfparams_and_state'
        endif
        call a%read_ctfparams_state_eo(fname)
    end subroutine binread_ctfparams_state_eo

    function binread_nlines( fname ) result( nl )
        character(len=*), intent(in) :: fname
        integer :: nl
        if( .not. file_exists(fname) )then
            write(*,*) 'file: ', trim(fname)
            stop 'does not exist in cwd; binoris_io :: binread_nlines'
        endif
        nl = nlines(fname)
    end function binread_nlines

    subroutine binwrite_oritab_1( fname, a, fromto )
        character(len=*), intent(in)    :: fname
        class(oris),      intent(inout) :: a
        integer,          intent(in)    :: fromto(2)
        call a%write(fname, fromto)
    end subroutine binwrite_oritab_1

end module simple_binoris_io
