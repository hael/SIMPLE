module simple_ctffind_iter
use simple_syscalls,  only: exec_cmdline
use simple_nrtxtfile, only: nrtxtfile
use simple_strings,   only: real2str
use simple_filehandling ! use all in there
implicit none

public :: ctffind_iter
private

type :: ctffind_iter
  contains
    procedure :: iterate
end type ctffind_iter

contains
    
    subroutine iterate( self, p, imovie, movie_counter, moviename_forctf, fname_ctrl, fname_output, os )
        use simple_params, only: params
        use simple_oris,   only: oris
        class(ctffind_iter),        intent(inout) :: self
        class(params),              intent(inout) :: p
        integer,                    intent(in)    :: imovie
        integer,                    intent(inout) :: movie_counter
        character(len=*),           intent(in)    :: moviename_forctf, fname_ctrl, fname_output
        class(oris),                intent(inout) :: os
        character(len=:), allocatable :: fname_diag
        real,             allocatable :: ctfparams(:,:)
        character(len=STDLEN)         :: cmd_str, fname_param
        type(nrtxtfile)               :: ctfparamfile
        integer                       :: funit, ndatlines, nrecs, j, file_stat
        if( .not. file_exists(moviename_forctf) )&
        & write(*,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        movie_counter = movie_counter + 1
        funit         = get_fileunit()
        if( p%dir_target .ne. '' )then
            fname_diag    = trim(p%dir_target)//'/'//add2fbody(moviename_forctf, p%ext, '_ctffind_diag')
            fname_param   = trim(p%dir_target)//'/'//fname_new_ext(fname_diag, 'txt')
        else
            fname_diag    = add2fbody(moviename_forctf, p%ext, '_ctffind_diag')
            fname_param   = fname_new_ext(fname_diag, 'txt')
        endif
        open(unit=funit, status='REPLACE', action='WRITE', file=fname_ctrl)
        write(funit,'(a)') trim(moviename_forctf)
        write(funit,'(a)') trim(fname_diag)
        write(funit,'(a)') real2str(p%smpd)
        write(funit,'(a)') real2str(p%kv)
        write(funit,'(a)') real2str(p%cs)
        write(funit,'(a)') real2str(p%fraca)
        write(funit,'(a)') real2str(real(p%pspecsz))
        write(funit,'(a)') real2str(p%hp)
        write(funit,'(a)') real2str(p%lp)
        write(funit,'(a)') real2str(1.0e4*p%dfmin)
        write(funit,'(a)') real2str(1.0e4*p%dfmax)
        write(funit,'(a)') real2str(1.0e4*p%astigstep)
        write(funit,'(a)') 'no'
        write(funit,'(a)') 'no'
        write(funit,'(a)') 'yes'
        write(funit,'(a)') real2str(1.0e4*p%expastig)
        write(funit,'(a)') trim(p%phaseplate)
        write(funit,'(a)') 'no';
        close(funit)
        cmd_str = 'cat ' // fname_ctrl//' | ctffind'
        call system(trim(cmd_str))
        call ctfparamfile%new(fname_param, 1)
        ndatlines = ctfparamfile%get_ndatalines()
        nrecs     = ctfparamfile%get_nrecs_per_line()
        allocate( ctfparams(ndatlines,nrecs) )
        do j=1,ndatlines
            call ctfparamfile%readNextDataLine(ctfparams(j,:))
        end do
        call os%set(movie_counter, 'kv',     p%kv                )
        call os%set(movie_counter, 'cs',     p%cs                )
        call os%set(movie_counter, 'fraca',  p%fraca             )
        call os%set(movie_counter, 'dfx',    ctfparams(1,2)/1.0e4)
        call os%set(movie_counter, 'dfy',    ctfparams(1,3)/1.0e4)
        call os%set(movie_counter, 'angast', ctfparams(1,4)      )
        deallocate(ctfparams)
        call ctfparamfile%kill
    end subroutine iterate

end module simple_ctffind_iter