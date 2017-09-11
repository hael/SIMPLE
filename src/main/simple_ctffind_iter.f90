! ctffind iterator
#include "simple_lib.f08"
module simple_ctffind_iter
use simple_defs
use simple_strings,   only: real2str
use simple_fileio,    only: fopen,fclose,fileio_errmsg,file_exists,add2fbody,fname_new_ext  ! use all in there
use simple_nrtxtfile, only: nrtxtfile
use simple_syslib,    only: exec_cmdline, alloc_errchk
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

        fname_diag    = add2fbody(moviename_forctf, p%ext, '_ctffind_diag')
        fname_param   = fname_new_ext(fname_diag, 'txt')

        call fopen(funit, status='REPLACE', action='WRITE', file=trim(fname_ctrl),iostat=file_stat)
        if(file_stat/=0) call fileio_errmsg("ctffind_iter:: iterate fopen failed "//trim(fname_ctrl),file_stat)
        write(funit,'(a)') trim(moviename_forctf)      ! integrated movie used for fitting
        write(funit,'(a)') trim(fname_diag)            ! diagnostic file
        write(funit,'(a)') real2str(p%smpd)            ! magnification dependent sampling distance
        write(funit,'(a)') real2str(p%kv)              ! acceleration voltage, default 300 kV (Titan)
        write(funit,'(a)') real2str(p%cs)              ! sperical aberration, default 2.7 mm (Titan/Arctica)
        write(funit,'(a)') real2str(p%fraca)           ! fraction of amplitude contrast, default 0.1 (10%)
        write(funit,'(a)') real2str(real(p%pspecsz))   ! size of power spectrum, default 1024 to avoid aliasing
        write(funit,'(a)') real2str(p%hp)              ! high-pass limit, default 30.0 A
        write(funit,'(a)') real2str(p%lp)              ! low-pass  limit, default 5.0 A
        write(funit,'(a)') real2str(1.0e4*p%dfmin)     ! minimum defocus, default 0.5 microns
        write(funit,'(a)') real2str(1.0e4*p%dfmax)     ! maximum defocus, default 5.0 microns
        write(funit,'(a)') real2str(1.0e4*p%dfstep)    ! defocus grid search step size, default 0.05 microns
        write(funit,'(a)') 'no'                        ! do you know what astigmatism is present?
        write(funit,'(a)') 'yes'                       ! slower, more exhaustive search
        write(funit,'(a)') 'yes'                       ! use a restraint on astigmatism
        write(funit,'(a)') real2str(1.0e4*p%astigtol)  ! defocus grid search step size, default 0.05 microns
        write(funit,'(a)') trim(p%phaseplate)          ! phase-plate or not (yes|no) {no}
        write(funit,'(a)') 'no';                       ! set expert options
        call fclose(funit,errmsg="ctffind_iter:: iterate fopen failed "//trim(fname_ctrl))
        cmd_str = 'cat ' // fname_ctrl//' | ctffind'
        call exec_cmdline(trim(cmd_str))
        call ctfparamfile%new(fname_param, 1)
        ndatlines = ctfparamfile%get_ndatalines()
        nrecs     = ctfparamfile%get_nrecs_per_line()
        allocate( ctfparams(ndatlines,nrecs) , stat=alloc_stat)
        if(alloc_stat /= 0) allocchk('In: iterate, module: simple_ctffind_iter ctfparams')
        do j=1,ndatlines
            call ctfparamfile%readNextDataLine(ctfparams(j,:))
        end do
        call os%set(movie_counter, 'kv',     p%kv                )
        call os%set(movie_counter, 'cs',     p%cs                )
        call os%set(movie_counter, 'fraca',  p%fraca             )
        call os%set(movie_counter, 'dfx',    ctfparams(1,2)/1.0e4)
        call os%set(movie_counter, 'dfy',    ctfparams(1,3)/1.0e4)
        call os%set(movie_counter, 'angast', ctfparams(1,4)      )
        call os%set(movie_counter, 'ctfres', ctfparams(1,7)      )
        deallocate(ctfparams)
        call ctfparamfile%kill
    end subroutine iterate

end module simple_ctffind_iter
