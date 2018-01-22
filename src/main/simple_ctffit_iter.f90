! ctffit iterator
module simple_ctffit_iter
#include "simple_lib.f08"
implicit none

public :: ctffit_iter
private

type :: ctffit_iter
  contains
    procedure :: iterate
end type ctffit_iter

contains

    subroutine iterate( self, p, imovie, movie_counter, moviename_forctf, os, dir_out )
        use simple_params, only: params
        use simple_oris,   only: oris
        use simple_image,  only: image
        use simple_ctffit  ! use all in there
        class(ctffit_iter), intent(inout)      :: self
        class(params),      intent(inout)      :: p
        integer,            intent(in)         :: imovie
        integer,            intent(inout)      :: movie_counter
        character(len=*),   intent(in)         :: moviename_forctf
        class(oris),        intent(inout)      :: os
        character(len=*), optional, intent(in) :: dir_out
        integer                        :: nframes, ldim(3), i, neven, nodd
        character(len=:), allocatable  :: fname_diag
        type(image)                    :: micrograph, pspec
        class(image),     allocatable  :: even_imgs(:), odd_imgs(:)
        real                           :: dfx, dfy, angast, phshift, cc, ccvalid
        if( .not. file_exists(moviename_forctf) )&
        & write(*,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        if( nframes /= 1 )then
            print *, 'nframes: ', nframes
            stop 'single frame input to ctffit assumed; simple_ctffit_iter :: iterate'
        endif
        ldim(3) = 1
        call micrograph%new(ldim, p%smpd)
        call micrograph%read(trim(adjustl(moviename_forctf)), 1)

        ! SHOULD HIGH-PASS FILTER BASED ON PSPECSZ AS IN EXTRACT

        pspec         = micrograph%mic2spec(p%pspecsz, 'sqrt')
        call micrograph%mic2eoimgs(p%pspecsz, even_imgs, odd_imgs)
        movie_counter = movie_counter + 1
        fname_diag    = add2fbody(moviename_forctf, p%ext, '_ctffit_diag')
        if( present(dir_out) )then
            fname_diag = remove_abspath(trim(fname_diag))
            fname_diag = trim(dir_out)//'/'//trim(fname_diag)
        endif
        call ctffit_init(pspec, p%smpd, p%kv, p%cs, p%fraca, [p%dfmin,p%dfmax],&
            &[p%hp,p%lp], p%phaseplate)
        call ctffit_srch(dfx, dfy, angast, phshift, cc, fname_diag)
        call ctffit_validate(even_imgs, odd_imgs, ccvalid)
        call ctffit_kill
        call os%set(movie_counter, 'kv',         p%kv   )
        call os%set(movie_counter, 'cs',         p%cs   )
        call os%set(movie_counter, 'fraca',      p%fraca)
        call os%set(movie_counter, 'dfx',        dfx    )
        call os%set(movie_counter, 'dfy',        dfy    )
        call os%set(movie_counter, 'angast',     angast )
        call os%set(movie_counter, 'phshift',    phshift)
        call os%set(movie_counter, 'ctffitcc',   cc     )
        call os%set(movie_counter, 'ctffiteocc', ccvalid)
        ! destruct (to ensure no memory leaks)
        neven = size(even_imgs)
        nodd  = size(odd_imgs)
        do i=1,neven
            call even_imgs(i)%kill
        end do
        do i=1,nodd
            call odd_imgs(i)%kill
        end do
        deallocate(even_imgs, odd_imgs)
    end subroutine iterate

end module simple_ctffit_iter
