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
        class(image),     allocatable  :: even_imgs(:), odd_imgs(:)
        type(image)                    :: micrograph, pspec_lower, pspec_upper, pspec_all
        real                           :: dfx, dfy, angast, phshift, cc, dferr
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
        ! filter out frequencies lower than the box can express to avoid aliasing
        call micrograph%bp(real(p%pspecsz) * p%smpd, 0.)
        ! extract even=lower and odd=upper powerspectra
        call micrograph%mic2eospecs(p%pspecsz, 'sqrt', even_imgs, odd_imgs)
        ! prepate lower/upper/all spectra
        call pspec_lower%new([p%pspecsz,p%pspecsz,1], p%smpd)
        call pspec_upper%new([p%pspecsz,p%pspecsz,1], p%smpd)
        do i=1,size(even_imgs)
            call pspec_lower%add(even_imgs(i))
        end do
        pspec_all = pspec_lower
        do i=1,size(odd_imgs)
            call pspec_upper%add(odd_imgs(i))
        end do
        call pspec_all%add(pspec_upper)
        ! deal with output
        movie_counter = movie_counter + 1
        fname_diag    = add2fbody(moviename_forctf, p%ext, '_ctffit_diag')
        if( present(dir_out) )then
            fname_diag = remove_abspath(trim(fname_diag))
            fname_diag = trim(dir_out)//'/'//trim(fname_diag)
        endif
        ! fitting
        call ctffit_init(pspec_all, pspec_lower, pspec_upper, p%smpd, p%kv,&
            &p%cs, p%fraca, [p%dfmin,p%dfmax], [p%hp,p%lp], p%phaseplate)
        call ctffit_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, fname_diag)
        call ctffit_kill
        ! reporting
        call os%set(movie_counter, 'kv',         p%kv   )
        call os%set(movie_counter, 'cs',         p%cs   )
        call os%set(movie_counter, 'fraca',      p%fraca)
        call os%set(movie_counter, 'dfx',        dfx    )
        call os%set(movie_counter, 'dfy',        dfy    )
        call os%set(movie_counter, 'angast',     angast )
        call os%set(movie_counter, 'phshift',    phshift)
        call os%set(movie_counter, 'ctffitcc',   cc     )
        call os%set(movie_counter, 'dferr',      dferr  )
    end subroutine iterate

end module simple_ctffit_iter
