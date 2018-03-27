! ctf_estimate iterator
module simple_ctf_estimate_iter
#include "simple_lib.f08"
implicit none

public :: ctf_estimate_iter
private

type :: ctf_estimate_iter
  contains
    procedure :: iterate
end type ctf_estimate_iter

contains

    subroutine iterate( self, p, moviename_forctf, orientation, dir_out )
        use simple_params, only: params
        use simple_ori,    only: ori
        use simple_image,  only: image
        use simple_ctf_estimate  ! use all in there
        class(ctf_estimate_iter), intent(inout)      :: self
        class(params),      intent(inout)      :: p
        character(len=*),   intent(in)         :: moviename_forctf
        class(ori),         intent(inout)      :: orientation
        character(len=*),   intent(in)         :: dir_out
        integer                       :: nframes, ldim(3), i
        character(len=:), allocatable :: fname_diag
        type(image)                   :: micrograph, pspec_lower, pspec_upper, pspec_all
        real                          :: dfx, dfy, angast, phshift, cc, dferr, ctfscore
        if( .not. file_exists(moviename_forctf) )&
        & write(*,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        if( nframes /= 1 )then
            print *, 'nframes: ', nframes
            stop 'single frame input to ctf_estimate assumed; simple_ctf_estimate_iter :: iterate'
        endif
        ldim(3) = 1
        call micrograph%new(ldim, p%smpd)
        call micrograph%read(trim(adjustl(moviename_forctf)), 1)
        ! filter out frequencies lower than the box can express to avoid aliasing
        call micrograph%bp(real(p%pspecsz) * p%smpd, 0.)
        ! extract powerspectra
        call pspec_lower%new([p%pspecsz,p%pspecsz,1], p%smpd)
        call pspec_upper%new([p%pspecsz,p%pspecsz,1], p%smpd)
        call pspec_all%new([p%pspecsz,p%pspecsz,1],   p%smpd)
        call micrograph%mic2eospecs(p%pspecsz, 'sqrt', pspec_lower, pspec_upper, pspec_all)
        ! deal with output
        fname_diag    = add2fbody(moviename_forctf, p%ext, '_ctf_estimate_diag')
        fname_diag = trim(dir_out)//'/'//remove_abspath(trim(fname_diag))
        ! fitting
        call ctf_estimate_init(pspec_all, pspec_lower, pspec_upper, p%smpd, p%kv,&
            &p%cs, p%fraca, [p%dfmin,p%dfmax], [p%hp,p%lp], p%astigtol, p%phaseplate)
        call ctf_estimate_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, ctfscore, fname_diag)
        call ctf_estimate_kill
        ! reporting
        call orientation%set('kv',         p%kv    )
        call orientation%set('cs',         p%cs    )
        call orientation%set('fraca',      p%fraca )
        call orientation%set('dfx',        dfx     )
        call orientation%set('dfy',        dfy     )
        call orientation%set('angast',     angast  )
        call orientation%set('phshift',    phshift )
        call orientation%set('ctf',           'yes')
        call orientation%set('ctf_estimatecc',   cc)
        call orientation%set('dferr',      dferr   )
        call orientation%set('ctfscore',   ctfscore)
        ! destruct (to avoid mem-leaks)
        call micrograph%kill
        call pspec_lower%kill
        call pspec_upper%kill
        call pspec_all%kill
    end subroutine iterate

end module simple_ctf_estimate_iter
