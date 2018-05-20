! ctf_estimate iterator
module simple_ctf_estimate_iter
include 'simple_lib.f08'
use simple_parameters, only: params_glob
implicit none

public :: ctf_estimate_iter
private

type :: ctf_estimate_iter
  contains
    procedure :: iterate
end type ctf_estimate_iter

contains

    subroutine iterate( self, ctfvars, moviename_forctf, orientation, dir_out )
        use simple_ori,    only: ori
        use simple_image,  only: image
        use simple_ctf_estimate
        class(ctf_estimate_iter), intent(inout) :: self
        type(ctfparams),          intent(in)    :: ctfvars
        character(len=*),         intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        character(len=*),         intent(in)    :: dir_out
        integer                       :: nframes, ldim(3)
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
        call micrograph%new(ldim, ctfvars%smpd)
        call micrograph%read(trim(adjustl(moviename_forctf)), 1)
        ! filter out frequencies lower than the box can express to avoid aliasing
        call micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.)
        ! extract powerspectra
        call pspec_lower%new([params_glob%pspecsz,params_glob%pspecsz,1], ctfvars%smpd)
        call pspec_upper%new([params_glob%pspecsz,params_glob%pspecsz,1], ctfvars%smpd)
        call pspec_all%new([params_glob%pspecsz,params_glob%pspecsz,1],   ctfvars%smpd)
        call micrograph%mic2eospecs(params_glob%pspecsz, 'sqrt', pspec_lower, pspec_upper, pspec_all)
        ! deal with output
        fname_diag = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
        fname_diag = swap_suffix(fname_diag, '_ctf_estimate_diag', FORCTF_SUFFIX)
        fname_diag = trim(dir_out)//'/'//trim(fname_diag)//trim(JPG_EXT)
        ! fitting
        call ctf_estimate_init(pspec_all, pspec_lower, pspec_upper, ctfvars%smpd, ctfvars%kv,&
            &ctfvars%cs, ctfvars%fraca, [params_glob%dfmin,params_glob%dfmax],&
            &[params_glob%hp,params_glob%lp], params_glob%astigtol, ctfvars%l_phaseplate)
        call ctf_estimate_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, ctfscore, fname_diag)
        call ctf_estimate_kill
        ! reporting
        call orientation%set('dfx',        dfx     )
        call orientation%set('dfy',        dfy     )
        call orientation%set('angast',     angast  )
        call orientation%set('phshift',    phshift )
        call orientation%set('ctf_estimatecc',   cc)
        call orientation%set('dferr',      dferr   )
        call orientation%set('ctfscore',   ctfscore)
        ! destruct
        call micrograph%kill
        call pspec_lower%kill
        call pspec_upper%kill
        call pspec_all%kill
    end subroutine iterate

end module simple_ctf_estimate_iter
