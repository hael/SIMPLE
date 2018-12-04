! ctf_estimate iterator
module simple_ctf_estimate_iter
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: ctf_estimate_iter
private
#include "simple_local_flags.inc"

type :: ctf_estimate_iter
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image) :: micrograph, pspec_lower, pspec_upper, pspec_all, thumbnail, img_jpg, pspec4viz
  contains
    procedure :: iterate
end type ctf_estimate_iter

contains

    subroutine iterate( self, ctfvars, moviename_forctf, orientation, dir_out, l_gen_thumb )
        use simple_ori, only: ori
        use simple_ctf_estimate
        class(ctf_estimate_iter), intent(inout) :: self
        type(ctfparams),          intent(in)    :: ctfvars
        character(len=*),         intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        character(len=*),         intent(in)    :: dir_out
        logical,                  intent(in)    :: l_gen_thumb
        character(len=:), allocatable :: fname_diag
        character(len=LONGSTRLEN)     :: moviename_thumb
        real                          :: dfx, dfy, angast, phshift, cc, dferr, ctfscore, cc90, scale
        integer                       :: nframes, ldim(3), ldim_thumb(3)
        if( .not. file_exists(moviename_forctf) )&
        & write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        if( nframes /= 1 )then
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('single frame input to ctf_estimate assumed; iterate')
        endif
        ldim(3) = 1
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(trim(adjustl(moviename_forctf)), 1)
        ! filter out frequencies lower than the box can express to avoid aliasing
        call self%micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.)
        ! extract powerspectra
        call self%pspec_lower%new([params_glob%pspecsz,params_glob%pspecsz,1], ctfvars%smpd)
        call self%pspec_upper%new([params_glob%pspecsz,params_glob%pspecsz,1], ctfvars%smpd)
        call self%pspec_all%new([params_glob%pspecsz,params_glob%pspecsz,1],   ctfvars%smpd)
        call self%micrograph%mic2eospecs(params_glob%pspecsz, 'sqrt', params_glob%hp, self%pspec_lower, self%pspec_upper, self%pspec_all)
        if( l_gen_thumb )then
            ! generate thumbnail
            scale         = real(params_glob%pspecsz)/real(ldim(1))
            ldim_thumb(1) = round2even(real(ldim(1))*scale)
            ldim_thumb(2) = round2even(real(ldim(2))*scale)
            ldim_thumb(3) = 1
            call self%thumbnail%new(ldim_thumb, ctfvars%smpd)
            call self%micrograph%fft()
            call self%micrograph%clip(self%thumbnail)
            call self%thumbnail%ifft()
            ! jpeg output
            call self%pspec4viz%copy(self%pspec_all)
            call self%pspec4viz%scale_pspec4viz
            call self%pspec4viz%collage(self%thumbnail, self%img_jpg)
            moviename_thumb = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
            moviename_thumb = swap_suffix(moviename_thumb, THUMBNAIL_SUFFIX, INTGMOV_SUFFIX)
            moviename_thumb = trim(dir_out)//trim(adjustl(moviename_thumb))//trim(JPG_EXT)
            call self%img_jpg%write_jpg(moviename_thumb, norm=.true., quality=90)
            call make_relativepath(CWD_GLOB,simple_abspath(moviename_thumb, errmsg='simple_ctf_estimate_iter'),moviename_thumb)
            call orientation%set('thumb', trim(moviename_thumb))
        endif
        ! deal with output
        fname_diag = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
        if( str_has_substr(fname_diag, FORCTF_SUFFIX) )then
            fname_diag = swap_suffix(fname_diag, '_ctf_estimate_diag', FORCTF_SUFFIX)
        else if( str_has_substr(fname_diag, INTGMOV_SUFFIX) )then
            fname_diag = swap_suffix(fname_diag, '_ctf_estimate_diag', INTGMOV_SUFFIX)
        endif
        fname_diag = filepath(trim(dir_out),trim(fname_diag)//trim(JPG_EXT))
        ! fitting
        call ctf_estimate_init(self%pspec_all, self%pspec_lower, self%pspec_upper, ctfvars%smpd, ctfvars%kv,&
            &ctfvars%cs, ctfvars%fraca, [params_glob%dfmin,params_glob%dfmax],&
            &[params_glob%hp,params_glob%lp], params_glob%astigtol, ctfvars%l_phaseplate, cc90)
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
        call orientation%set('cc90',       cc90)
        call orientation%set('ctfjpg',     fname_diag)
    end subroutine iterate

end module simple_ctf_estimate_iter
