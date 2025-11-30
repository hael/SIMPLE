! ctf_estimate iterator
module simple_ctf_estimate_iter
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_image,            only: image
use simple_ctf_estimate_fit, only: ctf_estimate_fit
implicit none

public :: ctf_estimate_iter
private
#include "simple_local_flags.inc"

type :: ctf_estimate_iter
    type(image) :: micrograph, pspec, thumbnail, img_jpg, pspec4viz
  contains
      procedure          :: iterate
      procedure, private :: gen_thumbnail
      procedure          :: kill
end type ctf_estimate_iter

contains

    subroutine iterate( self, ctfvars, moviename_forctf, orientation, dir_out, l_gen_thumb )
        class(ctf_estimate_iter), intent(inout) :: self
        class(ctfparams),         intent(inout) :: ctfvars
        class(string),            intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        class(string),            intent(in)    :: dir_out
        logical,                  intent(in)    :: l_gen_thumb
        type(ctf_estimate_fit) :: ctffit
        type(string) :: fname_diag, moviename_thumb, docname, tmpl_fname
        integer      :: nframes, ldim(3)
        if( .not. file_exists(moviename_forctf) )&
        & write(logfhandle,*) 'inputted micrograph does not exist: ', moviename_forctf%to_char()
        call find_ldim_nptcls(moviename_forctf, ldim, nframes)
        if( nframes /= 1 )then
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('single frame input to ctf_estimate assumed; iterate')
        endif
        ldim(3) = 1
        ! deal with output
        tmpl_fname = get_fbody(basename(moviename_forctf), params_glob%ext, separator=.false.)
        if( moviename_forctf%has_substr(FORCTF_SUFFIX) )then
            tmpl_fname = tmpl_fname%to_char([1,tmpl_fname%strlen_trim()-7])
        endif
        fname_diag      = filepath(dir_out, tmpl_fname)//'_ctf_estimate_diag'//JPG_EXT
        moviename_thumb = filepath(dir_out, tmpl_fname)//THUMBNAIL_SUFFIX//JPG_EXT
        docname         = filepath(dir_out, tmpl_fname)//'_ctf'//TXT_EXT !//C_NULL_CHAR ! for future use of star format?
        ! filter out frequencies lower than the box can express to avoid aliasing
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(moviename_forctf, 1)
        call self%micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.)
        call self%micrograph%ifft
        ! global fitting
        call ctffit%new(self%micrograph, params_glob%pspecsz, ctfvars, [params_glob%dfmin,params_glob%dfmax],&
            &[params_glob%hp,params_glob%lp], params_glob%astigtol)
        call ctffit%fit( ctfvars )
        if( l_gen_thumb )then
            call self%gen_thumbnail( ctffit )
            call self%img_jpg%write_jpg(moviename_thumb, norm=.true., quality=92)
            call orientation%set('thumb', simple_abspath(moviename_thumb))
        endif
        ! patch based fitting
        if( trim(params_glob%ctfpatch).eq.'yes' )then
            call ctffit%fit_patches
            call ctffit%write_doc(moviename_forctf, docname)
            call orientation%set('ctfdoc', simple_abspath(docname))
        endif
        ! diagnostic image
        call ctffit%write_diagnostic(fname_diag)
        ! reporting
        call orientation%set_dfx(              ctfvars%dfx)
        call orientation%set_dfy(              ctfvars%dfy)
        call orientation%set('df',             abs((ctfvars%dfx + ctfvars%dfy) / 2))
        call orientation%set('angast',         ctfvars%angast)
        call orientation%set('phshift',        ctfvars%phshift)
        call orientation%set('ctf_estimatecc', ctffit%get_ccfit())
        call orientation%set('ctfres',         ctffit%get_ctfres())
        call orientation%set('icefrac',        ctffit%get_icefrac())
        call orientation%set('astig',          ctffit%get_astig())
        call orientation%set('ctfjpg',         simple_abspath(fname_diag))
        ! clean
        call ctffit%kill
    end subroutine iterate

    ! generate thumbnail
    subroutine gen_thumbnail( self, ctffit )
        class(ctf_estimate_iter), intent(inout) :: self
        class(ctf_estimate_fit),  intent(inout) :: ctffit
        type(image) :: tmp
        real        :: scale, smpd
        integer     :: ldim(3), ldim_thumb(3)
        ! thumbnail
        smpd  = self%micrograph%get_smpd()
        ldim  = self%micrograph%get_ldim()
        scale = real(GUI_PSPECSZ)/real(maxval(ldim(1:2)))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_thumb(3)   = 1
        call self%thumbnail%new(ldim_thumb, smpd)
        call self%micrograph%fft()
        call self%micrograph%clip(self%thumbnail)
        call self%thumbnail%ifft
        ! spectrum
        call ctffit%get_pspec(self%pspec4viz)
        call self%pspec4viz%scale_pspec4viz
        if( params_glob%pspecsz > GUI_PSPECSZ )then
            call self%pspec4viz%fft
            call self%pspec4viz%clip_inplace([GUI_PSPECSZ,GUI_PSPECSZ,1])
        else if( params_glob%pspecsz < GUI_PSPECSZ )then
            tmp = self%pspec4viz
            call self%pspec4viz%zero_and_unflag_ft
            call tmp%fft
            call tmp%pad(self%pspec4viz)
            self%pspec4viz = tmp
        endif
        call self%pspec4viz%ifft
        ! assembly
        call self%pspec4viz%collage(self%thumbnail, self%img_jpg)
        ! cleanup
        call tmp%kill
    end subroutine gen_thumbnail

    subroutine kill( self )
        class(ctf_estimate_iter), intent(inout) :: self
        call self%micrograph%kill
        call self%pspec%kill
        call self%thumbnail%kill
        call self%img_jpg%kill
        call self%pspec4viz%kill
    end subroutine kill

end module simple_ctf_estimate_iter
