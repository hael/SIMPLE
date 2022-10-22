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
        character(len=*),         intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        character(len=*),         intent(in)    :: dir_out
        logical,                  intent(in)    :: l_gen_thumb
        type(ctf_estimate_fit)        :: ctffit
        character(len=:), allocatable :: fname_diag
        character(len=LONGSTRLEN)     :: moviename_thumb, rel_moviename_thumb, rel_fname, docname, tmpl_fname
        integer                       :: nframes, ldim(3)
        if( .not. file_exists(moviename_forctf) )&
        & write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        if( nframes /= 1 )then
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('single frame input to ctf_estimate assumed; iterate')
        endif
        ldim(3) = 1
        ! deal with output
        tmpl_fname = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
        if( str_has_substr(moviename_forctf, FORCTF_SUFFIX) )then
            tmpl_fname = tmpl_fname(1:len_trim(tmpl_fname)-7)
        endif
        fname_diag      = filepath(trim(dir_out),trim(adjustl(tmpl_fname))//'_ctf_estimate_diag'//trim(JPG_EXT))
        moviename_thumb = filepath(trim(dir_out),trim(adjustl(tmpl_fname))//trim(THUMBNAIL_SUFFIX)//trim(JPG_EXT))
        docname         = filepath(trim(dir_out),trim(adjustl(tmpl_fname))//'_ctf'//trim(TXT_EXT)) !//C_NULL_CHAR ! for future use of star format?
        ! filter out frequencies lower than the box can express to avoid aliasing
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(trim(adjustl(moviename_forctf)), 1)
        call self%micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.)
        call self%micrograph%ifft
        ! global fitting
        call ctffit%new(self%micrograph, params_glob%pspecsz, ctfvars, [params_glob%dfmin,params_glob%dfmax],&
            &[params_glob%hp,params_glob%lp], params_glob%astigtol)
        call ctffit%fit( ctfvars )
        if( l_gen_thumb )then
            call self%gen_thumbnail( ctffit )
            call self%img_jpg%write_jpg(moviename_thumb, norm=.true., quality=92)
            call make_relativepath(CWD_GLOB,moviename_thumb, rel_moviename_thumb)
            call orientation%set('thumb', trim(rel_moviename_thumb))
        endif
        ! patch based fitting
        if( trim(params_glob%ctfpatch).eq.'yes' )then
            call ctffit%fit_patches
            call ctffit%write_doc(moviename_forctf, docname)
            call make_relativepath(CWD_GLOB,docname,rel_fname)
            call orientation%set('ctfdoc', rel_fname)
        endif
        call ctffit%write_diagnostic(fname_diag)
        call make_relativepath(CWD_GLOB,fname_diag, rel_fname)
        ! reporting
        call orientation%set_dfx(              ctfvars%dfx)
        call orientation%set_dfy(              ctfvars%dfy)
        call orientation%set('angast',         ctfvars%angast)
        call orientation%set('phshift',        ctfvars%phshift)
        call orientation%set('ctf_estimatecc', ctffit%get_ccfit())
        call orientation%set('ctfres',         ctffit%get_ctfres())
        call orientation%set('ctfjpg',         rel_fname)
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
