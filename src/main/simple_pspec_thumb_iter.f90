! iterator for pspec_thumb (a program for motion correction, dose-weighting and frame-weighting of DDD movies)
module simple_pspec_thumb_iter
include 'simple_lib.f08'
use simple_image,      only: image
use simple_cmdline,    only: cmdline
use simple_parameters, only: params_glob
implicit none

public :: pspec_thumb_iter
private
#include "simple_local_flags.inc"

type :: pspec_thumb_iter
    private
    character(len=4) :: speckind = 'sqrt'
    type(string)     :: moviename_thumb, moviename_pspec
    type(image)      :: moviesum, pspec, thumbnail
  contains
    procedure :: iterate
end type pspec_thumb_iter

contains

    subroutine iterate( self, orientation, moviename_intg, dir_out )
        class(pspec_thumb_iter), intent(inout) :: self
        class(ori),              intent(inout) :: orientation
        class(string),           intent(in)    :: moviename_intg, dir_out
        type(string) :: fbody_here, ext
        type(image)  :: img_jpg
        integer      :: ldim(3), ldim_thumb(3), nframes
        real         :: scale, smpd
        ! check, increment counter & print
        if( .not. file_exists(moviename_intg) )then
            write(logfhandle,*) 'inputted integrated movie does not exist: ', moviename_intg%to_char()
        endif
        ! make filenames
        fbody_here = basename(moviename_intg)
        ext        = fname2ext(fbody_here)
        fbody_here = get_fbody(fbody_here, ext)
        self%moviename_pspec = dir_out%to_char()//fbody_here%to_char()//POWSPEC_SUFFIX//params_glob%ext%to_char()
        self%moviename_thumb = dir_out%to_char()//fbody_here%to_char()//THUMBNAIL_SUFFIX//JPG_EXT
        write(logfhandle,'(a,1x,a)') '>>> PROCESSING INTEGRATED MOVIE:', moviename_intg%to_char()
        call find_ldim_nptcls(moviename_intg, ldim, nframes)
        if( nframes /= 1 ) THROW_HARD('imported movie assumed to be integrated but nframes /= 1, aborting; simple_pspec_thumb_iter :: iterate')
        if( .not. orientation%isthere('smpd') ) THROW_HARD('smpd assumed to be set in input orientation, aborting; simple_pspec_thumb_iter :: iterate')
        smpd = orientation%get('smpd')
        call self%moviesum%new(ldim, smpd)
        call self%moviesum%read(moviename_intg)
        ! generate power-spectra
        call self%moviesum%mic2spec(params_glob%pspecsz, self%speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec)
        call self%pspec%write(self%moviename_pspec)
        ! generate thumbnail
        scale         = real(params_glob%pspecsz)/real(ldim(1))
        ldim_thumb(1) = round2even(real(ldim(1))*scale)
        ldim_thumb(2) = round2even(real(ldim(2))*scale)
        ldim_thumb(3) = 1
        call self%thumbnail%new(ldim_thumb, smpd)
        call self%moviesum%fft()
        call self%moviesum%clip(self%thumbnail)
        call self%thumbnail%ifft()
        ! jpeg output
        call self%pspec%collage(self%thumbnail, img_jpg)
        call img_jpg%write_jpg(self%moviename_thumb, norm=.true., quality=90)
        ! report to ori object
        call orientation%set('intg',    simple_abspath(moviename_intg))
        call orientation%set('thumb',   simple_abspath(self%moviename_thumb))
        call orientation%set('imgkind', 'intg')
        ! destruct
        call self%pspec%kill
        call self%moviesum%kill
        call img_jpg%kill
        call self%thumbnail%kill
    end subroutine iterate

end module  simple_pspec_thumb_iter
