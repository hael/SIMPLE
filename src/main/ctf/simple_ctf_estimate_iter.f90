!@descr: does one iteration of CTF estimation
module simple_ctf_estimate_iter
use simple_core_module_api
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_ctf_estimate_fit, only: ctf_estimate_fit
implicit none

public :: ctf_estimate_iter
private
#include "simple_local_flags.inc"

type :: ctf_estimate_iter
    type(image) :: micrograph, thumbnail, img_jpg, pspec4viz
  contains
      procedure          :: iterate
      procedure, private :: gen_thumbnail
      procedure          :: kill
end type ctf_estimate_iter

contains

    subroutine iterate( self, params, ctfvars, moviename_forctf, orientation, dir_out, l_gen_thumb )
        class(ctf_estimate_iter), intent(inout) :: self
        class(parameters),        intent(in)    :: params
        class(ctfparams),         intent(inout) :: ctfvars
        class(string),            intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        class(string),            intent(in)    :: dir_out
        logical,                  intent(in)    :: l_gen_thumb
        type(ctf_estimate_fit) :: ctffit
        type(string) :: fname_diag, moviename_thumb, docname, tmpl_fname
        integer      :: nframes, ldim(3)
        real         :: phlim_lo, phlim_hi
        character(len=STDLEN) :: phmsg
        if( .not. file_exists(moviename_forctf) )&
        & write(logfhandle,*) 'inputted micrograph does not exist: ', moviename_forctf%to_char()
        call find_ldim_nptcls(moviename_forctf, ldim, nframes)
        if( nframes /= 1 )then
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('single frame input to ctf_estimate assumed; iterate')
        endif
        ldim(3) = 1
        ! deal with output
        tmpl_fname = get_fbody(basename(moviename_forctf), MRC_EXT, separator=.false.)
        if( moviename_forctf%has_substr(FORCTF_SUFFIX) )then
            tmpl_fname = tmpl_fname%to_char([1,tmpl_fname%strlen_trim()-7])
        endif
        fname_diag      = filepath(dir_out, tmpl_fname)//'_ctf_estimate_diag'//JPG_EXT
        moviename_thumb = filepath(dir_out, tmpl_fname)//THUMBNAIL_SUFFIX//JPG_EXT
        docname         = filepath(dir_out, tmpl_fname)//'_ctf'//TXT_EXT !//C_NULL_CHAR ! for future use of star format?
        ! filter out frequencies lower than the box can express to avoid aliasing
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(moviename_forctf, 1)
        call self%micrograph%bp(real(params%pspecsz) * ctfvars%smpd, 0.)
        call self%micrograph%ifft
        ! Fitting policy belongs to the estimation request, not to acquisition
        ! metadata. Resolve it here so every caller of the shared iterator has
        ! identical behavior.
        select case(trim(params%fit_phshift))
            case('yes')
                ctfvars%l_fit_phshift = .true.
            case('no')
                ctfvars%l_fit_phshift = .false.
            case DEFAULT
                THROW_HARD('fit_phshift must be yes or no; iterate')
        end select
        ! global fitting
        call ctffit%new(self%micrograph, params%pspecsz, ctfvars, [params%dfmin,params%dfmax],&
            &[params%hp,params%lp], params%astigtol, &
            &[deg2rad(params%phshift_min),deg2rad(params%phshift_max)], deg2rad(params%phshift_step))
        call ctffit%fit(ctfvars, params%ctfresthreshold)
        ctfvars%phshift = canonical_phshift(ctfvars%phshift)
        ! The fitting objective is |CTF|, which has period pi, so the two ends of a
        ! pi-wide window score identically. A phase landing within a grid step of an
        ! edge cannot be told from its partner half a turn away, and micrographs will
        ! scatter between the two, yielding transfer functions of opposite sign that
        ! cancel in class averages. Report it: no window is safe for every phase
        ! plate, so this is the only way to catch a badly placed one on real data.
        if( ctfvars%l_fit_phshift )then
            phlim_lo = deg2rad(params%phshift_min)
            phlim_hi = deg2rad(params%phshift_max)
            if( phlim_hi - phlim_lo >= PI - deg2rad(params%phshift_step) )then
                if( min(ctfvars%phshift - phlim_lo, phlim_hi - ctfvars%phshift) < deg2rad(params%phshift_step) )then
                    write(phmsg,'(A,F7.2,A,F7.2,A,F7.2,A)') 'fitted CTF phase ', rad2deg(ctfvars%phshift), &
                        &' deg sits at the edge of the search window [', params%phshift_min, ',', &
                        &params%phshift_max, '] deg'
                    THROW_WARN(trim(phmsg)//'. Its 180-degree partner fits equally well, so micrographs may split &
                        &between the two and their CTFs will then have opposite sign. Re-centre the window on the &
                        &expected phase.')
                endif
            endif
        endif
        if( l_gen_thumb )then
            call self%gen_thumbnail( ctffit, params%pspecsz )
            call self%img_jpg%write_jpg(moviename_thumb, norm=.true., quality=92)
            call orientation%set('thumb', simple_abspath(moviename_thumb))
        endif
        ! patch based fitting
        if( trim(params%ctfpatch).eq.'yes' )then
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
    subroutine gen_thumbnail( self, ctffit, pspecsz )
        class(ctf_estimate_iter), intent(inout) :: self
        class(ctf_estimate_fit),  intent(inout) :: ctffit
        integer,                  intent(in)    :: pspecsz
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
        if( pspecsz > GUI_PSPECSZ )then
            call self%pspec4viz%fft
            call self%pspec4viz%clip_inplace([GUI_PSPECSZ,GUI_PSPECSZ,1])
        else if( pspecsz < GUI_PSPECSZ )then
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
        call self%thumbnail%kill
        call self%img_jpg%kill
        call self%pspec4viz%kill
    end subroutine kill

end module simple_ctf_estimate_iter
