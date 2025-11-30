! iterator for motion_correct (a program for motion correction, dose-weighting and frame-weighting of DDD movies)
module simple_motion_correct_iter
include 'simple_lib.f08'
use simple_image,        only: image
use simple_cmdline,      only: cmdline
use simple_parameters,   only: params_glob
use simple_stackops,     only: frameavg_stack
use simple_motion_correct
implicit none

public :: motion_correct_iter
private
#include "simple_local_flags.inc"

real,             parameter :: PATCH_FIT_THRESHOLD = 4.0 ! threshold for polynomial fitting in pixels
real,             parameter :: SMPD4VIZ_NANO = 0.9
character(len=*), parameter :: speckind = 'sqrt'
real                        :: effective_patch_fit_threshold = PATCH_FIT_THRESHOLD

type :: motion_correct_iter
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image)   :: moviesum, moviesum_corrected
    type(image)  :: moviesum_ctf, pspec_sum, pspec_ctf, img_jpg
    type(image)  :: pspec_half_n_half, thumbnail
    ! these strings are part of the instance for reporting purposes
    type(string) :: moviename, moviename_intg
    type(string) :: moviename_forctf, moviename_thumb
  contains
    procedure :: iterate
    procedure :: get_moviename
end type motion_correct_iter

contains

    subroutine iterate( self, cline, ctfvars, orientation, fbody, frame_counter, moviename, dir_out, gainref_fname, tseries )
        class(motion_correct_iter), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(ctfparams),            intent(inout) :: ctfvars
        class(ori),                 intent(inout) :: orientation
        integer,                    intent(inout) :: frame_counter
        class(string),              intent(in)    :: moviename, fbody, dir_out
        class(string),    optional, intent(in)    :: gainref_fname
        character(len=*), optional, intent(in)    :: tseries
        type(nrtxtfile)               :: boxfile
        real,             allocatable :: boxdata(:,:)
        type(string) :: fbody_here, ext, star_fname, poly_fname
        real         :: goodnessoffit(2), scale, bfac_here, bid
        integer      :: ldim(3), ldim_thumb(3), iptcl, nxpatch, nypatch
        logical      :: patch_success, l_tseries
        patch_success = .false.
        l_tseries     = .false.
        if( present(tseries) ) l_tseries = tseries.eq.'yes'
        ! check, increment counter & print
        if( .not. file_exists(moviename) )then
            write(logfhandle,*) 'inputted movie stack does not exist: ', moviename%to_char()
        endif
        ! sanity check
        if( fname2format(moviename) .eq. 'K' )then
            ! eer movie
            if( (.not.params_glob%l_dose_weight) .and. (.not.cline%defined('eer_fraction')) )then
                THROW_HARD('EER_FRACTION must be defined when TOTAL_DOSE is absent!')
            endif
        endif
        ! make filenames
        fbody_here = basename(moviename)
        ext        = fname2ext(fbody_here)
        if( l_tseries )then
            fbody_here = fbody
        else if( fbody.ne.'' )then
            fbody_here = get_fbody(fbody_here, ext)
            fbody_here = fbody%to_char()//'_'//fbody_here%to_char()
        else
            fbody_here = get_fbody(fbody_here, ext)
        endif
        ! shifts & star output
        patched_shift_fname   = dir_out%to_char()//fbody_here%to_char()//'_shifts.eps'
        star_fname            = dir_out%to_char()//fbody_here%to_char()//STAR_EXT
        poly_fname            = dir_out%to_char()//fbody_here%to_char()//'.poly'
        ! isotropic ones
        self%moviename_intg   = dir_out%to_char()//fbody_here%to_char()//INTGMOV_SUFFIX//params_glob%ext%to_char()
        self%moviename_forctf = dir_out%to_char()//fbody_here%to_char()//FORCTF_SUFFIX//params_glob%ext%to_char()
        self%moviename_thumb  = dir_out%to_char()//fbody_here%to_char()//THUMBNAIL_SUFFIX//JPG_EXT
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( params_glob%nframesgrp > 0 )then
            self%moviename = 'tmpnframesgrpmovie'//params_glob%ext%to_char()
            call frameavg_stack(moviename, self%moviename, params_glob%nframesgrp, ctfvars%smpd)
        else
            self%moviename = moviename
        endif
        ! determines whether to perform patch-based step & patch size
        call calc_npatches(self%moviename, ctfvars%smpd, cline, orientation)
        motion_correct_with_patched = (params_glob%mcpatch.eq.'yes') .and. (params_glob%nxpatch*params_glob%nypatch > 1)
        bid = 0.0
        ! ALIGNEMENT
        select case(params_glob%algorithm)
        case('wpatch','poly')
            write(logfhandle,'(a,1x,a)') '>>> PROCESSING MOVIE:', moviename%to_char()
            if( cline%defined('boxfile') )then
                if( file_exists(params_glob%boxfile) )then
                    if( nlines(params_glob%boxfile) > 0 )then
                        call boxfile%new(params_glob%boxfile, 1)
                        allocate( boxdata(boxfile%get_ndatalines(),boxfile%get_nrecs_per_line()))
                        do iptcl=1,boxfile%get_ndatalines()
                            call boxfile%readNextDataLine(boxdata(iptcl,:))
                        end do
                        call boxfile%kill
                    endif
                endif
                call motion_correct_dev(self%moviename, ctfvars, self%moviesum, self%moviesum_corrected,&
                    &self%moviesum_ctf, patch_success, goodnessoffit(1), gainref=gainref_fname, boxdata=boxdata)
            else
                call motion_correct_dev(self%moviename, ctfvars, self%moviesum, self%moviesum_corrected,&
                    &self%moviesum_ctf, patch_success, goodnessoffit(1), gainref=gainref_fname)
            endif
            ! STAR output
            if( .not. l_tseries ) call motion_correct_write2star(star_fname, self%moviename, patch_success, gainref_fname)
            call motion_correct_iso_kill
            call motion_correct_kill_common
            call motion_correct_patched_kill
            if( .not.patch_success )then
                THROW_WARN('Polynomial fitting to patch-determined shifts was of insufficient quality')
                THROW_WARN('Only isotropic/stage-drift correction will be used')
            endif
            call orientation%set('gof',   goodnessoffit(1))
            call motion_correct_mic2spec(self%moviesum, GUI_PSPECSZ, speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec_sum)
        case DEFAULT
            if( trim(params_glob%algorithm) .eq. 'iso' ) motion_correct_with_patched = .false.
            ! b-factors for alignment
            bfac_here = -1.
            if( cline%defined('bfac') ) bfac_here = params_glob%bfac
            ! check, increment counter & print
            write(logfhandle,'(a,1x,a)') '>>> PROCESSING MOVIE:', moviename%to_char()
            ! execute the motion_correction
            call motion_correct_iso(self%moviename, ctfvars, bfac_here, self%moviesum, gainref_fname=gainref_fname)
            call motion_correct_mic2spec(self%moviesum, GUI_PSPECSZ, speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec_sum)
            call self%moviesum%kill
            ! shifts frames accordingly
            call motion_correct_iso_shift_frames
            ! optionally calculate optimal weights
            call motion_correct_calc_opt_weights
            ! destruct before anisotropic correction
            call motion_correct_iso_kill
            ! Patch based approach
            if( motion_correct_with_patched ) then
                select case(trim(params_glob%mcconvention))
                    case('first','relion')
                        ! the threshold is slighly increased because goodness of fit is always higher
                        ! when calculated with reference to the first frame
                        effective_patch_fit_threshold = PATCH_FIT_THRESHOLD + 1.0
                    case DEFAULT
                        effective_patch_fit_threshold = PATCH_FIT_THRESHOLD
                end select
                nxpatch = params_glob%nxpatch
                nypatch = params_glob%nypatch
                call motion_correct_patched(bfac_here, effective_patch_fit_threshold, [nxpatch, nypatch], goodnessoffit)
                if( trim(params_glob%mcpatch_thres).eq.'no' )then
                    patch_success = .true. ! always accept patch solution
                    if( any(goodnessoffit >= effective_patch_fit_threshold) )then
                        THROW_WARN('Polynomial fitting to patch-determined shifts was unsatisfactory. The patch-based correction will however be used')
                    endif
                else
                    patch_success = all(goodnessoffit < effective_patch_fit_threshold)
                    if( patch_success )then
                        ! First pass of BIM correction was successful
                    else
                        THROW_WARN('Polynomial fitting to patch-determined shifts was unsatisfactory. Retrying with less patches')
                        nxpatch = max(1,nint(real(nxpatch)/2.))
                        nypatch = max(1,nint(real(nypatch)/2.))
                        if( nxpatch * nypatch >= 2 )then
                            patched_shift_fname = dir_out%to_char()//fbody_here%to_char()//'_shifts.eps'
                            call motion_correct_patched(bfac_here, effective_patch_fit_threshold, [nxpatch, nypatch], goodnessoffit)
                            patch_success = all(goodnessoffit < effective_patch_fit_threshold)
                        endif
                    endif
                endif
                ! generate sums
                if( patch_success )then
                    call motion_correct_patched_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
                else
                    call motion_correct_iso_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
                    THROW_WARN('Polynomial fitting to patch-determined shifts was unsatisfactory. Stage-drift correction will be used')
                endif
                call orientation%set('gofx',goodnessoffit(1))
                call orientation%set('gofy',goodnessoffit(2))
                ! cleanup
                call motion_correct_patched_kill
            else
                call motion_correct_iso_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
            endif
            ! STAR output
            if( .not. l_tseries )then
                call motion_correct_write_poly(poly_fname)
                call motion_correct_write2star(star_fname, self%moviename, patch_success, gainref_fname)
                call motion_correct_calc_bid(patch_success, bid)
            endif
        end select
        ! generate power-spectra
        call motion_correct_mic2spec(self%moviesum_ctf, GUI_PSPECSZ, speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec_ctf)
        call self%pspec_sum%before_after(self%pspec_ctf, self%pspec_half_n_half)
        if( l_tseries )then
            call self%pspec_half_n_half%scale_pspec4viz(rsmpd4viz=max(SMPD4VIZ_NANO,2.*ctfvars%smpd))
        else
            call self%pspec_half_n_half%scale_pspec4viz
        endif
        call self%pspec_sum%kill
        call self%pspec_ctf%kill
        ! write output
        call self%moviesum_corrected%write(self%moviename_intg)
        if( .not. l_tseries ) call self%moviesum_ctf%write(self%moviename_forctf)
        call self%moviesum_ctf%kill
        ! generate thumbnail
        ldim            = self%moviesum_corrected%get_ldim()
        scale           = real(GUI_PSPECSZ)/maxval(ldim(1:2))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_thumb(3)   = 1
        call orientation%set('smpd', ctfvars%smpd)
        call orientation%set('xdim', real(ldim(1)))
        call orientation%set('ydim', real(ldim(2)))
        call self%thumbnail%new(ldim_thumb, ctfvars%smpd)
        call self%moviesum_corrected%fft()
        call self%moviesum_corrected%clip(self%thumbnail)
        call self%thumbnail%ifft()
        call self%moviesum_corrected%kill
        ! jpeg output
        call self%pspec_half_n_half%collage(self%thumbnail, self%img_jpg)
        call self%img_jpg%write_jpg(self%moviename_thumb, norm=.true., quality=92)
        call self%pspec_half_n_half%kill
        call self%thumbnail%kill
        call self%img_jpg%kill
        ! report to ori object
        if( .not. l_tseries )then
            call orientation%set('movie',       simple_abspath(self%moviename))
            call orientation%set('forctf',      simple_abspath(self%moviename_forctf))
            call orientation%set('mc_starfile', simple_abspath(star_fname))
            call orientation%set('bid',         bid)
        endif
        call orientation%set('intg',    simple_abspath(self%moviename_intg))
        call orientation%set('thumb',   simple_abspath(self%moviename_thumb))
        call orientation%set('imgkind', 'mic')
        if( motion_correct_with_patched ) call orientation%set('mceps', simple_abspath(patched_shift_fname))
        call motion_correct_kill_common
    end subroutine iterate

    subroutine calc_npatches( moviename, smpd, cline, o )
        class(string),  intent(in)    :: moviename
        real,           intent(in)    :: smpd
        class(cmdline), intent(inout) :: cline
        class(ori),     intent(inout) :: o
        integer :: patchsz(2),pixsz(2),nx,ny
        real    :: physz(2)
        if( trim(params_glob%mcpatch).eq.'yes' )then
            pixsz = nint([o%get('xdim'),o%get('ydim')])     ! micrograph pixel dimensions
            physz = real(pixsz) * smpd                      ! stage size in Angs
            nx    = max(1, floor(physz(1) / MC_PATCHSZ))    ! # of patches along x
            ny    = max(1, floor(physz(2) / MC_PATCHSZ))    ! # of patches along y
            ! Micrograph scaled pixel dimensions
            pixsz = nint(real(pixsz) * params_glob%scale_movies)
            if( fname2format(moviename) .eq. 'K' ) pixsz = pixsz * params_glob%eer_upsampling
            ! Patch size in pixel
            patchsz(1) = nint(real(pixsz(1)) / real(nx))
            patchsz(2) = nint(real(pixsz(2)) / real(ny))
            ! Minimum patch size
            patchsz(1) = max(patchsz(1), MC_MINPATCHSZ)
            patchsz(2) = max(patchsz(2), MC_MINPATCHSZ)
            ! Effective number of patches
            nx = max(1, floor(real(pixsz(1)) / real(patchsz(1))))
            ny = max(1, floor(real(pixsz(2)) / real(patchsz(2))))
            if( .not.cline%defined('nxpatch') ) params_glob%nxpatch = nx
            if( .not.cline%defined('nypatch') ) params_glob%nypatch = ny
        else
            params_glob%nxpatch = 1
            params_glob%nypatch = 1
        endif
        call o%set('nxpatch', real(params_glob%nxpatch))
        call o%set('nypatch', real(params_glob%nypatch))
    end subroutine calc_npatches

    function get_moviename( self, which ) result( moviename )
        class(motion_correct_iter), intent(in) :: self
        character(len=*),           intent(in) :: which
        type(string) :: moviename
        select case( which )
            case('intg')
                moviename = self%moviename_intg
            case('forctf')
                moviename = self%moviename_forctf
            case('thumb')
                moviename = self%moviename_thumb
            case DEFAULT
                THROW_HARD('unsupported which flag; get_self%moviename')
        end select
    end function get_moviename

end module  simple_motion_correct_iter
