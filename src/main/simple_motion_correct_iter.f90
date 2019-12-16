! iterator for motion_correct (a program for motion correction, dose-weighting and frame-weighting of DDD movies)
module simple_motion_correct_iter
include 'simple_lib.f08'
use simple_image,        only: image
use simple_cmdline,      only: cmdline
use simple_parameters,   only: params_glob
use simple_ori,          only: ori
use simple_stackops,     only: frameavg_stack
use simple_motion_correct
implicit none

public :: motion_correct_iter
private
#include "simple_local_flags.inc"

real,             parameter :: PATCH_FIT_THRESHOLD = 4.0 ! threshold for polynomial fitting in pixels
character(len=*), parameter :: speckind = 'sqrt'
! benchmarking
logical                 :: L_BENCH = .false.
integer(timer_int_kind) :: t_postproc1
real(timer_int_kind)    :: rt_postproc1

type :: motion_correct_iter
    private
    ! these image objects are part of the instance to avoid excessive memory re-allocations
    type(image) :: moviesum, moviesum_corrected
    type(image) :: moviesum_ctf, pspec_sum, pspec_ctf, img_jpg
    type(image) :: pspec_half_n_half, thumbnail
    ! these strings are part of the instance for reporting purposes
    character(len=STDLEN) :: moviename, moviename_intg
    character(len=STDLEN) :: moviename_forctf, moviename_thumb
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
        character(len=*),           intent(in)    :: moviename, fbody, dir_out
        character(len=*), optional, intent(in)    :: gainref_fname
        character(len=*), optional, intent(in)    :: tseries
        character(len=:), allocatable :: fbody_here, ext, star_fname
        type(stats_struct)        :: shstats(2)
        character(len=LONGSTRLEN) :: rel_fname
        real    :: goodnessoffit(2), scale, var, bfac_here
        integer :: ldim(3), ldim_thumb(3)
        logical :: err, patch_success, l_tseries
        l_tseries = .false.
        if( present(tseries) ) l_tseries = tseries.eq.'yes'
        ! check, increment counter & print
        if( .not. file_exists(moviename) )then
            write(logfhandle,*) 'inputted movie stack does not exist: ', trim(moviename)
        endif
        ! make filenames
        fbody_here = basename(trim(moviename))
        ext        = fname2ext(trim(fbody_here))
        if( l_tseries )then
            fbody_here = trim(fbody)
        else if( fbody.ne.'' )then
            fbody_here = trim(fbody)//'_'//get_fbody(trim(fbody_here), trim(ext))
        else
            fbody_here = get_fbody(trim(fbody_here), trim(ext))
        endif
        ! star output
        star_fname = trim(dir_out)//trim(adjustl(fbody_here))//'.star'
        ! isotropic ones
        self%moviename_intg   = trim(dir_out)//trim(adjustl(fbody_here))//INTGMOV_SUFFIX//trim(params_glob%ext)
        self%moviename_forctf = trim(dir_out)//trim(adjustl(fbody_here))//FORCTF_SUFFIX//trim(params_glob%ext)
        self%moviename_thumb  = trim(dir_out)//trim(adjustl(fbody_here))//THUMBNAIL_SUFFIX//trim(JPG_EXT)
        ! determines whether to perform patch-based step & patch size
        if( params_glob%mcpatch.eq.'yes' )then
            if( .not.cline%defined('nxpatch') )then
                params_glob%nxpatch = min(MC_NPATCH, max(1,floor(orientation%get('xdim')*params_glob%scale/MC_PATCHSZ)) )
            endif
            if( .not.cline%defined('nypatch') )then
                params_glob%nypatch = min(MC_NPATCH, max(1,floor(orientation%get('ydim')*params_glob%scale/MC_PATCHSZ)) )
            endif
        endif
        motion_correct_with_patched = (params_glob%mcpatch.eq.'yes') .and. (params_glob%nxpatch*params_glob%nypatch > 1)
        if( params_glob%tomo.eq.'yes' ) motion_correct_with_patched = .false.
        ! b-factors for alignment
        bfac_here = -1.
        if( cline%defined('bfac') ) bfac_here = params_glob%bfac
        ! check, increment counter & print
        write(logfhandle,'(a,1x,a)') '>>> PROCESSING MOVIE:', trim(moviename)
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( params_glob%nframesgrp > 0 )then
            self%moviename = 'tmpnframesgrpmovie'//params_glob%ext
            call frameavg_stack(trim(moviename), trim(self%moviename), params_glob%nframesgrp, ctfvars%smpd)
        else
            self%moviename = trim(moviename)
        endif
        ! execute the motion_correction
        call motion_correct_iso(self%moviename, ctfvars, bfac_here, self%moviesum, gainref_fname=gainref_fname)
        call motion_correct_mic2spec(self%moviesum, params_glob%pspecsz, speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec_sum)
        call self%moviesum%kill
        ! shifts frames accordingly
        call motion_correct_iso_shift_frames
        ! optionally calculate optimal weights
        call motion_correct_calc_opt_weights
        ! destruct before anisotropic correction
        call motion_correct_iso_kill
        if( .not. l_tseries ) call write_iso2star(star_fname, self%moviename, gainref_fname)
        ! Patch based approach
        if( motion_correct_with_patched ) then
            patched_shift_fname = trim(dir_out)//trim(adjustl(fbody_here))//'_shifts.eps'
            call motion_correct_patched( bfac_here, goodnessoffit )
            if( .not. l_tseries ) call write_aniso2star
            patch_success = all(goodnessoffit < PATCH_FIT_THRESHOLD)
            ! generate sums
            if( patch_success )then
                call motion_correct_patched_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
            else
                if( params_glob%tomo .eq. 'yes' )then
                    call motion_correct_iso_calc_sums_tomo(frame_counter, params_glob%time_per_frame,&
                    &self%moviesum_corrected, self%moviesum_ctf)
                else
                    call motion_correct_iso_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
                endif
                THROW_WARN('Polynomial fitting to patch-determined shifts was of insufficient quality')
                THROW_WARN('Only isotropic/stage-drift correction will be used')
            endif
            call orientation%set('gofx',goodnessoffit(1))
            call orientation%set('gofy',goodnessoffit(2))
            ! cleanup
            call motion_correct_patched_kill
        else
            ! generate sums
            if( params_glob%tomo .eq. 'yes' )then
                call motion_correct_iso_calc_sums_tomo(frame_counter, params_glob%time_per_frame,&
                &self%moviesum_corrected, self%moviesum_ctf)
            else
                call motion_correct_iso_calc_sums(self%moviesum_corrected, self%moviesum_ctf)
            endif
        endif
        ! generate power-spectra
        if( L_BENCH ) t_postproc1 = tic()
        call motion_correct_mic2spec(self%moviesum_ctf, params_glob%pspecsz, speckind, LP_PSPEC_BACKGR_SUBTR, self%pspec_ctf)
        call self%pspec_sum%before_after(self%pspec_ctf, self%pspec_half_n_half)
        call self%pspec_half_n_half%scale_pspec4viz
        ! write output
        call self%moviesum_corrected%write(self%moviename_intg)
        if( .not. l_tseries ) call self%moviesum_ctf%write(self%moviename_forctf)
        if( L_BENCH ) rt_postproc1 = toc(t_postproc1)
        ! generate thumbnail
        ldim  = self%moviesum_corrected%get_ldim()
        scale = real(params_glob%pspecsz)/real(ldim(1))
        ldim_thumb(1:2) = round2even(real(ldim(1:2))*scale)
        ldim_thumb(3)   = 1
        call orientation%set('xdim', real(ldim(1)))
        call orientation%set('ydim', real(ldim(2)))
        call self%thumbnail%new(ldim_thumb, ctfvars%smpd)
        call self%moviesum_corrected%fft()
        call self%moviesum_corrected%clip(self%thumbnail)
        call self%thumbnail%ifft()
        ! jpeg output
        call self%pspec_half_n_half%collage(self%thumbnail, self%img_jpg)
        call self%img_jpg%write_jpg(self%moviename_thumb, norm=.true., quality=90)
        ! report to ori object
        call orientation%set('smpd',   ctfvars%smpd)
        if( .not. l_tseries )then
            call make_relativepath(CWD_GLOB,self%moviename,rel_fname)
            call orientation%set('movie',  trim(rel_fname))
            call make_relativepath(CWD_GLOB,self%moviename_forctf,rel_fname)
            call orientation%set('forctf', trim(rel_fname))
            call make_relativepath(CWD_GLOB,star_fname,rel_fname)
            call orientation%set("mc_starfile",rel_fname)
        endif
        call make_relativepath(CWD_GLOB,self%moviename_intg,rel_fname)
        call orientation%set('intg',   trim(rel_fname))
        call make_relativepath(CWD_GLOB,self%moviename_thumb,rel_fname)
        call orientation%set('thumb',  trim(rel_fname))
        call orientation%set('imgkind', 'mic')
        if( motion_correct_with_patched )then
            call make_relativepath(CWD_GLOB, patched_shift_fname, rel_fname)
            call orientation%set('mceps', rel_fname)
        endif
        call motion_correct_kill_common
        if( .not. l_tseries ) call close_starfile
        if( L_BENCH )then
            print *,'rt_postproc1: ',rt_postproc1
        endif
    end subroutine iterate

    function get_moviename( self, which ) result( moviename )
        class(motion_correct_iter), intent(in) :: self
        character(len=*),           intent(in) :: which
        character(len=:), allocatable  :: moviename
        select case( which )
            case('intg')
                allocate(moviename, source=trim(self%moviename_intg))
            case('forctf')
                allocate(moviename, source=trim(self%moviename_forctf))
            case('thumb')
                allocate(moviename, source=trim(self%moviename_thumb))
            case DEFAULT
                THROW_HARD('unsupported which flag; get_self%moviename')
        end select
    end function get_moviename

end module  simple_motion_correct_iter
