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

type :: motion_correct_iter
    private
    character(len=4)      :: speckind = 'sqrt'
    character(len=STDLEN) :: moviename, moviename_intg, moviename_intg_frames
    character(len=STDLEN) :: moviename_forctf, moviename_thumb, moviename_pspec
    type(image)           :: moviesum, moviesum_corrected, moviesum_corrected_frames
    type(image)           :: moviesum_ctf, pspec_sum, pspec_ctf
    type(image)           :: pspec_half_n_half, thumbnail
  contains
    procedure :: iterate
    procedure :: get_moviename
end type motion_correct_iter

contains

    subroutine iterate( self, cline, ctfvars, orientation, fbody, frame_counter, moviename, dir_out, gainref_fname )
        class(motion_correct_iter), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(ctfparams),            intent(inout) :: ctfvars
        class(ori),                 intent(inout) :: orientation
        integer,                    intent(inout) :: frame_counter
        character(len=*),           intent(in)    :: moviename, fbody, dir_out
        character(len=*), optional, intent(in)    :: gainref_fname
        character(len=:), allocatable :: fbody_here, ext, fname
        real,             allocatable :: shifts(:,:)
        type(stats_struct) :: shstats(2)
        type(image)        :: img_jpg
        integer            :: ldim(3), ldim_thumb(3)
        real               :: corr, scale, var
        logical            :: err
        ! check, increment counter & print
        if( .not. file_exists(moviename) )then
            write(*,*) 'inputted movie stack does not exist: ', moviename
        endif
        ! make filenames
        fbody_here = basename(trim(moviename))
        ext        = fname2ext(trim(fbody_here))
        if( fbody.ne.'' )then
            fbody_here = trim(fbody)//'_'//get_fbody(trim(fbody_here), trim(ext))
        else
            fbody_here = get_fbody(trim(fbody_here), trim(ext))
        endif
        self%moviename_intg   = trim(dir_out)//trim(adjustl(fbody_here))//INTGMOV_SUFFIX//trim(params_glob%ext)
        self%moviename_forctf = trim(dir_out)//trim(adjustl(fbody_here))//FORCTF_SUFFIX//trim(params_glob%ext)
        self%moviename_pspec  = trim(dir_out)//trim(adjustl(fbody_here))//POWSPEC_SUFFIX//trim(params_glob%ext)
        self%moviename_thumb  = trim(dir_out)//trim(adjustl(fbody_here))//THUMBNAIL_SUFFIX//trim(JPG_EXT)
        if( cline%defined('tof') )then
            self%moviename_intg_frames = trim(dir_out)//trim(adjustl(fbody_here))//'_frames'//int2str(params_glob%fromf)//'-'&
            &//int2str(params_glob%tof)//INTGMOV_SUFFIX//params_glob%ext
        endif
        ! check, increment counter & print
        write(*,'(a,1x,a)') '>>> PROCESSING MOVIE:', trim(moviename)
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( params_glob%nframesgrp > 0 )then
            self%moviename = 'tmpnframesgrpmovie'//params_glob%ext
            call frameavg_stack(trim(moviename), trim(self%moviename), params_glob%nframesgrp, ctfvars%smpd)
        else
            self%moviename = trim(moviename)
        endif
        ! execute the motion_correction
        call motion_correct_movie(self%moviename, ctfvars, corr, shifts, err, gainref_fname)
        if( err ) return
        ! return shift stats
        call moment(shifts(1,:), shstats(1)%avg, shstats(1)%sdev, var, err)
        call moment(shifts(2,:), shstats(2)%avg, shstats(2)%sdev, var, err)
        call orientation%set('xavg',  shstats(1)%avg)
        call orientation%set('xsdev', shstats(1)%sdev)
        call orientation%set('xmax',  maxval(shifts(1,:)))
        call orientation%set('xmin',  minval(shifts(1,:)))
        call orientation%set('yavg',  shstats(2)%avg)
        call orientation%set('ysdev', shstats(2)%sdev)
        call orientation%set('ymax',  maxval(shifts(2,:)))
        call orientation%set('ymin',  minval(shifts(2,:)))
        ! generate sums
        if( params_glob%tomo .eq. 'yes' )then
            call motion_correct_calc_sums_tomo(frame_counter, params_glob%time_per_frame,&
            &self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
        else
            if( cline%defined('tof') )then
                call motion_correct_calc_sums(self%moviesum_corrected_frames, [params_glob%fromf,params_glob%tof])
                call motion_correct_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            else
                call motion_correct_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            endif
            DebugPrint 'ldim(moviesum):           ', self%moviesum%get_ldim()
            DebugPrint 'ldim(moviesum_corrected): ', self%moviesum_corrected%get_ldim()
            DebugPrint 'ldim(moviesum_ctf):       ', self%moviesum_ctf%get_ldim()
        endif
        ! generate power-spectra and cleanup
        self%pspec_sum = self%moviesum%mic2spec(params_glob%pspecsz, self%speckind)
        self%pspec_ctf = self%moviesum_ctf%mic2spec(params_glob%pspecsz, self%speckind)
        call self%pspec_sum%dampen_central_cross
        call self%pspec_ctf%dampen_central_cross
        self%pspec_half_n_half = self%pspec_sum%before_after(self%pspec_ctf)
        ! write output
        if( cline%defined('tof') ) call self%moviesum_corrected_frames%write(self%moviename_intg_frames)
        call self%moviesum_corrected%write(self%moviename_intg)
        call self%moviesum_ctf%write(self%moviename_forctf)
        call self%pspec_ctf%write(self%moviename_pspec)
        ! generate thumbnail
        ldim          = self%moviesum_corrected%get_ldim()
        scale         = real(params_glob%pspecsz)/real(ldim(1))
        ldim_thumb(1) = round2even(real(ldim(1))*scale)
        ldim_thumb(2) = round2even(real(ldim(2))*scale)
        ldim_thumb(3) = 1
        call self%thumbnail%new(ldim_thumb, ctfvars%smpd)
        call self%moviesum_corrected%fft()
        call self%moviesum_corrected%clip(self%thumbnail)
        call self%thumbnail%ifft()
        ! jpeg output
        call self%pspec_half_n_half%collage(self%thumbnail, img_jpg)
        call img_jpg%write_jpg(self%moviename_thumb, quality=90)
        ! report to ori object
        call orientation%set('smpd',   ctfvars%smpd)
        fname = simple_abspath(moviename,  errmsg='simple_motion_correct_iter::iterate 1')
        call orientation%set('movie',  trim(fname))
        fname = simple_abspath(self%moviename_intg, errmsg='simple_motion_correct_iter::iterate 2')
        call orientation%set('intg',   trim(fname))
        fname = simple_abspath(self%moviename_forctf, errmsg='simple_motion_correct_iter::iterate 3')
        call orientation%set('forctf', trim(fname))
        fname = simple_abspath(self%moviename_thumb, errmsg='simple_motion_correct_iter::iterate 4')
        call orientation%set('thumb',  trim(fname))
        call orientation%set('imgkind', 'mic')
        if( cline%defined('tof') )then
            fname = simple_abspath(self%moviename_intg_frames, errmsg='simple_motion_correct_iter::iterate 5')
            call orientation%set('intg_frames', trim(self%moviename_intg_frames))
        endif
        ! destruct
        call self%pspec_sum%kill
        call self%pspec_ctf%kill
        call self%moviesum%kill
        call self%moviesum_corrected%kill
        call self%moviesum_corrected_frames%kill
        call self%moviesum_ctf%kill
        call self%pspec_half_n_half%kill
        call img_jpg%kill
        call self%thumbnail%kill
        deallocate(shifts)
    end subroutine iterate

    function get_moviename( self, which ) result( moviename )
        class(motion_correct_iter), intent(in) :: self
        character(len=*),   intent(in) :: which
        character(len=:), allocatable  :: moviename
        select case( which )
            case('intg')
                allocate(moviename, source=trim(self%moviename_intg))
            case('intg_frames')
                allocate(moviename, source=trim(self%moviename_intg_frames))
            case('forctf')
                allocate(moviename, source=trim(self%moviename_forctf))
            case('thumb')
                allocate(moviename, source=trim(self%moviename_thumb))
            case DEFAULT
                THROW_HARD('unsupported which flag; get_moviename')
        end select
    end function get_moviename

end module  simple_motion_correct_iter
