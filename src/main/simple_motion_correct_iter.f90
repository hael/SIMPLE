! iterator for motion_correct (a program for motion correction, dose-weighting and frame-weighting of DDD movies)
module simple_motion_correct_iter
include 'simple_lib.f08'
use simple_image,        only: image
use simple_cmdline,      only: cmdline
use simple_params,       only: params
use simple_ori,          only: ori
use simple_stackops,     only: frameavg_stack
use simple_motion_correct       ! use all in there
implicit none

public :: motion_correct_iter
private
#include "simple_local_flags.inc"

type :: motion_correct_iter
    private
    character(len=4)      :: speckind = 'sqrt'
    character(len=STDLEN) :: moviename, moviename_intg, moviename_intg_frames
    character(len=STDLEN) :: moviename_forctf, moviename_pspec
    character(len=STDLEN) :: moviename_thumb
    type(image)           :: moviesum, moviesum_corrected, moviesum_corrected_frames
    type(image)           :: moviesum_ctf, pspec_sum, pspec_ctf
    type(image)           :: pspec_half_n_half, thumbnail
  contains
    procedure :: iterate
    procedure :: get_moviename
end type motion_correct_iter

contains

    subroutine iterate( self, cline, p, orientation, fbody, frame_counter, moviename, smpd_out, dir_out )
        class(motion_correct_iter), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(params),              intent(inout) :: p
        class(ori),                 intent(inout) :: orientation
        integer,                    intent(inout) :: frame_counter
        character(len=*),           intent(in)    :: moviename, fbody
        real,                       intent(out)   :: smpd_out
        character(len=*),           intent(in)    :: dir_out
        character(len=:), allocatable :: fbody_here, ext, fname
        real,             allocatable :: shifts(:,:)
        type(image) :: img_jpg
        integer     :: ldim(3), ldim_thumb(3), status
        real        :: corr, scale
        logical     :: err
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
        self%moviename_intg   = trim(dir_out)//trim(adjustl(fbody_here))//'_intg'//trim(p%ext)
        self%moviename_forctf = trim(dir_out)//trim(adjustl(fbody_here))//'_forctf'//trim(p%ext)
        self%moviename_pspec  = trim(dir_out)//trim(adjustl(fbody_here))//'_pspec'//trim(JPG_EXT)
        self%moviename_thumb  = trim(dir_out)//trim(adjustl(fbody_here))//'_thumb'//trim(JPG_EXT)
        if( cline%defined('tof') )then
            self%moviename_intg_frames = trim(dir_out)//trim(adjustl(fbody_here))//'_frames'//int2str(p%fromf)//'-'&
            &//int2str(p%tof)//'_intg'//p%ext
        endif
        ! check, increment counter & print
        write(*,'(a,1x,a)') '>>> PROCESSING MOVIE:', trim(moviename)
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( p%nframesgrp > 0 )then
            self%moviename = 'tmpnframesgrpmovie'//p%ext
            call frameavg_stack(trim(moviename), trim(self%moviename), p%nframesgrp, p%smpd)
        else
            self%moviename = trim(moviename)
        endif
        ! execute the motion_correction
        call motion_correct_movie(self%moviename, p, corr, smpd_out, shifts, err)
        if( err ) return
        ! generate sums
        if( p%tomo .eq. 'yes' )then
            call motion_correct_calc_sums_tomo(frame_counter, p%time_per_frame,&
            &self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
        else
            if( cline%defined('tof') )then
                call motion_correct_calc_sums(self%moviesum_corrected_frames, [p%fromf,p%tof])
                call motion_correct_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            else
                call motion_correct_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            endif
            DebugPrint 'ldim(moviesum):           ', self%moviesum%get_ldim()
            DebugPrint 'ldim(moviesum_corrected): ', self%moviesum_corrected%get_ldim()
            DebugPrint 'ldim(moviesum_ctf):       ', self%moviesum_ctf%get_ldim()
        endif
        ! generate power-spectra and cleanup
        self%pspec_sum         = self%moviesum%mic2spec(p%pspecsz, self%speckind)
        self%pspec_ctf         = self%moviesum_ctf%mic2spec(p%pspecsz, self%speckind)
        self%pspec_half_n_half = self%pspec_sum%before_after(self%pspec_ctf)
        call self%pspec_sum%kill
        call self%pspec_ctf%kill
        ! write output
        if( cline%defined('tof') ) call self%moviesum_corrected_frames%write(self%moviename_intg_frames)
        call self%moviesum_corrected%write(self%moviename_intg)
        call self%moviesum_ctf%write(self%moviename_forctf)
        ! generate thumbnail
        ldim          = self%moviesum_corrected%get_ldim()
        scale         = real(p%pspecsz)/real(ldim(1))
        ldim_thumb(1) = round2even(real(ldim(1))*scale)
        ldim_thumb(2) = round2even(real(ldim(2))*scale)
        ldim_thumb(3) = 1
        call self%thumbnail%new(ldim_thumb, p%smpd)
        call self%moviesum_corrected%fft()
        call self%moviesum_corrected%clip(self%thumbnail)
        call self%thumbnail%ifft()
        ! jpeg output
        call self%pspec_half_n_half%collage(self%thumbnail, img_jpg)
        call img_jpg%write_jpg(self%moviename_thumb, quality=90, norm=.true.)
        ! report to ori object
        call orientation%set('smpd',   smpd_out)
        call simple_full_path(moviename, fname, 'simple_motion_correct_iter::iterate')
        call orientation%set('movie',  trim(fname))
        call simple_full_path(self%moviename_intg, fname,'simple_motion_correct_iter::iterate' )
        call orientation%set('intg',   trim(fname))
        call simple_full_path(self%moviename_forctf, fname, 'simple_motion_correct_iter::iterate')
        call orientation%set('forctf', trim(fname))
        call simple_full_path(self%moviename_pspec, fname,'simple_motion_correct_iter::iterate')
        call orientation%set('pspec',  trim(fname))
        call simple_full_path(self%moviename_thumb, fname, 'simple_motion_correct_iter::iterate')
        call orientation%set('thumb',  trim(fname))
        call orientation%set('imgkind', 'mic')
        if( cline%defined('tof') )then
            call simple_full_path(self%moviename_intg_frames, fname, 'simple_motion_correct_iter::iterate')
            call orientation%set('intg_frames', trim(self%moviename_intg_frames))
        endif
        ! destruct
        call self%moviesum%kill
        call self%moviesum_corrected%kill
        call self%moviesum_corrected_frames%kill
        call self%moviesum_ctf%kill
        call self%pspec_half_n_half%kill
        call img_jpg%kill
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
            case('pspec')
                allocate(moviename, source=trim(self%moviename_pspec))
            case('thumb')
                allocate(moviename, source=trim(self%moviename_thumb))
            case DEFAULT
                call simple_stop('unsupported which flag; simple_motion_correct_iter :: get_moviename')
        end select
    end function get_moviename

end module  simple_motion_correct_iter
