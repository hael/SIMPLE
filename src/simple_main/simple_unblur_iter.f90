module simple_unblur_iter
use simple_image,   only: image
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_strings, only: int2str_pad, int2str
use simple_math,    only: round2even
use simple_ori,     only: ori
use simple_filehandling ! use all in there
use simple_procimgfile  ! use all in there
use simple_unblur       ! use all in there
implicit none

public :: unblur_iter
private

logical, parameter :: DEBUG = .false.

type :: unblur_iter
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
end type unblur_iter

contains

    subroutine iterate( self, cline, p, orientation, imovie, movie_counter, frame_counter, moviename, smpd_out )
        class(unblur_iter), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(params),      intent(inout) :: p
        class(ori),         intent(inout) :: orientation
        integer,            intent(in)    :: imovie
        integer,            intent(inout) :: movie_counter, frame_counter
        character(len=*),   intent(in)    :: moviename
        real,               intent(out)   :: smpd_out
        integer :: ldim(3), ldim_thumb(3)
        real    :: corr, scale
        ! make names
        if( cline%defined('fbody') )then
            call orientation%set('movie', trim(moviename))
            self%moviename_intg   = trim(adjustl(p%fbody))//'_intg'//int2str_pad(imovie,p%numlen)//p%ext
            call orientation%set('intg', trim(self%moviename_intg))
            self%moviename_forctf = trim(adjustl(p%fbody))//'_forctf'//int2str_pad(imovie,p%numlen)//p%ext
            call orientation%set('forctf', trim(self%moviename_forctf))
            self%moviename_pspec  = trim(adjustl(p%fbody))//'_pspec'//int2str_pad(imovie,p%numlen)//p%ext
            call orientation%set('pspec', trim(self%moviename_pspec))
            self%moviename_thumb  = trim(adjustl(p%fbody))//'_thumb'//int2str_pad(imovie,p%numlen)//p%ext
            call orientation%set('thumb', trim(self%moviename_thumb))
            if( cline%defined('tof') )then
                self%moviename_intg_frames = trim(adjustl(p%fbody))//'_frames'//int2str(p%fromf)//'-'&
                &//int2str(p%tof)//'_intg'//int2str_pad(imovie,p%numlen)//p%ext
                call orientation%set('intg_frames', trim(self%moviename_intg_frames))
            endif
        else
            call orientation%set('movie', trim(moviename))
            self%moviename_intg   = int2str_pad(imovie,p%numlen)//'_intg'//p%ext
            call orientation%set('intg', trim(self%moviename_intg))
            self%moviename_forctf = int2str_pad(imovie,p%numlen)//'_forctf'//p%ext
            call orientation%set('forctf', trim(self%moviename_forctf))
            self%moviename_pspec  = int2str_pad(imovie,p%numlen)//'_pspec'//p%ext
            call orientation%set('pspec', trim(self%moviename_pspec))
            self%moviename_thumb  = int2str_pad(imovie,p%numlen)//'_thumb'//p%ext
            call orientation%set('thumb', trim(self%moviename_thumb))
            if( cline%defined('tof') )then
                self%moviename_intg_frames = int2str_pad(imovie,p%numlen)//'_frames'//int2str(p%fromf)//'-'&
                &//int2str(p%tof)//'_intg'//p%ext
                call orientation%set('intg_frames', trim(self%moviename_intg_frames))
            endif
        endif
        ! check, increment counter & print
        if( .not. file_exists(moviename) )then
            write(*,*) 'inputted movie stack does not exist: ', moviename
        endif
        movie_counter = movie_counter + 1
        write(*,'(a,1x,i5)') '>>> PROCESSING MOVIE:', imovie
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( p%nframesgrp > 0 )then
            self%moviename = 'tmpnframesgrpmovie'//p%ext
            call frameavg_imgfile(trim(moviename), trim(self%moviename), p%nframesgrp, p%smpd)
        else
            self%moviename = trim(moviename)
        endif
        ! execute the unblurring
        call unblur_movie(self%moviename, p, corr, smpd_out)
        call orientation%set('smpd', smpd_out)
        if( p%tomo .eq. 'yes' )then
            call unblur_calc_sums_tomo(frame_counter, p%time_per_frame,&
            &self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
        else
            if( cline%defined('tof') )then
                call unblur_calc_sums(self%moviesum_corrected_frames, [p%fromf,p%tof])
                call unblur_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            else
                call unblur_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
            endif
            if( DEBUG )then
                print *, 'ldim(moviesum):           ', self%moviesum%get_ldim()
                print *, 'ldim(moviesum_corrected): ', self%moviesum_corrected%get_ldim()
                print *, 'ldim(moviesum_ctf):       ', self%moviesum_ctf%get_ldim()
            endif
        endif
        ! generate power-spectra
        self%pspec_sum         = self%moviesum%mic2spec(p%pspecsz, self%speckind)
        self%pspec_ctf         = self%moviesum_ctf%mic2spec(p%pspecsz, self%speckind)
        self%pspec_half_n_half = self%pspec_sum%before_after(self%pspec_ctf)
        ! write output
        if( cline%defined('tof') ) call self%moviesum_corrected_frames%write(self%moviename_intg_frames)
        call self%moviesum_corrected%write(self%moviename_intg)
        call self%moviesum_ctf%write(self%moviename_forctf)
        call self%pspec_half_n_half%write(self%moviename_pspec)
        ! generate thumbnail
        ldim          = self%moviesum_corrected%get_ldim()
        scale         = real(p%pspecsz)/real(ldim(1))
        ldim_thumb(1) = round2even(real(ldim(1))*scale)
        ldim_thumb(2) = round2even(real(ldim(2))*scale)
        ldim_thumb(3) = 1
        call self%thumbnail%new(ldim_thumb, p%smpd)
        call self%moviesum_corrected%fwd_ft
        call self%moviesum_corrected%clip(self%thumbnail)
        call self%thumbnail%bwd_ft
        call self%thumbnail%write(self%moviename_thumb)
        ! destruct
        call self%moviesum%kill
        call self%moviesum_corrected%kill
        call self%moviesum_corrected_frames%kill
        call self%moviesum_ctf%kill
        call self%pspec_sum%kill
        call self%pspec_ctf%kill
        call self%pspec_half_n_half%kill
    end subroutine iterate

    function get_moviename( self, which ) result( moviename )
        class(unblur_iter), intent(in) :: self
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
                stop 'unsupported which flag; simple_unblur_iter :: get_moviename'
        end select
    end function get_moviename

end module  simple_unblur_iter