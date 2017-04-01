module simple_unblur_iter
use simple_image, only: image
use simple_filehandling ! use all in there
use simple_procimgfile  ! use all in there
use simple_unblur       ! use all in there
implicit none

public :: unblur_iter
private

logical, parameter :: DEBUG = .false.
real               :: native_smpd

type :: unblur_iter
    private
    character(len=4)      :: speckind = 'sqrt'
    character(len=STDLEN) :: moviename, moviename_intg
    character(len=STDLEN) :: moviename_forctf, moviename_pspec
    character(len=STDLEN) :: moviename_thumb
    type(image)           :: moviesum, moviesum_corrected
    type(image)           :: moviesum_ctf, pspec_sum, pspec_ctf
    type(image)           :: pspec_half_n_half, thumbnail
  contains
    procedure :: iterate
    procedure :: get_moviename
end type unblur_iter

contains

    subroutine iterate( self, cline, p, imovie, movie_counter, frame_counter, moviename )
        use simple_cmdline, only: cmdline
        use simple_params,  only: params
        use simple_strings, only: int2str_pad
        use simple_math,    only: round2even
        class(unblur_iter),         intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(params),              intent(inout) :: p
        integer,                    intent(in)    :: imovie
        integer,                    intent(inout) :: movie_counter, frame_counter
        character(len=*),           intent(in)    :: moviename
        integer :: ldim(3), ldim_thumb(3)
        real    :: corr, scale
        native_smpd = p%smpd
        if( cline%defined('scale') )then
            p%smpd = native_smpd
            call cline%set('scale', native_smpd)
        endif
        ! make names
        if( cline%defined('fbody') )then
            self%moviename_intg   = trim(adjustl(p%fbody))//'_intg'//int2str_pad(imovie,p%numlen)//p%ext
            self%moviename_forctf = trim(adjustl(p%fbody))//'_forctf'//int2str_pad(imovie,p%numlen)//p%ext
            self%moviename_pspec  = trim(adjustl(p%fbody))//'_pspec'//int2str_pad(imovie,p%numlen)//p%ext
            self%moviename_thumb  = trim(adjustl(p%fbody))//'_thumb'//int2str_pad(imovie,p%numlen)//p%ext
        else
            self%moviename_intg   = int2str_pad(imovie,p%numlen)//'_intg'//p%ext
            self%moviename_forctf = int2str_pad(imovie,p%numlen)//'_forctf'//p%ext
            self%moviename_pspec  = int2str_pad(imovie,p%numlen)//'_pspec'//p%ext
            self%moviename_thumb  = int2str_pad(imovie,p%numlen)//'_thumb'//p%ext
        endif
        ! check, increment counter & print
        if( .not. file_exists(moviename) )then
            write(*,*) 'inputted movie stack does not exist: ', moviename
        endif
        movie_counter = movie_counter + 1
        write(*,'(a,1x,i5)') '>>> PROCESSING MOVIE:', imovie
        ! averages frames as a pre-processing step (Falcon 3 with long exposures)
        if( p%frameavg > 0 )then
            self%moviename = 'tmpframeavgmovie'//p%ext
            call frameavg_imgfile(trim(moviename), trim(self%moviename), p%frameavg, p%smpd)
        else
            self%moviename = trim(moviename)
        endif
        ! execute the unblurring
        call unblur_movie(self%moviename, p, corr)
        if( p%tomo .eq. 'yes' )then
            call unblur_calc_sums_tomo(frame_counter, p%time_per_frame, self%moviesum, self%moviesum_corrected, self%moviesum_ctf)
        else
            if( cline%defined('fromf') )then
                if( .not. cline%defined('tof') ) stop 'need tof to be part of the command line in conjunction with fromf;&
                & simple_unblur_iter :: iterate'
                call unblur_calc_sums(self%moviesum, self%moviesum_corrected, self%moviesum_ctf, fromto=[p%fromf,p%tof])
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
        call self%moviesum_ctf%kill
        call self%pspec_sum%kill
        call self%pspec_ctf%kill
        call self%pspec_half_n_half%kill
        ! update sampling distance if scaling
        if( cline%defined('scale') )then
            p%smpd = p%smpd/p%scale
            call cline%set('scale', p%smpd)
        endif
    end subroutine iterate

    function get_moviename( self, which ) result( moviename )
        class(unblur_iter), intent(in) :: self
        character(len=*),   intent(in) :: which
        character(len=:), allocatable  :: moviename  
        select case( which )
            case('intg')
                allocate(moviename, source=trim(self%moviename_intg))
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