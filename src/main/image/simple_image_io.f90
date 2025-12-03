submodule (simple_image) simple_image_io
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
use simple_imgfile, only: imgfile
use gnufor2
implicit none
contains

    !===========================
    ! open
    !===========================
    module subroutine open( self, fname, ioimg, formatchar, readhead, rwaction )
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        class(imgfile),             intent(inout) :: ioimg
        character(len=1), optional, intent(in)    :: formatchar
        logical,          optional, intent(in)    :: readhead
        character(len=*), optional, intent(in)    :: rwaction
        character(len=1) :: form
        integer          :: mode
        if( self%existence )then
            if( .not. file_exists(fname) )then
                write(logfhandle,*) 'file: ', fname%to_char()
                THROW_HARD('The file you are trying to open does not exists; open')
            endif
            if( present(formatchar) )then
                form = formatchar
            else
                form = fname2format(fname)
            endif
            self%ft = .false.
            select case(form)
                case('M')
                    call ioimg%open(fname, self%ldim, self%smpd, formatchar=formatchar, readhead=readhead, rwaction=rwaction)
                    ! data type: 0 image: signed 8-bit bytes range -128 to 127
                    !            1 image: 16-bit halfwords
                    !            2 image: 32-bit reals (DEFAULT MODE)
                    !            3 transform: complex 16-bit integers
                    !            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
                    mode = ioimg%getMode()
                    if( mode == 3 .or. mode == 4 ) self%ft = .true.
                case('F','S','J','L')
                    call ioimg%open(fname, self%ldim, self%smpd, formatchar=formatchar, readhead=readhead, rwaction=rwaction)
            end select
        else
            THROW_HARD('image need to be constructed before read/write; open')
        endif
    end subroutine open

    !===========================
    ! read
    !===========================
    module subroutine read( self, fname, i, readhead )
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        integer,          optional, intent(in)    :: i
        logical,          optional, intent(in)    :: readhead
        type(imgfile)    :: ioimg
        character(len=1) :: form
        integer          :: ldim(3), first_slice
        integer          :: last_slice, ii
        real             :: smpd
        logical          :: isvol
        ldim  = self%ldim
        smpd  = self%smpd
        isvol = .true. ! assume volume by default
        ii    = 1      ! default location
        if( present(i) )then
            ! we are reading from a stack & in SIMPLE volumes are not allowed
            ! to be stacked so the image object must be 2D
            isvol = .false.
            ii = i ! replace default location
        endif
        form = fname2format(fname)
        select case(form)
            case('M', 'F', 'S', 'J', 'L')
                call self%open(fname, ioimg, form, readhead, rwaction='READ')
            case DEFAULT
                write(logfhandle,*) 'Trying to read from file: ', fname%to_char()
                THROW_HARD('unsupported file format; read')
        end select
        if( isvol )then
            first_slice = 1
            last_slice = ldim(3)
        else
            first_slice = ii
            last_slice = ii
        endif
        call ioimg%rSlices(first_slice,last_slice,self%rmat,is_mrc=form.eq.'M')
        call ioimg%close
    end subroutine read

    !===========================
    ! write
    !===========================
    module subroutine write( self, fname, i, del_if_exists)
        class(image),      intent(inout) :: self
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: i
        logical, optional, intent(in)    :: del_if_exists
        real             :: dev, mean
        type(imgfile)    :: ioimg
        character(len=1) :: form
        integer          :: first_slice, last_slice, iform, ii
        logical          :: isvol, die
        isvol = .false.
        die   = .false.
        isvol = self%is_3d()
        if( present(del_if_exists) ) die = del_if_exists
        ii = 1 ! default location
        if( present(i) )then
            ! we are writing to a stack & in SIMPLE volumes are not allowed
            ! to be stacked so the image object must be 2D
            if( isvol ) THROW_HARD('trying to write 3D image to stack ; write')
            ii = i ! replace default location
        endif
        ! work out the slice range
        if( isvol )then
            first_slice = 1
            last_slice = self%ldim(3)
        else
            first_slice = ii
            last_slice = ii
        endif
        form = fname2format(fname)
        select case(form)
            case('M','F')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die, formatchar=form, readhead=.false.)
                if( self%ft )then
                    call ioimg%setMode(4)
                else
                    call ioimg%setMode(2)
                endif
                call self%rmsd(dev, mean=mean)
                call ioimg%setRMSD(dev)
                call ioimg%setMean(mean)
                ! write slice(s) to disk & close
                call ioimg%wmrcSlices(first_slice,last_slice,self%rmat,self%ldim,self%ft)
                call ioimg%close
            case('S')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die, formatchar=form, readhead=.false.)
                ! iform file type specifier:
                !   1 = 2D image
                !   3 = 3D volume
                ! -11 = 2D Fourier odd
                ! -12 = 2D Fourier even
                ! -21 = 3D Fourier odd
                ! -22 = 3D Fourier even
                if( self%is_2d() )then
                    if( self%ft )then
                        if( self%even_dims() )then
                            iform = -12
                        else
                            iform = -11
                        endif
                    else
                        iform = 1
                    endif
                else
                    if( self%ft )then
                        if( self%even_dims() )then
                            iform = -22
                        else
                            iform = -21
                        endif
                    else
                        iform = 3
                    endif
                endif
                call ioimg%setIform(iform)
                ! write slice(s) to disk & close
                call ioimg%wSlices(first_slice,last_slice,self%rmat,self%ldim,self%ft,self%smpd)
                call ioimg%close
            case DEFAULT
                write(logfhandle,*) 'format descriptor: ', form
                THROW_HARD('unsupported file format; write')
        end select
    end subroutine write

    !===========================
    ! update_header_stats
    !===========================
    module subroutine update_header_stats( self, fname, stats)
        class(image),  intent(inout) :: self
        class(string), intent(in)    :: fname
        real,          intent(in)    :: stats(4) ! stats to update: min, max, mean, rms
        type(imgfile)    :: ioimg
        character(len=1) :: form
        if( self%existence )then
            if( self%ft )return
            form = fname2format(fname)
            select case(form)
            case('M','F')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=.false., formatchar=form, readhead=.true.)
                call ioimg%setMode(2) ! 32-bit reals (DEFAULT MODE)
                !  updates header
                call ioimg%update_MRC_stats(stats)
                ! writes header & close unit
                call ioimg%close
            case('S')
                ! spider stacks have one header each so we do nothing
                return
            case DEFAULT
                write(logfhandle,*) 'format descriptor: ', form
                THROW_HARD('unsupported file format; update_stats_header')
            end select
        else
            THROW_HARD('nonexisting image cannot be updated; update_header_stats')
        endif
    end subroutine update_header_stats

    !===========================
    ! write_jpg
    !===========================
    module subroutine write_jpg( self, fname, quality, colorspec, norm )
        use simple_jpg, only: jpg_img
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        integer,          optional, intent(in)    :: quality, colorspec
        logical,          optional, intent(in)    :: norm
        type(jpg_img)     :: jpg
        integer           :: status
        logical           :: norm_here
        norm_here = .false.
        if(present(norm))norm_here = norm
        if( norm_here )call self%norm4viz
        if( self%is_2d() )then
            status = jpg%writejpg(fname%to_char(), self%rmat(:self%ldim(1),:self%ldim(2),1),&
                &quality=quality, colorspec=colorspec)
        else
            status = jpg%writejpg(fname%to_char(), self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                &quality=quality, colorspec=colorspec)
        endif
    end subroutine write_jpg

end submodule simple_image_io
