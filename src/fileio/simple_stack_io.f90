module simple_stack_io
include 'simple_lib.f08'
use simple_image,   only: image
use simple_imgfile, only: imgfile
implicit none

public :: stack_io
private
#include "simple_local_flags.inc"

type stack_io
    private
    real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
    type(image)   :: buffer
    type(imgfile) :: ioimg
    integer       :: nptcls  = 0, fromp = 0, top = 0, bufsz = 0, n_in_buf = 0
    integer       :: ldim(3) = [0,0,0]
    real          :: smpd    = 0.
    logical       :: ft = .false., l_read = .false., is_open = .false.
contains
    procedure          :: open
    procedure          :: get_nptcls
    procedure          :: get_ldim
    procedure          :: read
    procedure          :: write
    procedure, private :: write_buffer
    procedure          :: close
end type stack_io

integer, parameter :: BUFSZ_DEFAULT = 1024

contains

    subroutine open( self, stkname, smpd, rwaction, is_ft, box, bufsz )
        class(stack_io),   intent(inout) :: self
        character(len=*),  intent(in)    :: stkname, rwaction
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: is_ft
        integer, optional, intent(in)    :: box
        integer, optional, intent(in)    :: bufsz
        character(len=1) :: form
        integer          :: mode
        call self%close
        ! extract info about the stack file and open it
        form        = fname2format(trim(stkname))
        if( form .ne. 'M') THROW_HARD('non MRC stacks unsupported')
        self%smpd   = smpd
        self%nptcls = 0.
        self%ft     = .false.
        select case(trim(rwaction))
            case('READ','read')
                if( .not. file_exists(trim(stkname)) ) THROW_HARD('input stack file does not exists')
                call find_ldim_nptcls(trim(stkname), self%ldim, self%nptcls)
                self%ldim(3) = 1
                call self%ioimg%open(trim(stkname), self%ldim, self%smpd, formatchar=form, readhead=.true., rwaction='READ')
                mode = self%ioimg%getMode()
                if( mode == 3 .or. mode == 4 ) self%ft = .true.
                self%l_read = .true.
            case DEFAULT
                if( present(box) )then
                    self%ldim = [box,box,1]
                    call self%ioimg%open(trim(stkname), self%ldim, self%smpd, formatchar=form, readhead=.false.)
                else
                    THROW_HARD('optional box dummy argument needed to write to stack')
                endif
                self%l_read = .false.
        end select
        if( present(is_ft) ) self%ft = is_ft
        ! allocate the buffer
        self%bufsz = BUFSZ_DEFAULT
        if( present(bufsz) )   self%bufsz = bufsz
        if( self%nptcls /= 0 ) self%bufsz = min(self%bufsz, self%nptcls)
        call self%buffer%new([self%ldim(1),self%ldim(2),self%bufsz], self%smpd)
        call self%buffer%get_rmat_ptr(self%rmat_ptr)
        self%is_open = .true.
    end subroutine open

    function get_nptcls( self ) result( nptcls )
        class(stack_io), intent(in) :: self
        integer :: nptcls
        nptcls = self%nptcls
    end function get_nptcls

    function get_ldim( self ) result( ldim )
        class(stack_io), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    subroutine read( self, i, img )
        class(stack_io), intent(inout) :: self
        integer,         intent(in)    :: i
        class(image),    intent(inout) :: img
        integer :: ind_in_buf, bufsz
        if( i < 1 .or. i > self%nptcls ) THROW_HARD('index i out of range')
        if( i >= self%fromp .and. i <= self%top )then
            ! the image is in buffer
        else
            ! read buffer
            if( self%fromp == 0 )then
                self%fromp = 1
                self%top   = self%bufsz
            else
                self%fromp = self%top + 1
                self%top   = min(self%fromp + self%bufsz - 1, self%nptcls)
                if( self%fromp == self%nptcls ) return
            endif
            bufsz = self%top - self%fromp + 1
            call self%ioimg%rSlices(self%fromp,self%top,self%rmat_ptr(:self%ldim(1),:self%ldim(2),:bufsz),is_mrc=.true.)
        endif
        ind_in_buf = i - self%fromp + 1
        call img%set_rmat(self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf), self%ft)
    end subroutine read

    subroutine write( self, i, img )
        class(stack_io), intent(inout) :: self
        integer,         intent(in)    :: i
        class(image),    intent(inout) :: img
        real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
        integer :: ind_in_buf
        if( self%fromp == 0 )then
            self%fromp = 1
            self%top   = self%bufsz
        endif
        if( i >= self%fromp .and. i <= self%top )then
            ! the buffer can be set
            call img%get_rmat_ptr(rmat_ptr)
            ind_in_buf = i - self%fromp + 1
            self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf) = rmat_ptr(:self%ldim(1),:self%ldim(2),1:1)
            self%n_in_buf  = self%n_in_buf + 1
        else
            call self%write_buffer
            ! update range
            self%fromp = self%top + 1
            self%top   = self%fromp + self%bufsz - 1
            if( self%fromp == self%nptcls ) return
            if( i >= self%fromp .and. i <= self%top )then
                ! the index is within the range, so the buffer can be set
                call img%get_rmat_ptr(rmat_ptr)
                ind_in_buf = i - self%fromp + 1
                self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf) = rmat_ptr(:self%ldim(1),:self%ldim(2),1:1)
                self%n_in_buf  = self%n_in_buf + 1
            else
                THROW_HARD('index i is out of range')
            endif
        endif
    end subroutine write

    subroutine write_buffer( self )
        class(stack_io), intent(inout) :: self
        integer :: top
        if( self%n_in_buf > 0 )then
            top = self%fromp + self%n_in_buf - 1
            call self%ioimg%wmrcSlices(self%fromp,top,self%rmat_ptr(:self%ldim(1),:self%ldim(2),:self%n_in_buf),self%ldim,is_ft=self%ft)
            self%n_in_buf = 0
        endif
    end subroutine write_buffer

    subroutine close( self )
        class(stack_io), intent(inout) :: self
        if( self%is_open )then
            if( .not. self%l_read ) call self%write_buffer
            call self%ioimg%close
            self%nptcls   = 0
            self%fromp    = 0
            self%top      = 0
            self%bufsz    = 0
            self%n_in_buf = 0
            self%ldim     = [0,0,0]
            self%smpd     = 0.
            self%ft       = .false.
            self%l_read   = .false.
            self%rmat_ptr => null()
            call self%buffer%kill
            self%is_open  = .false.
        endif
    end subroutine close

end module simple_stack_io
