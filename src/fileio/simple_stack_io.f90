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
    type(image)                 :: buffer
    type(imgfile)               :: ioimg
    integer                     :: nptcls  = 0, fromp = 0, top = 0, bufsz = 0
    integer                     :: ldim(3) = [0,0,0]
    real                        :: smpd    = 0.
    logical                     :: ft      = .false.
contains
    procedure          :: new
    procedure          :: read
    procedure, private :: read_buffer
    procedure          :: kill
end type stack_io

integer, parameter :: BUFSZ_DEFAULT = 1024

contains

    subroutine new( self, stkname, smpd, rwaction, bufsz )
        class(stack_io),   intent(inout) :: self
        character(len=*),  intent(in)    :: stkname, rwaction
        real,              intent(in)    :: smpd
        integer, optional, intent(in)    :: bufsz
        character(len=1) :: form
        integer          :: mode
        ! extract info about the stack file and open it
        if( .not. file_exists(trim(stkname)) ) THROW_HARD('input stack file does not exists')
        form = fname2format(trim(stkname))
        if( form .ne. 'M') THROW_HARD('non MRC stacks unsupported')
        call find_ldim_nptcls(trim(stkname), self%ldim, self%nptcls)
        self%ldim(3) = 1
        self%smpd    = smpd
        select case(trim(rwaction))
            case('READ','read')
                call self%ioimg%open(trim(stkname), self%ldim, self%smpd, formatchar=form, readhead=.true., rwaction='READ')
                mode = self%ioimg%getMode()
                self%ft = .false.
                if( mode == 3 .or. mode == 4 ) self%ft = .true.
            case('WRITE', 'write')
                THROW_HARD('write is not currently supported')
            case DEFAULT
                THROW_HARD('unsupported read/write action (rwaction)')
        end select
        ! allocate the buffer
        self%bufsz = BUFSZ_DEFAULT
        if( present(bufsz) ) self%bufsz = bufsz
        self%bufsz         = min(self%bufsz, self%nptcls)
        call self%buffer%new([self%ldim(1),self%ldim(2),self%bufsz], self%smpd)
        call self%buffer%get_rmat_ptr(self%rmat_ptr)
    end subroutine new

    subroutine read( self, i, img )
        class(stack_io), intent(inout) :: self
        integer,         intent(in)    :: i
        class(image),    intent(inout) :: img
        integer :: ind_in_buf
        if( i < 1 .or. i > self%nptcls ) THROW_HARD('index i out of range')
        if( i >= self%fromp .and. i <= self%top )then
            ! the image is in buffer
            ind_in_buf = i - self%fromp + 1
            call img%set_rmat(self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf), self%ft)
        else
            ! read buffer
            call self%read_buffer
        endif
    end subroutine read

    subroutine read_buffer( self )
        class(stack_io), intent(inout) :: self
        integer :: bufsz
        if( self%fromp == 0 )then
            self%fromp = 1
            self%top   = self%bufsz
        else
            self%fromp = self%top
            self%top   = min(self%top + self%bufsz - 1, self%nptcls)
            if( self%fromp == self%nptcls ) return
        endif
        bufsz = self%top - self%fromp + 1
        call self%ioimg%rSlices(self%fromp,self%top,self%rmat_ptr(:self%ldim(1),:self%ldim(2),:bufsz),is_mrc=.true.)
    end subroutine read_buffer

    subroutine kill( self )
        class(stack_io), intent(inout) :: self
        self%rmat_ptr => null()
        call self%buffer%kill
        call self%ioimg%close
    end subroutine kill

end module simple_stack_io
