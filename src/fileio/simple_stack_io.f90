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
    procedure :: new
    procedure :: read
    procedure :: write
    procedure :: kill
end type stack_io

integer, parameter :: BUFSZ_DEFAULT = 1024

contains

    subroutine new( self, stkname, smpd, rwaction, is_ft, box, bufsz )
        class(stack_io),   intent(inout) :: self
        character(len=*),  intent(in)    :: stkname, rwaction
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: is_ft
        integer, optional, intent(in)    :: box
        integer, optional, intent(in)    :: bufsz
        character(len=1) :: form
        integer          :: mode
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
            case DEFAULT
                if( present(box) )then
                    self%ldim = [box,box,1]
                    call self%ioimg%open(trim(stkname), self%ldim, self%smpd, formatchar=form, readhead=.true.)
                else
                    THROW_HARD('optional box dummy argument needed to write to stack')
                endif
        end select
        if( present(is_ft) ) self%ft = is_ft
        ! allocate the buffer
        self%bufsz = BUFSZ_DEFAULT
        if( present(bufsz) )   self%bufsz = bufsz
        if( self%nptcls /= 0 ) self%bufsz = min(self%bufsz, self%nptcls)
        call self%buffer%new([self%ldim(1),self%ldim(2),self%bufsz], self%smpd)
        call self%buffer%get_rmat_ptr(self%rmat_ptr)
    end subroutine new

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
        integer :: ind_in_buf, bufsz
        if( self%fromp == 0 )then
            self%fromp = 1
            self%top   = self%bufsz
        endif
        if( i >= self%fromp .and. i <= self%top )then
            ! the buffer can be set
            call img%get_rmat_ptr(rmat_ptr)
            ind_in_buf = i - self%fromp + 1
            self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf) = rmat_ptr(:self%ldim(1),:self%ldim(2),1:1)
        else
            ! write buffer
            bufsz = self%top - self%fromp + 1
            call self%ioimg%wmrcSlices(self%fromp,self%top,self%rmat_ptr(:self%ldim(1),:self%ldim(2),:bufsz),self%ldim,is_ft=self%ft)
            ! update range
            self%fromp = self%top + 1
            self%top   = min(self%fromp + self%bufsz - 1, self%nptcls)
            if( self%fromp == self%nptcls ) return
            if( i >= self%fromp .and. i <= self%top )then
                ! the index is within the range, so the buffer can be set
                call img%get_rmat_ptr(rmat_ptr)
                ind_in_buf = i - self%fromp + 1
                self%rmat_ptr(:self%ldim(1),:self%ldim(2),ind_in_buf:ind_in_buf) = rmat_ptr(:self%ldim(1),:self%ldim(2),1:1)
            else
                THROW_HARD('index i is out of range')
            endif
        endif
    end subroutine write

    subroutine kill( self )
        class(stack_io), intent(inout) :: self
        self%rmat_ptr => null()
        call self%buffer%kill
        call self%ioimg%close
    end subroutine kill

end module simple_stack_io
