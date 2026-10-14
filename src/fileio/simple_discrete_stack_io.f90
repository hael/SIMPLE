!@descr: for cached non-contiguous reading of image stacks
module simple_discrete_stack_io
use simple_core_module_api
use simple_image,   only: image
use simple_imgfile, only: imgfile
use simple_imghead, only: MRC_MODE_FLOAT16
implicit none

public :: dstack_io
private
#include "simple_local_flags.inc"

! one MRC stack open for reading; its header, read once by the open, decides dimensions and record offsets
type dstack_io
    private
    type(imgfile) :: ioimg
    type(string)  :: stkname
    real          :: smpd   = 0.
    integer       :: box    = 0, nptcls = 0, mode = -1
    logical       :: ft = .false., l_open = .false., l_threadsafe_reads = .true.
    logical       :: exists = .false.
contains
    procedure :: new
    procedure :: open => open_1
    procedure :: read
    procedure, private :: check_index
    procedure :: close   => close_stack
    procedure :: kill
end type dstack_io

contains

    subroutine new( self, smpd, box )
        class(dstack_io), intent(inout) :: self
        real,             intent(in)    :: smpd
        integer,          intent(in)    :: box
        call self%kill
        self%smpd               = smpd
        self%box                = box
        self%nptcls             = 0
        self%ft                 = .false.
        self%l_open             = .false.
        self%l_threadsafe_reads = .true.
        self%exists             = .true.
    end subroutine new

    !> One image; a stack other than the open one is opened first (serialized inside parallel regions)
    subroutine read( self, stkname, ind_in_stk, img )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer,          intent(in)    :: ind_in_stk
        class(image),     intent(inout) :: img
        if( .not. self%l_open .or. self%stkname.ne.stkname )then
            if( OMP_IN_PARALLEL() )then
                !$omp critical(dstack_io_open)
                call self%open(stkname)
                !$omp end critical(dstack_io_open)
            else
                call self%open(stkname)
            endif
        endif
        call self%check_index(ind_in_stk)
        call img%set_ft(self%ft)
        if( self%l_threadsafe_reads )then
            call img%read_single_mrc_image(self%ioimg, ind_in_stk)
        else
            !$omp critical(dstack_io_rslice_tmp)
            call img%read_single_mrc_image(self%ioimg, ind_in_stk)
            !$omp end critical(dstack_io_rslice_tmp)
        endif
    end subroutine read

    subroutine close_stack( self )
        class(dstack_io), intent(inout) :: self
        if( self%l_open ) call self%ioimg%close
        self%l_open = .false.
        self%nptcls = 0
        self%mode   = -1
    end subroutine close_stack

    subroutine check_index( self, ind_in_stk )
        class(dstack_io), intent(in) :: self
        integer,          intent(in) :: ind_in_stk
        if( ind_in_stk < 1 .or. ind_in_stk > self%nptcls )then
            write(logfhandle,*) 'stack: ', trim(self%stkname%to_char())
            THROW_HARD('index i out of range: '//int2str(ind_in_stk)//' / '//int2str(self%nptcls))
        endif
    end subroutine check_index

    !> Opens a stack read-only; the header the open reads gives its dimensions and image count
    subroutine open_1( self, stkname )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer :: ldim(3)
        call self%close
        self%stkname = stkname
        call self%ioimg%open(self%stkname, [self%box,self%box,1], self%smpd, formatchar='M', readhead=.true., rwaction='READ')
        ldim = self%ioimg%getDims()
        if( (ldim(1) /= self%box) .or. (ldim(2) /= self%box) )then
            write(logfhandle,*) 'ldim ', ldim
            write(logfhandle,*) 'box ', self%box
            write(logfhandle,*) 'stkname ', stkname%to_char()
            THROW_HARD('Incompatible dimensions!')
        endif
        self%nptcls             = ldim(3)
        self%mode               = self%ioimg%getMode()
        self%ft                 = (self%mode == 3) .or. (self%mode == 4)
        self%l_threadsafe_reads = any(self%mode == [0, 1, 2, 3, 4, 6, MRC_MODE_FLOAT16])
        self%l_open             = .true.
    end subroutine open_1

    subroutine kill( self )
        class(dstack_io), intent(inout) :: self
        if( self%exists )then
            call self%close
            call self%stkname%kill
            self%box    = 0
            self%smpd   = 0.
            self%exists = .false.
        endif
    end subroutine kill

end module simple_discrete_stack_io
