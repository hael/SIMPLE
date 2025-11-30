! This type is for single-threaded non-contiguous reading of image stacks
module simple_discrete_stack_io
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,   only: image, image_ptr
use simple_imgfile, only: imgfile
implicit none

public :: dstack_io
private
#include "simple_local_flags.inc"

type dstack_io
    private
    type(image_ptr), allocatable :: img_ptrs(:)
    type(imgfile),   allocatable :: ioimgs(:)
    type(string),    allocatable :: stknames(:)
    integer,         allocatable :: nptcls(:)
    logical,         allocatable :: fts(:), l_open(:)
    real                         :: smpd = 0.
    integer                      :: n    = 0, box = 0
    logical :: exists = .false.
contains
    procedure          :: new
    procedure, private :: open_1, open_2
    procedure          :: read
    generic,   private :: open => open_1, open_2
    procedure          :: does_exist
    procedure          :: kill
end type dstack_io

contains

    subroutine new( self, smpd, box )
        class(dstack_io),  intent(inout) :: self
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: box
        call self%kill
        self%n = 1
        allocate(self%stknames(self%n),self%ioimgs(self%n),self%img_ptrs(self%n),&
            &self%l_open(self%n),self%fts(self%n),self%nptcls(self%n))
        self%smpd        = smpd
        self%nptcls      = 0
        self%box         = box
        self%fts         = .false.
        self%l_open      = .false.
        self%exists      = .true.
    end subroutine new

    subroutine read( self, stkname, ind_in_stk, img )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer,          intent(in)    :: ind_in_stk
        class(image),     intent(inout) :: img
        integer, parameter :: ithr = 1
        ! multi-threaded
        ! ithr = omp_get_thread_num() + 1
        ! if( self%stknames(ithr)).ne.stkname ) call self%open(stkname, ithr)
        ! single-threaded
        if( self%stknames(ithr).ne.stkname ) call self%open(stkname)
        if( .not. self%l_open(ithr) ) THROW_HARD('stack not opened')
        if( ind_in_stk < 1 .or. ind_in_stk > self%nptcls(ithr) )then
            THROW_HARD('index i out of range: '//int2str(ind_in_stk)//' / '//int2str(self%nptcls(ithr)))
        endif
        call img%set_ft(self%fts(ithr))
        call img%get_rmat_ptr(self%img_ptrs(ithr)%rmat)
        call self%ioimgs(ithr)%rSlices(ind_in_stk,ind_in_stk,&
            &self%img_ptrs(ithr)%rmat(1:self%box,1:self%box,1:1),is_mrc=.true.)
    end subroutine read

    ! single-threaded
    subroutine open_1( self, stkname )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer :: ldim(3), mode ! FT or not in MRC file lingo
        self%stknames(1) = stkname
        call self%ioimgs(1)%close
        call find_ldim_nptcls(self%stknames(1), ldim, self%nptcls(1))
        if( (ldim(1)==self%box) .and. (ldim(2)==self%box) )then
            ldim(3) = 1
            call self%ioimgs(1)%open(self%stknames(1), ldim, self%smpd, formatchar='M', readhead=.true., rwaction='READ')
            mode = self%ioimgs(1)%getMode()
            self%fts(1)    = (mode == 3) .or. (mode == 4)
            self%l_open(1) = .true.
        else
            write(logfhandle,*) 'ldim ',ldim
            write(logfhandle,*) 'box ',self%box
            write(logfhandle,*) 'stkname ', stkname%to_char()
            THROW_HARD('Incompatible dimensions!')
        endif
    end subroutine open_1

    ! multi-threaded version
    subroutine open_2( self, stkname, ithr )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer,          intent(in)    :: ithr
        integer          :: ldim(3), mode ! FT or not in MRC file lingo
        self%stknames(ithr) = stkname
        !$omp critical
        call self%ioimgs(ithr)%close
        call find_ldim_nptcls(self%stknames(ithr), ldim, self%nptcls(ithr))
        ldim(3) = 1
        call self%ioimgs(ithr)%open(self%stknames(ithr), ldim, self%smpd, formatchar='M', readhead=.true., rwaction='READ')
        !$omp end critical
        if( (ldim(1)==self%box) .and. (ldim(2)==self%box) )then
            mode = self%ioimgs(ithr)%getMode()
            self%fts(ithr)    = (mode == 3) .or. (mode == 4)
            self%l_open(ithr) = .true.
        else
            THROW_HARD('Incompatible dimensions!')
        endif
    end subroutine open_2

    pure logical function does_exist(self)
        class(dstack_io), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    subroutine kill( self )
        class(dstack_io), intent(inout) :: self
        integer :: i
        if( self%exists )then
            do i = 1,self%n
                call self%ioimgs(i)%close
                nullify(self%img_ptrs(i)%rmat)
            enddo
            deallocate(self%stknames,self%nptcls,self%fts,self%l_open,self%ioimgs,self%img_ptrs)
            self%n      = 0
            self%box    = 0
            self%smpd   = 0.
            self%exists = .false.
        endif
    end subroutine kill

end module simple_discrete_stack_io