!@descr: for cached non-contiguous reading of image stacks
module simple_discrete_stack_io
use simple_core_module_api
use simple_image,   only: image
use simple_imgfile, only: imgfile
implicit none

public :: dstack_io
private
#include "simple_local_flags.inc"

type dstack_io
    private
    type(imgfile),   allocatable :: ioimgs(:)
    type(string),    allocatable :: stknames(:)
    type(string),    allocatable :: cache_stknames(:)
    integer,         allocatable :: nptcls(:)
    integer,         allocatable :: cache_ldims(:,:), cache_nptcls(:)
    logical,         allocatable :: fts(:), l_open(:), l_threadsafe_reads(:)
    real                         :: smpd = 0.
    integer                      :: n    = 0, box = 0, ncache = 0
    logical :: exists = .false.
contains
    procedure          :: new
    procedure          :: cache_stack_info
    procedure          :: open => open_1
    procedure, private :: get_stack_info, clear_stack_cache
    procedure          :: read
    procedure          :: does_exist
    procedure          :: kill
end type dstack_io

contains

    subroutine new( self, smpd, box )
        class(dstack_io), intent(inout) :: self
        real,             intent(in)    :: smpd
        integer,          intent(in)    :: box
        call self%kill
        self%n = 1
        allocate(self%stknames(self%n),self%ioimgs(self%n),&
            &self%l_open(self%n),self%fts(self%n),self%l_threadsafe_reads(self%n),self%nptcls(self%n))
        self%smpd               = smpd
        self%nptcls             = 0
        self%box                = box
        self%fts                = .false.
        self%l_open             = .false.
        self%l_threadsafe_reads = .true.
        self%exists             = .true.
    end subroutine new

    subroutine cache_stack_info( self, stkname, ldim, nptcls )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer,          intent(in)    :: ldim(3), nptcls
        call self%clear_stack_cache
        self%ncache = 1
        allocate(self%cache_stknames(self%ncache), self%cache_ldims(3,self%ncache), self%cache_nptcls(self%ncache))
        self%cache_stknames(1) = stkname
        self%cache_ldims(:,1)  = ldim
        self%cache_nptcls(1)   = nptcls
    end subroutine cache_stack_info

    subroutine read( self, stkname, ind_in_stk, img )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer,          intent(in)    :: ind_in_stk
        class(image),     intent(inout) :: img
        integer, parameter :: ithr = 1
        if( self%stknames(ithr).ne.stkname )then
            if( OMP_IN_PARALLEL() )then
                !$omp critical(dstack_io_open)
                if( self%stknames(ithr).ne.stkname ) call self%open(stkname)
                !$omp end critical(dstack_io_open)
            else
                call self%open(stkname)
            endif
        endif
        if( .not. self%l_open(ithr) ) THROW_HARD('stack not opened')
        if( ind_in_stk < 1 .or. ind_in_stk > self%nptcls(ithr) )then
            write(logfhandle,*) 'stack: ', trim(stkname%to_char())
            THROW_HARD('index i out of range: '//int2str(ind_in_stk)//' / '//int2str(self%nptcls(ithr)))
        endif
        call img%set_ft(self%fts(ithr))
        if( self%l_threadsafe_reads(ithr) )then
            call img%read_single_mrc_image(self%ioimgs(ithr), ind_in_stk)
        else
            !$omp critical(dstack_io_rslice_tmp)
            call img%read_single_mrc_image(self%ioimgs(ithr), ind_in_stk)
            !$omp end critical(dstack_io_rslice_tmp)
        endif
    end subroutine read

    ! single-threaded
    subroutine open_1( self, stkname )
        class(dstack_io), intent(inout) :: self
        class(string),    intent(in)    :: stkname
        integer :: ldim(3), mode ! FT or not in MRC file lingo
        self%stknames(1) = stkname
        call self%ioimgs(1)%close
        call self%get_stack_info(self%stknames(1), ldim, self%nptcls(1))
        if( (ldim(1)==self%box) .and. (ldim(2)==self%box) )then
            ldim(3) = 1
            call self%ioimgs(1)%open(self%stknames(1), ldim, self%smpd, formatchar='M', readhead=.true., rwaction='READ')
            mode = self%ioimgs(1)%getMode()
            self%fts(1)                = (mode == 3) .or. (mode == 4)
            self%l_threadsafe_reads(1) = (mode == 2) .or. (mode == 4)
            self%l_open(1)             = .true.
        else
            write(logfhandle,*) 'ldim ',ldim
            write(logfhandle,*) 'box ',self%box
            write(logfhandle,*) 'stkname ', stkname%to_char()
            THROW_HARD('Incompatible dimensions!')
        endif
    end subroutine open_1

    subroutine get_stack_info( self, stkname, ldim, nptcls )
        class(dstack_io), intent(in)  :: self
        class(string),    intent(in)  :: stkname
        integer,          intent(out) :: ldim(3), nptcls
        integer :: i
        do i = 1,self%ncache
            if( self%cache_stknames(i) .eq. stkname )then
                ldim   = self%cache_ldims(:,i)
                nptcls = self%cache_nptcls(i)
                return
            endif
        enddo
        call find_ldim_nptcls(stkname, ldim, nptcls)
    end subroutine get_stack_info

    subroutine clear_stack_cache( self )
        class(dstack_io), intent(inout) :: self
        if( allocated(self%cache_stknames) )then
            call self%cache_stknames(:)%kill
            deallocate(self%cache_stknames,self%cache_ldims,self%cache_nptcls)
        endif
        self%ncache = 0
    end subroutine clear_stack_cache

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
            enddo
            call self%stknames(:)%kill
            deallocate(self%stknames,self%nptcls,self%fts,self%l_open,self%l_threadsafe_reads,self%ioimgs)
            call self%clear_stack_cache
            self%n      = 0
            self%box    = 0
            self%smpd   = 0.
            self%ncache = 0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_discrete_stack_io
