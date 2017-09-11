! real hash data structure
#include "simple_lib.f08"
module simple_hash
use simple_defs    ! use all in there
use simple_strings, only: real2str, parse
use simple_syslib,  only: alloc_errchk
implicit none

public :: hash, test_hash
private

integer, parameter :: NMAX=100   !< maximum number of entries in hash table

!> hash stuct
type :: hash
    private
    character(len=32) :: keys(NMAX)
    real              :: vals(NMAX)
    integer           :: hash_index=0 !< index
  contains
    procedure, private :: copy
    procedure, private :: assign
    generic   :: assignment(=) => assign
    procedure :: push
    procedure :: size_of_hash
    procedure :: set
    procedure :: isthere
    procedure :: get
    procedure :: get_keys
    procedure :: get_vals
    procedure :: hash2str
    procedure :: print
    procedure :: write
    procedure :: read
end type hash

contains

    !>  \brief  copies a hash instance
    subroutine copy( self_out, self_in )
        class(hash), intent(inout) :: self_out
        class(hash), intent(in)    :: self_in
        integer :: i
        if( self_in%hash_index > 0 )then
            do i=1,self_in%hash_index
                self_out%keys(i) = self_in%keys(i)
                self_out%vals(i) = self_in%vals(i)
            end do
            self_out%hash_index = self_in%hash_index
        endif
    end subroutine copy

    !>  \brief  is a polymorphic assigner
    subroutine assign( self_out, self_in )
        class(hash), intent(inout) :: self_out
        class(hash), intent(in)    :: self_in
        call self_out%copy(self_in)
    end subroutine assign

    !>  \brief  pushes values to the hash
    subroutine push( self, key, val )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        self%hash_index = self%hash_index+1
        if( self%hash_index > NMAX )then
            write(*,*) 'Hash table full; push; simple_hash, nentries:', self%hash_index
            stop
        endif
        self%keys(self%hash_index) = trim(key)
        self%vals(self%hash_index) = val
    end subroutine push

    !>  \brief  returns size of hash
    function size_of_hash( self ) result( sz )
        class(hash), intent(in) :: self
        integer :: sz
        sz = self%hash_index
    end function size_of_hash

    !>  \brief  sets a value in the hash
    subroutine set( self, key, val )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        integer :: i
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( trim(self%keys(i)) .eq. trim(key) )then
                    self%vals(i) = val
                    return
                endif
            end do
        endif
        call self%push(key, val)
    end subroutine set

    !>  \brief  check for presence of key in the hash
    function isthere( self, key ) result( found )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: i
        logical :: found
        found = .false.
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( trim(self%keys(i)) .eq. trim(key) )then
                    found = .true.
                    return
                endif
            end do
        endif
    end function isthere

    !>  \brief  gets a value in the hash
    function get( self, key ) result( val )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        real    :: val
        integer :: i
        val = 0.
        do i=1,self%hash_index
            if( trim(self%keys(i)) .eq. trim(key) )then
                val = self%vals(i)
                return
            endif
        end do
    end function get

    !>  \brief  returns the keys of the hash
    function get_keys( self ) result( keys )
        class(hash), intent(inout) :: self
        character(len=32), allocatable :: keys(:)
        allocate(keys(self%hash_index), source=self%keys(:self%hash_index),stat=alloc_stat)
        if(alloc_stat /= 0) allocchk("In get_keys")
    end function get_keys

    !>  \brief  returns the values of the hash
    function get_vals( self ) result( vals )
        class(hash), intent(inout) :: self
        real(kind=4), allocatable  :: vals(:)
        allocate(vals(self%hash_index), source=self%vals(:self%hash_index),stat=alloc_stat)
        if(alloc_stat /= 0) allocchk("In get_vals")
    end function get_vals

    !>  \brief  convert hash to string
    function hash2str( self ) result( str )
        class(hash), intent(inout)    :: self
        character(len=:), allocatable :: str, str_moving
        integer :: i
        if( self%hash_index > 0 )then
            if( self%hash_index == 1 )then
                allocate(str, source=trim(self%keys(1))//'='//trim(real2str(self%vals(1))),stat=alloc_stat)
                if(alloc_stat /= 0) allocchk("In hash2str 1")
                return
            endif
            allocate(str_moving, source=trim(self%keys(1))//'='//trim(real2str(self%vals(1)))//' ',stat=alloc_stat)
            if(alloc_stat /= 0) allocchk("In hash2str 2")
            if( self%hash_index > 2 )then
                do i=2,self%hash_index-1
                    allocate(str, source=str_moving//trim(self%keys(i))//'='//trim(real2str(self%vals(i)))//' ',stat=alloc_stat)
                    if(alloc_stat /= 0) allocchk("In hash2str 2")
                    deallocate(str_moving)
                    allocate(str_moving,source=str,stat=alloc_stat)
                    if(alloc_stat /= 0) allocchk("In hash2str 4")
                    deallocate(str)
                end do
            endif
            allocate(str,source=trim(str_moving//trim(self%keys(self%hash_index))&
                &//'='//trim(real2str(self%vals(self%hash_index)))),stat=alloc_stat)
            if(alloc_stat /= 0) allocchk("In hash2str 5")
        endif
    end function hash2str

    !>  \brief  prints the hash
    subroutine print( self )
        class(hash), intent(inout) :: self
        integer :: i
        do i=1,self%hash_index-1
            write(*,"(1X,A,A)", advance="no") trim(self%keys(i)), '='
            write(*,"(A)", advance="no") trim(real2str(self%vals(i)))
        end do
        write(*,"(1X,A,A)", advance="no") trim(self%keys(self%hash_index)), '='
        write(*,"(A)") trim(real2str(self%vals(self%hash_index)))
    end subroutine print

    !>  \brief  writes the hash to file
    subroutine write( self, fnr )
        class(hash), intent(inout) :: self
        integer, intent(in)        :: fnr
        integer                    :: i
        do i=1,self%hash_index-1
            write(fnr,"(1X,A,A)", advance="no") trim(self%keys(i)), '='
            write(fnr,"(A)", advance="no") trim(real2str(self%vals(i)))
        end do
        write(fnr,"(1X,A,A)", advance="no") trim(self%keys(self%hash_index)), '='
        write(fnr,"(A)") trim(real2str(self%vals(self%hash_index)))
    end subroutine write

    !>  \brief  reads a row of a text-file into the inputted hash, assuming key=value pairs
    subroutine read( self, fnr )
        use simple_math, only: is_a_number
        class(hash), intent(inout) :: self
        integer,     intent(in)    :: fnr
        character(len=1024)        :: line
        character(len=64)          :: keyvalpairs(NMAX)
        character(len=64)          :: keyvalpair
        character(len=32)          :: key
        character(len=100)         :: io_message
        real    :: val
        integer :: j, pos1, nkeyvalpairs, io_stat
        read(fnr,fmt='(A)',advance='no',eor=100) line
        100 call parse(line, ' ', keyvalpairs, nkeyvalpairs)
        do j=1,nkeyvalpairs
            keyvalpair = keyvalpairs(j)
            pos1 = index(keyvalpair, '=') ! position of '='
            if( pos1 /= 0 )then
                key = keyvalpair(:pos1-1) ! key
                read(keyvalpair(pos1+1:),*,iostat=io_stat,iomsg=io_message) val
                ! Check the read was successful
                if( io_stat .ne. 0 )then
                    write(*,*) 'key = ', key
                    write(*,*) 'value = ', val
                    write(*,'(a,i0)') '**ERROR(simple_hash::read): I/O error ', io_stat
                    write(*,'(2a)') 'IO error message was: ', io_message
                    stop 'I/O error; read; simple_hash'
                endif
            endif
            if( .not. is_a_number(val) ) val = 0.
            call self%set(key, val)
        end do
    end subroutine read

    !>  \brief  is the hash unit test
    subroutine test_hash
        use simple_fileio      , only: fopen, fclose, fileio_errmsg
        type(hash) :: htab, htab2
        integer    :: fnr, file_stat
        write(*,'(a)') '**info(simple_hash_unit_test: testing all functionality'
        call htab%push('hans',  1.)
        call htab%push('donia', 2.)
        call htab%push('kenji', 3.)
        call htab%push('ralph', 4.)
        call htab%push('kyle',  5.)
        call htab%print

        call fopen(fnr, status='replace', action='write', file='hashtest.txt', iostat=file_stat)
        call fileio_errmsg("simple_hash_unit_test: test failed to open hastest.txt",file_stat)
        call htab%write(fnr)
        call fclose(fnr,errmsg="simple_hash_unit_test: test failed to close hastest.txt")

        call fopen(fnr, status='old', action='read', file='hashtest.txt', iostat=file_stat)
        call fileio_errmsg("simple_hash_unit_test: test failed to open 2nd  hastest.txt",file_stat)
        call htab2%read(fnr)
        call fclose(fnr,errmsg="simple_hash_unit_test: test failed to close 2nd hastest.txt")

        call htab2%print
        write(*,'(a)') 'SIMPLE_HASH_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_hash

end module simple_hash
