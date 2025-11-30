! real hash data structure
module simple_hash
use simple_error, only: simple_exception
use simple_string_utils
use simple_string
use simple_defs
implicit none

public :: hash
private
#include "simple_local_flags.inc"

integer, parameter :: BUFFSZ_DEFAULT = 10

!> hash stuct
type :: hash
    private
    type(string), allocatable :: keys(:)    !< hash keys
    real(dp),     allocatable :: values(:)  !< hash values
    integer :: buffsz     = BUFFSZ_DEFAULT  !< size of first buffer and subsequent allocation increments
    integer :: hash_index = 0               !< index
    logical :: exists     = .false.
  contains
    !< CONSTRUCTORS
    procedure, private :: new_1, new_2
    generic            :: new => new_1, new_2
    procedure, private :: alloc_hash
    procedure, private :: realloc_hash
    procedure, private :: copy
    generic            :: assignment(=) => copy
    !< SETTERS
    procedure, private :: push_1, push_2
    generic            :: push => push_1, push_2
    procedure, private :: set_1, set_2, set_3
    generic            :: set => set_1, set_2, set_3
    procedure          :: delete
    !< GETTERS
    procedure          :: isthere
    procedure          :: lookup
    procedure          :: get
    procedure          :: get_dp
    procedure, private :: getter_1, getter_2, getter_3
    generic            :: getter => getter_1, getter_2, getter_3
    procedure          :: get_keys
    procedure          :: get_key
    procedure          :: get_value_at
    procedure          :: hash2str
    procedure          :: hash_strlen
    procedure          :: size_of
    ! I/O
    procedure          :: print
    !< DESTRUCTORS
    procedure          :: kill
    procedure, private :: dealloc_hash
end type hash

interface hash
    module procedure constructor_1
    module procedure constructor_2
end interface hash

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor_1( ) result( self )
        type(hash) :: self
        call self%new_1
    end function constructor_1

    !>  \brief  is a constructor
    function constructor_2( sz, buffsz ) result( self )
        integer,           intent(in)    :: sz     !< initial size
        integer, optional, intent(in)    :: buffsz !< size of subsequent allocation increments
        type(hash) :: self
        call self%new_2(sz, buffsz)
    end function constructor_2

    !>  \brief  is a constructor
    subroutine new_1( self )
        class(hash), intent(inout) :: self
        call self%kill
        self%buffsz = BUFFSZ_DEFAULT
        call self%alloc_hash(self%buffsz)
        self%exists = .true.
    end subroutine new_1

    !>  \brief  is a constructor
    subroutine new_2( self, sz, buffsz )
        class(hash),       intent(inout) :: self   !< instance
        integer,           intent(in)    :: sz     !< initial size
        integer, optional, intent(in)    :: buffsz !< size of subsequent allocation increments
        call self%kill
        if( present(buffsz) )then
            self%buffsz = buffsz
        else
            self%buffsz = BUFFSZ_DEFAULT
        endif
        call self%alloc_hash(sz)
        self%exists = .true.
    end subroutine new_2

    !>  \brief  allocates keys and values without touching buffsz and hash_index
    subroutine alloc_hash( self, sz )
        class(hash), intent(inout) :: self !< instance
        integer,      intent(in)   :: sz   !< total size
        call self%dealloc_hash
        allocate(self%keys(sz), self%values(sz))
    end subroutine alloc_hash

    !>  \brief  re-allocates keys and values without touching buffsz and hash_index and preserving previous key-value pairs
    subroutine realloc_hash( self )
        class(hash), intent(inout)  :: self
        type(string),   allocatable :: keys_copy(:)
        real(dp),       allocatable :: values_copy(:)
        integer :: newsz, oldsz, i
        if( self%exists )then
            ! old/new size
            oldsz = size(self%keys)
            newsz = oldsz + self%buffsz
            ! copy key-value pairs
            allocate(keys_copy(oldsz), values_copy(oldsz))
            do i=1,oldsz
                if( self%keys(i)%is_allocated() ) keys_copy(i) = self%keys(i)
                values_copy(i) = self%values(i)
            enddo
            ! allocate extended size
            call self%alloc_hash(newsz)
            ! set back the copied key-value pairs
            do i=1,oldsz
                if( keys_copy(i)%is_allocated() ) self%keys(i) = keys_copy(i)
                self%values(i) = values_copy(i)
            enddo
        else
            THROW_HARD('cannot reallocate non-existent hash; realloc_hash')
        endif
    end subroutine realloc_hash

    !>  \brief  copies a hash instance
    subroutine copy( self_out, self_in )
        class(hash), intent(inout) :: self_out
        class(hash), intent(in)    :: self_in
        integer :: i, sz
        sz = size(self_in%keys)
        call self_out%alloc_hash(sz)
        self_out%buffsz     = self_in%buffsz
        self_out%hash_index = self_in%hash_index
        self_out%exists     = .true.
        if( self_in%hash_index > 0 )then
            do i=1,self_in%hash_index
                if( self_in%keys(i)%is_allocated() ) self_out%keys(i) = self_in%keys(i)
                self_out%values(i) = self_in%values(i)
            end do
        endif
    end subroutine copy

    ! SETTERS

    !>  \brief  pushes real values to the hash
    subroutine push_1( self, key, val )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: val
        integer :: sz
        if (.not. self%exists) call self%new
        ! increment index
        self%hash_index = self%hash_index + 1
        ! reallocate if needed
        sz = size(self%keys)
        if( self%hash_index > sz ) call self%realloc_hash
        ! set key
        self%keys(self%hash_index) = trim(adjustl(key))
        ! set value
        self%values(self%hash_index) = real(val,dp)
    end subroutine push_1

    !>  \brief  pushes real values to the hash
    subroutine push_2( self, key, val )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
        integer :: sz
        if (.not. self%exists) call self%new
        ! increment index
        self%hash_index = self%hash_index + 1
        ! reallocate if needed
        sz = size(self%keys)
        if( self%hash_index > sz ) call self%realloc_hash
        ! set key
        self%keys(self%hash_index) = trim(adjustl(key))
        ! set value
        self%values(self%hash_index) = real(val,dp)
    end subroutine push_2

    !>  \brief  sets a real value in the hash
    subroutine set_1( self, key, val )
        class(hash),     intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        integer :: i
        if (.not. self%exists) call self%new
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( self%keys(i) .eq. trim(key) )then
                    self%values(i) = real(val,dp)
                    return
                endif
            end do
        endif
        call self%push(key, real(val,dp))
    end subroutine set_1

    !>  \brief  sets an integer value in the hash
    subroutine set_2( self, key, ival )
        class(hash),     intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: ival
        integer :: i
        if (.not. self%exists) call self%new
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( self%keys(i) .eq. trim(key) )then
                    self%values(i) = real(ival,dp)
                    return
                endif
            end do
        endif
        call self%push(key, real(ival,dp))
    end subroutine set_2

    !>  \brief  sets a real(kind=8) value in the hash
    subroutine set_3( self, key, rval )
        class(hash),     intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: rval
        integer :: i
        if (.not. self%exists) call self%new
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( self%keys(i) .eq. trim(key) )then
                    self%values(i) = rval
                    return
                endif
            end do
        endif
        call self%push(key, rval)
    end subroutine set_3

    !>  \brief  deletes a value in the hash
    subroutine delete( self, key )
        class(hash),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: i, ind
        if( .not. self%exists ) return
        ind = self%lookup( key )
        if( ind==0 .or. ind > self%hash_index ) return
        do i=ind,self%hash_index - 1
            ! replace key
            self%keys(i)  = self%keys(i + 1)
            ! replace value
            self%values(i) = self%values(i + 1)
        enddo
        ! remove previous last entry
        if( self%keys(self%hash_index)%is_allocated() ) call self%keys(self%hash_index)%kill
        self%values( self%hash_index ) = 0.
        ! decrement index
        self%hash_index = self%hash_index - 1
    end subroutine delete

    ! GETTERS

    !>  \brief  check for presence of key in the hash
    pure logical function isthere( self, key )
        class(hash),      intent(in) :: self
        character(len=*), intent(in)  :: key
        integer :: i
        isthere = .false.
        if( self%hash_index >= 1 )then
            do i=1,self%hash_index
                if( self%keys(i) .eq. trim(key) )then
                    isthere = .true.
                    return
                endif
            end do
        endif
    end function isthere

    !>  \brief  gets the index to a key
    pure integer function lookup( self, key )
        class(hash),      intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: i
        lookup = 0
        do i=1,self%hash_index
            if( self%keys(i) .eq. trim(key) )then
                lookup = i
                return
            endif
        end do
    end function lookup

    !>  \brief  gets a value in the hash as real(kind=4)
    pure real(kind=sp) function get( self, key )
        class(hash),      intent(in) :: self
        character(len=*), intent(in) :: key
        integer  :: i
        get = 0.0
        do i=1,self%hash_index
            if( self%keys(i) .eq. trim(key) )then
                get = real(self%values(i))
                return
            endif
        end do
    end function get

    !>  \brief  gets a value in the hash as real(kind=4)
    pure real(kind=dp) function get_dp( self, key )
        class(hash),      intent(in) :: self
        character(len=*), intent(in) :: key
        integer  :: i
        get_dp = 0.0d0
        do i=1,self%hash_index
            if( self%keys(i) .eq. trim(key) )then
                get_dp = self%values(i)
                return
            endif
        end do
    end function get_dp

    !>  \brief  returns real(kind=4)
    pure subroutine getter_1( self, key, rvalsp )
        class(hash),      intent(in)  :: self
        character(len=*), intent(in)  :: key
        real(kind=sp),    intent(out) :: rvalsp
        rvalsp = real(self%get_dp(key))
    end subroutine getter_1

    !>  \brief  returns real(kind=8)
    pure subroutine getter_2( self, key, rvaldp )
        class(hash),      intent(in)  :: self
        character(len=*), intent(in)  :: key
        real(kind=dp),    intent(out) :: rvaldp
        rvaldp = self%get_dp(key)
    end subroutine getter_2

    !>  \brief  returns integer(kind=4)
    pure subroutine getter_3( self, key, ival )
        class(hash),      intent(in)  :: self
        character(len=*), intent(in)  :: key
        integer(kind=sp), intent(out) :: ival
        ival = nint(self%get_dp(key))
    end subroutine getter_3

    !>  \brief  gets a value in the hash
    function get_key( self, keyindx ) result( key )
        class(hash), intent(in) :: self
        integer,     intent(in) :: keyindx
        type(string) :: key
        if( keyindx > 0 .and.  keyindx .le. self%hash_index )then
            key = self%keys(keyindx)
        else
            key = ''
        endif
   end function get_key

    !>  \brief  gets a value in the hash
    pure real function get_value_at( self, keyindx )
        class(hash), intent(in) :: self
        integer,      intent(in) :: keyindx
        get_value_at = 0.0
        if(keyindx > 0 .and.  keyindx .le. self%hash_index)then
            get_value_at = real(self%values(keyindx))
        endif
    end function get_value_at

    function get_keys( self ) result( keys )
        class(hash), intent(in)   :: self
        type(string), allocatable :: keys(:)
        if(self%hash_index .gt. 0) then
            allocate(keys(self%hash_index), source=self%keys(:self%hash_index))
        else
            allocate(keys(0))
        end if
    end function get_keys

    !>  \brief  convert hash to string
    function hash2str( self ) result( str_out )
        class(hash),      intent(in)  :: self
        character(len=XLONGSTRLEN)    :: str_tmp
        character(len=:), allocatable :: str
        type(string) :: str_out
        integer :: cnt, i, n
        if( self%hash_index > 0 )then
            write(str_tmp,*)(self%keys(i)%to_char(),'=',self%values(i),'/', i=1,self%hash_index)
            n   = len_trim(str_tmp)
            str = repeat(' ',n)
            cnt = 0
            do i=1,n
                if( str_tmp(i:i) == ' ' ) cycle
                cnt = cnt+1
                str(cnt:cnt) = str_tmp(i:i)
            enddo
            do i=4,cnt
                if( str(i:i) == '/' ) str(i:i) = ' '
            enddo
            str_out = trim(adjustl(str))
            if( allocated(str) ) deallocate(str)
        endif
    end function hash2str

    !>  \brief  convert hash to string
    pure integer function hash_strlen( self )
        class(hash),      intent(in)  :: self
        character(len=XLONGSTRLEN)    :: str_tmp
        integer :: i
        hash_strlen = 0
        if( self%hash_index > 0 )then
            write(str_tmp,*)(self%keys(i)%to_char(), self%values(i), i=1,self%hash_index)
            do i=1,len_trim(str_tmp)
                if( str_tmp(i:i) == ' ' ) cycle
                hash_strlen = hash_strlen + 1
            enddo
            hash_strlen = hash_strlen + self%hash_index   ! for '=' separator
            hash_strlen = hash_strlen + self%hash_index-1 ! for ' ' separator
        endif
    end function hash_strlen

    !>  \brief  returns size of hash
    pure integer function size_of( self )
        class(hash), intent(in) :: self
        size_of = self%hash_index
    end function size_of

    ! I/O

    !>  \brief  prints the hash
    subroutine print( self )
        class(hash), intent(in) :: self
        integer :: i
        if( self%hash_index > 0 )then
            do i=1,self%hash_index-1,1
                write(logfhandle,"(1X,A,A)", advance="no") self%keys(i)%to_char(), '='
                write(logfhandle,"(A)", advance="no") trim(real2str(self%values(i)))
            end do
            write(logfhandle,"(1X,A,A)", advance="no") self%keys(self%hash_index)%to_char(), '='
            write(logfhandle,"(A)") trim(real2str(self%values(self%hash_index)))
        endif
    end subroutine print

    ! DESTRUCTORS

    subroutine kill( self )
        class(hash), intent(inout) :: self
        if( self%exists )then
            call self%dealloc_hash
            self%buffsz     = 0
            self%hash_index = 0
            self%exists     = .false.
        endif
    end subroutine kill

    subroutine dealloc_hash( self )
        class(hash), intent(inout) :: self
        integer :: i
        if( allocated(self%keys) )then
            do i=1,size(self%keys)
                call self%keys(i)%kill
            end do
            deallocate(self%keys)
        endif
        if(allocated(self%values)) deallocate(self%values)
    end subroutine dealloc_hash

end module simple_hash
