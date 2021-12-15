! character hash
module simple_chash
use simple_defs
use simple_error,   only: simple_exception
use simple_fileio,  only: fopen, fileiochk, fclose
use simple_strings, only: strisblank, striscomment, lexsort
use simple_syslib,  only: is_open
use simple_ansi_ctrls
implicit none

public :: chash
private
#include "simple_local_flags.inc"

integer, parameter :: BUFFSZ_DEFAULT = 5

!> Simple chash type
type :: chash
    private
    type(str4arr), allocatable :: keys(:)   !< chash keys
    type(str4arr), allocatable :: values(:) !< chash values
    integer :: buffsz      = BUFFSZ_DEFAULT !< size of first buffer and subsequent allocation increments
    integer :: chash_index = 0              !< current highest index in hash
    logical :: exists      = .false.        !< to indicate existence
  contains
    !< CONSTRUCTORS
    procedure, private :: new_1
    procedure, private :: new_2
    generic            :: new => new_1, new_2
    procedure, private :: alloc_chash
    procedure, private :: realloc_chash
    procedure, private :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: parse_cmdline
    procedure          :: gen_job_descr
    !< SETTERS
    procedure          :: push
    procedure          :: set
    procedure          :: delete
    !< GETTERS
    procedure          :: isthere
    procedure          :: lookup
    procedure          :: reverselookup
    procedure, private :: get_1
    procedure, private :: get_2
    generic            :: get => get_1, get_2
    procedure          :: get_static
    procedure          :: get_key
    procedure          :: chash2str
    procedure          :: chash_strlen
    procedure          :: size_of
    !< MODIFIERS
    procedure          :: sort
    !< PRINTERS
    procedure, private :: print_key_val_pair_1
    procedure, private :: print_key_val_pair_2
    procedure          :: print_key_val_pairs
    !< I/O
    procedure          :: read
    procedure          :: write
    !< DESTRUCTORS
    procedure          :: kill
    procedure, private :: dealloc_chash
end type chash

interface chash
    module procedure constructor_1
    module procedure constructor_2
end interface chash

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor_1( ) result( self )
        type(chash) :: self
        call self%new_1
    end function constructor_1

    !>  \brief  is a constructor
    function constructor_2( sz, buffsz ) result( self )
        integer,           intent(in)    :: sz     !< initial size
        integer, optional, intent(in)    :: buffsz !< size of subsequent allocation increments
        type(chash) :: self
        call self%new_2(sz, buffsz)
    end function constructor_2

    !>  \brief  is a constructor
    subroutine new_1( self )
        class(chash), intent(inout) :: self
        call self%kill
        self%buffsz = BUFFSZ_DEFAULT
        call self%alloc_chash(self%buffsz)
        self%exists = .true.
    end subroutine new_1

    !>  \brief  is a constructor
    subroutine new_2( self, sz, buffsz )
        class(chash),      intent(inout) :: self   !< instance
        integer,           intent(in)    :: sz     !< initial size
        integer, optional, intent(in)    :: buffsz !< size of subsequent allocation increments
        call self%kill
        if( present(buffsz) )then
            self%buffsz = buffsz
        else
            self%buffsz = BUFFSZ_DEFAULT
        endif
        call self%alloc_chash(sz)
        self%exists = .true.
    end subroutine new_2

    !>  \brief  allocates keys and values without touching buffsz and chash_index
    subroutine alloc_chash( self, sz )
        class(chash), intent(inout) :: self !< instance
        integer,      intent(in)    :: sz   !< total size
        call self%dealloc_chash
        allocate(self%keys(sz), self%values(sz))
    end subroutine alloc_chash

    !>  \brief  re-allocates keys and values without touching buffsz and chash_index and preserving previous key-value pairs
    subroutine realloc_chash( self )
        class(chash), intent(inout) :: self
        type(str4arr), allocatable :: keys_copy(:), values_copy(:)
        integer :: newsz, oldsz, i
        if( self%exists )then
            ! old/new size
            oldsz = size(self%keys)
            newsz = oldsz + self%buffsz
            ! copy key-value pairs
            allocate(keys_copy(oldsz), values_copy(oldsz))
            do i=1,oldsz
                if( allocated(self%keys(i)%str)   ) allocate(keys_copy(i)%str,   source=self%keys(i)%str  )
                if( allocated(self%values(i)%str) ) allocate(values_copy(i)%str, source=self%values(i)%str)
            enddo
            ! allocate extended size
            call self%alloc_chash(newsz)
            ! set back the copied key-value pairs
            do i=1,oldsz
                if( allocated(keys_copy(i)%str)   ) allocate(self%keys(i)%str,   source=keys_copy(i)%str  )
                if( allocated(values_copy(i)%str) ) allocate(self%values(i)%str, source=values_copy(i)%str)
            enddo
        else
            THROW_HARD('cannot reallocate non-existent chash; realloc_chash')
        endif
    end subroutine realloc_chash

    !>  \brief  copies a chash instance
    subroutine copy( self_out, self_in )
        class(chash), intent(inout) :: self_out
        class(chash), intent(in)    :: self_in
        integer :: i, sz
        sz = size(self_in%keys)
        call self_out%alloc_chash(sz)
        self_out%buffsz      = self_in%buffsz
        self_out%chash_index = self_in%chash_index
        self_out%exists      = .true.
        if( self_in%chash_index > 0 )then
            do i=1,self_in%chash_index
                if( allocated(self_in%keys(i)%str)   ) allocate(self_out%keys(i)%str,   source=self_in%keys(i)%str  )
                if( allocated(self_in%values(i)%str) ) allocate(self_out%values(i)%str, source=self_in%values(i)%str)
            end do
        endif
    end subroutine copy

    !>  \brief  is a polymorphic assigner
    subroutine assign( self_out, self_in )
        class(chash), intent(inout) :: self_out
        class(chash), intent(in)    :: self_in
        call self_out%copy(self_in)
    end subroutine assign

    !>  \brief  parse the ke-value pairs in the command line into chash
    subroutine parse_cmdline( self )
        class(chash), intent(inout) :: self
        character(len=STDLEN)       :: arg
        integer :: cmdstat, cmdlen, i, argcnt, pos1
        argcnt  = command_argument_count()
        if( .not. self%exists ) call self%new
        do i=1,argcnt
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(logfhandle,*) 'ERROR! while parsing the command line; simple_chash :: parse_cmdline'
                write(logfhandle,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(logfhandle,*) 'which likely exceeds the length limit STDLEN'
                write(logfhandle,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            pos1 = index(arg, '=') ! position of '='
            ! parse everyting containing '='
            if( pos1 /= 0 ) call self%push(arg(:pos1-1), trim(arg(pos1+1:)))
        end do
    end subroutine parse_cmdline

    !>  \brief  for generating a job description for distributed execution
    subroutine gen_job_descr( self, prg )
        class(chash),    intent(inout) :: self
        character(len=*), intent(in)    :: prg
        call self%parse_cmdline
        call self%set('prg', prg)
    end subroutine gen_job_descr

    ! SETTERS

    !>  \brief  pushes values to the chash
    subroutine push( self, key, val )
        class(chash),     intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        integer :: sz
        ! increment index
        self%chash_index = self%chash_index + 1
        ! reallocate if needed
        sz = size(self%keys)
        if( self%chash_index > sz ) call self%realloc_chash
        ! set key
        allocate(self%keys(self%chash_index)%str, source= trim(adjustl(key)))
        ! set value
        allocate(self%values(self%chash_index)%str, source=trim(adjustl(val)))
    end subroutine push

    !>  \brief  sets a value in the chash
    subroutine set( self, key, val )
        class(chash),     intent(inout) :: self
        character(len=*), intent(in)    :: key, val
        integer :: i
        if( self%chash_index >= 1 )then
            do i=1,self%chash_index
                ! replace if there (different from standard hash table)
                if( trim(self%keys(i)%str) .eq. trim(key) )then
                    if( allocated(self%values(i)%str) ) deallocate(self%values(i)%str)
                    allocate(self%values(i)%str, source=trim(adjustl(val)))
                    return
                endif
            end do
        endif
        ! otherwise, push
        call self%push(key, val)
    end subroutine set

    !>  \brief  deletes a value in the chash
    subroutine delete( self, key )
        class(chash),     intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: i, ind
        character(len=:), allocatable :: tmp
        ind = self%lookup( key )
        if( ind == 0 .or. ind > self%chash_index ) return
        do i=ind,self%chash_index - 1
            ! replace key
            if( allocated(self%keys(i)%str) ) deallocate(self%keys(i)%str)
            if(allocated(tmp)) deallocate(tmp)
            allocate(tmp, source=self%keys(i + 1)%str)
            self%keys(i)%str = tmp
            deallocate(tmp)
            ! replace value
            if( allocated(self%values(i)%str) ) deallocate(self%values(i)%str)
            allocate(tmp, source=self%values(i + 1)%str)
            self%values(i)%str = tmp
        enddo
        ! remove previous last entry
        if( allocated(self%keys(self%chash_index)%str)   ) deallocate(self%keys(self%chash_index)%str)
        if( allocated(self%values(self%chash_index)%str) ) deallocate(self%values(self%chash_index)%str)
        ! decrement index
        self%chash_index = self%chash_index - 1
    end subroutine delete

    ! GETTERS

    !>  \brief  check for presence of key in the chash
    pure function isthere( self, key ) result( found )
        class(chash),     intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: i
        logical :: found
        found = .false.
        if( self%chash_index >= 1 )then
            do i=1,self%chash_index
                if( trim(self%keys(i)%str) .eq. trim(key) )then
                    found = .true.
                    return
                endif
            end do
        endif
    end function isthere

    !>  \brief  looks up the index of a key value in the chash
    function lookup( self, key ) result( which )
        class(chash),     intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: ikey, which
        which = 0
        do ikey=1,self%chash_index
            if( trim(key) .eq. trim(self%keys(ikey)%str) )then
                which = ikey
                return
            endif
        end do
    end function lookup

    !>  \brief  looks up the index of a value in the chash
    !!* Warning* this will get the first occurance of the value
    !!
    function reverselookup( self, key ) result( which )
        class(chash),     intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: ikey, which
        which = 0
        do ikey=1,self%chash_index
            if( trim(key) .eq. trim(self%values(ikey)%str) )then
                which = ikey
                return
            endif
        end do
    end function reverselookup

    !>  \brief  gets a value in the chash
    function get_1( self, key ) result( val )
        class(chash),     intent(in)  :: self
        character(len=*), intent(in)  :: key
        character(len=:), allocatable :: val
        integer :: i
        do i=1,self%chash_index
            if( trim(self%keys(i)%str) .eq. trim(key) )then
                if( allocated(self%values(i)%str) ) allocate(val, source=trim(self%values(i)%str))
                return
            endif
        end do
    end function get_1

    !>  \brief  gets a value in the chash
    function get_2( self, ival ) result( val )
        class(chash), intent(in)      :: self
        integer,      intent(in)      :: ival
        character(len=:), allocatable :: val
        if( allocated(self%values(ival)%str) ) allocate(val, source=trim(self%values(ival)%str))
    end function get_2

    !>  \brief  returns fixed length value from the chash
    function get_static( self, key ) result( val )
        class(chash),     intent(in)  :: self
        character(len=*), intent(in)  :: key
        character(len=STDLEN)         :: val
        integer :: i
        do i=1,self%chash_index
            if( trim(self%keys(i)%str) .eq. trim(key) )then
                val = trim(self%values(i)%str)
                return
            endif
        end do
        val = ''
    end function get_static

    !>  \brief  gets a key in the chash
    function get_key( self, ikey ) result( val )
        class(chash), intent(in)      :: self
        integer,      intent(in)      :: ikey
        character(len=:), allocatable :: val
        if( allocated(self%keys(ikey)%str) ) allocate(val, source=trim(self%keys(ikey)%str))
    end function get_key

    !>  \brief  concatenates the chash into a string
    pure function chash2str( self ) result( str )
        class(chash), intent(in)  :: self
        character(len=XLONGSTRLEN) :: str
        integer :: i, len
        if( self%chash_index > 0 )then
            if( self%chash_index == 1 )then
                str = trim(self%keys(1)%str)//'='//trim(self%values(1)%str)
                return
            endif
            str = trim(self%keys(1)%str)//'='//trim(self%values(1)%str)//' '
            if( self%chash_index > 2 )then
                do i=2,self%chash_index-1
                    len = len_trim(str)
                    str = str(1:len+1)//trim(self%keys(i)%str)//'='//trim(self%values(i)%str)//' '
                end do
            endif
            len = len_trim(str)
            str = trim(str(1:len+1)//trim(self%keys(self%chash_index)%str)//'='//&
                &trim(self%values(self%chash_index)%str))
        endif
    end function chash2str

    !>  \brief  convert hash to string
    pure integer function chash_strlen( self )
        class(chash),intent(in)    :: self
        integer :: i
        chash_strlen = 0
        if( self%chash_index == 0 )then
            return
        elseif( self%chash_index == 1 )then
            chash_strlen = len_trim(self%keys(1)%str)+len_trim(self%values(1)%str)
            chash_strlen = chash_strlen + 1 ! for '=' separator
        else
            do i=1,self%chash_index
                chash_strlen = chash_strlen+len_trim(self%keys(i)%str)+len_trim(self%values(i)%str)
            end do
            chash_strlen = chash_strlen + self%chash_index ! for '=' separator
            chash_strlen = chash_strlen + self%chash_index-1 ! for ' ' separator
        endif
    end function chash_strlen

    !>  \brief  returns size of chash
    pure integer function size_of( self )
        class(chash), intent(in) :: self
        size_of = self%chash_index
    end function size_of

    ! MODIFIERS

    !>  \brief  sort chash lexographically based on the keys
    subroutine sort( self )
        class(chash),          intent(inout) :: self
        character(len=KEYLEN), allocatable   :: keys_sorted(:)
        character(len=:),      allocatable   :: val
        type(chash) :: tmp
        integer     :: ikey, sz
        if( self%chash_index > 1 )then
            ! make conforming temporary chash
            sz = size(self%keys)
            call tmp%alloc_chash(sz)
            tmp%buffsz      = self%buffsz
            tmp%chash_index = self%chash_index
            tmp%exists      = .true.
            ! fill in keys
            allocate(keys_sorted(self%chash_index))
            do ikey=1,self%chash_index
                keys_sorted(ikey) = self%keys(ikey)%str
            end do
            ! sort keys
            call lexSort(keys_sorted)
            ! generate the sorted hash
            do ikey=1,self%chash_index
                val = self%get(keys_sorted(ikey))
                call tmp%push(keys_sorted(ikey), val)
                deallocate(val)
            end do
            self = tmp
            call tmp%kill
        endif
    end subroutine sort

    ! PRINTERS

    !>  \brief  pretty printing of a key-value pair
    subroutine print_key_val_pair_1( self, key, maxlen, fhandle, ansi_flag, latex_bold )
        class(chash),               intent(in) :: self
        character(len=*),           intent(in) :: key
        integer,                    intent(in) :: maxlen
        integer,                    intent(in) :: fhandle
        character(len=*), optional, intent(in) :: ansi_flag
        logical,          optional, intent(in) :: latex_bold
        character(len=:), allocatable :: aansi_flag, str2print
        character(len=maxlen) :: key_padded
        integer               :: which
        logical               :: latex_format
        if( present(ansi_flag) )then
            allocate(aansi_flag, source=trim(ansi_flag))
        else
            allocate(aansi_flag, source=C_WHITE)
        endif
        latex_format = present(latex_bold)
        which = self%lookup(key)
        if( which > 0 )then
            ! print table (for easy to read instructions)
            key_padded = trim(self%keys(which)%str)
            do while(len(key_padded) < maxlen)
                key_padded = trim(key_padded)//' '
            end do
            if( latex_format )then
                if( latex_bold )then
                    allocate(str2print, source='+textbf['//key_padded//' = '//trim(self%values(which)%str)//']')
                else
                    allocate(str2print, source=key_padded//' = '//trim(self%values(which)%str))
                endif
            else
                allocate(str2print, source=format_str(key_padded//' = '//trim(self%values(which)%str), aansi_flag))
            endif
            write(fhandle, '(a,1x,a)') str2print
        else
            THROW_HARD('key: '//trim(key)//' does not exist in the chash; print_key_val_pair_1')
        endif
    end subroutine print_key_val_pair_1

    !>  \brief  pretty printing of a key-value pair
    subroutine print_key_val_pair_2( self, ikey, maxlen, fhandle, ansi_flag, latex_bold )
        class(chash), intent(in)      :: self
        integer,                    intent(in) :: ikey
        integer,                    intent(in) :: maxlen
        integer,                    intent(in) :: fhandle
        character(len=*), optional, intent(in) :: ansi_flag
        logical,          optional, intent(in) :: latex_bold
        character(len=maxlen)         :: key_padded
        character(len=:), allocatable :: aansi_flag, str2print
        logical :: latex_format
        if( present(ansi_flag) )then
            allocate(aansi_flag, source=trim(ansi_flag))
        else
            allocate(aansi_flag, source=C_WHITE)
        endif
        latex_format = present(latex_bold)
        key_padded   = trim(self%keys(ikey)%str)
        do while(len(key_padded) < maxlen)
            key_padded = trim(key_padded)//' '
        end do
        if( latex_format )then
            if( latex_bold )then
                allocate(str2print, source='+textbf['//key_padded//' = '//trim(self%values(ikey)%str)//']')
            else
                allocate(str2print, source=key_padded//' = '//trim(self%values(ikey)%str))
            endif
        else
            allocate(str2print, source=format_str(key_padded//' = '//trim(self%values(ikey)%str), aansi_flag))
        endif
        write(fhandle, '(a,1x,a)') str2print
    end subroutine print_key_val_pair_2

    !>  \brief  pretty printing of all key-value  pairs
    subroutine print_key_val_pairs( self, fhandle, keys2print, mask, latex )
        class(chash),               intent(in) :: self
        integer,                    intent(in) :: fhandle
        character(len=*), optional, intent(in) :: keys2print(:)
        logical,          optional, intent(in) :: mask(:) ! mask true -> C_BOLD
        logical,          optional, intent(in) :: latex
        logical, allocatable :: mmask(:)
        logical :: present_mask, llatex
        integer :: maxlen, ikey, sz, keylen
        maxlen  = 0
        ikey    = 0
        sz      = 0
        keylen  = 0
        present_mask = present(mask)
        llatex  = .false.
        if( present(latex) ) llatex = latex
        if( present(keys2print) )then
            sz = size(keys2print)
            if( present_mask )then
                if( size(mask) /= sz ) THROW_HARD('Nonconforming sizes keys2print / mask; print_key_val_pairs')
                allocate(mmask(sz), source=mask)
            else
                allocate(mmask(sz), source=.false.)
            endif
            ! find max key len in keys2print
            maxlen = 0
            do ikey=1,sz
                keylen = len_trim(keys2print(ikey))
                if( keylen > maxlen ) maxlen = keylen
            end do
            ! print selected key-value pairs
            do ikey=1,sz
                if( llatex )then
                    call self%print_key_val_pair_1(keys2print(ikey), maxlen, fhandle, latex_bold=mmask(ikey))
                else
                    if( mmask(ikey) )then
                        call self%print_key_val_pair_1(keys2print(ikey), maxlen, fhandle, ansi_flag=C_BOLD)
                    else
                        call self%print_key_val_pair_1(keys2print(ikey), maxlen, fhandle)
                    endif
                endif
            end do
        else
            if( self%chash_index > 0 )then
                if( present_mask )then
                    if( size(mask) /= self%chash_index ) THROW_HARD('Nonconforming sizes self%chash_index / size(mask); print_key_val_pairs')
                    allocate(mmask(self%chash_index), source=mask)
                else
                    allocate(mmask(self%chash_index), source=.false.)
                endif
                ! find max key len among all keys in chash
                maxlen = 0
                do ikey=1,self%chash_index
                    keylen = len_trim(self%keys(ikey)%str)
                    if( keylen > maxlen ) maxlen = keylen
                end do
                ! print all key-value pairs
                do ikey=1,self%chash_index
                    if( llatex )then
                        call self%print_key_val_pair_2(ikey, maxlen, fhandle, latex_bold=mmask(ikey))
                    else
                        if( mmask(ikey) )then
                            call self%print_key_val_pair_2(ikey, maxlen, fhandle, ansi_flag=C_BOLD)
                        else
                            call self%print_key_val_pair_2(ikey, maxlen, fhandle)
                        endif
                    endif
                end do
            endif
        endif
    end subroutine print_key_val_pairs

    ! I/O

    !>  \brief  for reading key-vals from a text file
    subroutine read( self, fname )
        class(chash), intent(inout)  :: self   !< instance
        character(len=*), intent(in) :: fname  !< name of file
        character(len=LINE_MAX_LEN)  :: buffer !< will hold a line from the file
        character(len=KEYLEN)        :: key
        character(len=STDLEN)        :: val
        integer :: ios, pos, funit
        call fopen(funit, file=fname, iostat=ios, status='old')
        call fileiochk('simple_chash :: read; Error when opening file for reading: '//trim(fname), ios)
        do
            read(funit, '(a)', iostat=ios) buffer
            if ( ios == 0 ) then
                if( strIsComment(buffer) .or. strIsBlank(buffer) )then
                    ! don't do anything
                else
                    pos = index(buffer, '=') ! position of '='
                    key = trim(adjustl(buffer(:pos-1)))
                    val = trim(adjustl(buffer(pos+1:)))
                    if( .not. strIsBlank(val) ) call self%set(key, val)
                endif
            else
                exit ! we found the end...
            endif
        enddo
        call fclose(funit)
    end subroutine read

    !>  \brief  for writing key-vals to a text file
    subroutine write( self, fname )
        class(chash), intent(inout)  :: self   !< instance
        character(len=*), intent(in) :: fname  !< name of file
        integer :: ios, funit
        call fopen(funit, fname, iostat=ios, status='replace')
        call fileiochk('simple_chash :: write; Error when opening file for writing: '&
            //trim(fname), ios)
        call self%print_key_val_pairs(funit)
        call fclose(funit)
    end subroutine write

    ! DESTRUCTORS

    subroutine kill( self )
        class(chash), intent(inout) :: self
        if( self%exists )then
            call self%dealloc_chash
            self%buffsz      = 0
            self%chash_index = 0
            self%exists      = .false.
        endif
    end subroutine kill

    subroutine dealloc_chash( self )
        class(chash), intent(inout) :: self
        integer :: i
        if( allocated(self%keys) )then
            do i=1,size(self%keys)
                if( allocated(self%keys(i)%str) ) deallocate(self%keys(i)%str)
            end do
            deallocate(self%keys)
        endif
        if( allocated(self%values) )then
            do i=1,size(self%values)
                if( allocated(self%values(i)%str) ) deallocate(self%values(i)%str)
            end do
            deallocate(self%values)
        endif
    end subroutine dealloc_chash

end module simple_chash
