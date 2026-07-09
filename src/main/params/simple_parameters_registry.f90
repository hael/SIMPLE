!@descr: typed registry used to bind command-line arguments onto the flat SIMPLE parameters object
module simple_parameters_registry
use simple_core_module_api
use simple_cmdline, only: cmdline
implicit none
private
#include "simple_local_flags.inc"

public :: param_registry

type :: char_binding
    character(len=KEYLEN) :: key = ''
    character(len=:), pointer :: var => null()
end type char_binding

type :: string_binding
    character(len=KEYLEN) :: key = ''
    type(string), pointer :: var => null()
end type string_binding

type :: int_binding
    character(len=KEYLEN) :: key = ''
    integer, pointer :: var => null()
end type int_binding

type :: real_binding
    character(len=KEYLEN) :: key = ''
    real, pointer :: var => null()
end type real_binding

type :: file_binding
    character(len=KEYLEN) :: key = ''
    type(string), pointer :: var => null()
    character(len=1)      :: allowed1 = ' '
    character(len=1)      :: allowed2 = ' '
    character(len=1)      :: allowed3 = ' '
    character(len=1)      :: not_allowed = ' '
    logical               :: has_allowed1 = .false.
    logical               :: has_allowed2 = .false.
    logical               :: has_allowed3 = .false.
    logical               :: has_not_allowed = .false.
end type file_binding

type :: dir_binding
    character(len=KEYLEN) :: key = ''
    type(string), pointer :: var => null()
end type dir_binding

type :: param_registry
    private
    integer :: nchars = 0
    integer :: nstrings = 0
    integer :: nints = 0
    integer :: nreals = 0
    integer :: nfiles = 0
    integer :: ndirs = 0
    type(char_binding),   allocatable :: char_bindings(:)
    type(string_binding), allocatable :: string_bindings(:)
    type(int_binding),    allocatable :: int_bindings(:)
    type(real_binding),   allocatable :: real_bindings(:)
    type(file_binding),   allocatable :: file_bindings(:)
    type(dir_binding),    allocatable :: dir_bindings(:)
contains
    procedure :: clear
    procedure, private :: add_char_fixed
    procedure, private :: add_char_string
    procedure, private :: ensure_char_capacity
    procedure, private :: ensure_string_capacity
    procedure, private :: ensure_int_capacity
    procedure, private :: ensure_real_capacity
    procedure, private :: ensure_file_capacity
    procedure, private :: ensure_dir_capacity
    generic   :: add_char => add_char_fixed, add_char_string
    procedure :: add_int
    procedure :: add_real
    procedure :: add_file
    procedure :: add_dir
    procedure :: apply_cmdline
end type param_registry

contains

    subroutine clear(self)
        class(param_registry), intent(inout) :: self
        if( allocated(self%char_bindings)   ) deallocate(self%char_bindings)
        if( allocated(self%string_bindings) ) deallocate(self%string_bindings)
        if( allocated(self%int_bindings)    ) deallocate(self%int_bindings)
        if( allocated(self%real_bindings)   ) deallocate(self%real_bindings)
        if( allocated(self%file_bindings)   ) deallocate(self%file_bindings)
        if( allocated(self%dir_bindings)    ) deallocate(self%dir_bindings)
        self%nchars   = 0
        self%nstrings = 0
        self%nints    = 0
        self%nreals   = 0
        self%nfiles   = 0
        self%ndirs    = 0
    end subroutine clear

    subroutine add_char_fixed(self, key, var)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        character(len=*), target, intent(inout) :: var
        self%nchars = self%nchars + 1
        call self%ensure_char_capacity(self%nchars)
        self%char_bindings(self%nchars)%key = trim(key)
        self%char_bindings(self%nchars)%var => var
    end subroutine add_char_fixed

    subroutine add_char_string(self, key, var)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        type(string), target,  intent(inout) :: var
        self%nstrings = self%nstrings + 1
        call self%ensure_string_capacity(self%nstrings)
        self%string_bindings(self%nstrings)%key = trim(key)
        self%string_bindings(self%nstrings)%var => var
    end subroutine add_char_string

    subroutine add_int(self, key, var)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        integer, target,       intent(inout) :: var
        self%nints = self%nints + 1
        call self%ensure_int_capacity(self%nints)
        self%int_bindings(self%nints)%key = trim(key)
        self%int_bindings(self%nints)%var => var
    end subroutine add_int

    subroutine add_real(self, key, var)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        real, target,          intent(inout) :: var
        self%nreals = self%nreals + 1
        call self%ensure_real_capacity(self%nreals)
        self%real_bindings(self%nreals)%key = trim(key)
        self%real_bindings(self%nreals)%var => var
    end subroutine add_real

    subroutine add_file(self, key, var, allowed1, allowed2, allowed3, notAllowed)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        type(string), target,  intent(inout) :: var
        character(len=1), optional, intent(in) :: allowed1, allowed2, allowed3, notAllowed
        self%nfiles = self%nfiles + 1
        call self%ensure_file_capacity(self%nfiles)
        self%file_bindings(self%nfiles)%key = trim(key)
        self%file_bindings(self%nfiles)%var => var
        if( present(allowed1) )then
            self%file_bindings(self%nfiles)%allowed1 = allowed1
            self%file_bindings(self%nfiles)%has_allowed1 = .true.
        endif
        if( present(allowed2) )then
            self%file_bindings(self%nfiles)%allowed2 = allowed2
            self%file_bindings(self%nfiles)%has_allowed2 = .true.
        endif
        if( present(allowed3) )then
            self%file_bindings(self%nfiles)%allowed3 = allowed3
            self%file_bindings(self%nfiles)%has_allowed3 = .true.
        endif
        if( present(notAllowed) )then
            self%file_bindings(self%nfiles)%not_allowed = notAllowed
            self%file_bindings(self%nfiles)%has_not_allowed = .true.
        endif
    end subroutine add_file

    subroutine add_dir(self, key, var)
        class(param_registry), intent(inout) :: self
        character(len=*),      intent(in)    :: key
        type(string), target,  intent(inout) :: var
        self%ndirs = self%ndirs + 1
        call self%ensure_dir_capacity(self%ndirs)
        self%dir_bindings(self%ndirs)%key = trim(key)
        self%dir_bindings(self%ndirs)%var => var
    end subroutine add_dir

    subroutine ensure_char_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(char_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%char_bindings) )then
            old_size = size(self%char_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%char_bindings
        call move_alloc(tmp, self%char_bindings)
    end subroutine ensure_char_capacity

    subroutine ensure_string_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(string_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%string_bindings) )then
            old_size = size(self%string_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%string_bindings
        call move_alloc(tmp, self%string_bindings)
    end subroutine ensure_string_capacity

    subroutine ensure_int_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(int_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%int_bindings) )then
            old_size = size(self%int_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%int_bindings
        call move_alloc(tmp, self%int_bindings)
    end subroutine ensure_int_capacity

    subroutine ensure_real_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(real_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%real_bindings) )then
            old_size = size(self%real_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%real_bindings
        call move_alloc(tmp, self%real_bindings)
    end subroutine ensure_real_capacity

    subroutine ensure_file_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(file_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%file_bindings) )then
            old_size = size(self%file_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%file_bindings
        call move_alloc(tmp, self%file_bindings)
    end subroutine ensure_file_capacity

    subroutine ensure_dir_capacity(self, needed)
        class(param_registry), intent(inout) :: self
        integer,               intent(in)    :: needed
        type(dir_binding), allocatable :: tmp(:)
        integer :: new_size, old_size
        if( allocated(self%dir_bindings) )then
            old_size = size(self%dir_bindings)
        else
            old_size = 0
        endif
        if( needed <= old_size ) return
        new_size = max(needed, max(16, 2*old_size))
        allocate(tmp(new_size))
        if( old_size > 0 ) tmp(1:old_size) = self%dir_bindings
        call move_alloc(tmp, self%dir_bindings)
    end subroutine ensure_dir_capacity

    subroutine apply_cmdline(self, cline, cntfile, checkupfile)
        class(param_registry), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        integer, optional,     intent(inout) :: cntfile
        character(len=1), optional, intent(inout) :: checkupfile(:)
        type(string)     :: str, abspath_name
        character(len=1) :: file_descr
        logical          :: raise_exception
        integer          :: i
        do i=1,self%nchars
            if( cline%defined(self%char_bindings(i)%key) )then
                str = cline%get_carg(self%char_bindings(i)%key)
                self%char_bindings(i)%var = str%to_char()
            endif
        end do
        do i=1,self%nstrings
            if( cline%defined(self%string_bindings(i)%key) )then
                self%string_bindings(i)%var = cline%get_carg(self%string_bindings(i)%key)
            endif
        end do
        do i=1,self%nints
            if( cline%defined(self%int_bindings(i)%key) )then
                self%int_bindings(i)%var = cline%get_iarg(self%int_bindings(i)%key)
            endif
        end do
        do i=1,self%nreals
            if( cline%defined(self%real_bindings(i)%key) )then
                self%real_bindings(i)%var = cline%get_rarg(self%real_bindings(i)%key)
            endif
        end do
        do i=1,self%nfiles
            if( .not. present(cntfile) .or. .not. present(checkupfile) )then
                THROW_HARD('File bindings require cntfile and checkupfile; simple_parameters_registry')
            endif
            if( cline%defined(self%file_bindings(i)%key) )then
                self%file_bindings(i)%var = cline%get_carg(self%file_bindings(i)%key)
                file_descr      = fname2format(self%file_bindings(i)%var)
                raise_exception = .false.
                if( self%file_bindings(i)%has_allowed1 )then
                    if( self%file_bindings(i)%has_allowed2 .or. self%file_bindings(i)%has_allowed3 )then
                        if( self%file_bindings(i)%allowed1 == file_descr .or. self%file_bindings(i)%allowed2 == file_descr .or. &
                            &self%file_bindings(i)%allowed3 == file_descr )then
                        else
                            raise_exception = .true.
                        endif
                    else
                        if( self%file_bindings(i)%allowed1 /= file_descr ) raise_exception = .true.
                    endif
                endif
                if( self%file_bindings(i)%has_not_allowed )then
                    if( self%file_bindings(i)%not_allowed == file_descr ) raise_exception = .true.
                endif
                if( raise_exception )then
                    write(logfhandle,*) 'This format: ', file_descr, ' is not allowed for this file: ', self%file_bindings(i)%var%to_char()
                    write(logfhandle,*) 'flag:', trim(self%file_bindings(i)%key)
                    stop
                endif
                select case(file_descr)
                    case('I')
                        THROW_HARD('Support for IMAGIC files is not implemented!')
                    case('M')
                        cntfile = cntfile + 1
                        checkupfile(cntfile) = 'M'
                    case('S')
                        cntfile = cntfile + 1
                        checkupfile(cntfile) = 'S'
                    case('J')
                        cntfile = cntfile + 1
                        checkupfile(cntfile) = 'J'
                    case('L')
                        cntfile = cntfile + 1
                        checkupfile(cntfile) = 'L'
                    case('K')
                        cntfile = cntfile + 1
                        checkupfile(cntfile) = 'K'
                    case('T','B','P','O','R')
                    case('N')
                        write(logfhandle,*) 'file: ', self%file_bindings(i)%var%to_char()
                        THROW_HARD('This file format is not supported by SIMPLE')
                    case DEFAULT
                        write(logfhandle,*) 'file: ', self%file_bindings(i)%var%to_char()
                        THROW_HARD('This file format is not supported by SIMPLE')
                end select
                if( file_exists(self%file_bindings(i)%var) )then
                    abspath_name = simple_abspath(self%file_bindings(i)%var, check_exists=.false.)
                    self%file_bindings(i)%var = abspath_name%to_char()
                    call cline%set(self%file_bindings(i)%key, self%file_bindings(i)%var)
                    call abspath_name%kill
                endif
            endif
        end do
        do i=1,self%ndirs
            if( cline%defined(self%dir_bindings(i)%key) )then
                self%dir_bindings(i)%var = cline%get_carg(self%dir_bindings(i)%key)
                if( file_exists(self%dir_bindings(i)%var) )then
                    abspath_name = simple_abspath(self%dir_bindings(i)%var, check_exists=.false.)
                    self%dir_bindings(i)%var = abspath_name%to_char()
                    call cline%set(self%dir_bindings(i)%key, self%dir_bindings(i)%var)
                    call abspath_name%kill
                endif
            endif
        end do
    end subroutine apply_cmdline

end module simple_parameters_registry
