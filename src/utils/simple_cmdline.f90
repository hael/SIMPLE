! command line parser
module simple_cmdline
include 'simple_lib.f08'
use simple_user_interface ! use all in there
use simple_private_prgs   ! use all in there
use simple_args, only: args
implicit none
private
public :: cmdline, cmdline_err
#include "simple_local_flags.inc"

integer, parameter :: MAX_CMDARGS = 100
logical, parameter :: DEBUG_HERE  = .false.

!> cmdarg key/value pair command-line type
type cmdarg
    character(len=:), allocatable :: key, carg
    logical                       :: defined=.false.
    real(dp)                      :: rarg=0.d0
end type cmdarg

type cmdline
    private
    type(cmdarg)          :: cmds(MAX_CMDARGS)
    character(len=KEYLEN) :: checker(MAX_CMDARGS)
    character(len=8192)   :: entire_line        !! PATH_MAX is 4096 https://github.com/torvalds/linux/blob/master/include/uapi/linux/limits.h
    integer               :: argcnt=0, ncheck=0 !! max system args `getconf ARG_MAX`
  contains
    procedure          :: parse
    procedure          :: parse_private
    procedure          :: parse_oldschool
    procedure, private :: parse_command_line_value
    procedure, private :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: set_1, set_2, set_3, set_4
    generic            :: set => set_1, set_2, set_3, set_4
    procedure, private :: lookup
    procedure          :: get_keys
    procedure          :: delete
    procedure          :: checkvar
    procedure          :: check
    procedure          :: printline
    procedure          :: readline
    procedure          :: writeline
    procedure          :: defined
    procedure          :: get_rarg
    procedure          :: get_iarg
    procedure          :: get_carg
    procedure          :: gen_job_descr
    procedure          :: kill
end type cmdline

contains

    !> \brief for parsing the command line arguments passed as key=val
    subroutine parse( self )
        class(cmdline), intent(inout) :: self
        type(str4arr), allocatable    :: keys_required(:)
        type(args)                    :: allowed_args
        type(simple_program), pointer :: ptr2prg => null()
        character(len=LONGSTRLEN)     :: arg
        character(len=STDLEN)         :: exec_cmd
        character(len=:), allocatable :: prgname, exec_cmd_ui
        integer :: i, cmdstat, cmdlen, ikey, pos, nargs_required, sz_keys_req
        ! parse command line
        self%argcnt = command_argument_count()
        call get_command(self%entire_line)
        cmdline_glob = trim(self%entire_line)
        if( .not. str_has_substr(self%entire_line,'prg=') ) THROW_HARD('prg=flag must be set on command line')
        call get_command_argument(0,exec_cmd)
        if( DEBUG_HERE ) print *, 'exec_cmd from get_command_argument in cmdline class: ', trim(exec_cmd)
        ! parse program name
        call get_command_argument(1, arg, cmdlen, cmdstat)
        pos = index(arg, '=') ! position of '='
        call cmdline_err(cmdstat, cmdlen, arg, pos)
        allocate(prgname, source=arg(pos+1:))
        if( DEBUG_HERE ) print *, 'prgname from command-line in cmdline class: ', trim(prgname)
        if( str_has_substr(prgname, 'report_selection') ) prgname = 'selection' ! FIX4NOW
        ! obtain pointer to the program in the simple_user_interface specification
        call get_prg_ptr(prgname, ptr2prg)
        if( .not. associated(ptr2prg) ) THROW_HARD(trim(prgname)//' is not part of any executable in the code base')
        if( DEBUG_HERE ) print *, 'prgname from UI in cmdline class: ', trim(prgname)
        exec_cmd_ui = ptr2prg%get_executable()
        if( DEBUG_HERE ) print *, 'exec_cmd from UI in cmdline class: ', exec_cmd_ui
        if( trim(exec_cmd_ui) .ne. 'all' )then
            if( basename(trim(exec_cmd)) .ne. trim(exec_cmd_ui) )then
                THROW_HARD(trim(prgname)//' is not executed by '//trim(exec_cmd)//' but by '//exec_cmd_ui)
            endif
        endif
        ! list programs if so instructed
        if( str_has_substr(self%entire_line, 'prg=list') )then
            select case(trim(exec_cmd))
                case('simple_exec')
                    call list_simple_prgs_in_ui
                    stop
                case('single_exec')
                    call list_single_prgs_in_ui
                    stop
                case('simple_stream')
                    call list_stream_prgs_in_ui
                    stop
                case DEFAULT
                    THROW_HARD('Program '//prgname//' not supported; parse')
            end select
        endif
        ! describe program if so instructed
        if( str_has_substr(self%entire_line, 'describe=yes') )then
            call ptr2prg%print_prg_descr_long()
            stop
        endif
        if( .not. associated(ptr2prg) )then
            THROW_HARD('inputted prg not supported, use prg=list to list all available programs')
        endif
        ! get required keys
        sz_keys_req = ptr2prg%get_nrequired_keys()
        if( sz_keys_req > 0 )then
            keys_required  = ptr2prg%get_required_keys()
            nargs_required = sz_keys_req + 1 ! +1 because prg part of command line
        else
            nargs_required = 1
        endif
        if( self%argcnt < nargs_required )then
            call ptr2prg%print_cmdline()
            stop
        else
            ! indicate which variables are required
            if( sz_keys_req > 0 )then
                do ikey=1,sz_keys_req
                    if( allocated(keys_required(ikey)%str) ) call self%checkvar(keys_required(ikey)%str, ikey)
                end do
            endif
        endif
        allowed_args = args()
        do i=1,self%argcnt
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(logfhandle,*) 'ERROR! while parsing the command line: simple_cmdline :: parse'
                write(logfhandle,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(logfhandle,*) 'which likely exceeds the length limit LONGSTRLEN'
                write(logfhandle,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            call self%parse_command_line_value(i, arg, allowed_args)
        end do
        if( sz_keys_req > 0 ) call self%check
    end subroutine parse

    !> \brief for parsing the command line arguments passed as key=val
    subroutine parse_private( self )
        class(cmdline), intent(inout) :: self
        type(str4arr), allocatable    :: keys_required(:)
        type(args)                    :: allowed_args
        character(len=LONGSTRLEN)     :: arg
        type(simple_program), pointer :: ptr2prg => null()
        integer :: i, ikey, pos, cmdstat, cmdlen, sz_keys_req, nargs_required
        ! parse command line
        self%argcnt = command_argument_count()
        call get_command(self%entire_line)
        cmdline_glob = trim(self%entire_line)
        if( .not. str_has_substr(self%entire_line,'prg=') ) THROW_HARD('prg= flag must be set on command line')
        ! parse program name
        call get_command_argument(1, arg, cmdlen, cmdstat)
        pos = index(arg, '=') ! position of '='
        call cmdline_err(cmdstat, cmdlen, arg, pos)
        ! obtain pointer to the program in the simple_user_interface specification
        call get_prg_ptr(arg(pos+1:), ptr2prg)
        ! decide wether to print the command line instructions or not
        select case(trim(arg(pos+1:)))
            case('write_ui_json','print_ui_json','print_ui_latex','print_sym_subgroups','print_ui_stream')
                ! no command line arguments
                return
            case DEFAULT
                if( associated(ptr2prg) )then
                    ! get required keys
                    sz_keys_req = ptr2prg%get_nrequired_keys()
                    if( sz_keys_req > 0 )then
                        keys_required  = ptr2prg%get_required_keys()
                        nargs_required = sz_keys_req + 1 ! +1 because prg part of command line
                    else
                        nargs_required = 1
                    endif
                    if( self%argcnt < nargs_required )then
                        call ptr2prg%print_cmdline()
                        stop
                    else
                        ! indicate which variables are required
                        if( sz_keys_req > 0 )then
                            do ikey=1,sz_keys_req
                                if( allocated(keys_required(ikey)%str) ) call self%checkvar(keys_required(ikey)%str, ikey)
                            end do
                        endif
                    endif
                else
                    ! get required keys
                    sz_keys_req = get_n_private_keys_required(trim(arg(pos+1:)))
                    if( sz_keys_req > 0 )then
                        keys_required  = get_private_keys_required(trim(arg(pos+1:)))
                        nargs_required = sz_keys_req + 1 ! +1 because prg part of command line
                    else
                        nargs_required = 1
                    endif
                    if( self%argcnt < nargs_required .or. (sz_keys_req == 0 .and. self%argcnt == 1) )then
                        call print_private_cmdline(arg(pos+1:))
                        stop
                    else
                        ! indicate which variables are required
                        if( sz_keys_req > 0 )then
                            do ikey=1,sz_keys_req
                                if( allocated(keys_required(ikey)%str) ) call self%checkvar(keys_required(ikey)%str, ikey)
                            end do
                        endif
                    endif
                endif
        end select
        allowed_args = args()
        do i=1,self%argcnt
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(logfhandle,*) 'ERROR! while parsing the command line: simple_cmdline :: parse_private'
                write(logfhandle,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(logfhandle,*) 'which likely exceeds the length limit LONGSTRLEN'
                write(logfhandle,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            call self%parse_command_line_value(i, arg, allowed_args)
        end do
        if( sz_keys_req > 0 ) call self%check
    end subroutine parse_private

    !> \brief for parsing the command line arguments passed as key=val
    subroutine parse_oldschool( self, keys_required, keys_optional )
        use simple_args, only: args
        class(cmdline),             intent(inout) :: self
        character(len=*), optional, intent(in)    :: keys_required(:), keys_optional(:)
        character(len=STDLEN)     :: exec_name
        character(len=LONGSTRLEN) :: arg
        type(args)                :: allowed_args
        integer                   :: i, cmdstat, cmdlen, ikey
        integer                   :: nreq, cmdargcnt
        logical                   :: distr_exec
        call get_command_argument(0,exec_name)
        distr_exec = str_has_substr(exec_name,'distr')
        cmdargcnt = command_argument_count()
        call get_command(self%entire_line)
        cmdline_glob = trim(self%entire_line)
        if( present(keys_required) )then
            if( str_has_substr(self%entire_line,'prg=') )then
                nreq = size(keys_required) + 1 ! +1 because prg part of command line
            else
                nreq = size(keys_required)
            endif
            if( cmdargcnt < nreq )then
                call print_cmdline_oldschool(keys_required, keys_optional, distr=distr_exec)
                stop
            else
                ! indicate which variables are required
                do ikey=1,size(keys_required)
                    call self%checkvar(keys_required(ikey), ikey)
                end do
            endif
        else
            if( cmdargcnt < 2 )then
                call print_cmdline_oldschool(keys_required, keys_optional, distr=distr_exec)
                stop
            endif
        endif
        allowed_args = args()
        self%argcnt  = command_argument_count()
        do i=1,self%argcnt
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(logfhandle,*) 'ERROR! while parsing the command line: simple_cmdline :: parse_oldschool'
                write(logfhandle,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(logfhandle,*) 'which likely exceeds the length limit LONGSTRLEN'
                write(logfhandle,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            call self%parse_command_line_value(i, arg, allowed_args)
        end do
        if( present(keys_required) ) call self%check
    end subroutine parse_oldschool

    subroutine parse_command_line_value( self, i, arg, allowed_args )
        class(cmdline),   intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: arg
        class(args),      intent(in)    :: allowed_args
        character(len=:), allocatable   :: form
        real(kind=dp) :: rval
        integer       :: pos1, io_stat, ival
        pos1 = index(arg, '=') ! position of '='
        ! parse everyting containing '='
        if( pos1 /= 0 )then
            self%cmds(i)%key = arg(:pos1-1) ! KEY
            if( .not. allowed_args%is_present(self%cmds(i)%key) )then
                write(logfhandle,'(a,a)') trim(self%cmds(i)%key), ' argument is not allowed'
                write(logfhandle,'(a)') 'Perhaps you have misspelled?'
                stop
            endif
            call str2format(arg(pos1+1:), form, rval, ival)
            select case(form)
                case('file', 'dir', 'char')
                    self%cmds(i)%carg = adjustl(arg(pos1+1:))
                case('real')
                    self%cmds(i)%rarg = rval
                case('int')
                    ival = str2int(adjustl(arg(pos1+1:)), io_stat)
                    if( io_stat == 0 )then
                        self%cmds(i)%rarg = real(ival,kind=dp)
                    else
                        self%cmds(i)%carg = adjustl(arg(pos1+1:))
                    endif
            end select
            self%cmds(i)%defined = .true.
        endif
    end subroutine parse_command_line_value

    elemental subroutine copy( self, self2copy )
        class(cmdline), intent(inout) :: self
        class(cmdline), intent(in) :: self2copy
        integer :: icmd
        ! copy cmds
        do icmd=1,MAX_CMDARGS
            if( allocated(self2copy%cmds(icmd)%key) )then
                self%cmds(icmd)%key = self2copy%cmds(icmd)%key
            else
                if( allocated(self%cmds(icmd)%key) ) deallocate(self%cmds(icmd)%key)
            endif
            if( allocated(self2copy%cmds(icmd)%carg) )then
                self%cmds(icmd)%carg = self2copy%cmds(icmd)%carg
            else
                if( allocated(self%cmds(icmd)%carg) ) deallocate(self%cmds(icmd)%carg)
            endif
        end do
        self%cmds(:)%defined = self2copy%cmds(:)%defined
        self%cmds(:)%rarg    = self2copy%cmds(:)%rarg
        ! copy rest
        self%checker(:)  = self2copy%checker(:)
        self%entire_line = trim(self2copy%entire_line)
        self%argcnt      = self2copy%argcnt
        self%ncheck      = self2copy%ncheck
    end subroutine copy

    !>  \brief  polymorphic assignment (=)
    elemental subroutine assign( self, self2copy )
        class(cmdline), intent(inout) :: self
        class(cmdline), intent(in) :: self2copy
        call self%copy(self2copy)
    end subroutine assign

    !> \brief  for setting a real(kind=8) valued command line argument
    subroutine set_1( self, key, rarg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(kind=dp),    intent(in)    :: rarg
        integer :: which
        which = self%lookup(key)
        if( which == 0 )then
            self%argcnt = self%argcnt + 1
            if( self%argcnt > MAX_CMDARGS )then
                write(logfhandle,*) 'self%argcnt: ', self%argcnt
                write(logfhandle,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                THROW_HARD('stack overflow; set_1')
            endif
            self%cmds(self%argcnt)%key     = trim(key)
            self%cmds(self%argcnt)%rarg    = rarg
            self%cmds(self%argcnt)%defined = .true.
        else
            self%cmds(which)%rarg = rarg
        endif
    end subroutine set_1

    !> \brief  for setting a real valued command line argument
    subroutine set_2( self, key, rarg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: rarg
        call self%set_1(key, real(rarg,dp))
    end subroutine set_2

    !> \brief  for setting a char valued command line argument
    subroutine set_3( self, key, carg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: carg
        integer :: which
        which = self%lookup(key)
        if( which == 0 )then
            self%argcnt = self%argcnt + 1
            if( self%argcnt > MAX_CMDARGS )then
                write(logfhandle,*) 'self%argcnt: ', self%argcnt
                write(logfhandle,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                THROW_HARD('stack overflow; set_2')
            endif
            self%cmds(self%argcnt)%key = trim(key)
            self%cmds(self%argcnt)%carg = trim(carg)
            self%cmds(self%argcnt)%defined = .true.
        else
            self%cmds(which)%carg = trim(carg)
        endif
    end subroutine set_3

    !> \brief  for setting an integer valued command line argument as a real
    subroutine set_4( self, key, iarg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: iarg
        call self%set_1(key, real(iarg,dp))
    end subroutine set_4

    !> \brief for removing a command line argument
    subroutine delete( self, key )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: which
        which = self%lookup(key)
        if( which > 0 )then
            self%cmds( which:self%argcnt-1 ) = self%cmds( which+1:self%argcnt )
            if( allocated(self%cmds(self%argcnt)%key)  ) deallocate(self%cmds(self%argcnt)%key)
            if( allocated(self%cmds(self%argcnt)%carg) ) deallocate(self%cmds(self%argcnt)%carg)
            self%cmds(self%argcnt)%rarg = 0.d0
            self%cmds(self%argcnt)%defined = .true.
            self%argcnt = self%argcnt - 1
        endif
    end subroutine delete

    !> \brief  to set variables to be checked
    subroutine checkvar( self, str, nr )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: str
        integer,          intent(in)    :: nr
        self%checker(nr) = str
        self%ncheck = nr
    end subroutine checkvar

    !> \brief  for checking that the passed array of keys are present in the struct
    subroutine check( self )
        class(cmdline), intent(inout) :: self
        logical, allocatable  :: cmderr(:)
        integer               :: i, nstates
        character(len=STDLEN) :: str
        logical               :: vol_defined
        allocate( cmderr(self%ncheck) )
        cmderr = .false.
        do i=1,self%ncheck
           cmderr(i) = .not. self%defined(self%checker(i))
           if( cmderr(i) ) write(logfhandle,'(a,a)') 'Missing key on command line: ', trim(self%checker(i))
        end do
        ! to take care of the case where the first state has vanished (eg vol1 key absent)
        if( self%defined('nstates') )then
            nstates = nint( self%get_rarg('nstates') )
            if( nstates>1 .and. .not.self%defined('vol1') )then
                ! if 'vol1' absent . . .
                vol_defined = .false.
                do i=1,nstates
                    str = 'vol'//int2str(i)
                    if( self%lookup(str) > 0) then
                        ! and another volX key is defined
                        vol_defined = .true.
                    endif
                end do
                ! ...then there is no error
                if( vol_defined )then
                    do i=1,self%ncheck
                        str = self%checker(i)
                        if( self%checker(i) .eq. 'vol1' )cmderr(i) = .false.
                    end do
                endif
            endif
        endif
        ! output
        if( any(cmderr) )then
            write(logfhandle,'(a)') 'ERROR, not enough input variables defined!'
            stop
        endif
        deallocate( cmderr )
    end subroutine check

    !> \brief  for writing the command line
    subroutine printline( self, tag )
        class(cmdline),             intent(inout) :: self
        character(len=*), optional, intent(in)    :: tag
        integer :: i
        if( present(tag) )then
            do i=1,self%argcnt
                if( self%cmds(i)%defined .and. allocated(self%cmds(i)%carg) )then
                    write(logfhandle,*) trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg), ' ', trim(tag)
                else if( self%cmds(i)%defined )then
                    write(logfhandle,*) trim(self%cmds(i)%key), ' ', self%cmds(i)%rarg, ' ',trim(tag)
                endif
            end do
        else
            do i=1,self%argcnt
                if( self%cmds(i)%defined .and. allocated(self%cmds(i)%carg) )then
                    write(logfhandle,*) trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg)
                else if( self%cmds(i)%defined )then
                    write(logfhandle,*) trim(self%cmds(i)%key), ' ', self%cmds(i)%rarg
                endif
            end do
        endif
    end subroutine printline

    !> \brief  for reading the command line from file
    subroutine readline( self, filename )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: filename
        character(len=:), allocatable   :: key, arg
        character(len=XLONGSTRLEN)      :: line, tmp
        real(kind=dp) :: rarg
        integer       :: i, funit, ok, nl
        call self%kill
        if( .not.file_exists(adjustl(filename))) THROW_HARD('Could not find file: '//trim(filename))
        nl = nlines(filename)
        if( nl > MAX_CMDARGS ) THROW_HARD('Number of entries > maximum!')
        call fopen(funit, file=trim(adjustl(filename)), status='old', iostat=ok)
        call fileiochk('readline fopen failed '//trim(filename),ok)
        do i = 1,nl
            read(funit,'(A)',iostat=ok) line
            call fileiochk('readline reading failed'//trim(filename),iostat=ok)
            line = adjustl(trim(line))
            call split_str(line, ' ', tmp)
            key = adjustl(trim(tmp))
            arg = adjustl(trim(line))
            read(arg,*,iostat=ok) rarg
            if( ok == 0 )then
                call self%set(key, rarg)
            else
                call self%set(key, arg)
            endif
        enddo
        call fclose(funit)
    end subroutine readline

    !> \brief  for writing the command line to file
    subroutine writeline( self, filename, tag )
        class(cmdline),             intent(inout) :: self
        character(len=*),           intent(in)    :: filename
        character(len=*), optional, intent(in)    :: tag
        integer :: i, parameters_fhandle, ok
        if( file_exists(trim(adjustl(filename)))) call del_file(trim(adjustl(filename)))
        call fopen(parameters_fhandle, file=trim(adjustl(filename)), status='new', iostat=ok)
        if( present(tag) )then
            do i=1,self%argcnt
                if( self%cmds(i)%defined .and. allocated(self%cmds(i)%carg) )then
                    write(parameters_fhandle, '(RD,A,A,A,A,A)') trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg), ' ', trim(tag)
                else if( self%cmds(i)%defined )then
                    write(parameters_fhandle, '(RD,A,F14.4,A)') trim(self%cmds(i)%key), self%cmds(i)%rarg, trim(tag)
                endif
            end do
        else
            do i=1,self%argcnt
                if( self%cmds(i)%defined .and. allocated(self%cmds(i)%carg) )then
                    write(parameters_fhandle, '(RD,A,A,A)') trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg)
                else if( self%cmds(i)%defined )then
                    write(parameters_fhandle, '(RD,A,F14.4,A)') trim(self%cmds(i)%key), self%cmds(i)%rarg
                endif
            end do
        endif
        call fclose(parameters_fhandle)
    end subroutine writeline
    
    !> \brief  for checking the existence of of arg
    function defined( self, key ) result( def )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: def
        integer :: i
        def = .false.
        do i=1,self%argcnt
            if( trim(self%cmds(i)%key) .eq. trim(key) )then
                def = .true.
                return
            endif
        end do
    end function defined

    !> \brief for looking up a key
    function lookup( self, key ) result( which )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: i, which
        which = 0
        do i=1,self%argcnt
            if( allocated(self%cmds(i)%key) )then
                if( trim(self%cmds(i)%key) .eq. trim(key) )then
                    which = i
                    return
                endif
            endif
        end do
    end function lookup

    !> \brief for getting the keys
    function get_keys( self ) result( keys )
        class(cmdline),   intent(in) :: self
        type(str4arr), allocatable :: keys(:)
        integer :: i
        if( self%argcnt > 0 )then
            allocate(keys(self%argcnt))
            do i=1,self%argcnt
                if( allocated(self%cmds(i)%key) ) allocate(keys(i)%str, source=self%cmds(i)%key)
            end do
        endif
    end function get_keys

    !> \brief for getting real args
    elemental function get_rarg( self, key ) result( rval )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        real :: rval
        integer :: i
        rval = 0.
        do i=1,self%argcnt
            if( allocated(self%cmds(i)%key) )then
                if( trim(self%cmds(i)%key) .eq. trim(key) )then
                    rval = real(self%cmds(i)%rarg)
                    return
                endif
            endif
        end do
    end function get_rarg

    !> \brief for getting real args
    elemental integer function get_iarg( self, key )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        integer :: i
        get_iarg = 0
        do i=1,self%argcnt
            if( allocated(self%cmds(i)%key) )then
                if( trim(self%cmds(i)%key) .eq. trim(key) )then
                    get_iarg = nint(self%cmds(i)%rarg)
                    return
                endif
            endif
        end do
    end function get_iarg

    !> \brief for getting char args
    pure function get_carg( self, key ) result( cval )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: cval
        integer :: i
        do i=1,self%argcnt
            if( allocated(self%cmds(i)%key) )then
                if( trim(self%cmds(i)%key) .eq. trim(key) )then
                    cval = trim(self%cmds(i)%carg)
                    return
                endif
            endif
        end do
    end function get_carg

    !> \brief for generating a job description (prg overrides prg in cline)
    subroutine gen_job_descr( self, hash, prg )
        class(cmdline),             intent(in)    :: self
        class(chash),               intent(inout) :: hash
        character(len=*), optional, intent(in)    :: prg
        integer :: i
        call hash%new(MAX_CMDARGS)
        do i=1,self%argcnt
            if( .not. allocated(self%cmds(i)%carg) )then
                ! value is real
                call hash%push(trim(self%cmds(i)%key), real2str(self%cmds(i)%rarg))
            else
                ! value is char
                call hash%push(trim(self%cmds(i)%key), trim(self%cmds(i)%carg))
            endif
        end do
        if( present(prg) ) call hash%set('prg', trim(prg))
    end subroutine gen_job_descr

    elemental subroutine kill( self )
        class(cmdline), intent(inout) :: self
        integer :: icmd
        do icmd=1,MAX_CMDARGS
            if( allocated(self%cmds(icmd)%key) ) deallocate(self%cmds(icmd)%key)
            if( allocated(self%cmds(icmd)%carg) ) deallocate(self%cmds(icmd)%carg)
        end do
        self%argcnt = 0
        self%ncheck = 0
    end subroutine kill

    ! EXCEPTION-HANDLING

    !> \brief  is for raising command line exception
    subroutine cmdline_err( cmdstat, cmdlen, arg, pos )
        integer, intent(in)          :: cmdstat, cmdlen, pos
        character(len=*), intent(in) :: arg
        if( cmdstat == -1 )then
            write(logfhandle,*) 'ERROR! while parsing the command line'
            write(logfhandle,*) 'The string length of argument: ', arg, 'is: ', cmdlen
            write(logfhandle,*) 'which likely exeeds the lenght limit LONGSTRLEN'
            write(logfhandle,*) 'Create a symbolic link with shorter name in the cwd'
            stop
        endif
        if( arg(:pos-1) .ne. 'prg' )then
            write(logfhandle,'(a)') 'ERROR!'
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'prg=program required to be first on command line'
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'Please, execute with prg=list to print all available programs'
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'Please, execute with prg flag and describe=yes to obtain a short description'
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'Executing with prg flag only prints all available command-line options'
            stop
        endif
    end subroutine cmdline_err

end module simple_cmdline
