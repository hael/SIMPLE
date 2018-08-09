! command line parser
module simple_cmdline
include 'simple_lib.f08'
use simple_cmd_dict ! use all in there
implicit none

public :: cmdarg, cmdline, cmdline_err
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_CMDARGS = 60

!> cmdarg key/value pair command-line type
type cmdarg
    character(len=:), allocatable :: key, carg
    logical                       :: defined=.false.
    real                          :: rarg=0.
contains
    procedure :: new
end type cmdarg
interface cmdarg
    module procedure constructor
end interface

type cmdline
    type(cmdarg)          :: cmds(MAX_CMDARGS)
    character(len=KEYLEN) :: checker(MAX_CMDARGS)
    character(len=8192)   :: entire_line
    integer               :: argcnt=0, ncheck=0
contains
    procedure          :: parse
    procedure          :: parse_oldschool
    procedure          :: parse_command_line_value
    procedure, private :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: set_1
    procedure          :: set_2
    procedure          :: set_3
    procedure          :: set_4
    procedure, private :: lookup
    procedure          :: get_keys
    generic            :: set => set_1, set_2,set_3, set_4
    procedure          :: delete
    procedure          :: checkvar
    procedure          :: check
    procedure          :: printline
    procedure          :: defined
    procedure          :: get_rarg
    procedure          :: get_carg
    procedure          :: gen_job_descr
end type cmdline

contains
    function constructor( key, carg, rarg ) result( self )
        character(len=*),  intent(in) :: key
        character(len=*),  intent(in), optional :: carg
        real,              intent(in), optional :: rarg
        type(cmdarg) :: self
        call self%new( key, carg, rarg )
    end function constructor

    subroutine new( self, key, carg, rarg )
        class(cmdarg) :: self
        character(len=*),  intent(in) :: key
        character(len=*),  intent(in), optional :: carg
        real,              intent(in), optional :: rarg
        allocate(self%key, source=key)
        if(present(carg)) allocate(self%carg, source=carg)
        if(present(rarg)) self%rarg = rarg
    end subroutine new

    !> \brief for parsing the command line arguments passed as key=val
    subroutine parse( self )
        use simple_args, only: args
        use simple_user_interface ! use all in there
        class(cmdline), intent(inout) :: self
        type(str4arr), allocatable    :: keys_required(:)
        type(args)                    :: allowed_args
        type(simple_program), pointer :: ptr2prg => null()
        character(len=LONGSTRLEN)     :: arg
        integer :: i, cmdstat, cmdlen, ikey, pos, nargs_required, sz_keys_req
        ! parse command line
        self%argcnt  = command_argument_count()
        call get_command(self%entire_line)
        cmdline_glob = trim(self%entire_line)
        ! parse program name
        call get_command_argument(1, arg, cmdlen, cmdstat)
        pos = index(arg, '=') ! position of '='
        call cmdline_err(cmdstat, cmdlen, arg, pos)
        if( str_has_substr(arg(pos+1:), 'simple_') ) stop 'giving program names with simple_* prefix is depreciated'
        ! obtain pointer to the program in the simple_user_interface specification
        call get_prg_ptr(arg(pos+1:), ptr2prg)
        ! list programs if so instructed
        if( str_has_substr(self%entire_line, 'prg=list') )then
            if( ptr2prg%is_distr() )then
                call list_distr_prgs_in_ui
                stop
            else
                call list_shmem_prgs_in_ui
                stop
            endif
        endif
        ! describe program if so instructed
        if( str_has_substr(self%entire_line, 'describe=yes') )then
            call ptr2prg%print_prg_descr_long()
            stop
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
                write(*,*) 'ERROR! while parsing the command line: simple_cmdline :: parse'
                write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(*,*) 'which likely exceeds the length limit LONGSTRLEN'
                write(*,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            call self%parse_command_line_value(i, arg, allowed_args)
        end do
        if( sz_keys_req > 0 ) call self%check
    end subroutine parse

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
                call print_cmdline(keys_required, keys_optional, distr=distr_exec)
                stop
            else
                ! indicate which variables are required
                do ikey=1,size(keys_required)
                    call self%checkvar(keys_required(ikey), ikey)
                end do
            endif
        else
            if( cmdargcnt < 2 )then
                call print_cmdline(keys_required, keys_optional, distr=distr_exec)
                stop
            endif
        endif
        allowed_args = args()
        self%argcnt  = command_argument_count()
        do i=1,self%argcnt
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(*,*) 'ERROR! while parsing the command line: simple_cmdline :: parse_oldschool'
                write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(*,*) 'which likely exceeds the length limit LONGSTRLEN'
                write(*,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            call self%parse_command_line_value(i, arg, allowed_args)
        end do
        if( present(keys_required) ) call self%check
    end subroutine parse_oldschool

    subroutine parse_command_line_value( self, i, arg, allowed_args )
        use simple_args, only: args
        class(cmdline),   intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: arg
        class(args),      intent(in)    :: allowed_args
        character(len=:), allocatable   :: form
        integer :: pos1, io_stat, ri
        pos1 = index(arg, '=') ! position of '='
        ! parse everyting containing '='
        if( pos1 /= 0 )then
            self%cmds(i)%key = arg(:pos1-1) ! KEY
            if( .not. allowed_args%is_present(self%cmds(i)%key) )then
                write(*,'(a,a)') trim(self%cmds(i)%key), ' argument is not allowed'
                write(*,'(a)') 'Perhaps you have misspelled?'
                stop
            endif
            form = str2format(arg(pos1+1:))
            select case(form)
                case('file', 'dir', 'char')
                    self%cmds(i)%carg = adjustl(arg(pos1+1:))
                case('real')
                    self%cmds(i)%rarg = str2real(adjustl(arg(pos1+1:)))
                case('int')
                    call str2int(adjustl(arg(pos1+1:)), io_stat, ri )
                    if( io_stat == 0 )then
                        self%cmds(i)%rarg = real(ri)
                    else
                        self%cmds(i)%carg = adjustl(arg(pos1+1:))
                    endif
            end select
            self%cmds(i)%defined = .true.
        endif
    end subroutine parse_command_line_value

    !>  \brief  private copier
    subroutine copy( self, self2copy )
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
    subroutine assign( self, self2copy )
        class(cmdline), intent(inout) :: self
        class(cmdline), intent(in) :: self2copy
        call self%copy(self2copy)
    end subroutine assign

    !> \brief  for setting a real valued command line argument
    subroutine set_1( self, key, rarg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: rarg
        integer :: which
        which = self%lookup(key)
        if( which == 0 )then
            self%argcnt = self%argcnt + 1
            if( self%argcnt > MAX_CMDARGS )then
                write(*,*) 'self%argcnt: ', self%argcnt
                write(*,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                stop 'ERROR! stack overflow; simple_cmdline :: set_1'
            endif
            self%cmds(self%argcnt)%key = trim(key)
            self%cmds(self%argcnt)%rarg = rarg
            self%cmds(self%argcnt)%defined = .true.
        else
            self%cmds(which)%rarg = rarg
        endif
    end subroutine set_1

    !> \brief  for setting a char valued command line argument
    subroutine set_2( self, key, carg )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: carg
        integer :: which
        which = self%lookup(key)
        if( which == 0 )then
            self%argcnt = self%argcnt + 1
            if( self%argcnt > MAX_CMDARGS )then
                write(*,*) 'self%argcnt: ', self%argcnt
                write(*,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                stop 'ERROR! stack overflow; simple_cmdline :: set_2'
            endif
            self%cmds(self%argcnt)%key = trim(key)
            self%cmds(self%argcnt)%carg = trim(carg)
            self%cmds(self%argcnt)%defined = .true.
        else
            self%cmds(which)%carg = trim(carg)
        endif
    end subroutine set_2

   !> \brief  for setting a command line argument with cmdarg object
    subroutine set_3( self, cmarg )
        class(cmdline),   intent(inout) :: self
        type(cmdarg), intent(in)    :: cmarg
        integer :: which
     !   do n=1,size(cmarg)
            which = self%lookup(cmarg%key)
            if( which == 0 )then
                self%argcnt = self%argcnt + 1
                if( self%argcnt > MAX_CMDARGS )then
                    write(*,*) 'self%argcnt: ', self%argcnt
                    write(*,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                    stop 'ERROR! stack overflow; simple_cmdline :: set_3'
                endif
                self%cmds(self%argcnt)%key = trim(cmarg%key)
                if(allocated(cmarg%carg))then
                    self%cmds(self%argcnt)%carg = trim(cmarg%carg)
                else
                    self%cmds(self%argcnt)%rarg = cmarg%rarg
                endif
                self%cmds(self%argcnt)%defined = .true.
            else
                if(allocated(cmarg%carg))then
                    self%cmds(which)%carg = trim(cmarg%carg)
                else
                    self%cmds(self%argcnt)%rarg = cmarg%rarg
                endif
            endif
      !  end do
    end subroutine set_3

    !> \brief  for setting a command line argument with a list of cmdarg objects
    subroutine set_4( self, cmarg )
        class(cmdline),   intent(inout) :: self
        type(cmdarg), intent(in)    :: cmarg(:)
        integer :: which, n
        do n=1,size(cmarg)
            which = self%lookup(cmarg(n)%key)
            if( which == 0 )then
                self%argcnt = self%argcnt + 1
                if( self%argcnt > MAX_CMDARGS )then
                    write(*,*) 'self%argcnt: ', self%argcnt
                    write(*,*) 'MAX_CMDARGS: ', MAX_CMDARGS
                    stop 'ERROR! stack overflow; simple_cmdline :: set_3'
                endif
                self%cmds(self%argcnt)%key = trim(cmarg(n)%key)
                if(allocated(cmarg(n)%carg))then
                    self%cmds(self%argcnt)%carg = trim(cmarg(n)%carg)
                else
                    self%cmds(self%argcnt)%rarg = cmarg(n)%rarg
                endif
                self%cmds(self%argcnt)%defined = .true.
            else
                if(allocated(cmarg(n)%carg))then
                    self%cmds(which)%carg = trim(cmarg(n)%carg)
                else
                    self%cmds(self%argcnt)%rarg = cmarg(n)%rarg
                endif
            endif
        end do
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
            self%cmds(self%argcnt)%rarg = 0.
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
        allocate( cmderr(self%ncheck), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('check; simple_cmdline',alloc_stat)
        cmderr = .false.
        do i=1,self%ncheck
           cmderr(i) = .not. self%defined(self%checker(i))
           if( cmderr(i) ) write(*,'(a,a)') 'Missing key on command line: ', trim(self%checker(i))
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
            write(*,'(a)') 'ERROR, not enough input variables defined!'
            stop
        endif
        deallocate( cmderr )
    end subroutine check

    !> \brief  for writing the command line
    subroutine printline( self )
        class(cmdline), intent(inout) :: self
        integer :: i
        do i=1,self%argcnt
            if( self%cmds(i)%defined .and. allocated(self%cmds(i)%carg) )then
                write(*,*) trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg)
            else if( self%cmds(i)%defined )then
                write(*,*) trim(self%cmds(i)%key), ' ', self%cmds(i)%rarg
            endif
        end do
    end subroutine printline

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
    pure function get_rarg( self, key ) result( rval )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        real :: rval
        integer :: i
        rval = 0.
        do i=1,self%argcnt
            if( allocated(self%cmds(i)%key) )then
                if( trim(self%cmds(i)%key) .eq. trim(key) )then
                    rval = self%cmds(i)%rarg
                    return
                endif
            endif
        end do
    end function get_rarg

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

    ! EXCEPTION-HANDLING

    !> \brief  is for raising command line exception
    subroutine cmdline_err( cmdstat, cmdlen, arg, pos )
        integer, intent(in)          :: cmdstat, cmdlen, pos
        character(len=*), intent(in) :: arg
        if( cmdstat == -1 )then
            write(*,*) 'ERROR! while parsing the command line'
            write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
            write(*,*) 'which likely exeeds the lenght limit LONGSTRLEN'
            write(*,*) 'Create a symbolic link with shorter name in the cwd'
            stop
        endif
        if( arg(:pos-1) .ne. 'prg' )then
            write(*,'(a)') 'ERROR!'
            write(*,'(a)') 'prg=simple_program required to be first on command line'
            write(*,'(a)') 'Please, refer to the manual for a comprehensive '
            write(*,'(a)') 'list of all programs and their specific documentation'
            stop
        endif
    end subroutine cmdline_err

end module simple_cmdline
