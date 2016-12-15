!>  \brief  SIMPLE command line parsing singleton
!! The code is distributed with the hope that it will be useful, but WITHOUT ANY WARRANTY. 
!! Redistribution or modification is regulated by the GNU General Public License.
!! Hans Elmlund, 2011-08-17.
module simple_cmdline
use simple_chash, only: chash
use simple_defs     ! singleton
use simple_cmd_dict ! singleton
use simple_yaml_output
use simple_yaml_strings
use simple_file_utils
use simple_file_defs
use simple_err_defs
use simple_eglossary
use simple_eglossary_lowlev
use simple_dynamic_memory

implicit none

public :: cmdline
private

logical, parameter :: debug=.false.

type cmdarg
    character(len=32)     :: key=''
    character(len=STDLEN) :: carg=''
    logical               :: defined=.false.
    real                  :: rarg=0.
end type cmdarg

type cmdline
    integer               :: NMAX=60
    type(cmdarg)          :: cmds(60)
    character(len=32)     :: checker(60)
    character(len=STDLEN) :: entire_line=''
    integer               :: argcnt=0, totlen=0, ncheck=0
  contains
    procedure             :: parse
    procedure, private    :: copy
    procedure, private    :: assign
    generic               :: assignment(=) => assign
    procedure, private    :: set_1
    procedure, private    :: set_2
    procedure, private    :: lookup
    generic               :: set => set_1, set_2
    procedure             :: delete
    procedure             :: checkvar
    procedure             :: check
    procedure, private    :: write
    procedure             :: print
    procedure             :: defined
    procedure             :: get_rarg
    procedure             :: get_carg
    procedure             :: gen_job_descr
end type cmdline

contains

    !> \brief  for parsing the command line arguments passed as key=val 
    subroutine parse( self, keys_required, keys_optional )
        use simple_args, only: args
        use simple_jiffys
        class(cmdline),             intent(inout) :: self
        character(len=*), optional, intent(in)    :: keys_required(:), keys_optional(:)
        character(len=STDLEN) :: arg
        type(args)            :: allowed_args
        integer               :: i, ri, ierr, cmdstat, cmdlen, cntbin, ikey
        integer               :: cnttxt, n, fnr, file_stat, io_stat, nreq, cmdargcnt
        logical               :: here
        cmdargcnt = command_argument_count()
        call get_command(self%entire_line)
        if( debug ) print *, 'DEBUG(simple_cmdline, L61), command_argument_count: ', cmdargcnt
        if( present(keys_required) )then
            if( str_has_substr(self%entire_line,'prg=') )then
                nreq = size(keys_required)+1 ! +1 because prg part of command line
            else
                nreq = size(keys_required)
            endif
            if( debug )then
                print *, 'DEBUG(simple_cmdline, L69), no keys required:       ', nreq
                print *, 'DEBUG(simple_cmdline, L70), command_argument_count: ', cmdargcnt
            endif
            if( cmdargcnt < nreq )then
                call print_cmdline(keys_required, keys_optional)
                stop
            else
                ! indicate which variables are required
                do ikey=1,size(keys_required)
                    call self%checkvar(keys_required(ikey), ikey)
                end do
            endif
        else
            if( cmdargcnt < 2 )then
                call print_cmdline(keys_required, keys_optional)
                stop
            endif
        endif    
        allowed_args = args()
        self%totlen  = 0
        self%argcnt  = command_argument_count()
        cntbin = 0
        cnttxt = 0
        do i=1,self%argcnt 
            call get_command_argument(i, arg, cmdlen, cmdstat)
            if( cmdstat == -1 )then
                write(*,*) 'ERROR! while parsing the command line; simple_cmdline :: parse'
                write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
                write(*,*) 'which likely exeeds the length limit STDLEN'
                write(*,*) 'Create a symbolic link with shorter name in the cwd'
                stop
            endif
            self%totlen = self%totlen+cmdlen+1
            call parse_local(arg)
        end do
        if( present(keys_required) )then
            if( debug ) print *, 'DEBUG(simple_cmdline, L93), reached checker'
            call self%check
        endif
        call self%write
        
        contains
        
            subroutine parse_local( arg ) 
                character(len=*) :: arg
                integer          :: pos1
                pos1 = index(arg, '=') ! position of '='
                ! parse everyting containing '='
                if( pos1 /= 0 )then
                    self%cmds(i)%key = arg(:pos1-1) ! KEY
                    if( .not. allowed_args%is_present(self%cmds(i)%key) )then
                        write(*,'(a,a)') trim(self%cmds(i)%key), ' argument is not allowed'
                        write(*,'(a)') 'Perhaps you have misspelled?'
                        stop
                    endif
                    if( index(arg(pos1+1:), '.spi') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.mrc') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.map') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.ctf') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.bin') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.txt') /= 0 )then
                        cnttxt = cnttxt+1
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.asc') /= 0 )then
                        cnttxt = cnttxt+1
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.log') /= 0 )then
                         stop '.log files are obsolete!'
                    else if( index(arg(pos1+1:), '.dat') /= 0 )then
                         stop '.dat files are obsolete! Use .txt files exclusively!'
                    else if( index(arg(pos1+1:), '.box') /= 0 )then
                        self%cmds(i)%carg =arg(pos1+1:)
                    else if( index(arg(pos1+1:), '.') /= 0 )then
                        call str2real(arg(pos1+1:), io_stat, self%cmds(i)%rarg)
                        if( io_stat .ne. 0 )then
                            stop 'ERROR(I/O), str2real; simple_cmdline :: parse_local'
                        endif
                    else
                        call str2int(arg(pos1+1:), io_stat, ri )
                        if( io_stat==0 )then 
                            self%cmds(i)%rarg = real(ri)
                        else
                            self%cmds(i)%carg = arg(pos1+1:)
                        endif
                    endif
                    self%cmds(i)%defined = .true.
                    if( debug )then
                        write(*,*) 'DEBUG(simple_cmdline, L161), key/cval/rval/defined: ',&
                        trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg), self%cmds(i)%rarg, self%cmds(i)%defined
                    endif
                endif
            end subroutine
            
    end subroutine parse
    
    !>  \brief  private copier
    subroutine copy( self, self2copy )
        class(cmdline), intent(inout) :: self
        class(cmdline), intent(in)    :: self2copy
        integer :: i
        ! copy cmds
        self%cmds(:)%key     = self2copy%cmds(:)%key
        self%cmds(:)%carg    = self2copy%cmds(:)%carg
        self%cmds(:)%defined = self2copy%cmds(:)%defined
        self%cmds(:)%rarg    = self2copy%cmds(:)%rarg
        ! copy rest
        self%checker(:)      = self2copy%checker(:)
        self%entire_line     = self2copy%entire_line
        self%argcnt          = self2copy%argcnt
        self%totlen          = self2copy%totlen
        self%ncheck          = self2copy%ncheck
    end subroutine copy
    
    !>  \brief  polymorphic assignment (=)
    subroutine assign( self, self2copy )
        class(cmdline), intent(inout) :: self
        class(cmdline), intent(in)    :: self2copy
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
            self%argcnt                    = self%argcnt + 1
            self%cmds(self%argcnt)%key     = key
            self%cmds(self%argcnt)%rarg    = rarg
            self%cmds(self%argcnt)%defined = .true. 
        else
            self%cmds(which)%rarg    = rarg
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
            self%argcnt                    = self%argcnt + 1
            self%cmds(self%argcnt)%key     = key
            self%cmds(self%argcnt)%carg    = carg
            self%cmds(self%argcnt)%defined = .true. 
        else
            self%cmds(which)%carg = carg
        endif
    end subroutine set_2

    !> \brief  for removing a command line argument
    subroutine delete( self, key )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: which
        which = self%lookup(key)
        if( which > 0 )then
            self%cmds( which:self%argcnt-1 ) = self%cmds( which+1:self%argcnt )
            self%cmds( self%argcnt )%key     = ''
            self%cmds( self%argcnt )%carg    = ''
            self%cmds( self%argcnt )%rarg    = 0.
            self%cmds( self%argcnt )%defined = .true. 
            self%argcnt                      = self%argcnt - 1
        endif
    end subroutine delete

    
    !> \brief  to set variables to be checked
    subroutine checkvar( self, str, nr )
        class(cmdline),   intent(inout) :: self
        character(len=*), intent(in)    :: str
        integer,          intent(in)    :: nr
        self%checker(nr) = str
        self%ncheck      = nr
    end subroutine checkvar
    
    !> \brief  for checking that the passed array of keys are present in the struct
    subroutine check( self )
        use simple_jiffys, only: alloc_err
        class(cmdline), intent(inout) :: self
        logical, allocatable :: cmderr(:)
        integer :: i, alloc_stat 
        allocate( cmderr(self%ncheck), stat=alloc_stat )
        call alloc_err('check; simple_cmdline', alloc_stat)
        cmderr = .false.
        do i=1,self%ncheck
            cmderr(i) = .not. self%defined(self%checker(i))
        end do
        if( any(cmderr) )then
            write(*,'(a)') 'ERROR, not enough input variables defined!'
            stop
        endif
        deallocate( cmderr )
    end subroutine check
    
    !> \brief  for writing the command line
    subroutine write( self )
        use simple_jiffys, only: get_fileunit
        class(cmdline), intent(inout) :: self
        integer           :: funit, file_stat, i
        integer           :: yunit
        integer           :: err ! error code
        !filename string
        character(len=3)  :: char_out
        character(len=80) :: tmr_name
        !function calls
        integer           :: convert_int2char_indexed_c
        integer           :: strlen
        yunit = 1
        err = convert_int2char_indexed_c(char_out,yunit,1,1)
        tmr_name = 'cmdline'
        tmr_name = tmr_name(1:strlen(tmr_name))//char_out(1:strlen(char_out))
        tmr_name = tmr_name(1:strlen(tmr_name))//".yaml"
        call file_open(tmr_name,yunit,'unknown','asis','readwrite')
        write(funit,*) trim(self%entire_line)
        call yaml_comment('The comand line details')
        call yaml_map('full command line',trim(self%entire_line))
        do i=1,self%argcnt
            if( self%cmds(i)%defined .and. self%cmds(i)%carg .eq. '' )then
                call yaml_map(trim(self%cmds(i)%key),self%cmds(i)%rarg)
             else if( self%cmds(i)%defined)then
                call yaml_map(trim(self%cmds(i)%key),trim(self%cmds(i)%carg))
            endif
        end do
        close(yunit)
        return
    end subroutine write

    !> \brief  for writing the command line
    subroutine print( self )
        class(cmdline), intent(inout) :: self
        integer :: i
        do i=1,self%argcnt
            if( self%cmds(i)%defined .and. self%cmds(i)%carg .eq. '' )then
                write(*,*) trim(self%cmds(i)%key), ' ', self%cmds(i)%rarg
            else if( self%cmds(i)%defined)then
                write(*,*) trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg)
            endif
        end do
    end subroutine print

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

    !> \brief  for looking up a key
    function lookup( self, key ) result( which )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in)  :: key
        integer :: i, which
        which = 0
        do i=1,self%argcnt
            if( trim(self%cmds(i)%key) .eq. trim(key) )then
                which = i
                return
            endif
        end do
    end function lookup
    
    !> \brief  for getting real args
    pure function get_rarg( self, key ) result( rval )
        class(cmdline),   intent(in) :: self
        character(len=*), intent(in) :: key
        real :: rval
        integer :: i
        rval = 0.
        do i=1,self%argcnt
            if( trim(self%cmds(i)%key) .eq. trim(key) )then
                rval = self%cmds(i)%rarg
                return
            endif
        end do
    end function get_rarg
    
    !> \brief  for getting char args
    pure function get_carg( self, key ) result( cval )
        class(cmdline),   intent(in)  :: self
        character(len=*), intent(in)  :: key
        character(len=STDLEN) :: cval
        integer :: i
        cval = ''
        do i=1,self%argcnt
            if( trim(self%cmds(i)%key) .eq. trim(key) )then
                cval = self%cmds(i)%carg
                return
            endif
        end do
    end function get_carg

    !> \brief  for generating a job description (prg overrides prg in cline)
    subroutine gen_job_descr( self, hash, prg )
        use simple_jiffys, only: real2str
        class(cmdline),             intent(in)    :: self
        class(chash),               intent(inout) :: hash
        character(len=*), optional, intent(in)    :: prg
        character(len=STDLEN) :: str
        integer               :: i
        call hash%new(self%NMAX)
        do i=1,self%argcnt
            if( self%cmds(i)%carg .eq. '' )then
                ! value is real
                call real2str(self%cmds(i)%rarg, str)
                call hash%push(trim(self%cmds(i)%key), trim(str))
            else
                ! value is char
                call hash%push(trim(self%cmds(i)%key), trim(self%cmds(i)%carg))
            endif
        end do
        if( present(prg) )then
            call hash%set('prg', trim(prg))
        endif
    end subroutine gen_job_descr

end module simple_cmdline
