!@descr: coarray launcher for SIMPLE distributed partition scripts
program simple_coarray_exec
use simple_defs,         only: XLONGSTRLEN
use simple_error,        only: simple_exception
use simple_string_utils, only: int2str, int2str_pad
implicit none
#include "simple_local_flags.inc"

character(len=XLONGSTRLEN)    :: script_prefix, script_name, arg, cmd, cmdmsg
character(len=:), allocatable :: errmsg
integer                       :: first_part, last_part, numlen
integer                       :: ipart
integer                       :: parse_stat, cmdstat, exitstat

call parse_arguments

do ipart = first_part + this_image() - 1, last_part, num_images()
    script_name = trim(script_prefix)//int2str_pad(ipart, numlen)
    cmd         = 'bash '//trim(script_name)
    cmdstat     = 0
    exitstat    = 0
    cmdmsg      = ''
    call execute_command_line(trim(cmd), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
    if( cmdstat /= 0 .or. exitstat /= 0 )then
        if( len_trim(cmdmsg) > 0 ) THROW_WARN(trim(cmdmsg))
        errmsg = 'simple_coarray_exec image '//int2str(this_image())//' failed part '//int2str(ipart)
        errmsg = errmsg//' exitstat='//int2str(exitstat)//' cmdstat='//int2str(cmdstat)
        THROW_HARD(errmsg)
    endif
end do

sync all

contains

    subroutine parse_arguments
        if( command_argument_count() /= 4 )then
            THROW_HARD('Usage: simple_coarray_exec <script_prefix> <first_part> <last_part> <numlen>')
        endif
        call get_command_argument(1, script_prefix)
        call parse_int_arg(2, first_part)
        call parse_int_arg(3, last_part)
        call parse_int_arg(4, numlen)
        if( first_part < 1 .or. last_part < first_part .or. numlen < 1 )then
            THROW_HARD('Invalid coarray partition range or numlen')
        endif
    end subroutine parse_arguments

    subroutine parse_int_arg( iarg, val )
        integer, intent(in)  :: iarg
        integer, intent(out) :: val
        call get_command_argument(iarg, arg)
        read(arg,*,iostat=parse_stat) val
        if( parse_stat /= 0 )then
            THROW_HARD('Invalid integer argument '//int2str(iarg)//': '//trim(arg))
        endif
    end subroutine parse_int_arg

end program simple_coarray_exec
