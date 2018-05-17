! SAURON: SIMPLE Attempt to a Unified Resources and Orientations Notebook
! provides support for one-line per one particle input/output
module simple_sauron
use simple_defs
use simple_strings, only: str2format, str2real, str2int, parsestr, split_str
use simple_error,   only: allocchk
use simple_hash,    only: hash
use simple_chash,   only: chash
implicit none

contains

    subroutine sauron_line_parser( line, htab, chtab )
        character(len=*), intent(inout)    :: line
        class(hash),      intent(inout)    :: htab
        class(chash),     intent(inout)    :: chtab
        character(len=32),     allocatable :: keys(:)
        character(len=STDLEN), allocatable :: vals(:)
        character(len=:),      allocatable :: line_trimmed
        character(len=STDLEN) :: args(128), args_pair(2)
        integer :: nargs, iarg, ival, io_stat
        allocate(line_trimmed, source=trim(line))
        call parsestr(line_trimmed,' ', args, nargs)
        allocate(keys(nargs), vals(nargs), stat=alloc_stat)
        if(alloc_stat /= 0)call allocchk("simple_sauron::sauron_line_parser ",alloc_stat)
        do iarg=1,nargs
            args_pair(2) = args(iarg)
            call split_str(args_pair(2), '=', args_pair(1))
            keys(iarg) = trim(args_pair(1))
            vals(iarg) = args_pair(2)
            select case(str2format(vals(iarg)))
                case( 'real' )
                    call htab%set(keys(iarg), str2real(trim(vals(iarg))))
                case( 'file', 'dir', 'char' )
                    call chtab%set(trim(keys(iarg)), trim(vals(iarg)))
                case( 'int'  )
                    call str2int(trim(vals(iarg)), io_stat, ival)
                    call htab%set(keys(iarg), real(ival))
            end select
        end do
        deallocate(keys, vals, line_trimmed)
    end subroutine sauron_line_parser

end module simple_sauron
