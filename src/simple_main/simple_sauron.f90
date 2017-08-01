!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> SAURON: SIMPLE Attempt to a Unified Resources and Orientations Notebook
!> is a module that provides support for one-line per one particle input/output for the SIMPLE suite
module simple_sauron
use simple_strings
implicit none

contains

    subroutine sauron_line_parser( line, htab, chtab )
        use simple_hash,    only: hash
        use simple_chash,   only: chash
        character(len=*), intent(inout)    :: line
        class(hash),      intent(inout)    :: htab
        class(chash),     intent(inout)    :: chtab
        character(len=32),     allocatable :: keys(:)
        character(len=STDLEN), allocatable :: vals(:)
        character(len=:),      allocatable :: line_trimmed
        character(len=STDLEN) :: args(100), args_pair(5), format
        integer :: nargs, iarg, nargs_pair, ival, io_stat
        real    :: rval
        allocate(line_trimmed, source=trim(line))
        call parse(line_trimmed,' ',args,nargs)
        allocate(keys(nargs), vals(nargs))
        do iarg=1,nargs
            call parse(args(iarg),'=',args_pair,nargs_pair)
            if( nargs_pair > 2 ) write(*,'(a)')&
            &'WARNING! nr of args in key-val pair > 2; simple_strings :: simple_line_parser'
            if( nargs_pair < 1 ) write(*,'(a)')&
            &'WARNING! nr of args in key-val pair < 1; simple_strings :: simple_line_parser'
            keys(iarg) = args_pair(1)
            vals(iarg) = args_pair(2)
            select case(str2format(vals(iarg)))
                case( 'file' )
                    call chtab%set(trim(keys(iarg)), trim(vals(iarg)))
                case( 'real' )
                    rval = str2real(trim(vals(iarg)))
                    call htab%set(trim(keys(iarg)), rval)
                case( 'int'  )
                    call str2int(trim(vals(iarg)), io_stat, ival)
                    rval = real(ival)
                    call htab%set(trim(keys(iarg)), rval)
                case( 'char' )
                    call chtab%set(trim(keys(iarg)), trim(vals(iarg)))
            end select
        end do
        deallocate(keys, vals, line_trimmed)
    end subroutine sauron_line_parser

end module simple_sauron
