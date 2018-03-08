! SAURON: SIMPLE Attempt to a Unified Resources and Orientations Notebook
! provides support for one-line per one particle input/output
module simple_sauron
include 'simple_lib.f08'
implicit none

contains

    subroutine sauron_line_parser( line, htab, chtab )
        character(len=*), intent(inout)    :: line
        class(hash),      intent(inout)    :: htab
        class(chash),     intent(inout)    :: chtab
        character(len=32),     allocatable :: keys(:)
        character(len=STDLEN), allocatable :: vals(:)
        character(len=:),      allocatable :: line_trimmed
        character(len=STDLEN) :: args(128), args_pair(2), format
        integer :: nargs, iarg, nargs_pair, ival, io_stat
        allocate(line_trimmed, source=trim(line))
        call parsestr(line_trimmed,' ',args,nargs)
        allocate(keys(nargs), vals(nargs), stat=alloc_stat)
        if(alloc_stat /= 0)call allocchk("simple_sauron::sauron_line_parser ",alloc_stat)
        do iarg=1,nargs
            args_pair(2) = args(iarg)
            call split_str(args_pair(2), '=', args_pair(1))
            if( index(args_pair(1),' ')==0 ) write(*,'(a,a)')&
                &'WARNING! Invalid key; simple_sauron :: sauron_line_parser:', trim(args_pair(1))
            if( index(args_pair(2),' ')==0 ) write(*,'(a)')&
                &'WARNING! Invalid arg; simple_sauron :: sauron_line_parser:', trim(args_pair(2))
            if(len(args_pair(1)) == 0 ) write(*,'(a)')&
                &'WARNING! Invalid key; simple_sauron :: sauron_line_parser'
            if(len(args_pair(2)) == 0 ) write(*,'(a)')&
                &'WARNING! Invalid arg; simple_sauron :: sauron_line_parser'
            keys(iarg) = args_pair(1)
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
