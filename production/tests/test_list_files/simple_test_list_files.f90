program simple_test_list_files
use simple_defs
implicit none
character(len=STDLEN), allocatable :: sp_files(:)
character(len=STDLEN) :: cwd
integer :: iostat_projfile, nsp_files, i
call simple_list_files('*.simple', sp_files)
if( allocated(sp_files) )then
    nsp_files = size(sp_files)
else
    nsp_files = 0
endif
if( nsp_files > 0 )then
    do i=1,nsp_files
        print *, trim(sp_files(i))
    end do
endif

contains

    subroutine simple_list_files( pattern, list )
        character(len=*), intent(in) :: pattern
        character(len=STDLEN) :: cmd, tmpfile
        character(len=STDLEN), allocatable, intent(inout) :: list(:)
        character(len=1) :: junk
        integer :: sz, funit, ios, i, nlines
        tmpfile = '__simple_filelist__'
        cmd     = 'ls '//trim(pattern)//' > '//trim(tmpfile)
        call exec_cmdline( cmd, suppress_errors=.true.)
        inquire(file=trim(tmpfile), size=sz)
        if( allocated(list) ) deallocate(list)
        if( sz > 0 )then
            open(newunit=funit, file=trim(tmpfile))
            nlines = 0
            do
                read(funit,*,iostat=ios) junk
                if(ios /= 0)then
                    exit
                else
                    nlines = nlines + 1
                endif
            end do
            rewind(funit)
            allocate( list(nlines) )
            do i=1,nlines
                read(funit, '(a)') list(i)
            enddo
            close(funit, status='delete')
        else
            open(newunit=funit, file=trim(tmpfile))
            close(funit, status='delete')
        endif
    end subroutine simple_list_files

    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline, waitflag, suppress_errors)
        character(len=*),  intent(in) :: cmdline
        logical, optional, intent(in) :: waitflag, suppress_errors
        character(len=:), allocatable :: cmdmsg, tmp,  cmsg
        character(11) :: suppress_msg='2>/dev/null'
        integer ::  cstat, exec_stat
        logical :: l_doprint, wwait, l_suppress_errors
        l_doprint = .false.
        wwait     = .true.
        if( present(waitflag)        ) wwait = waitflag
        if( present(suppress_errors) ) l_suppress_errors = suppress_errors
        allocate(cmsg, source=trim(adjustl(cmdline)))
        if( l_suppress_errors )then
            allocate(tmp, source=cmsg//' '//suppress_msg)
            cmsg = tmp
        endif
        call execute_command_line(cmsg, wait=wwait, exitstat=exec_stat, cmdstat=cstat, cmdmsg=cmdmsg)
        if( .not. l_suppress_errors ) call raise_sys_error( cmsg, exec_stat, cstat, cmdmsg )
        if( l_doprint )then
            write(logfhandle,*) 'command            : ', cmsg
            write(logfhandle,*) 'status of execution: ', exec_stat
        endif
    end subroutine exec_cmdline

    !>  Handles error from system call
    subroutine raise_sys_error( cmd, exitstat, cmdstat, cmdmsg )
        integer,               intent(in) :: exitstat, cmdstat
        character(len=*),      intent(in) :: cmd
        character(len=STDLEN), intent(in) :: cmdmsg
        logical :: err
        err = .false.
        if( exitstat /= 0 )then
            write(logfhandle,*)'System error', exitstat,' for command: ', trim(adjustl(cmd))
            err = .true.
        endif
        if( cmdstat /= 0 )then
            ! call simple_error_check(cmdstat,' command could not be executed: '//trim(adjustl(cmd)))
            write(logfhandle,*)'cmdstat = ',cmdstat,' command could not be executed: ', trim(adjustl(cmd))
            err = .true.
        endif
        ! if( err ) write(logfhandle,*) trim(adjustl(cmdmsg))
    end subroutine raise_sys_error

end program simple_test_list_files
