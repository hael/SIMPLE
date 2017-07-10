!> Simple system call module
!
!! simple_syscalls is a little module for using basic system functions e.g. CPU_time.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-10-01.
!
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-10-01
!
module simple_syscalls
use simple_jiffys       ! use all in there
use simple_defs         ! use all in there
use simple_filehandling ! use all in there
implicit none

private :: raise_sys_error

contains

    !> is the fortran 90 variant of the classic dtime
    real function dtime( time )
        real                  :: time(2)
        double precision,save :: last_time = 0
        double precision      :: this_time
        intrinsic cpu_time
        call cpu_time(this_time)
        time(1) = real(this_time-last_time)
        time(2) = 0.
        dtime = time(1)
        last_time = this_time
    end function dtime

    !> is the fortran 90 variant of the classic etime
    real function etime( time )
        real :: time(2)
        call cpu_time(etime)
        time(1) = etime
        time(2) = 0
    end function etime

    !> is for getting the actual cputime
    function getabscpu( lprint ) result( actual )
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: actual
        actual = etime( tarray )
        if( lprint )then
            write(*,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        endif
    end function getabscpu

    !> is for getting the relative cpu-time
    function getdiffcpu( lprint ) result( delta )
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: delta
        delta = dtime( tarray )
        if( lprint )then
            write(*,'(A,F9.2)') 'Relative cpu-time:', delta
        endif
    end function getdiffcpu

    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline, wait )
#if defined(INTEL)
        use ifport
#endif
        character(len=*),  intent(in) :: cmdline
        logical, optional, intent(in) :: wait
        character(len=STDLEN) :: cmsg
        integer               :: estat, cstat, exec_stat
        logical               :: doprint = .true., wwait = .true.
#if defined(PGI)
        call system(trim(adjustl(cmdline)))
#elif defined(INTEL)
        exec_stat = system(trim(adjustl(cmdline)))
        if( doprint )then
            write(*,*) 'command: ', trim(adjustl(cmdline))
            write(*,*) 'status of execution: ', exec_stat
        endif
#else
        wwait = .true.
        if( present(wait) ) wwait = wait
        call execute_command_line( trim(adjustl(cmdline)), wait=wwait, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
        call raise_sys_error( cmdline, estat, cstat, cmsg )
#endif
    end subroutine exec_cmdline

    !>  Handles error from system call
    subroutine raise_sys_error( cmd, exitstat, cmdstat, cmdmsg )
        integer,               intent(in) :: exitstat, cmdstat
        character(len=*),      intent(in) :: cmd
        character(len=STDLEN), intent(in) :: cmdmsg
        logical :: err
        err = .false.
        if( exitstat /= 0 )then
            write(*,*)'System error', exitstat,' for command: ', trim(adjustl(cmd))
            err = .true.
        endif
        if( cmdstat /= 0 )then
            write(*,*)'cmdstat = ',cmdstat,' command could not be executed: ', trim(adjustl(cmd))
            err = .true.
        endif
        ! if( err ) write(*,*) 'cmdmsg',trim(adjustl(cmdmsg))
    end subroutine raise_sys_error

    function sys_get_env_var( name ) result( varval )
        character(len=*), intent(in)  :: name
        character(len=STDLEN)         :: value
        character(len=:), allocatable :: varval
        integer :: length, status
        call get_environment_variable( trim(name), value=value, length=length, status=status)
        if( status == -1 ) write(*,*) 'value string too short; simple_syscalls :: get_env_var'
        if( status ==  1 ) write(*,*) 'environment variable: ', trim(name), ' is not defined; simple_syscalls :: get_env_var'
        if( status ==  2 ) write(*,*) 'environment variables not supported by system; simple_syscalls :: get_env_var'
        if( length ==  0 .or. status /= 0 ) return
        allocate(varval, source=trim(value))
    end function sys_get_env_var

    subroutine sys_gen_mrcfiletab( dir, filetabname )
        character(len=*),intent(in)  :: dir, filetabname
        character(len=STDLEN) :: cmd
        cmd = 'ls -tr '//trim(dir)//'/*.mrc*'//' > '//trim(filetabname)
        call exec_cmdline(cmd)
    end subroutine sys_gen_mrcfiletab

    subroutine sys_gen_filetab( fbody, ext, filetabname )
        character(len=*), intent(in)  :: fbody, ext, filetabname
        character(len=STDLEN) :: cmd
        cmd = 'ls -tr '//trim(fbody)//'*'//trim(ext)//' > '//trim(filetabname)
        call exec_cmdline(cmd)
    end subroutine sys_gen_filetab

    subroutine sys_del_files( fbody, ext )
        character(len=*),      intent(in)  :: fbody, ext
        character(len=STDLEN), allocatable :: fnames(:)
        character(len=STDLEN), parameter   :: FTAB = 'ftab_from_sys_del_files.txt'
        integer :: i, last
        call sys_gen_filetab(fbody, ext, FTAB) ! filetable written to disc
        call read_filetable(FTAB, fnames)      ! filetable read back in
        last = size(fnames)
        do i=1,last
            call del_file(fnames(i))
        end do
        call del_file(FTAB)
        deallocate(fnames)
    end subroutine sys_del_files

    function sys_get_last_fname( fbody, ext ) result( fname )
        character(len=*),      intent(in)  :: fbody, ext
        character(len=STDLEN), allocatable :: fnames(:)
        character(len=STDLEN), parameter   :: FTAB = 'ftab_from_sys_find_last_fname.txt'
        character(len=STDLEN) :: fname
        integer :: last
        call sys_gen_filetab(fbody, ext, FTAB) ! filetable written to disc
        call read_filetable(FTAB, fnames)      ! filetable read back in
        last = size(fnames)
        fname = fnames(last)
        call del_file(FTAB)
        deallocate(fnames)
    end function sys_get_last_fname

    subroutine sys_merge_docs( docnames, fname_merged )
        character(len=STDLEN), intent(in) :: docnames(:)
        character(len=*),      intent(in) :: fname_merged
        character(len=STDLEN) :: cmd
        integer :: ndocs, idoc
        call del_file(fname_merged)
        ndocs = size(docnames)
        do idoc=1,ndocs
            cmd = 'cat '//trim(docnames(idoc))//' >> '//trim(fname_merged)
            call exec_cmdline(cmd)
        end do
    end subroutine sys_merge_docs

    subroutine simple_sleep( secs )
#if defined(INTEL)
        use ifport
#endif
        integer, intent(in) :: secs
        integer  :: msecs
        msecs = 1000*secs
#if defined(INTEL)
        call sleepqq(msecs)
#else
        call sleep(secs)
#endif
    end subroutine simple_sleep

    !>  Wrapper for system call
    subroutine sys_stat( filename, buffer, status )
#if defined(INTEL)
        use ifport
#endif
        character(len=*),  intent(in) :: filename
        integer,dimension(13), intent(inout) :: buffer
        integer, intent(inout) :: status
        character(len=STDLEN) :: cmsg
        integer               :: estat, cstat, exec_stat
        logical               :: doprint = .true.
        logical exists
#if defined(INTEL)
        cmsg = ' failed to find '//trim(adjustl(filename))
        inquire(FILE = trim(adjustl(filename)), IOMSG=cmsg, EXIST = exists )
        if(exists)then
            status = 0
        else
            status = -1
        end if
#else
        call stat(trim(adjustl(filename)), buffer, status)
        if( doprint )then
            write(*,*) 'command: stat ', trim(adjustl(filename))
            write(*,*) 'status of execution: ', status
        endif
        status = exec_stat
#endif
    end subroutine sys_stat

    

    function waitforfileclose( funit )result( all_good )
        integer, intent(inout) :: funit
        integer :: rec_prev, rec, iostat, waits, maxwaits = 60
        logical :: all_good, opened
        all_good = .false.
        waits    = 0
        inquire(unit=funit, opened=opened, recl=rec_prev)
        do while(.not.all_good)
            call simple_sleep(1)
            inquire(unit=funit, opened=opened, recl=rec, iostat=iostat)
            all_good = .not.opened
            all_good = all_good .and. (rec == rec_prev)
            all_good = all_good .and. (iostat == 0)
            rec_prev = rec
            waits    = waits + 1
            if(waits > 10)then
                write(*,'(A,I3)')'This is taking a long time for file unit:', funit
                if(waits > maxwaits)then
                    write(*,'(A,I3)')'This is taking too long a time for file unit:', funit
                    stop
                endif
            endif
        enddo
    end function waitforfileclose

end module simple_syscalls
