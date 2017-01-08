module simple_generic_exec
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_jiffys          ! singleton
use simple_defs            ! singleton
implicit none

contains

    subroutine restarted_exec( cmder, cline, dirbody, nrestarts )
        use simple_syscalls, only: sys_mkdir
        class(commander_base), intent(inout) :: cmder
        class(cmdline),        intent(inout) :: cline
        character(len=*),      intent(in)    :: dirbody
        integer,               intent(in)    :: nrestarts
        character(len=STDLEN), allocatable   :: dirnames(:)
        integer :: iexec
        dirnames = make_numbered_names(dirbody, nrestarts)
        do iexec=1,nrestarts
            call sys_mkdir(trim(dirnames(iexec)))
            call chdir(trim(dirnames(iexec)))
            call cmder%execute(cline)
            call chdir('../')
        end do
        deallocate(dirnames)
    end subroutine restarted_exec

end module simple_generic_exec
