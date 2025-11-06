! unix_wait.F90
!
! Author:  Philipp Engel
! Licence: ISC
module unix_wait
    use, intrinsic :: iso_c_binding
    use :: unix_types
    implicit none
    private

    public :: c_wait
    public :: c_waitpid
    
    interface
        ! pid_t wait(int *stat_loc)
        function c_wait(stat_loc) bind(c, name='wait')
            import :: c_int, c_pid_t
            implicit none
            integer(kind=c_int), intent(out) :: stat_loc
            integer(kind=c_pid_t)            :: c_wait
        end function c_wait

        ! pid_t waitpid(pid_t pid, int *stat_loc, int options)
        function c_waitpid(pid, stat_loc, options) bind(c, name='waitpid')
            import :: c_int, c_pid_t
            implicit none
            integer(kind=c_pid_t), intent(in), value :: pid
            integer(kind=c_int),   intent(in), value :: options
            integer(kind=c_int),   intent(out)       :: stat_loc
            integer(kind=c_pid_t)                    :: c_waitpid
        end function c_waitpid
    end interface
end module unix_wait
