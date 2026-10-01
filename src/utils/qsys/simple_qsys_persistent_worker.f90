!@descr: SIMPLE worker backend submission/environment manager module
! Writes job-description keys as header directives into per-job scripts. The TCP server lives in
! simple_persistent_worker_server; its shared handle persistent_worker%server is allocated by
! simple_qsys_env and never touched here.
module simple_qsys_persistent_worker
    use simple_core_module_api
    use simple_persistent_worker_server, only: persistent_worker_server
    use simple_qsys_base,                only: qsys_base
    implicit none

    public :: qsys_persistent_worker
    private

    ! ------------------------------------------------------------------
    ! Module-level constants
    ! ------------------------------------------------------------------
    integer, parameter :: MAXENVITEMS = 100   !< chash capacity for env vars

    ! ------------------------------------------------------------------
    ! qsys_persistent_worker type definition
    ! ------------------------------------------------------------------

    !> Submission/environment manager for the worker backend, extends qsys_base.
    !> Holds a key-value store mapping job-description keys to script header
    !> directives.  Initialise with new(); write headers with write_instr().
    type, extends(qsys_base) :: qsys_persistent_worker
        private
        type(chash) :: env !< submission environment key-value store
    contains
        procedure :: new               => new_worker_env
        procedure :: submit_cmd        => get_worker_submit_cmd
        procedure :: write_instr       => write_worker_header
        procedure :: write_array_instr => write_worker_array_header
        procedure :: kill              => kill_worker_env
    end type qsys_persistent_worker

contains

    ! ------------------------------------------------------------------
    ! qsys_persistent_worker lifecycle
    ! ------------------------------------------------------------------

    !> Constructor: initialise the environment key-value store and create
    !> the standard error/output staging directory.
    subroutine new_worker_env( self )
        class(qsys_persistent_worker), intent(inout) :: self
        character(len=STDLEN) :: stderrout
        call self%env%new(MAXENVITEMS)
        stderrout = PATH_HERE // trim(STDERROUT_DIR)
        call simple_mkdir(trim(stderrout))
    end subroutine new_worker_env

    !> Return the submission command string stored in the environment hash.
    function get_worker_submit_cmd( self ) result( cmd )
        class(qsys_persistent_worker), intent(in) :: self
        type(string) :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_worker_submit_cmd

    ! ------------------------------------------------------------------
    ! Script header I/O
    ! ------------------------------------------------------------------

    !> Write per-job script header lines derived from the job description hash.
    !> Writes to \p fhandle when present, otherwise to logfhandle.
    subroutine write_worker_header( self, q_descr, fhandle )
        class(qsys_persistent_worker), intent(in) :: self
        class(chash),       intent(in) :: q_descr
        integer, optional,  intent(in) :: fhandle
        type(string) :: key, sbatch_cmd, sbatch_val
        integer      :: i, which
        logical      :: write2file
        write2file = present(fhandle)
        do i = 1, q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key%to_char())
            if( which > 0 ) then
                sbatch_cmd = self%env%get(which)
                sbatch_val = q_descr%get(i)
                if( write2file ) then
                    write(fhandle,   '(a)') sbatch_cmd%to_char() // '=' // sbatch_val%to_char()
                else
                    write(logfhandle,'(a)') sbatch_cmd%to_char() // '=' // sbatch_val%to_char()
                end if
                call sbatch_cmd%kill()
                call sbatch_val%kill()
            end if
            call key%kill()
        end do
        if( write2file ) then
            write(fhandle,   '(a)') '#WRITTEN FOR EXECUTION BY SIMPLE_WORKER'
        else
            write(logfhandle,'(a)') '#WRITTEN FOR EXECUTION BY SIMPLE_WORKER'
        end if
    end subroutine write_worker_header

    !> Write per-array-job script header lines (qsys_base interface requirement).
    !> The implementation is identical to write_worker_header; array-specific
    !> parameters \p parts_fromto and \p nactive are required by the base
    !> interface but are not used by this backend.
    subroutine write_worker_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_persistent_worker), intent(in) :: self
        class(chash),       intent(in) :: q_descr
        integer,            intent(in) :: parts_fromto(2) !< reserved: array partition range
        integer, optional,  intent(in) :: fhandle
        integer, optional,  intent(in) :: nactive         !< reserved: active array element count
        type(string) :: key, sbatch_cmd, sbatch_val
        integer      :: i, which
        logical      :: write2file
        write2file = present(fhandle)
        do i = 1, q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key%to_char())
            if( which > 0 ) then
                sbatch_cmd = self%env%get(which)
                sbatch_val = q_descr%get(i)
                if( write2file ) then
                    write(fhandle,   '(a)') sbatch_cmd%to_char() // '=' // sbatch_val%to_char()
                else
                    write(logfhandle,'(a)') sbatch_cmd%to_char() // '=' // sbatch_val%to_char()
                end if
                call sbatch_cmd%kill()
                call sbatch_val%kill()
            end if
            call key%kill()
        end do
        if( write2file ) then
            write(fhandle,   '(a)') '#WRITTEN FOR EXECUTION BY SIMPLE_WORKER'
        else
            write(logfhandle,'(a)') '#WRITTEN FOR EXECUTION BY SIMPLE_WORKER'
        end if
    end subroutine write_worker_array_header

    !> Destructor: release the environment key-value store.
    !> persistent_worker%server belongs to simple_qsys_env; do not tear it down here.
    subroutine kill_worker_env( self )
        class(qsys_persistent_worker), intent(inout) :: self
        call self%env%kill()
    end subroutine kill_worker_env

end module simple_qsys_persistent_worker
