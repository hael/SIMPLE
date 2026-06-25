!@descr: batch-processing manager - Fortran coarray launcher backend
module simple_qsys_coarray
use simple_core_module_api
use simple_qsys_base, only: qsys_base
implicit none

public :: qsys_coarray
private
#include "simple_local_flags.inc"

integer, parameter :: MAXENVITEMS = 100

type, extends(qsys_base) :: qsys_coarray
    private
    type(chash) :: env !< defines the coarray launcher environment
  contains
    procedure :: new               => new_coarray_env
    procedure :: submit_cmd        => get_coarray_submit_cmd
    procedure :: write_instr       => write_coarray_header
    procedure :: write_array_instr => write_coarray_array_header
    procedure :: kill              => kill_coarray_env
end type qsys_coarray

contains

    !> Constructor.  SIMPLE_COARRAY_SUBMIT_CMD may override the default
    !! OpenCoarrays launcher prefix, for example "cafrun -n" or "mpiexec -n".
    subroutine new_coarray_env( self )
        class(qsys_coarray), intent(inout) :: self
        character(len=XLONGSTRLEN) :: submit_cmd
        integer :: envlen, envstat
        call self%env%new(MAXENVITEMS)
        call get_environment_variable('SIMPLE_COARRAY_SUBMIT_CMD', value=submit_cmd, length=envlen, status=envstat)
        if( envstat == 0 .and. envlen > 0 )then
            call self%env%push('qsys_submit_cmd', trim(adjustl(submit_cmd(:envlen))))
        else
            call self%env%push('qsys_submit_cmd', 'cafrun -np')
        endif
    end subroutine new_coarray_env

    !> Return launcher prefix.  qsys_ctrl appends the image count and the
    !! simple_private_exec --coarray invocation.
    function get_coarray_submit_cmd( self ) result( cmd )
        class(qsys_coarray), intent(in) :: self
        type(string) :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_coarray_submit_cmd

    !> Coarray jobs are launched directly through simple_private_exec.  Generating
    !! scheduler script headers for this backend indicates a wrong execution path.
    subroutine write_coarray_header( self, q_descr, fhandle )
        class(qsys_coarray), intent(in) :: self
        class(chash),        intent(in) :: q_descr
        integer, optional,   intent(in) :: fhandle
        THROW_HARD('coarray backend does not write scheduler scripts; use simple_private_exec --coarray dispatch')
    end subroutine write_coarray_header

    !> Array scripts are not used by the coarray backend.
    subroutine write_coarray_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_coarray), intent(in) :: self
        class(chash),        intent(in) :: q_descr
        integer,             intent(in) :: parts_fromto(2)
        integer, optional,   intent(in) :: fhandle, nactive
        THROW_HARD('coarray backend does not write array scripts; use simple_private_exec --coarray dispatch')
    end subroutine write_coarray_array_header

    !> Destructor.
    subroutine kill_coarray_env( self )
        class(qsys_coarray), intent(inout) :: self
        call self%env%kill
    end subroutine kill_coarray_env

end module simple_qsys_coarray
