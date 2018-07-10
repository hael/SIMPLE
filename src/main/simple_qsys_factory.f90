! batch-processing manager - Factory class
module simple_qsys_factory
use simple_qsys_base,  only: qsys_base
use simple_qsys_local, only: qsys_local
use simple_qsys_slurm, only: qsys_slurm
use simple_qsys_sge,   only: qsys_sge
use simple_qsys_pbs,   only: qsys_pbs
implicit none

public :: qsys_factory
private

type :: qsys_factory
    private
    class(qsys_base), allocatable :: qsys_base_type
  contains
    procedure :: new
    procedure :: kill
end type qsys_factory

contains

    !> \brief  is a constructor
    subroutine new( self, which, ptr )
        class(qsys_factory), target, intent(inout) :: self  !< instance
        character(len=*),            intent(in)    :: which !< which qsys
        class(qsys_base),    pointer               :: ptr   !< pointer to constructed object
        ptr => null()
        call self%kill
        select case(which)
            case('local')
                allocate(qsys_local :: self%qsys_base_type)
            case('slurm')
                allocate(qsys_slurm :: self%qsys_base_type)
            case('sge')
                allocate(qsys_sge   :: self%qsys_base_type)
            case('pbs')
                allocate(qsys_pbs   :: self%qsys_base_type)
            case DEFAULT
                write(*,*) 'class:', which
                stop 'unsupported in qsys_factory constructor'
        end select
        call self%qsys_base_type%new
        ptr => self%qsys_base_type
    end subroutine new

    !> \brief  is a destructor
    subroutine kill( self )
        class(qsys_factory), intent(inout) :: self
        if( allocated(self%qsys_base_type) )then
            call self%qsys_base_type%kill
            deallocate(self%qsys_base_type)
        endif
    end subroutine kill

end module simple_qsys_factory
