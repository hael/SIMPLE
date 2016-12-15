module simple_opt_factory
use simple_optimizer,          only: optimizer
use simple_opt_spec,           only: opt_spec
use simple_bfgs_opt,           only: bfgs_opt
use simple_powell_opt,         only: powell_opt
use simple_simplex_opt,        only: simplex_opt
use simple_oasis_opt,          only: oasis_opt
use simple_bforce_opt,         only: bforce_opt
use simple_particle_swarm_opt, only: particle_swarm_opt
use simple_de_opt,             only: de_opt
use simple_bforce_opt,         only: bforce_opt
implicit none

public :: opt_factory
private

type :: opt_factory
    private
    class(optimizer), allocatable :: optimizer_type
  contains
    procedure :: new                   
    procedure :: kill
end type

contains

    !> \brief  is a constructor
    subroutine new( self, spec, ptr )
        class(opt_factory), intent(inout), target :: self !< instance 
        class(opt_spec), intent(inout)            :: spec !< specification
        class(optimizer), pointer                 :: ptr  !< pointer to constructed object
        ptr => null()
        call self%kill
        select case(spec%str_opt)
            case('bfgs')
                allocate(bfgs_opt           :: self%optimizer_type)
            case('powell')
                allocate(powell_opt         :: self%optimizer_type)
            case('simplex')
                allocate(simplex_opt        :: self%optimizer_type)
            case('oasis')
                allocate(oasis_opt          :: self%optimizer_type)
            case('bforce')
                allocate(bforce_opt         :: self%optimizer_type)
            case('pso')
                allocate(particle_swarm_opt :: self%optimizer_type)
            case('de')
                allocate(de_opt             :: self%optimizer_type)
            case DEFAULT
                write(*,*) 'class:', spec%str_opt
                stop 'unsupported in opt_factory constructor'
        end select
        call self%optimizer_type%new(spec)
        ptr => self%optimizer_type
    end subroutine
    
    !> \brief  is a destructor
    subroutine kill( self )
        class(opt_factory), intent(inout) :: self
        if( allocated(self%optimizer_type) )then
            call self%optimizer_type%kill
            deallocate(self%optimizer_type)
        endif
    end subroutine

end module simple_opt_factory
