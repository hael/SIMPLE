! factory pattern class for the SIMPLE optimisers
module simple_opt_factory
use simple_optimizer,          only: optimizer
use simple_opt_spec,           only: opt_spec
use simple_bfgs_opt,           only: bfgs_opt
use simple_bfgs2_opt,          only: bfgs2_opt
use simple_fr_cg_opt,          only: fr_cg_opt
use simple_pr_cg_opt,          only: pr_cg_opt
use simple_stde_opt,           only: stde_opt
use simple_powell_opt,         only: powell_opt
use simple_simplex_opt,        only: simplex_opt
use simple_bforce_opt,         only: bforce_opt
use simple_particle_swarm_opt, only: particle_swarm_opt
use simple_de_opt,             only: de_opt
use simple_bforce_opt,         only: bforce_opt
implicit none

public :: opt_factory
private
!> abstract optimisation factory type
type :: opt_factory
    private
    class(optimizer), pointer :: optimizer_type
  contains
    procedure :: new
end type opt_factory

contains

    !> \brief  is a constructor
    subroutine new( self, spec, ptr )
        class(opt_factory), target, intent(inout) :: self !< instance
        class(opt_spec),            intent(inout) :: spec !< specification
        class(optimizer), pointer                 :: ptr  !< pointer to constructed object
        select case(spec%str_opt)
            case('bfgs')
                allocate(bfgs_opt           :: self%optimizer_type)
            case('powell')
                allocate(powell_opt         :: self%optimizer_type)
            case('simplex')
                allocate(simplex_opt        :: self%optimizer_type)
            case('bforce')
                allocate(bforce_opt         :: self%optimizer_type)
            case('pso')
                allocate(particle_swarm_opt :: self%optimizer_type)
            case('de')
                allocate(de_opt             :: self%optimizer_type)
            case('fr_cg')
                allocate(fr_cg_opt          :: self%optimizer_type)
            case('pr_cg')
                allocate(pr_cg_opt          :: self%optimizer_type)
            case('bfgs2')
                allocate(bfgs2_opt          :: self%optimizer_type)
            case('stde')
                allocate(stde_opt           :: self%optimizer_type)                
            case DEFAULT
                write(*,*) 'class:', spec%str_opt
                stop 'unsupported in opt_factory constructor'
        end select
        call self%optimizer_type%new(spec)
        ptr => self%optimizer_type
    end subroutine

end module simple_opt_factory
