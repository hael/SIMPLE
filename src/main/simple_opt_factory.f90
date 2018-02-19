! factory pattern class for the SIMPLE optimisers
module simple_opt_factory
use simple_optimizer,          only: optimizer
use simple_opt_spec,           only: opt_spec
use simple_opt_bfgs,           only: opt_bfgs
use simple_opt_bfgs2,          only: opt_bfgs2
use simple_opt_lbfgsb,         only: opt_lbfgsb
use simple_opt_fr_cg,          only: opt_fr_cg
use simple_opt_pr_cg,          only: opt_pr_cg
use simple_opt_stde,           only: opt_stde
use simple_opt_powell,         only: opt_powell
use simple_opt_simplex,        only: opt_simplex
use simple_opt_bforce,         only: opt_bforce
use simple_opt_particle_swarm, only: opt_particle_swarm
use simple_opt_de,             only: opt_de
use simple_opt_bforce,         only: opt_bforce
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
                allocate(opt_bfgs           :: self%optimizer_type)
            case('powell')
                allocate(opt_powell         :: self%optimizer_type)
            case('simplex')
                allocate(opt_simplex        :: self%optimizer_type)
            case('bforce')
                allocate(opt_bforce         :: self%optimizer_type)
            case('pso')
                allocate(opt_particle_swarm :: self%optimizer_type)
            case('de')
                allocate(opt_de             :: self%optimizer_type)
            case('fr_cg')
                allocate(opt_fr_cg          :: self%optimizer_type)
            case('pr_cg')
                allocate(opt_pr_cg          :: self%optimizer_type)
            case('bfgs2')
                allocate(opt_bfgs2          :: self%optimizer_type)
            case('lbfgsb')
                allocate(opt_lbfgsb         :: self%optimizer_type)
            case('stde')
                allocate(opt_stde           :: self%optimizer_type)                
            case DEFAULT
                write(*,*) 'class:', spec%str_opt
                stop 'unsupported in opt_factory constructor'
        end select
        call self%optimizer_type%new(spec)
        ptr => self%optimizer_type
    end subroutine

end module simple_opt_factory
