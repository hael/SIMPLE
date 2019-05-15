! factory pattern class for the SIMPLE optimisers
module simple_opt_factory
include 'simple_lib.f08'
use simple_optimizer,   only: optimizer
use simple_opt_spec,    only: opt_spec
use simple_opt_lbfgsb,  only: opt_lbfgsb
use simple_opt_bfgs2,   only: opt_bfgs2
use simple_opt_simplex, only: opt_simplex
use simple_opt_bforce,  only: opt_bforce
use simple_opt_de,      only: opt_de
use simple_opt_stde,    only: opt_stde
implicit none

public :: opt_factory
private
#include "simple_local_flags.inc"

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
            case('simplex')
                allocate(opt_simplex        :: self%optimizer_type)
            case('bforce')
                allocate(opt_bforce         :: self%optimizer_type)
            case('de')
                allocate(opt_de             :: self%optimizer_type)
            case('lbfgsb')
                allocate(opt_lbfgsb         :: self%optimizer_type)
            case('bfgs')
                allocate(opt_bfgs2          :: self%optimizer_type)
            case('stde')
                allocate(opt_stde           :: self%optimizer_type)
            case DEFAULT
                THROW_HARD('class: '//trim(spec%str_opt)//' unsupported in opt_factory constructor')
        end select
        call self%optimizer_type%new(spec)
        ptr => self%optimizer_type
    end subroutine

end module simple_opt_factory
