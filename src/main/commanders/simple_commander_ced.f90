! concrete commander: ced for coherence-enhancing diffusion filter
module simple_commander_ced

include 'simple_lib.f08'
use simple_parameters,     only: parameters, params_glob
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image

public :: ced_2D_filter_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: ced_2D_filter_commander
  contains
    procedure :: execute      => exec_ced_2D_filter
end type ced_2D_filter_commander

contains

    ! Applies the ced filter to an image
    subroutine exec_ced_2D_filter(self, cline)
        use simple_ced_filter, only: ced_filter_2D
        class(ced_2D_filter_commander), intent(inout)   :: self
        class(cmdline),                 intent(inout)   :: cline
        type(parameters)    :: params
        type(image)         :: img
        integer             :: nptcls, iptcl

        ! Setup
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        call img%new(params%ldim, params%smpd)

        ! Apply filter to each particle
        do iptcl = 1, nptcls
            write(*, *) 'Particle # ', iptcl
            call img%read(params%stk, iptcl)
            call ced_filter_2D(img, params%sigma)
            call img%write('ced_2D_filter.mrc', iptcl)  ! Might want to include sigma in file name
        enddo
    end subroutine exec_ced_2D_filter

end module simple_commander_ced