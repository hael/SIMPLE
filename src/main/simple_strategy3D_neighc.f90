! concrete strategy3D: continuous stochastic neighborhood refinement
module simple_strategy3D_neighc
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: params_glob
use simple_builder,         only: build_glob
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
use simple_cartft_corrcalc, only: cartftcc_glob
use simple_ori,             only: ori
implicit none

public :: strategy3D_neighc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_neighc
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_neighc
    procedure          :: srch        => srch_neighc
    procedure          :: oris_assign => oris_assign_neighc
    procedure          :: kill        => kill_neighc
end type strategy3D_neighc

contains

    subroutine new_neighc( self, spec )
        class(strategy3D_neighc), intent(inout) :: self
        class(strategy3D_spec),   intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_neighc

    subroutine srch_neighc( self, ithr )
        class(strategy3D_neighc), intent(inout) :: self
        integer,                  intent(in)    :: ithr
        integer   :: isample
        type(ori) :: o, osym
        real      :: corr, euldist, dist_inpl
        ! continuous sochastic search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4_cont_srch
            o = self%s%o_prev
            ! zero shifts because particle is shifted to its previous origin
            call o%set('x', 0.)
            call o%set('y', 0.)
            ! init counter
            self%s%nrefs_eval = 0
            do isample=1,self%s%nsample
                ! make a random rotation matrix within the assymetric unit
                call build_glob%pgrpsyms%rnd_euler(self%s%o_prev, self%s%athres, o)
                ! calculate Cartesian corr
                call cartftcc_glob%project_and_correlate(self%s%iptcl, o, corr)
                ! keep track of how many references we are evaluating
                self%s%nrefs_eval = self%s%nrefs_eval + 1
                ! exit condition
                if( corr > self%s%prev_corr )then
                    call build_glob%pgrpsyms%sym_dists(self%s%o_prev, o, osym, euldist, dist_inpl)
                    call o%set('dist',      euldist)
                    call o%set('dist_inpl', dist_inpl)
                    call o%set('corr',      corr)
                    call build_glob%spproj_field%set_ori(self%s%iptcl, o)
                    exit
                endif
            end do
            !!!!!!!!!!!!!!!!!!! CARTESIAN SHIFT SRCH TO BE IMPLEMENTED
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_neighc

    subroutine oris_assign_neighc( self )
        class(strategy3D_neighc), intent(inout) :: self
    end subroutine oris_assign_neighc

    subroutine kill_neighc( self )
        class(strategy3D_neighc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_neighc

end module simple_strategy3D_neighc
