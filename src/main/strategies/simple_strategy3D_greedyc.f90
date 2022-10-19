! concrete strategy3D: continuous stochastic neighborhood refinement
module simple_strategy3D_greedyc
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

public :: strategy3D_greedyc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_greedyc
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_greedyc
    procedure          :: srch        => srch_greedyc
    procedure          :: oris_assign => oris_assign_greedyc
    procedure          :: kill        => kill_greedyc
end type strategy3D_greedyc

contains

    subroutine new_greedyc( self, spec )
        class(strategy3D_greedyc), intent(inout) :: self
        class(strategy3D_spec),   intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_greedyc

    subroutine srch_greedyc( self, ithr )
        class(strategy3D_greedyc), intent(inout) :: self
        integer,                  intent(in)    :: ithr
        integer   :: iproj, irot
        type(ori) :: o, osym, obest
        real      :: cc(3), euldist, dist_inpl, euls(3), cc_best(3)
        ! continuous sochastic search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4_cont_srch
            o     = self%s%o_prev
            obest = self%s%o_prev
            ! zero shifts because particle is shifted to its previous origin
            call o%set('x', 0.)
            call o%set('y', 0.)
            cc_best = -huge(cc(1))
            ! init counter
            self%s%nrefs_eval = 0
            do iproj=1,self%s%nprojs
                euls = build_glob%eulspace%get_euler(iproj)
                do irot = 1,self%s%nrots
                    euls(3) = build_glob%inpl_rots(irot)
                    call o%set_euler(euls)
                    ! calculate Cartesian corr
                    cc = cartftcc_glob%project_and_correlate(self%s%iptcl, o)
                    if( cc(1) > cc_best(1) )then
                        cc_best = cc
                        obest = o
                    endif
                end do
            end do
            self%s%nrefs_eval = self%s%nprojs * self%s%nrots

            !!!!!!!!!!!!!!!!!!! CARTESIAN SHIFT SRCH TO BE IMPLEMENTED

            call build_glob%pgrpsyms%sym_dists(self%s%o_prev, obest, osym, euldist, dist_inpl)
            call obest%set('dist',      euldist)
            call obest%set('dist_inpl', dist_inpl)
            call obest%set('corr',      norm_corr(cc))
            call obest%set('cc_unnorm', cc(1))
            call obest%set('frac',      100.0)
            call build_glob%spproj_field%set_ori(self%s%iptcl, obest)
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_greedyc

    subroutine oris_assign_greedyc( self )
        class(strategy3D_greedyc), intent(inout) :: self
    end subroutine oris_assign_greedyc

    subroutine kill_greedyc( self )
        class(strategy3D_greedyc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_greedyc

end module simple_strategy3D_greedyc
