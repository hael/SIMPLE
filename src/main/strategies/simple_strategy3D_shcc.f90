! concrete strategy3D: continuous stochastic refinement
module simple_strategy3D_shcc
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: params_glob
use simple_builder,         only: build_glob
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
use simple_cartft_corrcalc, only: cartftcc_glob
implicit none

public :: strategy3D_shcc
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_shcc
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
contains
    procedure          :: new         => new_shcc
    procedure          :: srch        => srch_shcc
    procedure          :: oris_assign => oris_assign_shcc
    procedure          :: kill        => kill_shcc
end type strategy3D_shcc

contains

    subroutine new_shcc( self, spec )
        class(strategy3D_shcc), intent(inout) :: self
        class(strategy3D_spec), intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
    end subroutine new_shcc

    subroutine srch_shcc( self, ithr )
        class(strategy3D_shcc), intent(inout) :: self
        integer,                 intent(in)   :: ithr
        integer   :: isample, irot
        type(ori) :: o, osym, obest
        real      :: corr, euldist, dist_inpl, corr_best, frac, corr_inpl, e3
        logical   :: got_better
        ! continuous sochastic search
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4_cont_srch
            ! transfer critical per-particle params
            o = self%s%o_prev
            obest = self%s%o_prev
            ! zero shifts because particle is shifted to its previous origin
            call o%set('x', 0.)
            call o%set('y', 0.)
            ! currently the best correlation is the previous one
            corr_best  = self%s%prev_corr
            got_better = .false.
            do isample=1,self%s%nsample
                ! make a random rotation matrix within the assymetric unit
                call build_glob%pgrpsyms%rnd_euler(o)

                ! greedy optimization over in-plane angle
                ! corr_inpl = -1.
                ! do irot = 1,self%s%nrots
                !     call o%e3set(build_glob%inpl_rots(irot))
                !     ! calculate Cartesian corr
                !     corr = cartftcc_glob%project_and_correlate(self%s%iptcl, o)
                !     if( corr > corr_inpl )then
                !         corr_inpl = corr
                !         e3 = build_glob%inpl_rots(irot)
                !     endif
                ! end do
                ! call o%e3set(e3)

                corr = cartftcc_glob%project_and_correlate(self%s%iptcl, o)
                if( corr > corr_best )then
                    corr_best  = corr
                    obest      = o
                    got_better = .true.
                endif
            end do
            if( got_better )then
                call build_glob%pgrpsyms%sym_dists(self%s%o_prev, obest, osym, euldist, dist_inpl)
                call obest%set('dist',      euldist)
                call obest%set('dist_inpl', dist_inpl)
                call obest%set('corr',      corr_best)
                call obest%set('frac',      100.0)
                call build_glob%spproj_field%set_ori(self%s%iptcl, obest)
            endif
            ! local refinement step
            got_better = .false.
            do isample=1,self%s%nsample
                ! make a random rotation matrix neighboring the previous best within the assymetric unit
                call build_glob%pgrpsyms%rnd_euler(obest, self%s%athres, o)
                ! calculate Cartesian corr
                corr = cartftcc_glob%project_and_correlate(self%s%iptcl, o)
                if( corr > corr_best )then
                    corr_best  = corr
                    obest      = o
                    got_better = .true.
                endif
            end do
            if( got_better )then
                call build_glob%pgrpsyms%sym_dists(self%s%o_prev, obest, osym, euldist, dist_inpl)
                call obest%set('dist',      euldist)
                call obest%set('dist_inpl', dist_inpl)
                call obest%set('corr',      corr_best)
                call obest%set('frac',      100.0)
                call build_glob%spproj_field%set_ori(self%s%iptcl, obest)
            endif


            !!!!!!!!!!!!!!!!!!! CARTESIAN SHIFT SRCH TO BE IMPLEMENTED

        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    end subroutine srch_shcc

    subroutine oris_assign_shcc( self )
        class(strategy3D_shcc), intent(inout) :: self
    end subroutine oris_assign_shcc

    subroutine kill_shcc( self )
        class(strategy3D_shcc),   intent(inout) :: self
        call self%s%kill
    end subroutine kill_shcc

end module simple_strategy3D_shcc
