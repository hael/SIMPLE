! concrete strategy3D: 5-degree derivative-free optimization using NLOpt
module simple_strategy3D_nlopt
include 'simple_lib.f08'
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: params_glob
use simple_builder,         only: build_glob
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch, strategy3D_spec
use simple_cartft_corrcalc, only: cftcc_glob
use nlopt_wrap,             only: nlopt_opt, nlopt_func, create, destroy
use nlopt_enum,             only: NLOPT_SUCCESS, algorithm_from_string
implicit none

public :: strategy3D_nlopt
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_nlopt
    type(strategy3D_srch) :: s
    type(strategy3D_spec) :: spec
    type(nlopt_opt)       :: opt
contains
    procedure          :: new         => new_nlopt
    procedure          :: srch        => srch_nlopt
    procedure          :: oris_assign => oris_assign_nlopt
    procedure          :: kill        => kill_nlopt
end type strategy3D_nlopt

contains

    subroutine new_nlopt( self, spec )
        class(strategy3D_nlopt), intent(inout) :: self
        class(strategy3D_spec),  intent(inout) :: spec
        call self%s%new(spec)
        self%spec = spec
        call create(self%opt, algorithm_from_string('LN_COBYLA'), 5)
    end subroutine new_nlopt

    subroutine srch_nlopt( self, ithr )
        class(strategy3D_nlopt), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        integer,                 parameter     :: wp  = kind(0.0d0)
        real(wp),                parameter     :: TOL = 0.001_wp     ! tolerance for success
        type(ori) :: o, osym
        real      :: corr, euldist, dist_inpl, corr_inpl, shvec(2)
        real(wp)  :: x_nlopt(5), best_corr
        integer   :: stat
        if( build_glob%spproj_field%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4_cont_srch
            ! transfer critical per-particle params
            o    = self%s%o_prev
            corr = self%s%prev_corr
            ! zero shifts because particle is shifted to its previous origin
            shvec(1) = 0.
            shvec(2) = 0.
            associate(f => nlopt_func(nloptf_myfunc))
                call self%opt%set_min_objective(f)
                call self%opt%set_ftol_rel(TOL)
                x_nlopt(1) = self%s%o_prev%e1get()
                x_nlopt(2) = self%s%o_prev%e2get()
                x_nlopt(3) = self%s%o_prev%e3get()
                x_nlopt(4) = shvec(1)
                x_nlopt(5) = shvec(2)
                call self%opt%optimize(x_nlopt, best_corr, stat)
            end associate
            call o%e1set(real(x_nlopt(1)))
            call o%e2set(real(x_nlopt(2)))
            call o%e3set(real(x_nlopt(3)))
            shvec = x_nlopt(4:5)
            if( best_corr > corr )then
                call build_glob%pgrpsyms%sym_dists(self%s%o_prev, o, osym, euldist, dist_inpl)
                call o%set('dist',      euldist)
                call o%set('dist_inpl', dist_inpl)
                call o%set('corr',      real(best_corr))
                call o%set('frac',      100.0)
                call build_glob%spproj_field%set_ori(  self%s%iptcl, o)
                call build_glob%spproj_field%set_shift(self%s%iptcl, shvec)
            endif
        else
            call build_glob%spproj_field%reject(self%s%iptcl)
        endif
    contains
        function nloptf_myfunc(x_in, gradient, func_data) result(f)
            real(wp), intent(in)              :: x_in(:)
            real(wp), intent(inout), optional :: gradient(:)
            class(*), intent(in),    optional :: func_data
            real(wp)                          :: f
            real                              :: f_sp
            call o%e1set(real(x_in(1)))
            call o%e2set(real(x_in(2)))
            call o%e3set(real(x_in(3)))
            f_sp = cftcc_glob%project_and_correlate( self%s%iptcl, o, real(x_in(4:5)))
            f = f_sp
        end function nloptf_myfunc
    end subroutine srch_nlopt

    subroutine oris_assign_nlopt( self )
        class(strategy3D_nlopt), intent(inout) :: self
    end subroutine oris_assign_nlopt

    subroutine kill_nlopt( self )
        class(strategy3D_nlopt),   intent(inout) :: self
        call self%s%kill
        call destroy(self%opt)
    end subroutine kill_nlopt

end module simple_strategy3D_nlopt
