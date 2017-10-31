! Cartesian volume shift alignment based on band-pass limited cross-correlation
module simple_vol_shsrch
use simple_image,       only: image
use simple_opt_spec,    only: opt_spec
use simple_simplex_opt, only: simplex_opt
use simple_defs         ! use all in there
implicit none

public :: vol_shsrch_init, vol_shsrch_minimize
private

type(opt_spec)       :: ospec           !< optimizer specification object
type(simplex_opt)    :: nlopt           !< optimizer object
integer              :: nrestarts = 3   !< simplex restarts (randomized bounds)
real                 :: hp, lp, trs     !< corr/srch ctrl params
logical              :: shbarr = .true. !< shift barrier constraint or not
type(image), pointer :: vref  => null() !< pointer to reference
type(image), pointer :: vtarg => null() !< pointer to target

contains

    subroutine vol_shsrch_init( vol_ref, vol_target, hp_in, lp_in, trs_in, shbarrier, nrestarts_in )
        class(image),     target,   intent(in) :: vol_ref, vol_target
        real,                       intent(in) :: hp_in, lp_in, trs_in
        character(len=*), optional, intent(in) :: shbarrier
        integer,          optional, intent(in) :: nrestarts_in
        real :: lims(3,2)
        ! associate pointers
        vref  => vol_ref
        vtarg => vol_target
        ! set corr/srch ctrl params
        hp  = hp_in
        lp  = lp_in
        trs = trs_in
        ! flag the barrier constraint
        shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) shbarr = .false.
        endif
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! make optimizer spec
        lims(:,1) = -trs
        lims(:,2) =  trs
        call ospec%specify('simplex', 3, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=30)
        call ospec%set_costfun(vol_shsrch_costfun)
        ! generate the simplex optimizer object 
        call nlopt%new(ospec)
    end subroutine vol_shsrch_init

    function vol_shsrch_costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: vec_here(3), cost
        vec_here = vec
        where( abs(vec) < 1.e-6 ) vec_here = 0.
        if( shbarr )then
            if( any(vec_here(:) < ospec%limits(:,1)) .or.&
               &any(vec_here(:) > ospec%limits(:,2)) )then
                cost = 1.
                return
            endif
        endif
        cost = -vref%corr_shifted(vtarg, vec, lp_dyn=lp, hp_dyn=hp) 
    end function vol_shsrch_costfun

    function vol_shsrch_minimize( ) result( cxyz )
        real              :: cost_init, cost, cxyz(4)
        class(*), pointer :: fun_self => null()
        ospec%x = 0.0
        cost_init = vol_shsrch_costfun(fun_self, ospec%x, ospec%ndim)
        call nlopt%minimize(ospec, fun_self, cost)
        if( cost < cost_init )then
            cxyz(1)  = -cost   ! correlation
            cxyz(2:) = ospec%x ! shift
        else
            cxyz(1)  = -cost_init ! correlation
            cxyz(2:) = 0.0        ! shift
        endif
    end function vol_shsrch_minimize

end module simple_vol_shsrch
